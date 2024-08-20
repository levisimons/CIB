rm(list=ls())
require(dplyr)
require(plyr)
require(tidyr)
require(phyloseq)
require(decontam)
require(data.table)
require(stringr)
require(biomformat)

wd <- ""

setwd(wd)

#Set lists for filtering taxa.
Invertebrate_Phyla <- c("Platyhelminthes","Nemertea","Rotifera","Gastrotricha","Acanthocephala","Nematoda","Nematomorpha","Priapulida","Kinorhyncha","Loricifera","Entoprocta","Cycliophora","Gnathostomulida","Micrognathozoa","Chaetognatha","Hemichordata","Bryozoa","Brachiopoda","Phoronida","Annelida","Mollusca","Arthropoda")
Fungal_Phyla <- c("Ascomycota", "Basidiomycota", "Entorrhizomycetes", "Blastocladiomycota", "Chytridiomycota", "Cryptomycota", "Microsporidia", "Mucoromycota", "Nephridiophaga", "Olpidiomycota", "Sanchytriomycota", "Zoopagomycota")

#Select a primer
Primers <- c("12S_MiFish_U","16S_Bacteria","18S_Euk","CO1_Metazoa","ITS1_Fungi","ITS2_Plants","vert12S")
Primer <- "ITS1_Fungi"

##Read in and process Tronko-assign eDNA data.
#Get project directories
Project_Directories <- list.dirs(path=paste(wd,"CALeDNA",sep="/"),recursive=F)

#Aggregate tronko-assign results to a particular taxonomic level.
TaxonomicRanks <- c("superkingdom","phylum","class","order","family","genus","species")
project_stats <- c()
#Select a tronko-assign mismatch cutoff.
Mismatch <- 25
i=1
for(Project_Directory in Project_Directories){
  ProjectID <- basename(Project_Directory)
  #Get processed project metadata
  Metadata_Processed_Path <- paste(wd,"CALeDNA",ProjectID,"terradactyl","metabarcoding_metadata_terradactyl.csv",sep="/")
  Metadata_Processed <- read.table(Metadata_Processed_Path, header=TRUE, sep=",",as.is=T,skip=0,fill=TRUE,check.names=FALSE,quote = "\"", encoding = "UTF-8",na.strings=c("n/a","NA"))
  #Store all extracted metadata.
  Metadata_Processed$projectid <- ProjectID
  #Get uploaded project metadata.
  Metadata_Path <- paste(wd,"CALeDNA",ProjectID,"terradactyl","metabarcoding_metadata_original.csv",sep="/")
  MetadataTable <- read.table(Metadata_Path, header=TRUE, sep=",",as.is=T,skip=0,fill=TRUE,check.names=FALSE,quote = "\"", encoding = "UTF-8",na.strings=c("n/a","NA"))
  MetadataTable <- MetadataTable[,c("Sample ID","Sample Type","Latitude","Longitude","Sample Date","Fastq Forward Reads Filename")]
  #Merge in cleaned dates
  Metadata_Processed$Sample_Date <- as.Date(Metadata_Processed$Sample_Date)
  MetadataTable <- dplyr::left_join(MetadataTable, Metadata_Processed[,c("name","Sample_Date")],by=c("Sample ID"="name"),multiple="all")
  MetadataTable$`Sample Date` <- NULL
  names(MetadataTable)[names(MetadataTable) == 'Sample_Date'] <- 'Sample Date'
  #Add in eDNA Explorer project ID
  MetadataTable$projectid <- basename(Project_Directory)
  #Designate a fastqid.
  MetadataTable$fastqid <- gsub("(_L001_R1_001.fastq.gz|_R1_001.fastq.gz)$", "", MetadataTable$`Fastq Forward Reads Filename`)
  MetadataTable$fastqid <- gsub("_","-",MetadataTable$fastqid)
  MetadataTable$fastqid <- paste(Primer,MetadataTable$fastqid,sep="_")
  #Deduplicate metadata
  MetadataTable <- MetadataTable[!duplicated(MetadataTable),]
  #Add sample/control variable
  MetadataTable$Sample_or_Control <- ifelse(trimws(MetadataTable$`Sample Type`)!="Sample","Control Sample","True Sample")
  
  #Get source metadata
  Metadata_Source_Path <- paste(wd,"CALeDNA",ProjectID,"terradactyl","metabarcoding_metadata_original.csv",sep="/")
  Metadata_Source <- read.table(Metadata_Source_Path, header=TRUE, sep=",",as.is=T,skip=0,fill=TRUE,check.names=FALSE,quote = "\"", encoding = "UTF-8",na.strings=c("n/a","NA"))
  #Merge in substrate
  MetadataTable <- dplyr::left_join(MetadataTable,Metadata_Source[,c("Sample ID","Substrate")])
  MetadataTable <- MetadataTable[!duplicated(MetadataTable),]
  
  #Get tronko-assign data.
  TronkoFile <- paste(Primer,"_Max",Mismatch,".txt",sep="")
  TronkoTable_Path <- paste(wd,"CALeDNA",ProjectID,"tronko",Primer,TronkoFile,sep="/")
  #Only save tronko-assign data and metadata if the tronko-assign files exist.
  if(file.exists(TronkoTable_Path)){
    TronkoTable <- read.table(TronkoTable_Path, header=TRUE, sep="\t",as.is=T,skip=0,fill=TRUE,check.names=FALSE,quote = "\"", encoding = "UTF-8")
    #Create otu phyloseq object
    otu_mat <- TronkoTable[!duplicated(TronkoTable),]
    #colnames(otu_mat) <- gsub("-L001$", "", colnames(otu_mat))
    rownames(otu_mat) <- otu_mat[,paste(Primer,"assignment",sep="_")]
    #Create taxa table phyloseq object
    tax_mat <- as.data.frame(otu_mat[,paste(Primer,"assignment",sep="_")])
    colnames(tax_mat) <- "Taxon"
    rownames(tax_mat) <- tax_mat$Taxon
    tax_mat <- separate(data = tax_mat, col = "Taxon", sep=";",into = TaxonomicRanks)
    if(nrow(tax_mat)>1){
      otu_mat[,paste(Primer,"assignment",sep="_")] <- NULL
      otu_mat[is.na(otu_mat)] <- 0
      otu_mat <- as.matrix(otu_mat)
      OTU <- otu_table(otu_mat,taxa_are_rows=T)
      TAX <- tax_table(as.matrix(tax_mat))
      #Create sample phyloseq object
      sample_mat <- MetadataTable
      sample_mat <- MetadataTable[MetadataTable$`Sample ID` %in% colnames(otu_mat),]
      rownames(sample_mat) <- sample_mat$`Sample ID`
      Sample <- sample_data(sample_mat)
      #Create full phyloseq object
      physeq_initial <- phyloseq(OTU,TAX,Sample)
      #Add in check on sample/controls
      sample_data(physeq_initial)$is.neg <- sample_data(physeq_initial)$Sample_or_Control == "Control Sample"
      #Create a prevalence filter phyloseq object
      contamdf.prev_threshold <- isContaminant(physeq_initial, method="prevalence", neg="is.neg", threshold=0.1)
      #Filter on threshold
      contamdf.prev_threshold_filtered <-prune_taxa(!contamdf.prev_threshold$contaminant, physeq_initial)
      #Filter on taxa prevalence
      contamdf.prev_no_blanks = subset_samples(contamdf.prev_threshold_filtered, Sample_or_Control == "True Sample")
      contamdf.prev_min_threshold_a = prune_taxa(taxa_sums(contamdf.prev_no_blanks) > 2, contamdf.prev_no_blanks)
      contamdf.prev_min_threshold_p = prune_samples(sample_sums(contamdf.prev_min_threshold_a)>=2, contamdf.prev_min_threshold_a)
      #Save stats on the number samples and taxa per project before and after decontamination.
      project_stat <- data.frame(matrix(nrow=1,ncol=5))
      colnames(project_stat) <- c("ProjectID","Number_samples_inital","Number_taxa_initial","Number_samples_decontaminated","Number_taxa_decontaminated")
      project_stat$ProjectID <- ProjectID
      project_stat$Number_samples_inital <- nsamples(physeq_initial)
      project_stat$Number_taxa_initial <- ntaxa(physeq_initial)
      project_stat$Number_samples_decontaminated <- nsamples(contamdf.prev_min_threshold_p)
      project_stat$Number_taxa_decontaminated <- ntaxa(contamdf.prev_min_threshold_p)
      project_stats[[i]] <- project_stat
      #Store phyloseq objects
      if(!exists("physeq_total")){
        physeq_total <- contamdf.prev_min_threshold_p
      } else {
        physeq_total <- merge_phyloseq(physeq_total,contamdf.prev_min_threshold_p)
      }
      #Save filtered metadata
      Metadata_Filtered <- as.data.frame(sample_data(contamdf.prev_min_threshold_p))
      colnames(Metadata_Filtered) <- gsub("\\."," ",colnames(Metadata_Filtered))
      names(Metadata_Filtered)[names(Metadata_Filtered) == 'is neg'] <- 'is.neg'
    }
  }
  i=i+1
  print(paste(i,TronkoTable_Path))
}
#Store project stats
project_stats <- rbind.fill(project_stats)

#Export total taxonomy table, using Tronko-assign taxonomy
if(Primer=="CO1_Metazoa"){taxa_retain <- taxa_names(physeq_total)[tax_table(physeq_total)[, "phylum"] %in% Invertebrate_Phyla]}
if(Primer=="ITS1_Fungi"){taxa_retain <- taxa_names(physeq_total)[tax_table(physeq_total)[, "phylum"] %in% Fungal_Phyla]}
physeq_retain <- prune_taxa(taxa_retain, physeq_total)
#Remove samples without location and date data.
physeq_retain <- subset_samples(physeq_retain, !is.na(Latitude))
physeq_retain <- subset_samples(physeq_retain, !is.na(Longitude))
physeq_retain <- subset_samples(physeq_retain, !is.na(Sample.Date))

#Create presence/absence taxa by sample table.
Project_OTU_Table <- as.data.frame(otu_table(physeq_retain))
Project_OTU_Table <- Project_OTU_Table %>% mutate_all(~ ifelse(. > 0, 1, 0))
Project_OTU_Table$Taxonomic_Path <- rownames(Project_OTU_Table)
Project_OTU_Table <- separate(data = Project_OTU_Table, col = "Taxonomic_Path", sep=";",into = TaxonomicRanks)
Project_OTU_Table[Project_OTU_Table == "NA"] <- NA
#Select taxonomic group to calculate richness
selected_taxon <- "Pseudobaeospora"
selected_rank <- "genus"
#Obtain coordinates for selected taxon
Selected_OTU_Table <- Project_OTU_Table[!is.na(Project_OTU_Table[,selected_rank]) & Project_OTU_Table[,selected_rank]==selected_taxon,]
Selected_OTU_Table <- Selected_OTU_Table[,!colnames(Selected_OTU_Table) %in% TaxonomicRanks]
site_list <- colnames(Selected_OTU_Table[colSums(Selected_OTU_Table)>0])
Project_Metadata <- as.data.frame(sample_data(physeq_retain))
if(sum(Project_Metadata$Sample.ID %in% site_list)>0){
  Occurrence_Export <- as.data.frame(Project_Metadata[Project_Metadata$Sample.ID %in% site_list,c("Latitude","Longitude")])
  Occurrence_Export <- as.data.frame(Occurrence_Export@.Data)
  colnames(Occurrence_Export) <- c("Latitude","Longitude")
  #Add in taxonomic information,presence, and source information
  Occurrence_Export$Taxon <- selected_taxon
  Occurrence_Export$taxonRank <- selected_rank
  Occurrence_Export$scientificName <- selected_taxon
  Occurrence_Export$Source <- "eDNA"
  Occurrence_Export$Present <- 1
  Occurrence_Export <- Occurrence_Export[,c("Taxon","taxonRank","scientificName","Latitude","Longitude","Source","Present")]
} else{
  Occurrence_Export <- data.frame()
}

#Read in randomly generated background point locations.
Background <- read.table("BackgroundPoints.csv", header=TRUE, sep=",",as.is=T,skip=0,fill=TRUE,check.names=FALSE,quote = "\"", encoding = "UTF-8",na.strings=c("n/a","NA"))
#Add in taxonomic information,presence, and source information
Background$Taxon <- selected_taxon
Background$taxonRank <- selected_rank
Background$scientificName <- selected_taxon
Background$Source <- "ArcGIS Pro"
Background$Present <- 0
Background <- Background[,c("Taxon","taxonRank","scientificName","Latitude","Longitude","Source","Present")]

#Obtain occurrences found in GBIF in California
GBIF_Export <- as.data.frame(occ_data(scientificName = selected_taxon, publishingCountry = "US",stateProvince = "California")$data)
GBIF_Export <- GBIF_Export[!is.na(GBIF_Export$decimalLatitude),]
GBIF_Export <- GBIF_Export[GBIF_Export$taxonomicStatus=="ACCEPTED",]
GBIF_Export$Taxon <- selected_taxon
GBIF_Export$taxonRank <- selected_rank
GBIF_Export <- GBIF_Export %>%
  mutate(scientificName = coalesce(species,genus,family,order,class,phylum,kingdom))
GBIF_Export$Source <- "GBIF"
GBIF_Export$Present <- 1
GBIF_Export$Longitude <- GBIF_Export$decimalLongitude
GBIF_Export$Latitude <- GBIF_Export$decimalLatitude
GBIF_Export <- GBIF_Export[,c("Taxon","taxonRank","scientificName","Latitude","Longitude","Source","Present")]

set.seed(1)
#Randomly select an equal number of background points to presence points.
#Rationale: https://doi.org/10.1016/j.ecoinf.2022.101914
#Add with presence data.
#Get the number of presence points
n_occurrences <- nrow(Occurrence_Export)+nrow(GBIF_Export)
#Combine presence and randomly selected background points
Occurrences_Export <- rbind(Occurrence_Export,GBIF_Export,Background[sample(nrow(Background), n_occurrences), ])
#Export
write.table(Occurrences_Export,paste(Primer,selected_taxon,selected_rank,"Occurrence.csv",sep="_"), row.names=FALSE, sep=",",quote = FALSE)
