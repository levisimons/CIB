rm(list=ls())
require(dplyr)
require(plyr)
require(tidyr)
require(phyloseq)
require(decontam)
require(dismo)
require(randomForest)
require(caret)

wd <- "/Users/levisimons/Desktop/CIB"

setwd(wd)

#Get project directories
Project_Directories <- list.dirs(path=wd,recursive=F)
Project_Directories <- Project_Directories[Project_Directories != paste(wd,"sratoolkit.3.0.10-mac-x86_64",sep="/")]
Project_Directories <- Project_Directories[Project_Directories != paste(wd,"ca-state-boundary",sep="/")]

#Select a primer
Primers <- c("12S_MiFish_U","16S_Bacteria","18S_Euk","CO1_Metazoa","ITS1_Fungi","ITS2_Plants","vert12S")
Primer <- "ITS1_Fungi"

#Aggregate tronko-assign results to a particular taxonomic level.
TaxonomicRanks <- c("superkingdom","phylum","class","order","family","genus","species")
TaxonomicRank <- "species"
Taxon_selected <- "Tremella mesenterica"

#Select a tronko-assign mismatch cutoff.
Mismatch <- 25
#Get tronko-assign tables, and sample metadata, across projects
TronkoTables <- c()
MetadataTables <- c()
#Get the number of samples before and after decontamination.
Project_Samples <-c()
#Store all environmental metadata.
Environmental_Metadata <- c()
i=1
for(Project_Directory in Project_Directories){
  ProjectID <- basename(Project_Directory)
  #Get processed project metadata
  Metadata_Processed_Path <- paste(wd,ProjectID,"terradactyl","metabarcoding_metadata_terradactyl.csv",sep="/")
  Metadata_Processed <- read.table(Metadata_Processed_Path, header=TRUE, sep=",",as.is=T,skip=0,fill=TRUE,check.names=FALSE,quote = "\"", encoding = "UTF-8",na.strings=c("n/a","NA"))
  #Store all extracted metadata.
  Metadata_Processed$projectid <- ProjectID
  Environmental_Metadata[[i]] <- Metadata_Processed
  #Get uploaded project metadata.
  Metadata_Path <- paste(wd,ProjectID,"terradactyl","metabarcoding_metadata_original.csv",sep="/")
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
  Metadata_Source_Path <- paste(wd,ProjectID,"terradactyl","metabarcoding_metadata_original.csv",sep="/")
  Metadata_Source <- read.table(Metadata_Source_Path, header=TRUE, sep=",",as.is=T,skip=0,fill=TRUE,check.names=FALSE,quote = "\"", encoding = "UTF-8",na.strings=c("n/a","NA"))
  #Merge in substrate
  MetadataTable <- dplyr::left_join(MetadataTable,Metadata_Source[,c("Sample ID","Substrate")])
  MetadataTable <- MetadataTable[!duplicated(MetadataTable),]
  
  #Get tronko-assign data.
  TronkoFile <- paste(Primer,"_Max",Mismatch,".txt",sep="")
  TronkoTable_Path <- paste(wd,ProjectID,"tronko",Primer,TronkoFile,sep="/")
  #Only save tronko-assign data and metadata if the tronko-assign files exist.
  if(file.exists(TronkoTable_Path)){
    TronkoTable <- read.table(TronkoTable_Path, header=TRUE, sep="\t",as.is=T,skip=0,fill=TRUE,check.names=FALSE,quote = "\"", encoding = "UTF-8")
    #Create otu phyloseq object
    otu_mat <- TronkoTable[!duplicated(TronkoTable),]
    colnames(otu_mat) <- gsub("-L001$", "", colnames(otu_mat))
    rownames(otu_mat) <- otu_mat[,paste(Primer,"seq_number",sep="_")]
    #Create taxa table phyloseq object
    tax_mat <- as.data.frame(otu_mat[,paste(Primer,"seq_number",sep="_")])
    colnames(tax_mat) <- "Taxon"
    rownames(tax_mat) <- tax_mat$Taxon
    tax_mat <- separate(data = tax_mat, col = "Taxon", sep=";",into = TaxonomicRanks)
    if(nrow(tax_mat)>1){
      otu_mat[,paste(Primer,"seq_number",sep="_")] <- NULL
      otu_mat[is.na(otu_mat)] <- 0
      otu_mat <- as.matrix(otu_mat)
      OTU <- otu_table(otu_mat,taxa_are_rows=T)
      TAX <- tax_table(as.matrix(tax_mat))
      #Create sample phyloseq object
      sample_mat <- MetadataTable
      sample_mat <- MetadataTable[MetadataTable$fastqid %in% colnames(otu_mat),]
      rownames(sample_mat) <- sample_mat$fastqid
      Sample <- sample_data(sample_mat)
      #Create full phyloseq object
      physeq_initial <- phyloseq(OTU,TAX,Sample)
      #Add in check on sample/controls
      sample_data(physeq_initial)$is.neg <- sample_data(physeq_initial)$Sample_or_Control == "Control Sample"
      #Create a prevalence filter phyloseq object
      contamdf.prev_threshold <- isContaminant(physeq_initial, method="prevalence", neg="is.neg", threshold=0.5)
      #Filter on threshold
      contamdf.prev_threshold_filtered <-prune_taxa(!contamdf.prev_threshold$contaminant, physeq_initial)
      #Filter on taxa prevalence
      contamdf.prev_prevalence_filtered = filter_taxa(contamdf.prev_threshold_filtered,function(x) sum(x > 2) > 1, TRUE)
      contamdf.prev_no_blanks = subset_samples(contamdf.prev_prevalence_filtered, Sample_or_Control == "True Sample")
      contamdf.prev_min_threshold_a = prune_taxa(taxa_sums(contamdf.prev_no_blanks) > 19, contamdf.prev_no_blanks)
      contamdf.prev_min_threshold_p = prune_samples(sample_sums(contamdf.prev_min_threshold_a)>=20, contamdf.prev_min_threshold_a)
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
      MetadataTables[[i]] <- Metadata_Filtered
      #Get the number of samples retained per project.
      Project_Sample <- data.frame(matrix(nrow=1,ncol=3))
      Project_Sample$projectid <- ProjectID
      Project_Sample$Num_Samples <- nrow(MetadataTable)
      Project_Sample$Num_Samples_Filtered <- nrow(Metadata_Filtered)
      Project_Samples[[i]] <- Project_Sample
      #Store tronko results by project.
      Tronko_Filtered <- as.data.frame(otu_table(contamdf.prev_min_threshold_p))
      Tronko_Filtered$Taxon <- rownames(Tronko_Filtered)
      rownames(Tronko_Filtered) <- NULL
      TronkoTables[[i]] <- Tronko_Filtered
    }
  }
  i=i+1
  print(paste(i,TronkoTable_Path))
}

#Rarefy phyloseq object.
physeq_total <- rarefy_even_depth(physeq_total, rngseed=1, sample.size=1000, replace=F)

#Aggregate tronko-assign results into a single table.
TronkoTables <- rbind.fill(TronkoTables)
TronkoTables[is.na(TronkoTables)] <- 0

#Aggregate summary statistics
Project_Samples <- rbind.fill(Project_Samples)

#Aggregate project metadata results into a single table.
MetadataTables <- rbind.fill(MetadataTables)
#Standardize substrate types for metadata
MetadataTables$Substrate <- tolower(MetadataTables$Substrate)
MetadataTables$Substrate <- ifelse(MetadataTables$Substrate %in% c("sand", "organic_sand"), "soil", MetadataTables$Substrate)
MetadataTables$Substrate <- ifelse(MetadataTables$Substrate %in% c("silt"), "sediment", MetadataTables$Substrate)
MetadataTables$Substrate <- ifelse(MetadataTables$Substrate %in% c("", "NA"), NA, MetadataTables$Substrate)

#Aggreate environmental metadata
Environmental_Metadata <- rbind.fill(Environmental_Metadata)
#Set date and NA values
Environmental_Metadata$Sample_Date <- as.Date(Environmental_Metadata$Sample_Date)
Environmental_Metadata[Environmental_Metadata==-999999] <- NA
Environmental_Metadata[Environmental_Metadata==""] <- NA
Environmental_Metadata[Environmental_Metadata==-32768] <- NA
#Remove empty columns
Environmental_Metadata <- Filter(function(x)!all(is.na(x)), Environmental_Metadata)

#Aggregate tronko-assign results by a selected taxonomic rank
TronkoTables <- separate(data = TronkoTables, col = "Taxon", sep=";",into = TaxonomicRanks)
SampleIDs <- colnames(TronkoTables)[!(colnames(TronkoTables) %in% TaxonomicRanks)]
#Filter to only contain data within a taxon
TronkoTables <- TronkoTables[TronkoTables[,TaxonomicRank]==Taxon_selected & !is.na(TronkoTables[,TaxonomicRank]),]
#Sum read counts by taxon across samples
TronkoTables <- TronkoTables %>% group_by(superkingdom, phylum, class, order, family, genus, species) %>% summarise_all(sum)
#Convert read counts to presence/absence
TronkoTables <- TronkoTables %>% mutate_if(is.numeric, ~as.numeric(. > 0))
#Get sample IDs for where taxon is found.
SampleIDs_Present <- names(TronkoTables)[colSums(TronkoTables[,!(colnames(TronkoTables) %in% TaxonomicRanks)]) != 0]

#Create occurrence table for tronko-assign results
TronkoExport <- MetadataTables[MetadataTables$fastqid %in% SampleIDs_Present,c("Latitude","Longitude")]
TronkoExport$Taxon <- Taxon_selected
TronkoExport$taxonRank <- TaxonomicRank
TronkoExport$scientificName <- Taxon_selected
TronkoExport$Source <- "eDNA"
#Add in presence column
TronkoExport$Present <- 1
TronkoExport <- TronkoExport[,c("Taxon","taxonRank","scientificName","Latitude","Longitude","Source","Present")]

#GBIF.org (08 April 2024) GBIF Occurrence Download  https://doi.org/10.15468/dl.x5kfsx
#All observations with less than 1km uncertainty in California taken using
#Human observation, Machine observation, Observation, or Preserved specimen
GBIF_ranks <- c("kingdom","phylum","class","order","family","genus","species")
CA_GBIF_Input <- fread(input="CA_GBIF.csv",select=c(GBIF_ranks,"taxonRank","scientificName","decimalLatitude","decimalLongitude"),quote="")
#Filter out to only include selected taxa.
CA_GBIF_Input <- CA_GBIF_Input[CA_GBIF_Input[,get(TaxonomicRank)]==Taxon_selected,]
#Clean and standardize taxonomic data.
CA_GBIF_Input$taxonRank <- tolower(CA_GBIF_Input$taxonRank)
CA_GBIF_Input[CA_GBIF_Input==""] <- NA
CA_GBIF_Input <- CA_GBIF_Input[CA_GBIF_Input$taxonRank %in% GBIF_ranks,]
#Filter out to only include observations within coordinates.
CA_GBIF_Input <- CA_GBIF_Input[!is.na(CA_GBIF_Input$decimalLatitude),]
#Get unique taxa occurrences within a selected taxon
CA_GBIF_Occurrences <- c()
for(i in 1:nrow(CA_GBIF_Input)){
  tmp <- CA_GBIF_Input[i,]
  tmp <- tmp %>% rowwise() %>% mutate(Taxon = get(taxonRank))
  tmp <- tmp[,c("Taxon","taxonRank","scientificName","decimalLatitude","decimalLongitude")]
  CA_GBIF_Occurrences[[i]] <- tmp
  print(i)
}
CA_GBIF_Occurrences <- rbind.fill(CA_GBIF_Occurrences)
names(CA_GBIF_Occurrences)[names(CA_GBIF_Occurrences) == "decimalLatitude"] <- "Latitude"
names(CA_GBIF_Occurrences)[names(CA_GBIF_Occurrences) == "decimalLongitude"] <- "Longitude"
CA_GBIF_Occurrences <- CA_GBIF_Occurrences[!duplicated(CA_GBIF_Occurrences),]
#Add in source column
CA_GBIF_Occurrences$Source <- "GBIF"
#Add in presence column
CA_GBIF_Occurrences$Present <- 1

#Read in randomly generated background point locations.
Background <- read.table("BackgroundPoints.csv", header=TRUE, sep=",",as.is=T,skip=0,fill=TRUE,check.names=FALSE,quote = "\"", encoding = "UTF-8",na.strings=c("n/a","NA"))
Background$Taxon <- Taxon_selected
Background$taxonRank <- TaxonomicRank
Background$scientificName <- Taxon_selected
Background$Source <- "ArcGIS Pro"
#Add in presence column
Background$Present <- 0
Background <- Background[,c("Taxon","taxonRank","scientificName","Latitude","Longitude","Source","Present")]

#Combine presence and background points
Occurrences_Export <- rbind(TronkoExport,CA_GBIF_Occurrences,Background)

#Export presence/background table.
write.table(Occurrences_Export,paste(gsub(" ","_",Taxon_selected),"_PresenceBackground.csv",sep=""),quote=FALSE,sep=",",row.names = FALSE)
