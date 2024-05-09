rm(list=ls())
require(dplyr)
require(plyr)
require(tidyr)
require(phyloseq)
require(decontam)
require(data.table)

wd <- "/Users/levisimons/Desktop/CIB"

setwd(wd)

#Select a primer
Primers <- c("12S_MiFish_U","16S_Bacteria","18S_Euk","CO1_Metazoa","ITS1_Fungi","ITS2_Plants","vert12S")
Primer <- "ITS1_Fungi"

#Read in GBIF data sets.
#GBIF.org (03 May 2024) GBIF Occurrence Download  https://doi.org/10.15468/dl.j6wbeq
#All taxa in California
CA_GBIF_All <- read.table("CA_GBIF_Taxa_Full.csv", header=TRUE, sep="\t",as.is=T,skip=0,fill=TRUE,check.names=FALSE,quote = "\"", encoding = "UTF-8")

#GBIF.org (03 May 2024) GBIF Occurrence Download  https://doi.org/10.15468/dl.2jwkqg
#All taxa in California.  Observation categories: Preserved Specimen, Fossil Specimen, Living Specimen.
CA_GBIF_Collections <- read.table("CA_GBIF_Taxa_Collections.csv", header=TRUE, sep="\t",as.is=T,skip=0,fill=TRUE,check.names=FALSE,quote = "\"", encoding = "UTF-8")

#GBIF.org (03 May 2024) GBIF Occurrence Download  https://doi.org/10.15468/dl.djd56f
#All taxa in California. Spatial uncertainty under 1 km.
CA_GBIF_1km <- read.table("CA_GBIF_Taxa_1km.csv", header=TRUE, sep="\t",as.is=T,skip=0,fill=TRUE,check.names=FALSE,quote = "\"", encoding = "UTF-8")

#Read in iNaturalist taxon ID to GBIF ID backbone.
iNaturalist_backbone <- fread(input="iNaturalist_taxa.csv",sep=",",select=c("id","kingdom","phylum","class","order","family","genus","scientificName","taxonRank"))

#Filter GBIF occurrences, and read in the appropriate genera-level iNaturalist observations.
if(Primer=="ITS1_Fungi"){
  CA_GBIF_All <- CA_GBIF_All[CA_GBIF_All$kingdom=="Fungi",]
  CA_GBIF_Collections <- CA_GBIF_Collections[CA_GBIF_Collections$kingdom=="Fungi",]
  CA_GBIF_1km <- CA_GBIF_1km[CA_GBIF_1km$kingdom=="Fungi",]
  iNaturalist_Locations <- fread("iNaturalist_Fungi.csv",sep=",",select=c("taxon_id","scientific_name","latitude","longitude","positional_accuracy"))
  iNaturalist_Locations <- dplyr::left_join(iNaturalist_Locations,iNaturalist_backbone,by=c("taxon_id"="id"))
  #Get unique iNaturalist taxa
  iNaturalist_taxa <- iNaturalist_Locations[,c("taxon_id","kingdom","phylum","class","order","family","genus","scientificName","taxonRank")]
  iNaturalist_taxa <- iNaturalist_taxa[!duplicated(iNaturalist_taxa),]
  #Get unique iNaturalist occurrences with a 1 km spatial uncertainty filter.
  iNaturalist_occurrences <- iNaturalist_Locations[iNaturalist_Locations$positional_accuracy < 1000 & !is.na(iNaturalist_Locations$positional_accuracy),]
  #Assign superkingdoms
  iNaturalist_occurrences <- iNaturalist_occurrences %>%
    mutate(superkingdom = case_when(
      kingdom %in% c("Animalia", "Plantae", "Fungi", "Protista") ~ "Eukaryota",
      kingdom == "Bacteria" ~ "Bacteria",
      kingdom == "Archaea" ~ "Archaea",
      TRUE ~ NA_character_
    ))
  iNaturalist_occurrences <- iNaturalist_occurrences %>% dplyr::count(superkingdom,phylum,class,order,family,genus)
  names(iNaturalist_occurrences)[names(iNaturalist_occurrences) == "n"] <- "iNaturalist_prevalence"
}
if(Primer=="CO1_Metazoa"){
  CA_GBIF_All <- CA_GBIF_All[CA_GBIF_All$phylum %in% c("Platyhelminthes","Turbellaria","Trematoda","Cestoda","Nemertea","Nemertea","Rotifera","Gastrotricha","Acanthocephala","Nematoda","Nematomorpha","Priapulida","Kinorhyncha","Loricifera","Entoprocta","Cycliophora","Gnathostomulida","Micrognathozoa","Chaetognatha","Hemichordata","Bryozoa","Brachiopoda","Phoronida","Ectoprocta","Annelida","Polychaeta","Clitellata","Mollusca","Gastropoda","Bivalvia","Arthropoda","Insecta","Arachnida","Crustacea","Myriapoda"),]
  CA_GBIF_Collections <- CA_GBIF_Collections[CA_GBIF_Collections$phylum %in% c("Platyhelminthes","Turbellaria","Trematoda","Cestoda","Nemertea","Nemertea","Rotifera","Gastrotricha","Acanthocephala","Nematoda","Nematomorpha","Priapulida","Kinorhyncha","Loricifera","Entoprocta","Cycliophora","Gnathostomulida","Micrognathozoa","Chaetognatha","Hemichordata","Bryozoa","Brachiopoda","Phoronida","Ectoprocta","Annelida","Polychaeta","Clitellata","Mollusca","Gastropoda","Bivalvia","Arthropoda","Insecta","Arachnida","Crustacea","Myriapoda"),]
}

#GBIF.org (23 April 2024) GBIF Occurrence Download https://doi.org/10.15468/dl.3skmb8
#All observations with less than 1km uncertainty in California taken using
#All observation categories.
GBIF_ranks <- c("kingdom","phylum","class","order","family","genus","species")
CA_GBIF_Input <- fread(input="CA_GBIF_1km.csv",select=c(GBIF_ranks,"taxonRank","scientificName","decimalLatitude","decimalLongitude"),quote="")
#Clean and standardize taxonomic data.
CA_GBIF_Input$taxonRank <- tolower(CA_GBIF_Input$taxonRank)
CA_GBIF_Input[CA_GBIF_Input==""] <- NA
#CA_GBIF_Input <- CA_GBIF_Input[CA_GBIF_Input$taxonRank %in% GBIF_ranks,]
#Filter out to only include observations with coordinates.
CA_GBIF_Input <- CA_GBIF_Input[!is.na(CA_GBIF_Input$decimalLatitude),]
#Assign superkingdoms
CA_GBIF_Input <- CA_GBIF_Input %>%
  mutate(superkingdom = case_when(
    kingdom %in% c("Animalia", "Plantae", "Fungi", "Protista") ~ "Eukaryota",
    kingdom == "Bacteria" ~ "Bacteria",
    kingdom == "Archaea" ~ "Archaea",
    TRUE ~ NA_character_
  ))
#Filter GBIF occurrences
CA_GBIF_Input <- as.data.frame(CA_GBIF_Input)
if(Primer=="ITS1_Fungi"){
  CA_GBIF_Input <- CA_GBIF_Input[CA_GBIF_Input$kingdom=="Fungi",]
}
if(Primer=="CO1_Metazoa"){
  CA_GBIF_Input <- CA_GBIF_Input[CA_GBIF_Input$phylum %in% c("Platyhelminthes","Turbellaria","Trematoda","Cestoda","Nemertea","Nemertea","Rotifera","Gastrotricha","Acanthocephala","Nematoda","Nematomorpha","Priapulida","Kinorhyncha","Loricifera","Entoprocta","Cycliophora","Gnathostomulida","Micrognathozoa","Chaetognatha","Hemichordata","Bryozoa","Brachiopoda","Phoronida","Ectoprocta","Annelida","Polychaeta","Clitellata","Mollusca","Gastropoda","Bivalvia","Arthropoda","Insecta","Arachnida","Crustacea","Myriapoda"),]
}
#Set columns to match eDNA prevalence table
CA_GBIF_Input <- CA_GBIF_Input[,TaxonomicRanks]
#Aggregate GBIF occurrences by unique taxon.
CA_GBIF_Input <- CA_GBIF_Input %>% dplyr::count(superkingdom,phylum,class,order,family,genus,species)
names(CA_GBIF_Input)[names(CA_GBIF_Input) == "n"] <- "GBIF_prevalence"

##Read in and process Tronko-assign eDNA data.
#Get project directories
Project_Directories <- list.dirs(path=wd,recursive=F)
Project_Directories <- Project_Directories[Project_Directories != paste(wd,"sratoolkit.3.0.10-mac-x86_64",sep="/")]
Project_Directories <- Project_Directories[Project_Directories != paste(wd,"ca-state-boundary",sep="/")]

#Aggregate tronko-assign results to a particular taxonomic level.
TaxonomicRanks <- c("superkingdom","phylum","class","order","family","genus","species")

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

#Aggregate environmental metadata
Environmental_Metadata <- rbind.fill(Environmental_Metadata)
#Set date and NA values
Environmental_Metadata$Sample_Date <- as.Date(Environmental_Metadata$Sample_Date)
Environmental_Metadata[Environmental_Metadata==-999999] <- NA
Environmental_Metadata[Environmental_Metadata==""] <- NA
Environmental_Metadata[Environmental_Metadata==-32768] <- NA
#Remove empty columns
Environmental_Metadata <- Filter(function(x)!all(is.na(x)), Environmental_Metadata)

#Aggregate tronko-assign reads by unique taxon.
TronkoTables <- aggregate(. ~ Taxon, data = TronkoTables, FUN = sum)
#Convert read counts to presence/absence.
TronkoTables[,!(colnames(TronkoTables) %in% "Taxon")] <- TronkoTables[,!(colnames(TronkoTables) %in% "Taxon")] %>% mutate_all(~ ifelse(. > 0, 1, 0))
#Get sample prevalence
TronkoTables$eDNA_prevalence <- rowSums(TronkoTables[,!(colnames(TronkoTables) %in% "Taxon")])

#Generate sample prevalence table for tronko-assign results
TronkoExport <- TronkoTables[,c("Taxon","eDNA_prevalence")]
TronkoExport <- separate(data = TronkoExport, col = "Taxon", sep=";",into = TaxonomicRanks)
#Clean NA values.
TronkoExport <- TronkoExport %>% mutate_all(~ifelse(. == "NA", NA, .))
#Filter taxa to match GBIF filtering.
if(Primer=="CO1_Metazoa"){
  TronkoExport <- TronkoExport[TronkoExport$phylum %in% c("Platyhelminthes","Turbellaria","Trematoda","Cestoda","Nemertea","Nemertea","Rotifera","Gastrotricha","Acanthocephala","Nematoda","Nematomorpha","Priapulida","Kinorhyncha","Loricifera","Entoprocta","Cycliophora","Gnathostomulida","Micrognathozoa","Chaetognatha","Hemichordata","Bryozoa","Brachiopoda","Phoronida","Ectoprocta","Annelida","Polychaeta","Clitellata","Mollusca","Gastropoda","Bivalvia","Arthropoda","Insecta","Arachnida","Crustacea","Myriapoda"),]
}
#Aggregate prevalence by duplicate taxa.
TronkoExport <- TronkoExport %>% mutate_all(~replace_na(., "NA"))
TronkoExport <- aggregate(. ~superkingdom+phylum+class+order+family+genus+species, data=TronkoExport, sum, na.rm=TRUE)
TronkoExport <- TronkoExport %>% mutate_all(~ifelse(. == "NA", NA, .))

#Add in column showing if eDNA taxa show up in California GBIF data (No spatial filter).
tmp <- CA_GBIF_All[,c("kingdom","phylum","class","order","family","genus","species")]
#Assign superkingdoms
tmp <- tmp %>%
  mutate(superkingdom = case_when(
    kingdom %in% c("Animalia", "Plantae", "Fungi", "Protista") ~ "Eukaryota",
    kingdom == "Bacteria" ~ "Bacteria",
    kingdom == "Archaea" ~ "Archaea",
    TRUE ~ NA_character_
  ))
tmp$kingdom <- NULL
tmp$GBIF_in_CA <- 1
tmp <- tmp[!duplicated(tmp),]
Prevalence_Export <- dplyr::left_join(TronkoExport,tmp)
Prevalence_Export$GBIF_in_CA <- ifelse(is.na(Prevalence_Export$GBIF_in_CA),0,1)

#Add in column showing if eDNA taxa show up in California GBIF collections.
tmp <- CA_GBIF_Collections[,c("kingdom","phylum","class","order","family","genus","species")]
#Assign superkingdoms
tmp <- tmp %>%
  mutate(superkingdom = case_when(
    kingdom %in% c("Animalia", "Plantae", "Fungi", "Protista") ~ "Eukaryota",
    kingdom == "Bacteria" ~ "Bacteria",
    kingdom == "Archaea" ~ "Archaea",
    TRUE ~ NA_character_
  ))
tmp$kingdom <- NULL
tmp$GBIF_in_CA_Collections <- 1
tmp <- tmp[!duplicated(tmp),]
Prevalence_Export <- dplyr::left_join(Prevalence_Export,tmp)
Prevalence_Export$GBIF_in_CA_Collections <- ifelse(is.na(Prevalence_Export$GBIF_in_CA_Collections),0,1)

#Add in column showing if eDNA taxa show up in California iNaturalist data.
tmp <- iNaturalist_taxa[,c("kingdom","phylum","class","order","family","genus"),]
#Assign superkingdoms
tmp <- tmp %>%
  mutate(superkingdom = case_when(
    kingdom %in% c("Animalia", "Plantae", "Fungi", "Protista") ~ "Eukaryota",
    kingdom == "Bacteria" ~ "Bacteria",
    kingdom == "Archaea" ~ "Archaea",
    TRUE ~ NA_character_
  ))
tmp$kingdom <- NULL
tmp$iNaturalist_in_CA <- 1
tmp <- tmp[!duplicated(tmp),]
Prevalence_Export <- dplyr::left_join(Prevalence_Export,tmp)
Prevalence_Export$iNaturalist_in_CA <- ifelse(is.na(Prevalence_Export$iNaturalist_in_CA),0,1)

#Combine eDNA, iNaturalist, GBIF prevalence tables.
#Note that GBIF and iNaturalist prevalence tables are filtered to under 1 km spatial uncertainty.
Prevalence_Export <- dplyr::full_join(Prevalence_Export,CA_GBIF_Input)
Prevalence_Export <- Prevalence_Export[!duplicated(Prevalence_Export),]
Prevalence_Export <- dplyr::full_join(Prevalence_Export,iNaturalist_occurrences,by=c("superkingdom","phylum","class","order","family","genus"))
Prevalence_Export <- Prevalence_Export[!duplicated(Prevalence_Export),]
Prevalence_Export$Total_prevalence <- rowSums(Prevalence_Export[,c("GBIF_prevalence","eDNA_prevalence","iNaturalist_prevalence")], na.rm = TRUE)
Prevalence_Export <- Prevalence_Export %>% mutate_if(is.numeric, ~replace(., is.na(.), 0))

#Export prevalence tables
write.table(Prevalence_Export,paste(Primer,"_CA_Taxa_Prevalence.csv",sep=""),quote=FALSE,sep=",",row.names = FALSE)
write.table(Prevalence_Export[Prevalence_Export$Total_prevalence>=30,],paste(Primer,"_CA_Taxa_Prevalence_Filtered.csv",sep=""),quote=FALSE,sep=",",row.names = FALSE)
