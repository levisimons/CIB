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
      #contamdf.prev_prevalence_filtered = filter_taxa(contamdf.prev_threshold_filtered,function(x) sum(x > 2) > 1, TRUE)
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
#
project_stats <- rbind.fill(project_stats)

#Export total taxonomy table, using NCBI taxonomy
if(Primer=="CO1_Metazoa"){taxa_retain <- taxa_names(physeq_total)[tax_table(physeq_total)[, "phylum"] %in% Invertebrate_Phyla]}
if(Primer=="ITS1_Fungi"){taxa_retain <- taxa_names(physeq_total)[tax_table(physeq_total)[, "phylum"] %in% Fungal_Phyla]}
physeq_retain <- prune_taxa(taxa_retain, physeq_total)
OTUTable_Tronko <- as.data.frame(otu_table(physeq_retain))
#Convert read counts to presence/absence.
OTUTable_Tronko <- OTUTable_Tronko %>% mutate_all(~ ifelse(. > 0, 1, 0))
#Get taxa by ranks
TaxonomyExport_Tronko <- OTUTable_Tronko
TaxonomyExport_Tronko$Taxon <- rownames(TaxonomyExport_Tronko)
TaxonomyExport_Tronko <- separate(data = TaxonomyExport_Tronko, col = "Taxon", sep=";",into = paste0("tronko_",TaxonomicRanks,sep=""))
TaxonomyExport_Tronko$eDNA_prevalence_Tronko <- rowSums(TaxonomyExport_Tronko[,!(colnames(TaxonomyExport_Tronko) %in% paste0("tronko_",TaxonomicRanks,sep=""))])
TaxonomyExport_Tronko <- TaxonomyExport_Tronko[,c(paste0("tronko_",TaxonomicRanks,sep=""),"eDNA_prevalence_Tronko")]
TaxonomyExport_Tronko <- TaxonomyExport_Tronko %>% mutate_all(~ifelse(. == "NA", NA, .))
#Assign taxon and rank explicitly.
TaxonomyExport_Tronko <- TaxonomyExport_Tronko %>% mutate(tronko_name = coalesce(tronko_species, tronko_genus, tronko_family, tronko_order, tronko_class, tronko_phylum, tronko_superkingdom))
TaxonomyExport_Tronko <- TaxonomyExport_Tronko %>% mutate(tronko_rank = case_when(!is.na(tronko_species) ~ "species",!is.na(tronko_genus) ~ "genus",!is.na(tronko_family) ~ "family",!is.na(tronko_order) ~ "order",!is.na(tronko_class) ~ "class",!is.na(tronko_phylum) ~ "phylum",!is.na(tronko_superkingdom) ~ "superkingdom",TRUE ~ NA_character_))

#Read in GBIF-NCBI backbone
GBIF_NCBI <- fread(input="~/Desktop/Archive/backbone/GBIF_NCBI_export.tsv",sep="\t")
#Read in NCBI backbone and determine the most resolved taxon and rank per entry.
ncbi_taxa_full <- fread(input="~/Desktop/Archive/backbone/ncbi_taxa_full.tsv",sep="\t")
#Get the names of the most resolved NCBI taxa
ncbi_taxa_full <- ncbi_taxa_full %>% mutate(ncbi_name = coalesce(ncbi_species, ncbi_genus, ncbi_family, ncbi_order, ncbi_class, ncbi_phylum, ncbi_kingdom, ncbi_superkingdom))
NCBI_ranks <- c("ncbi_species", "ncbi_genus", "ncbi_family", "ncbi_order", "ncbi_class", "ncbi_phylum", "ncbi_kingdom", "ncbi_superkingdom")
#Get ranks of the most resolved NCBI taxa
tmp <- ncbi_taxa_full[,..NCBI_ranks]
first_non_na_index <- function(row) {
  which(!is.na(row))[1]
}
ncbi_taxa_full$ncbi_rank <- NCBI_ranks[apply(tmp[,..NCBI_ranks], 1, first_non_na_index)]
ncbi_taxa_full$ncbi_rank <- gsub("ncbi_","",ncbi_taxa_full$ncbi_rank)
ncbi_taxa_full <- ncbi_taxa_full[!duplicated(ncbi_taxa_full),]

#Check for where Tronko outputs don't fully align with NCBI taxonomies.
TaxonomyExport_NCBI <- TaxonomyExport_Tronko
tmp <- ncbi_taxa_full[,..NCBI_ranks]
tmp$in_ncbi <- 1
tmp <- tmp[!duplicated(tmp),]
TaxonomyExport_NCBI <- dplyr::left_join(TaxonomyExport_NCBI,tmp[,.(ncbi_species,ncbi_genus,ncbi_family,ncbi_order,ncbi_class,ncbi_phylum,ncbi_superkingdom,in_ncbi)],
                                        by=c("tronko_species"="ncbi_species","tronko_genus"="ncbi_genus","tronko_family"="ncbi_family","tronko_order"="ncbi_order","tronko_class"="ncbi_class","tronko_phylum"="ncbi_phylum","tronko_superkingdom"="ncbi_superkingdom"))
TaxonomyExport_NCBI_Mismatched <- TaxonomyExport_NCBI[is.na(TaxonomyExport_NCBI$in_ncbi),]
#Set aside entries where Tronko outputs fully align with NCBI taxonomies.
TaxonomyExport_NCBI <- TaxonomyExport_NCBI[!is.na(TaxonomyExport_NCBI$in_ncbi),]
TaxonomyExport_NCBI <- dplyr::left_join(TaxonomyExport_NCBI,ncbi_taxa_full[,.(ncbi_name,ncbi_rank,ncbi_id,ncbi_species,ncbi_genus,ncbi_family,ncbi_order,ncbi_class,ncbi_phylum,ncbi_superkingdom)],by=c("tronko_species"="ncbi_species","tronko_genus"="ncbi_genus","tronko_family"="ncbi_family","tronko_order"="ncbi_order","tronko_class"="ncbi_class","tronko_phylum"="ncbi_phylum","tronko_superkingdom"="ncbi_superkingdom"))
for(TaxonomicRank in TaxonomicRanks){
  TaxonomyExport_NCBI[,paste("ncbi_",TaxonomicRank,sep="")] <- TaxonomyExport_NCBI[,paste("tronko_",TaxonomicRank,sep="")]
}
#Find the most resolved NCBI taxonomies to match Tronko outputs
Tronko_ranks <- paste0("tronko_",TaxonomicRanks,sep="")
TaxonomyExport_NCBI_corrected <- c()
for(i in 1:nrow(TaxonomyExport_NCBI_Mismatched)){
  row_check <- TaxonomyExport_NCBI_Mismatched[i,]
  ncbi_check <- ncbi_taxa_full[ncbi_taxa_full$ncbi_name==row_check$tronko_name & ncbi_taxa_full$ncbi_rank==row_check$tronko_rank,]
  if(nrow(ncbi_check)==1){
    TaxonomyExport_NCBI_corrected[[i]] <- dplyr::left_join(row_check,ncbi_check[,.(ncbi_name,ncbi_rank,ncbi_id,ncbi_species,ncbi_genus,ncbi_family,ncbi_order,ncbi_class,ncbi_phylum,ncbi_superkingdom)],by=c("tronko_name"="ncbi_name","tronko_rank"="ncbi_rank"))
  }
  if(nrow(ncbi_check)==2){
    #If there are multiple matches with a given taxonomic name and rank for the most resolved Tronko taxon,
    # choose the one with the largest number of matches across all of the taxonomic ranks.
    ncbi_check$match_count <- as.integer(ncbi_check$ncbi_species %in% row_check$tronko_species)+
      as.integer(ncbi_check$ncbi_genus %in% row_check$tronko_genus)+
      as.integer(ncbi_check$ncbi_family %in% row_check$tronko_family)+
      as.integer(ncbi_check$ncbi_order %in% row_check$tronko_order)+
      as.integer(ncbi_check$ncbi_class %in% row_check$tronko_class)+
      as.integer(ncbi_check$ncbi_phylum %in% row_check$tronko_phylum)+
      as.integer(ncbi_check$ncbi_superkingdom %in% row_check$tronko_superkingdom)
    ncbi_check <- ncbi_check[ncbi_check$match_count==max(ncbi_check$match_count),]
    ncbi_check$match_count <- NULL
    TaxonomyExport_NCBI_corrected[[i]] <- dplyr::left_join(row_check,ncbi_check[,.(ncbi_name,ncbi_rank,ncbi_id,ncbi_species,ncbi_genus,ncbi_family,ncbi_order,ncbi_class,ncbi_phylum,ncbi_superkingdom)],by=c("tronko_name"="ncbi_name","tronko_rank"="ncbi_rank"))
  }
  #If no match occurrs between a Tronko-assigned taxon and an NCBI entry at the most
  #resolved taxonomic level, the keep moving up ranks to find the most resolved NCBI match.
  if(nrow(ncbi_check)==0){
    rank_num <- which(row_check[,Tronko_ranks]==row_check$tronko_name)
    while(nrow(ncbi_check)==0){
      
      new_tronko_name <- row_check[,Tronko_ranks[rank_num-1]]
      new_tronko_rank <- gsub("tronko_","",Tronko_ranks[rank_num-1])
      ncbi_check <- ncbi_taxa_full[ncbi_taxa_full$ncbi_name==new_tronko_name & ncbi_taxa_full$ncbi_rank==new_tronko_rank,]
      rank_num <- rank_num-1
      print(paste(new_tronko_name,new_tronko_rank))
    }
    row_check$tronko_name_tmp <- new_tronko_name
    row_check$tronko_rank_tmp <- new_tronko_rank
    tmp <- dplyr::left_join(row_check,ncbi_check[,.(ncbi_name,ncbi_rank,ncbi_id,ncbi_species,ncbi_genus,ncbi_family,ncbi_order,ncbi_class,ncbi_phylum,ncbi_superkingdom)],by=c("tronko_name_tmp"="ncbi_name","tronko_rank_tmp"="ncbi_rank"))
    tmp$tronko_name_tmp <- NULL
    tmp$tronko_rank_tmp <- NULL
    TaxonomyExport_NCBI_corrected[[i]] <- tmp
  }
}
TaxonomyExport_NCBI_corrected <- rbind.fill(TaxonomyExport_NCBI_corrected)
TaxonomyExport_NCBI_corrected <- dplyr::left_join(TaxonomyExport_NCBI_corrected,ncbi_taxa_full[,.(ncbi_id,ncbi_name,ncbi_rank)])
#Merge back into the full Tronko-assign output
TaxonomyExport_NCBI <- rbind.fill(TaxonomyExport_NCBI,TaxonomyExport_NCBI_corrected)
TaxonomyExport_NCBI$in_ncbi <- NULL

# Function to find non-ambiguous NCBI to GBIF taxon matches on a rank-by-rank basis
#given only NCBI names as input.
#Create a new data table where if there are more than one entry for TaxonomicRank exists
#for a given unique value of ncbi_TaxonomicRank, then preferentially retain the rows
#where the values for TaxonomicRank and ncbi_TaxonomicRank match.
#TaxonomicRank is a general placeholder for species, genus, etc.
# List of taxonomic ranks
taxonomic_ranks <- c("species", "genus", "family", "order", "class", "phylum")
retain_preferential_rows <- function(data, taxonomic_ranks) {
  # Identify the columns with ncbi_ prefix
  ncbi_cols <- grep("^ncbi_", colnames(data), value = TRUE)
  
  # Process each taxonomic rank
  results <- lapply(taxonomic_ranks, function(rank) {
    ncbi_col <- paste0("ncbi_", rank)
    
    if (ncbi_col %in% ncbi_cols) {
      tmp_filtered <- data %>%
        filter(!is.na(!!sym(ncbi_col))) %>%
        group_by(!!sym(ncbi_col)) %>%
        filter(n() == 1 | !!sym(ncbi_col) == !!sym(rank)) %>%
        ungroup()
      
      return(tmp_filtered)
    }
    return(NULL)
  })
  
  # Combine the results
  final_result <- bind_rows(results)
  return(final_result)
}

#Add in the corresponding GBIF taxonomies, on a rank-by-rank basis to the NCBI inputs.
TaxonomyExport_GBIF <- TaxonomyExport_NCBI
#TaxonomyExport_GBIF <- dplyr::left_join(TaxonomyExport_GBIF,GBIF_NCBI[,c("ncbi_id","species","genus","family","order","class","phylum","kingdom")],by=c("ncbi_id"),na_matches="na")
for(TaxonomicRank in c("species","genus","family","order","class","phylum")){
  select_cols <- c(paste("ncbi_",TaxonomicRank,sep=""),TaxonomicRank)
  tmp <- GBIF_NCBI[,..select_cols]
  tmp <- tmp[!duplicated(tmp),]
  tmp <- retain_preferential_rows(tmp, taxonomic_ranks)
  TaxonomyExport_GBIF <- dplyr::left_join(TaxonomyExport_GBIF,tmp[,c(paste("ncbi_",TaxonomicRank,sep=""),TaxonomicRank)])
}
#Add in GBIF kingdoms
tmp <- GBIF_NCBI[,.(phylum,kingdom)]
tmp <- tmp[complete.cases(tmp),]
tmp <- tmp[!duplicated(tmp),]
TaxonomyExport_GBIF <- dplyr::left_join(TaxonomyExport_GBIF,tmp)

#Determine which eDNA occurrences are found in GBIF in California
#Read in GBIF data sets.
#GBIF.org (25 June 2024) GBIF Occurrence Download https://doi.org/10.15468/dl.u8mz2u
#All taxa in California
#All entries in the initial zip file (0093168-240506114902167.zip) are converted to a
#tab-delimited file of the unique taxonomic ranks (CA_GBIF_All_Taxa.tsv) using the
#bash script CA_GBIF_All.sh: https://github.com/levisimons/CIB/blob/main/CA_GBIF_All.sh
CA_GBIF_All <- fread(input="CA_GBIF_All_Taxa.tsv",sep="\t",na.strings=c("NA",""))

#Assign rank and most resolved name.
GBIF_ranks <- c("species", "genus", "family", "order", "class", "phylum", "kingdom")
tmp <- TaxonomyExport_GBIF[,GBIF_ranks]
TaxonomyExport_GBIF$gbif_rank <- GBIF_ranks[apply(tmp, 1, first_non_na_index)]
TaxonomyExport_GBIF <- TaxonomyExport_GBIF %>%
  mutate(gbif_name = coalesce(species,genus,family,order,class,phylum,kingdom))
#Check if eDNA data shows up in GBIF California data.
TaxonomyExport_GBIF <- setDT(TaxonomyExport_GBIF)
TaxonomyExport_GBIF[, GBIF_in_CA := gbif_name %in% CA_GBIF_All[[gbif_rank]], by = gbif_rank]

#Determine which eDNA occurrences are found in GBIF in California collections
#GBIF.org (25 June 2024) GBIF Occurrence Download https://doi.org/10.15468/dl.e3kttn
#All entries in the initial zip file (0093569-240506114902167.zip) are converted to a
#tab-delimited file of the unique taxonomic ranks ( CA_GBIF_Taxa_Collections.tsv) using the
#bash script CA_GBIF_Collections.sh: https://github.com/levisimons/CIB/blob/main/CA_GBIF_Collections.sh
CA_GBIF_Collections <- fread(input="CA_GBIF_Taxa_Collections.tsv",sep="\t",na.strings=c("NA",""))
TaxonomyExport_GBIF[, GBIF_in_CA_Collections := gbif_name %in% CA_GBIF_Collections[[gbif_rank]], by = gbif_rank]

#Determine which eDNA occurrences are found in iNaturalist in California
#Read in iNaturalist taxon ID to GBIF ID backbone.
iNaturalist_backbone <- fread(input="iNaturalist_backbone.csv",sep=",")
indx <- which(sapply(iNaturalist_backbone, is.character)) 
for (j in indx) set(iNaturalist_backbone, i = grep("^$|^ $", iNaturalist_backbone[[j]]), j = j, value = NA_character_)
#Get corresponding GBIF taxonomies for iNaturalist taxon IDs
if(Primer=="CO1_Metazoa"){iNaturalist_taxa <- fread(input="iNaturalist_Invertebrates.csv",sep=",")}
if(Primer=="ITS1_Fungi"){iNaturalist_taxa <- fread(input="iNaturalist_Fungi.csv",sep=",")}
iNaturalist_backbone <- iNaturalist_backbone[iNaturalist_backbone$id %in% iNaturalist_taxa$taxon_id,]
#Get unique taxa, genus level up, for iNaturalist California data.
iNaturalist_backbone <- iNaturalist_backbone[,.(genus,family,order,class,phylum,kingdom)]
iNaturalist_backbone <- iNaturalist_backbone[!duplicated(iNaturalist_backbone),]
#Determine which eDNA occurrences are found in iNaturalist in California
TaxonomyExport_GBIF$gbif_name_tmp <- ifelse(TaxonomyExport_GBIF$gbif_rank=="species",TaxonomyExport_GBIF$genus,TaxonomyExport_GBIF$gbif_name)
TaxonomyExport_GBIF$gbif_rank_tmp <- ifelse(TaxonomyExport_GBIF$gbif_rank=="species","genus",TaxonomyExport_GBIF$gbif_rank)
TaxonomyExport_GBIF[, iNaturalist_in_CA := gbif_name_tmp %in% iNaturalist_backbone[[gbif_rank_tmp]], by = gbif_rank_tmp]
#Remove temporary columns
TaxonomyExport_GBIF$gbif_name_tmp <- NULL
TaxonomyExport_GBIF$gbif_rank_tmp <- NULL

#Determine which eDNA occurrences are found in BOLD in California
#Tables are generated here: https://github.com/levisimons/CIB/blob/main/CIB_BOLD.R
if(Primer=="CO1_Metazoa"){BOLD_Input <- fread(input="BOLD_Invertebrates.tsv",sep="\t")}
if(Primer=="ITS1_Fungi"){BOLD_Input <- fread(input="BOLD_Fungi.tsv",sep="\t")}
TaxonomyExport_GBIF[, BOLD_in_CA := gbif_name %in% BOLD_Input[[gbif_rank]], by = gbif_rank]

#Set eDNA to observation/collection comparison columns to numeric
TaxonomyExport_GBIF$GBIF_in_CA <- as.integer(TaxonomyExport_GBIF$GBIF_in_CA)
TaxonomyExport_GBIF$GBIF_in_CA_Collections <- as.integer(TaxonomyExport_GBIF$GBIF_in_CA_Collections)
TaxonomyExport_GBIF$iNaturalist_in_CA <- as.integer(TaxonomyExport_GBIF$iNaturalist_in_CA)
TaxonomyExport_GBIF$BOLD_in_CA <- as.integer(TaxonomyExport_GBIF$BOLD_in_CA)

#Add in check if taxa are detected via eDNA, but not yet found in a GBIF collection or the BOLD database.
TaxonomyExport_GBIF$Uncollected <- ifelse(TaxonomyExport_GBIF$GBIF_in_CA_Collections==0 & TaxonomyExport_GBIF$BOLD_in_CA==0,1,0)

#Retain only fungal or invertebrate results
if(Primer=="CO1_Metazoa"){TaxonomyExport_GBIF <- TaxonomyExport_GBIF[TaxonomyExport_GBIF$kingdom=="Animalia",]}
if(Primer=="ITS1_Fungi"){TaxonomyExport_GBIF <- TaxonomyExport_GBIF[TaxonomyExport_GBIF$kingdom=="Fungi",]}

#Export table.
write.table(TaxonomyExport_GBIF,paste(Primer,"_eDNAPrevalence_GBIF.tsv",sep=""),quote=FALSE,sep="\t",row.names = FALSE)
