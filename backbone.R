rm(list=ls())
require(dplyr)
require(plyr)
require(tidyr)
require(phyloseq)
require(decontam)
require(data.table)
require(stringr)

wd <- ""

setwd(wd)

##Match eDNA data with NCBI taxonomies to the closest match in the GBIF backbone.
#Read in GBIF - NCBI backbone generated via https://github.com/eDNA-Explorer/eDNAExplorer/blob/main/Backbone_generator.R
backbone <- fread(input="GBIF_NCBI.csv",sep=",")
names(backbone)[names(backbone) == 'kingdom'] <- 'gbif_kingdom'
names(backbone)[names(backbone) == 'phylum'] <- 'gbif_phylum'
names(backbone)[names(backbone) == 'class'] <- 'gbif_class'
names(backbone)[names(backbone) == 'order'] <- 'gbif_order'
names(backbone)[names(backbone) == 'family'] <- 'gbif_family'
names(backbone)[names(backbone) == 'genus'] <- 'gbif_genus'
names(backbone)[names(backbone) == 'species'] <- 'gbif_species'
backbone <- backbone %>% filter(!is.na(ncbi_kingdom) | !is.na(ncbi_phylum) | !is.na(ncbi_class) | !is.na(ncbi_order) | !is.na(ncbi_family) | !is.na(ncbi_genus) | !is.na(ncbi_species))
backbone <- backbone %>% distinct(taxonID, .keep_all = TRUE) %>% distinct(ncbi_id, .keep_all = TRUE)
backbone <- as.data.frame(backbone)

#Set lists for filtering taxa.
Invertebrate_Phyla <- c("Platyhelminthes","Turbellaria","Trematoda","Cestoda","Nemertea","Nemertea","Rotifera","Gastrotricha","Acanthocephala","Nematoda","Nematomorpha","Priapulida","Kinorhyncha","Loricifera","Entoprocta","Cycliophora","Gnathostomulida","Micrognathozoa","Chaetognatha","Hemichordata","Bryozoa","Brachiopoda","Phoronida","Ectoprocta","Annelida","Polychaeta","Clitellata","Mollusca","Gastropoda","Bivalvia","Arthropoda","Insecta","Arachnida","Crustacea","Myriapoda")
Fungal_Phyla <- c("Ascomycota","Basidiomycota","Blastocladiomycota","Chytridiomycota","Entorrhizomycota","Glomeromycota","Mucoromycota","Neocallimastigomycota","Sanchytriomycota","Zoopagomycota","Zygomycota")

#Select a primer
Primers <- c("12S_MiFish_U","16S_Bacteria","18S_Euk","CO1_Metazoa","ITS1_Fungi","ITS2_Plants","vert12S")
Primer <- "CO1_Metazoa"

##Read in and process Tronko-assign eDNA data.
#Get project directories
Project_Directories <- list.dirs(path=paste(wd,"CALeDNA",sep="/"),recursive=F)

#Aggregate tronko-assign results to a particular taxonomic level.
TaxonomicRanks <- c("superkingdom","phylum","class","order","family","genus","species")

#Select a tronko-assign mismatch cutoff.
Mismatch <- 25
#Get the number of unique taxa before and after decontamination.
TaxaBySamples_Unfiltered_NCBI <- c()
TaxaBySamples_Filtered_NCBI <- c()
TaxaBySamples_Unfiltered_GBIF <- c()
TaxaBySamples_Filtered_GBIF <- c()
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
      contamdf.prev_threshold <- isContaminant(physeq_initial, method="prevalence", neg="is.neg", threshold=0.1)
      #Filter on threshold
      contamdf.prev_threshold_filtered <-prune_taxa(!contamdf.prev_threshold$contaminant, physeq_initial)
      #Filter on taxa prevalence
      #contamdf.prev_prevalence_filtered = filter_taxa(contamdf.prev_threshold_filtered,function(x) sum(x > 2) > 1, TRUE)
      contamdf.prev_no_blanks = subset_samples(contamdf.prev_threshold_filtered, Sample_or_Control == "True Sample")
      contamdf.prev_min_threshold_a = prune_taxa(taxa_sums(contamdf.prev_no_blanks) > 2, contamdf.prev_no_blanks)
      contamdf.prev_min_threshold_p = prune_samples(sample_sums(contamdf.prev_min_threshold_a)>=2, contamdf.prev_min_threshold_a)
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
      
      #Find the number of unique taxa per project, with no filtering, using NCBI taxonomy
      Taxonomy_Unfiltered <- as.data.frame(tax_table(physeq_initial))
      Taxonomy_Unfiltered <-Taxonomy_Unfiltered %>% mutate_all(~ifelse(. == "NA", NA, .))
      if(Primer=="CO1_Metazoa"){
        Taxonomy_Unfiltered <- Taxonomy_Unfiltered[Taxonomy_Unfiltered$phylum %in% Invertebrate_Phyla,]
      }
      if(Primer=="ITS1_Fungi"){
        Taxonomy_Unfiltered <- Taxonomy_Unfiltered[Taxonomy_Unfiltered$phylum %in% Fungal_Phyla,]
      }
      row_unfiltered <- data.frame(matrix(nrow=1,ncol=7))
      colnames(row_unfiltered) <- c("unique_superkingdom","unique_phylum","unique_class","unique_order","unique_family","unique_genus","unique_species")
      for(TaxonomicRank in TaxonomicRanks){
        row_unfiltered[,paste("unique_",TaxonomicRank,sep="")] <- length(na.omit(unique(Taxonomy_Unfiltered[,TaxonomicRank])))
      }
      row_unfiltered$ProjectID <- ProjectID
      row_unfiltered$Num_Samples <- nrow(MetadataTable)
      row_unfiltered$Num_Samples_Filtered <- nrow(Metadata_Filtered)
      row_unfiltered$Filter_Status <- "Unfiltered"
      row_unfiltered$Taxonomy <- "NCBI"
      #Find the number of unique taxa per project, with filtering, using NCBI taxonomy
      Taxonomy_Filtered <- as.data.frame(tax_table(contamdf.prev_min_threshold_p))
      Taxonomy_Filtered <-Taxonomy_Filtered %>% mutate_all(~ifelse(. == "NA", NA, .))
      if(Primer=="CO1_Metazoa"){
        Taxonomy_Filtered <- Taxonomy_Filtered[Taxonomy_Filtered$phylum %in% Invertebrate_Phyla,]
      }
      if(Primer=="ITS1_Fungi"){
        Taxonomy_Filtered <- Taxonomy_Filtered[Taxonomy_Filtered$phylum %in% Fungal_Phyla,]
      }
      row_filtered <- data.frame(matrix(nrow=1,ncol=7))
      colnames(row_filtered) <- c("unique_superkingdom","unique_phylum","unique_class","unique_order","unique_family","unique_genus","unique_species")
      for(TaxonomicRank in TaxonomicRanks){
        row_filtered[,paste("unique_",TaxonomicRank,sep="")] <- length(na.omit(unique(Taxonomy_Filtered[,TaxonomicRank])))
      }
      row_filtered$ProjectID <- ProjectID
      row_filtered$Num_Samples <- nrow(MetadataTable)
      row_filtered$Num_Samples_Filtered <- nrow(Metadata_Filtered)
      row_filtered$Filter_Status <- "Filtered"
      row_filtered$Taxonomy <- "NCBI"
      #Store results.
      TaxaBySamples_Unfiltered_NCBI[[i]] <- row_unfiltered
      TaxaBySamples_Filtered_NCBI[[i]] <- row_filtered
      #Find the number of unique taxa per project, with no filtering, using GBIF taxonomy
      Taxonomy_Unfiltered <- as.data.frame(tax_table(physeq_initial))
      Taxonomy_Unfiltered <-Taxonomy_Unfiltered %>% mutate_all(~ifelse(. == "NA", NA, .))
      if(Primer=="CO1_Metazoa"){
        Taxonomy_Unfiltered <- Taxonomy_Unfiltered[Taxonomy_Unfiltered$phylum %in% Invertebrate_Phyla,]
      }
      if(Primer=="ITS1_Fungi"){
        Taxonomy_Unfiltered <- Taxonomy_Unfiltered[Taxonomy_Unfiltered$phylum %in% Fungal_Phyla,]
      }
      #Add in GBIF taxonomies linked to NCBI taxonomies in eDNA data.
      Taxonomy_Unfiltered <- dplyr::left_join(Taxonomy_Unfiltered,backbone[backbone$ncbi_rank=="species",c("ncbi_species","gbif_species")],by=c("species"="ncbi_species"))
      Taxonomy_Unfiltered <- dplyr::left_join(Taxonomy_Unfiltered,backbone[backbone$ncbi_rank=="genus",c("ncbi_genus","gbif_genus")],by=c("genus"="ncbi_genus"))
      Taxonomy_Unfiltered <- dplyr::left_join(Taxonomy_Unfiltered,backbone[backbone$ncbi_rank=="family",c("ncbi_family","gbif_family")],by=c("family"="ncbi_family"))
      Taxonomy_Unfiltered <- dplyr::left_join(Taxonomy_Unfiltered,backbone[backbone$ncbi_rank=="order",c("ncbi_order","gbif_order")],by=c("order"="ncbi_order"))
      Taxonomy_Unfiltered <- dplyr::left_join(Taxonomy_Unfiltered,backbone[backbone$ncbi_rank=="class",c("ncbi_class","gbif_class")],by=c("class"="ncbi_class"))
      Taxonomy_Unfiltered <- dplyr::left_join(Taxonomy_Unfiltered,backbone[backbone$ncbi_rank=="phylum",c("ncbi_phylum","gbif_phylum")],by=c("phylum"="ncbi_phylum"))
      row_unfiltered <- data.frame(matrix(nrow=1,ncol=7))
      colnames(row_unfiltered) <- c("unique_superkingdom","unique_phylum","unique_class","unique_order","unique_family","unique_genus","unique_species")
      for(TaxonomicRank in TaxonomicRanks){
        if(TaxonomicRank!="superkingdom"){row_unfiltered[,paste("unique_",TaxonomicRank,sep="")] <- length(na.omit(unique(Taxonomy_Unfiltered[,paste("gbif_",TaxonomicRank,sep="")])))}
        if(TaxonomicRank=="superkingdom"){row_unfiltered[,paste("unique_",TaxonomicRank,sep="")] <- NA}
      }
      row_unfiltered$ProjectID <- ProjectID
      row_unfiltered$Num_Samples <- nrow(MetadataTable)
      row_unfiltered$Num_Samples_Filtered <- nrow(Metadata_Filtered)
      row_unfiltered$Filter_Status <- "Unfiltered"
      row_unfiltered$Taxonomy <- "GBIF"
      #Find the number of unique taxa per project, with no filtering, using GBIF taxonomy
      Taxonomy_Filtered <- as.data.frame(tax_table(contamdf.prev_min_threshold_p))
      Taxonomy_Filtered <-Taxonomy_Filtered %>% mutate_all(~ifelse(. == "NA", NA, .))
      if(Primer=="CO1_Metazoa"){
        Taxonomy_Filtered <- Taxonomy_Filtered[Taxonomy_Filtered$phylum %in% Invertebrate_Phyla,]
      }
      if(Primer=="ITS1_Fungi"){
        Taxonomy_Filtered <- Taxonomy_Filtered[Taxonomy_Filtered$phylum %in% Fungal_Phyla,]
      }
      #Add in GBIF taxonomies linked to NCBI taxonomies in eDNA data.
      Taxonomy_Filtered <- dplyr::left_join(Taxonomy_Filtered,backbone[backbone$ncbi_rank=="species",c("ncbi_species","gbif_species")],by=c("species"="ncbi_species"))
      Taxonomy_Filtered <- dplyr::left_join(Taxonomy_Filtered,backbone[backbone$ncbi_rank=="genus",c("ncbi_genus","gbif_genus")],by=c("genus"="ncbi_genus"))
      Taxonomy_Filtered <- dplyr::left_join(Taxonomy_Filtered,backbone[backbone$ncbi_rank=="family",c("ncbi_family","gbif_family")],by=c("family"="ncbi_family"))
      Taxonomy_Filtered <- dplyr::left_join(Taxonomy_Filtered,backbone[backbone$ncbi_rank=="order",c("ncbi_order","gbif_order")],by=c("order"="ncbi_order"))
      Taxonomy_Filtered <- dplyr::left_join(Taxonomy_Filtered,backbone[backbone$ncbi_rank=="class",c("ncbi_class","gbif_class")],by=c("class"="ncbi_class"))
      Taxonomy_Filtered <- dplyr::left_join(Taxonomy_Filtered,backbone[backbone$ncbi_rank=="phylum",c("ncbi_phylum","gbif_phylum")],by=c("phylum"="ncbi_phylum"))
      row_filtered <- data.frame(matrix(nrow=1,ncol=7))
      colnames(row_filtered) <- c("unique_superkingdom","unique_phylum","unique_class","unique_order","unique_family","unique_genus","unique_species")
      for(TaxonomicRank in TaxonomicRanks){
        if(TaxonomicRank!="superkingdom"){row_filtered[,paste("unique_",TaxonomicRank,sep="")] <- length(na.omit(unique(Taxonomy_Filtered[,paste("gbif_",TaxonomicRank,sep="")])))}
        if(TaxonomicRank=="superkingdom"){row_filtered[,paste("unique_",TaxonomicRank,sep="")] <- NA}
      }
      row_filtered$ProjectID <- ProjectID
      row_filtered$Num_Samples <- nrow(MetadataTable)
      row_filtered$Num_Samples_Filtered <- nrow(Metadata_Filtered)
      row_filtered$Filter_Status <- "Filtered"
      row_filtered$Taxonomy <- "GBIF"
      #Store results
      TaxaBySamples_Unfiltered_GBIF[[i]] <- row_unfiltered
      TaxaBySamples_Filtered_GBIF[[i]] <- row_filtered
    }
  }
  i=i+1
  print(paste(i,TronkoTable_Path))
}

#Aggregate summary statistics
TaxaBySamples_Unfiltered_NCBI <- rbind.fill(TaxaBySamples_Unfiltered_NCBI)
TaxaBySamples_Filtered_NCBI <- rbind.fill(TaxaBySamples_Filtered_NCBI)
TaxaBySamples_Unfiltered_GBIF <- rbind.fill(TaxaBySamples_Unfiltered_GBIF)
TaxaBySamples_Filtered_GBIF <- rbind.fill(TaxaBySamples_Filtered_GBIF)
TaxaBySamples <- rbind(TaxaBySamples_Unfiltered_NCBI,TaxaBySamples_Filtered_NCBI,TaxaBySamples_Unfiltered_GBIF,TaxaBySamples_Filtered_GBIF)
#write.table(TaxaBySamples,paste(Primer,"_UniqueTaxa.csv",sep=""),quote=FALSE,sep=",",row.names = FALSE)

#Visualize unique taxa by projects
require(ggplot2)
ggplot(TaxaBySamples,aes(ProjectID,unique_species))+geom_point()+
  scale_y_continuous(trans='log10')+
  facet_grid(Filter_Status ~ Taxonomy)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

#Export total taxonomy table, using NCBI taxonomy
if(Primer=="CO1_Metazoa"){taxa_retain <- taxa_names(physeq_total)[tax_table(physeq_total)[, "phylum"] %in% Invertebrate_Phyla]}
if(Primer=="ITS1_Fungi"){taxa_retain <- taxa_names(physeq_total)[tax_table(physeq_total)[, "phylum"] %in% Fungal_Phyla]}
physeq_retain <- prune_taxa(taxa_retain, physeq_total)
OTUTable_NCBI <- as.data.frame(otu_table(physeq_retain))
#Convert read counts to presence/absence.
OTUTable_NCBI <- OTUTable_NCBI %>% mutate_all(~ ifelse(. > 0, 1, 0))
#Calculate eDNA taxa prevalence
OTUTable_NCBI$eDNA_prevalence <- rowSums(OTUTable_NCBI)
#Get taxa by ranks
TaxonomyExport_NCBI <- OTUTable_NCBI
TaxonomyExport_NCBI$Taxon <- rownames(TaxonomyExport_NCBI)
TaxonomyExport_NCBI <- separate(data = TaxonomyExport_NCBI, col = "Taxon", sep=";",into = TaxonomicRanks)
TaxonomyExport_NCBI <- TaxonomyExport_NCBI[,c(TaxonomicRanks,"eDNA_prevalence")]
TaxonomyExport_NCBI <- TaxonomyExport_NCBI %>% mutate_all(~ifelse(. == "NA", NA, .))
#write.table(TaxonomyExport_NCBI,paste(Primer,"_eDNAPrevalence_NCBI.csv",sep=""),quote=FALSE,sep=",",row.names = FALSE)
#Export total taxonomy table, using GBIF taxonomy
TaxonomyExport_GBIF <- TaxonomyExport_NCBI
#Add in GBIF taxonomies linked to NCBI taxonomies in eDNA data.
TaxonomyExport_GBIF <- dplyr::left_join(TaxonomyExport_GBIF,backbone[backbone$ncbi_rank=="species",c("ncbi_species","gbif_species","taxonID")],by=c("species"="ncbi_species"))
names(TaxonomyExport_GBIF)[names(TaxonomyExport_GBIF) == "taxonID"] <- "species_id" 
TaxonomyExport_GBIF <- dplyr::left_join(TaxonomyExport_GBIF,backbone[backbone$ncbi_rank=="genus",c("ncbi_genus","gbif_genus","taxonID")],by=c("genus"="ncbi_genus"))
names(TaxonomyExport_GBIF)[names(TaxonomyExport_GBIF) == "taxonID"] <- "genus_id" 
TaxonomyExport_GBIF <- dplyr::left_join(TaxonomyExport_GBIF,backbone[backbone$ncbi_rank=="family",c("ncbi_family","gbif_family","taxonID")],by=c("family"="ncbi_family"))
names(TaxonomyExport_GBIF)[names(TaxonomyExport_GBIF) == "taxonID"] <- "family_id"
TaxonomyExport_GBIF <- dplyr::left_join(TaxonomyExport_GBIF,backbone[backbone$ncbi_rank=="order",c("ncbi_order","gbif_order","taxonID")],by=c("order"="ncbi_order"))
names(TaxonomyExport_GBIF)[names(TaxonomyExport_GBIF) == "taxonID"] <- "order_id"
TaxonomyExport_GBIF <- dplyr::left_join(TaxonomyExport_GBIF,backbone[backbone$ncbi_rank=="class",c("ncbi_class","gbif_class","taxonID")],by=c("class"="ncbi_class"))
names(TaxonomyExport_GBIF)[names(TaxonomyExport_GBIF) == "taxonID"] <- "class_id"
TaxonomyExport_GBIF <- dplyr::left_join(TaxonomyExport_GBIF,backbone[backbone$ncbi_rank=="phylum",c("ncbi_phylum","gbif_phylum","taxonID")],by=c("phylum"="ncbi_phylum"))
names(TaxonomyExport_GBIF)[names(TaxonomyExport_GBIF) == "taxonID"] <- "phylum_id"
TaxonomyExport_GBIF <- TaxonomyExport_GBIF %>%
  mutate(taxonID = coalesce(species_id, genus_id, family_id, order_id, class_id, phylum_id))

#Relabel NCBI taxonomic ranks
for(TaxonomicRank in TaxonomicRanks){
  names(TaxonomyExport_GBIF)[names(TaxonomyExport_GBIF) == TaxonomicRank] <- paste("ncbi_",TaxonomicRank,sep="") 
}
#Determine which eDNA occurrences are found in GBIF in California
#Read in GBIF data sets.
#GBIF.org (03 May 2024) GBIF Occurrence Download  https://doi.org/10.15468/dl.j6wbeq
#All taxa in California
CA_GBIF_All <- fread(input="CA_GBIF_Taxa_Full.csv",sep="\t")
CA_GBIF_All <- CA_GBIF_All %>% mutate_all(~ifelse(. == "NA", NA, .))
CA_GBIF_All <- CA_GBIF_All %>% mutate_all(~ifelse(. == "", NA, .))
#Check if taxa, according to their GBIF taxonomy, were found in GBIF records for California
GBIF_all_taxa <- na.omit(unique(unlist(CA_GBIF_All[,c("acceptedTaxonKey","phylumKey","classKey","orderKey","familyKey","genusKey","speciesKey")])))
TaxonomyExport_GBIF$GBIF_in_CA <- ifelse(TaxonomyExport_GBIF$taxonID %in% GBIF_all_taxa,1,0)

#Determine which eDNA occurrences are found in GBIF in California collections
#GBIF.org (03 May 2024) GBIF Occurrence Download  https://doi.org/10.15468/dl.2jwkqg
#All taxa in California.  Observation categories: Preserved Specimen, Fossil Specimen, Living Specimen.
CA_GBIF_Collections <- fread(input="CA_GBIF_Taxa_Collections.csv",sep="\t")
CA_GBIF_Collections <- CA_GBIF_Collections %>% mutate_all(~ifelse(. == "NA", NA, .))
CA_GBIF_Collections <- CA_GBIF_Collections %>% mutate_all(~ifelse(. == "", NA, .))
#Check if taxa, according to their GBIF taxonomy, were found in GBIF collection records for California
GBIF_collections_taxa <- na.omit(unique(unlist(CA_GBIF_Collections[,c("acceptedTaxonKey","phylumKey","classKey","orderKey","familyKey","genusKey","speciesKey")])))
TaxonomyExport_GBIF$GBIF_in_CA_Collections <- ifelse(TaxonomyExport_GBIF$taxonID %in% GBIF_collections_taxa,1,0)

#Determine which eDNA occurrences are found in iNaturalist in California
#Read in iNaturalist taxon ID to GBIF ID backbone.
iNaturalist_backbone <- fread(input="iNaturalist_taxa.csv",sep=",")
indx <- which(sapply(iNaturalist_backbone, is.character)) 
for (j in indx) set(iNaturalist_backbone, i = grep("^$|^ $", iNaturalist_backbone[[j]]), j = j, value = NA_character_)
#Clean names
iNaturalist_backbone <- iNaturalist_backbone %>% dplyr::mutate(specificEpithet = gsub("[^[:alnum:]]", "", specificEpithet))
iNaturalist_backbone <- iNaturalist_backbone %>% dplyr::mutate(specificEpithet = gsub("[^[:alnum:]]", "", specificEpithet))
#Split iNaturalist taxa into entries with standard and non-standard taxonomic ranks
iNaturalist_backbone_standard <- iNaturalist_backbone[iNaturalist_backbone$taxonRank %in% TaxonomicRanks,]
iNaturalist_backbone_nonstandard <- iNaturalist_backbone[!(iNaturalist_backbone$taxonRank %in% TaxonomicRanks),]

#Read in Open Tree of Life backbone https://tree.opentreeoflife.org/about/taxonomy-version/ott3.5
OTL_taxonomy <- fread(input="otl_taxonomy.tsv",sep="\t")
#Filter by matching taxonomic names with GBIF entries at the rank of genus and above.
OTL_taxonomy_subset <- OTL_taxonomy[OTL_taxonomy$name %in% na.omit(unique(iNaturalist_backbone_nonstandard$scientificName)),c("name","rank","sourceinfo")]
OTL_bridge <- OTL_taxonomy_subset %>% filter(grepl("gbif",sourceinfo))
OTL_bridge <- OTL_bridge %>% mutate(taxonID = str_extract(sourceinfo, "(?<=gbif:)\\d+"))
OTL_bridge$taxonID <- as.numeric(OTL_bridge$taxonID)
OTL_bridge <- OTL_bridge[,c("name","rank","taxonID")]
OTL_bridge <- OTL_bridge[!duplicated(OTL_bridge),]
OTL_bridge <- OTL_bridge[OTL_bridge$taxonID %in% backbone$taxonID & OTL_bridge$rank %in% c("species","genus","family","order","class","phylum"),]
#Add in GBIF IDs to match taxonomic names and ranks .
iNaturalist_backbone_nonstandard$taxonID <- NULL
iNaturalist_backbone_nonstandard <- dplyr::left_join(iNaturalist_backbone_nonstandard,OTL_bridge,by=c("scientificName"="name"))

#Add in GBIF IDs to match taxonomic names and ranks.
tmp <- backbone[backbone$taxonRank=="species",c("gbif_species","speciesKey")]
tmp <- tmp[!duplicated(tmp),]
iNaturalist_backbone_standard <- dplyr::left_join(iNaturalist_backbone_standard,tmp,by=c("scientificName"="gbif_species"))
tmp <- backbone[backbone$taxonRank=="genus",c("gbif_genus","genusKey")]
tmp <- tmp[!duplicated(tmp),]
iNaturalist_backbone_standard <- dplyr::left_join(iNaturalist_backbone_standard,tmp,by=c("genus"="gbif_genus"))
tmp <- backbone[backbone$taxonRank=="family",c("gbif_family","familyKey")]
tmp <- tmp[!duplicated(tmp),]
iNaturalist_backbone_standard <- dplyr::left_join(iNaturalist_backbone_standard,tmp,by=c("family"="gbif_family"))
tmp <- backbone[backbone$taxonRank=="order",c("gbif_order","orderKey")]
tmp <- tmp[!duplicated(tmp),]
iNaturalist_backbone_standard <- dplyr::left_join(iNaturalist_backbone_standard,tmp,by=c("order"="gbif_order"))
tmp <- backbone[backbone$taxonRank=="class",c("gbif_class","classKey")]
tmp <- tmp[!duplicated(tmp),]
iNaturalist_backbone_standard <- dplyr::left_join(iNaturalist_backbone_standard,tmp,by=c("class"="gbif_class"))
tmp <- backbone[backbone$taxonRank=="phylum",c("gbif_phylum","phylumKey")]
tmp <- tmp[!duplicated(tmp),]
iNaturalist_backbone_standard <- dplyr::left_join(iNaturalist_backbone_standard,tmp,by=c("phylum"="gbif_phylum"))
iNaturalist_backbone_standard <- as.data.frame(iNaturalist_backbone_standard)
#Check if taxa, according to their GBIF taxonomy, were found in iNaturalist's records within California
iNaturalist_taxa <- c(na.omit(unique(unlist(iNaturalist_backbone_standard[,c("phylumKey","classKey","orderKey","familyKey","genusKey","speciesKey")]))),na.omit(unique(iNaturalist_backbone_nonstandard$taxonID)))
TaxonomyExport_tmp <- TaxonomyExport_GBIF[,c("gbif_species","gbif_genus","gbif_family","gbif_order","gbif_class","gbif_phylum")]
for(TaxonomicRank in c("species","genus","family","order","class","phylum")){
  tmp <- backbone[backbone$taxonRank==TaxonomicRank,c(paste("gbif_",TaxonomicRank,sep=""),paste(TaxonomicRank,"Key",sep=""))]
  tmp <- tmp[!duplicated(tmp),]
  TaxonomyExport_tmp <- dplyr::left_join(TaxonomyExport_tmp,tmp)
}
TaxonomyExport_tmp <- TaxonomyExport_tmp %>%
  mutate(taxonID = coalesce(speciesKey,genusKey,familyKey,orderKey,classKey,phylumKey))
TaxonomyExport_tmp$iNaturalist_in_CA <- ifelse(TaxonomyExport_tmp$taxonID %in% iNaturalist_taxa,1,0)
TaxonomyExport_tmp <- TaxonomyExport_tmp[,c("taxonID","iNaturalist_in_CA")]
TaxonomyExport_tmp <- TaxonomyExport_tmp[!duplicated(TaxonomyExport_tmp),]
TaxonomyExport_GBIF <- dplyr::left_join(TaxonomyExport_GBIF,TaxonomyExport_tmp[,c("taxonID","iNaturalist_in_CA")])
#Check for any taxa which may be in iNaturalist, but may initially be mismatched due to taxonomic synonym issues.
TaxonomyExport_GBIF_recheck <- TaxonomyExport_GBIF[is.na(TaxonomyExport_GBIF$iNaturalist_in_CA) | TaxonomyExport_GBIF$iNaturalist_in_CA==0,]
TaxonomyExport_GBIF_recheck$iNaturalist_in_CA <- NULL
TaxonomyExport_GBIF_recheck <- TaxonomyExport_GBIF_recheck %>%
  mutate(Taxon = coalesce(gbif_genus,gbif_family,gbif_order,gbif_class,gbif_phylum))
TaxonomyExport_GBIF_recheck$iNaturalist_in_CA <- ifelse(TaxonomyExport_GBIF_recheck$Taxon %in% na.omit(unique(iNaturalist_backbone$scientificName)),1,0)
TaxonomyExport_GBIF_recheck <- TaxonomyExport_GBIF_recheck[,colnames(TaxonomyExport_GBIF)]
TaxonomyExport_GBIF <- TaxonomyExport_GBIF[!is.na(TaxonomyExport_GBIF$iNaturalist_in_CA) & TaxonomyExport_GBIF$iNaturalist_in_CA==1,]
TaxonomyExport_GBIF <- rbind(TaxonomyExport_GBIF,TaxonomyExport_GBIF_recheck)
TaxonomyExport_GBIF <- TaxonomyExport_GBIF[,c("ncbi_species","ncbi_genus","ncbi_family","ncbi_order","ncbi_class","ncbi_phylum","taxonID","gbif_species","gbif_genus","gbif_family","gbif_order","gbif_class","gbif_phylum","eDNA_prevalence","GBIF_in_CA","GBIF_in_CA_Collections","iNaturalist_in_CA")]
#Export table
write.table(TaxonomyExport_GBIF,paste(Primer,"_eDNAPrevalence_GBIF.csv",sep=""),quote=FALSE,sep=",",row.names = FALSE)
