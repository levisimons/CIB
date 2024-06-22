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
Invertebrate_Phyla <- c("Platyhelminthes","Turbellaria","Trematoda","Cestoda","Nemertea","Rotifera","Gastrotricha","Acanthocephala","Nematoda","Nematomorpha","Priapulida","Kinorhyncha","Loricifera","Entoprocta","Cycliophora","Gnathostomulida","Micrognathozoa","Chaetognatha","Hemichordata","Bryozoa","Brachiopoda","Phoronida","Annelida","Polychaeta","Clitellata","Mollusca","Gastropoda","Bivalvia","Arthropoda","Insecta","Arachnida","Crustacea","Myriapoda")
Fungal_Phyla <- c("Ascomycota", "Basidiomycota", "Entorrhizomycetes", "Blastocladiomycota", "Chytridiomycota", "Cryptomycota", "Microsporidia", "Mucoromycota", "Nephridiophaga", "Olpidiomycota", "Sanchytriomycota", "Zoopagomycota")

#Select a primer
Primers <- c("12S_MiFish_U","16S_Bacteria","18S_Euk","CO1_Metazoa","ITS1_Fungi","ITS2_Plants","vert12S")
Primer <- "CO1_Metazoa"

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
      #Save stats on the number samples and taxa per project before and after decontamination.
      project_stat <- data.frame(matrix(nrow=1,ncol=5))
      colnames(project_stat) <- c("ProjectID","Number_samples_inital","Number_taxa_initial","Number_samples_decontaminated","Number_taxa_decontaminated")
      project_stat$ProjectID <- ProjectID
      project_stat$Number_samples_inital <- nsamples(physeq_initial)
      project_stat$Number_taxa_initial <- ntaxa(physeq_initial)
      project_stat$Number_samples_decontaminated <- nsamples(contamdf.prev_min_threshold_p)
      project_stat$Number_taxa_decontaminated <- ntaxa(contamdf.prev_min_threshold_p)
      project_stats[[i]] <- project_stat
      #Export data on projects before and after decontamination.
      saveRDS(physeq_initial, paste(ProjectID,"_",Primer,"_initial_phyloseq.rds",sep=""))
      saveRDS(contamdf.prev_min_threshold_p, paste(ProjectID,"_",Primer,"_decontaminated_phyloseq.rds",sep=""))
      write.table(data.frame("sum.taxonomy"=rownames(otu_table(physeq_initial)),otu_table(physeq_initial)),paste(ProjectID,"_",Primer,"_initial_tax_table.txt",sep=""), row.names=FALSE, sep="\t",quote = FALSE)
      write.table(data.frame("sum.taxonomy"=rownames(otu_table(contamdf.prev_min_threshold_p)),otu_table(contamdf.prev_min_threshold_p)),paste(ProjectID,"_",Primer,"_decontaminated_tax_table.txt",sep=""), row.names=FALSE, sep="\t",quote = FALSE)
      write.table(as.data.frame(sample_data(physeq_initial)),paste(ProjectID,"_",Primer,"_initial_metadata.txt",sep=""), row.names=FALSE, sep="\t",quote = FALSE)
      write.table(as.data.frame(sample_data(contamdf.prev_min_threshold_p)),paste(ProjectID,"_",Primer,"_decontaminated_metadata.txt",sep=""), row.names=FALSE, sep="\t",quote = FALSE)
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
#Save merged phyloseq results
if(Primer=="CO1_Metazoa"){
  saveRDS(physeq_retain, "AllInvertebrates.rds")
  write.table(data.frame("sum.taxonomy"=rownames(otu_table(physeq_retain)),otu_table(physeq_retain)),"AllInvertebrates_tax_table.txt", row.names=FALSE, sep="\t",quote = FALSE)
  write.table(as.data.frame(sample_data(physeq_retain)),"AllInvertebrates_metadata.txt", row.names=FALSE, sep="\t",quote = FALSE)
}
if(Primer=="ITS1_Fungi"){
  saveRDS(physeq_retain, "AllFungi.rds")
  write.table(data.frame("sum.taxonomy"=rownames(otu_table(physeq_retain)),otu_table(physeq_retain)),"AllFungi_tax_table.txt", row.names=FALSE, sep="\t",quote = FALSE)
  write.table(as.data.frame(sample_data(physeq_retain)),"AllFungi_metadata.txt", row.names=FALSE, sep="\t",quote = FALSE)
}
OTUTable_NCBI <- as.data.frame(otu_table(physeq_retain))
#Convert read counts to presence/absence.
OTUTable_NCBI <- OTUTable_NCBI %>% mutate_all(~ ifelse(. > 0, 1, 0))
#Calculate eDNA taxa prevalence
OTUTable_NCBI$eDNA_prevalence <- rowSums(OTUTable_NCBI)
#Get taxa by ranks
TaxonomyExport_NCBI <- OTUTable_NCBI
TaxonomyExport_NCBI$Taxon <- rownames(TaxonomyExport_NCBI)
TaxonomyExport_NCBI <- separate(data = TaxonomyExport_NCBI, col = "Taxon", sep=";",into = paste0("ncbi_",TaxonomicRanks,sep=""))
TaxonomyExport_NCBI$eDNA_prevalence_NCBI <- rowSums(TaxonomyExport_NCBI[,!(colnames(TaxonomyExport_NCBI) %in% paste0("ncbi_",TaxonomicRanks,sep=""))])
TaxonomyExport_NCBI <- TaxonomyExport_NCBI[,c(paste0("ncbi_",TaxonomicRanks,sep=""),"eDNA_prevalence_NCBI")]
TaxonomyExport_NCBI <- TaxonomyExport_NCBI %>% mutate_all(~ifelse(. == "NA", NA, .))
#Assign taxon and rank explicitly.
TaxonomyExport_NCBI <- TaxonomyExport_NCBI %>% mutate(ncbi_name = coalesce(ncbi_species, ncbi_genus, ncbi_family, ncbi_order, ncbi_class, ncbi_phylum, ncbi_superkingdom))
TaxonomyExport_NCBI <- TaxonomyExport_NCBI %>% mutate(ncbi_rank = case_when(!is.na(ncbi_species) ~ "species",!is.na(ncbi_genus) ~ "genus",!is.na(ncbi_family) ~ "family",!is.na(ncbi_order) ~ "order",!is.na(ncbi_class) ~ "class",!is.na(ncbi_phylum) ~ "phylum",!is.na(ncbi_superkingdom) ~ "superkingdom",TRUE ~ NA_character_))

#Read in GBIF-NCBI backbone
GBIF_NCBI <- fread(input="~/Desktop/Archive/backbone/GBIF_NCBI_export.tsv",sep="\t")
#Export total taxonomy table, using GBIF taxonomy
TaxonomyExport_GBIF <- TaxonomyExport_NCBI

# Function to find non-ambiguous NCBI to GBIF taxon matches on a rank-by-rank basis
#given only NCBI names as input.
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
for(TaxonomicRank in c("species","genus","family","order","class","phylum")){
  #select_cols <- c(paste("ncbi_",TaxonomicRank,sep=""),paste("ncbi_",TaxonomicRank,"_id",sep=""),TaxonomicRank)
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

#Determine taxa prevalence using GBIF taxonomy.
TaxonomyExport_GBIF <- TaxonomyExport_GBIF %>%
  dplyr::group_by(species, genus, family, order, class, phylum, kingdom) %>%
  dplyr::mutate(count = dplyr::n(),
         eDNA_prevalence_GBIF = if_else(count > 1, sum(eDNA_prevalence_NCBI, na.rm = TRUE), eDNA_prevalence_NCBI)) %>%
  dplyr::ungroup() %>%
  select(-count)
