rm(list=ls())
require(dplyr)
require(plyr)
require(tidyr)
require(phyloseq)
require(decontam)
require(data.table)
require(stringr)
require(biomformat)
require(sf)
require(ggpubr)
require(ggVennDiagram)
require(vegan)

wd <- ""

setwd(wd)

#Read in Metropolitan Statistical Area (MSA) shapefiles
#https://www.census.gov/cgi-bin/geo/shapefiles/index.php?year=2020&layergroup=Urban+Areas
MSA_Boundaries <- sf::st_read(paste(wd,"/tl_2020_us_uac20/tl_2020_us_uac20.shp",sep=""))
#Reproject
MSA_Boundaries <- sf::st_transform(MSA_Boundaries,crs=4326)

#Set lists for filtering taxa.
Invertebrate_Phyla <- c("Platyhelminthes","Nemertea","Rotifera","Gastrotricha","Acanthocephala","Nematoda","Nematomorpha","Priapulida","Kinorhyncha","Loricifera","Entoprocta","Cycliophora","Gnathostomulida","Micrognathozoa","Chaetognatha","Hemichordata","Bryozoa","Brachiopoda","Phoronida","Annelida","Mollusca","Arthropoda")
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
#
project_stats <- rbind.fill(project_stats)

#Export total taxonomy table, using Tronko-assign taxonomy
if(Primer=="CO1_Metazoa"){taxa_retain <- taxa_names(physeq_total)[tax_table(physeq_total)[, "phylum"] %in% Invertebrate_Phyla]}
if(Primer=="ITS1_Fungi"){taxa_retain <- taxa_names(physeq_total)[tax_table(physeq_total)[, "phylum"] %in% Fungal_Phyla]}
physeq_retain <- prune_taxa(taxa_retain, physeq_total)
#Remove empty samples for beta diversity calculations
physeq_retain <- prune_samples(sample_sums(physeq_retain) > 0, physeq_retain)
#Remove samples without location and date data.
physeq_retain <- subset_samples(physeq_retain, !is.na(Latitude))
physeq_retain <- subset_samples(physeq_retain, !is.na(Longitude))
physeq_retain <- subset_samples(physeq_retain, !is.na(Sample.Date))

#Get sample coordinates
Sample_Points <- data.frame(sample_data(physeq_retain)[,c("Sample.ID","Latitude","Longitude","Substrate")])
Spatial_Points <- sf::st_as_sf(Sample_Points, coords = c("Longitude","Latitude"),crs=4326)

#Check which samples are found within a MSA.
MSA_Points <- sf::st_join(Spatial_Points,MSA_Boundaries,suffix=c("_Samples","_Urban"))
MSA_Points_df <- as.data.frame(MSA_Points[,c("Sample.ID","NAME20","Substrate")])
MSA_Points_df$geometry <- NULL
MSA_Points_df$Area <- ifelse(is.na(MSA_Points_df$NAME20),"Non-Urban","Urban")
Urban_Samples <- MSA_Points_df[!is.na(MSA_Points_df$NAME20),"Sample.ID"]
NonUrban_Samples <- MSA_Points_df[is.na(MSA_Points_df$NAME20),"Sample.ID"]

#Add in urban / non-urban categorization to phyloseq data
sample_data(physeq_retain)$Area <- MSA_Points_df[,"Area"]
#Extract retained sample metadata.
metadata_total <- data.frame(sample_data(physeq_retain))

#Test alpha diversity differences between areas and substrates
plot_richness(physeq_retain,x="Area",measures=c("Chao1"))
richness_data <- plot_richness(physeq_retain,x="Area",measures=c("Chao1"))$data
wilcox.test(value ~ Area,alternative = "two.sided",data=richness_data)
mean(richness_data[richness_data$Area %in% c("Urban"),"value"])
sd(richness_data[richness_data$Area %in% c("Urban"),"value"])/sqrt(nrow(richness_data[richness_data$Area %in% c("Urban"),]))
mean(richness_data[richness_data$Area %in% c("Non-Urban"),"value"])
sd(richness_data[richness_data$Area %in% c("Non-Urban"),"value"])/sqrt(nrow(richness_data[richness_data$Area %in% c("Non-Urban"),]))
anova(lm(value ~ Area+Substrate,data=richness_data))

#Get unique sampling substrates
unique_substrates <- na.omit(unique(sample_data(physeq_retain)$Substrate))

#Check if alpha diversity distributions are significantly different between urban and non-urban areas, split by sampling substrate.
Alpha_Tests <- c()
i=1
for(unique_substrate in unique_substrates){
  physeq_retain_subset <- subset_samples(physeq_retain,Substrate==unique_substrate)
  richness_data <- plot_richness(physeq_retain_subset,x="Area",measures=c("Chao1"))$data
  
  Alpha_Test <- data.frame(matrix(nrow=1,ncol=6))
  colnames(Alpha_Test) <- c("Substrate","Mean urban diversity","SE on urban diversity","Mean non-urban diversity","SE on non-urban diversity","p")
  Alpha_Test$Substrate <- unique_substrate
  Alpha_Test$`Mean urban diversity` <- mean(richness_data[richness_data$Area %in% c("Urban"),"value"])
  Alpha_Test$`SE on urban diversity` <- sd(richness_data[richness_data$Area %in% c("Urban"),"value"])/sqrt(nrow(richness_data[richness_data$Area %in% c("Urban"),]))
  Alpha_Test$`Mean non-urban diversity`<- mean(richness_data[richness_data$Area %in% c("Non-Urban"),"value"])
  Alpha_Test$`SE on non-urban diversity` <- sd(richness_data[richness_data$Area %in% c("Non-Urban"),"value"])/sqrt(nrow(richness_data[richness_data$Area %in% c("Non-Urban"),]))
  if(is.na(Alpha_Test$`Mean urban diversity`) | is.na(Alpha_Test$`Mean non-urban diversity`)){
    Alpha_Test$p <- NA
  } else{
    Alpha_Test$p <- wilcox.test(x=richness_data[richness_data$Area %in% c("Urban"),"value"],y=richness_data[richness_data$Area %in% c("Non-Urban"),"value"],alternative="two.sided")$p.value
  }
  Alpha_Tests[[i]] <- Alpha_Test
  i=i+1
}
Alpha_Tests <- rbindlist(Alpha_Tests)

#Check if beta diversity distributions are significantly different between urban and non-urban areas
beta_dist <- phyloseq::distance(physeq_retain, "chao")
tmp <- as.data.frame(TukeyHSD(betadisper(beta_dist,group=metadata_total$Area))$group)

Beta_Tests <- c()
i=1
#Get unique sampling substrates
unique_substrates <- na.omit(unique(sample_data(physeq_retain)$Substrate))
#Check if beta diversity distributions are significantly different between urban and non-urban areas, split by sampling substrate.
for(unique_substrate in unique_substrates){
  physeq_retain_subset <- subset_samples(physeq_retain,Substrate==unique_substrate)
  
  Beta_Test <- as.data.frame(matrix(nrow=1,ncol=3))
  colnames(Beta_Test) <- c("Substrate","NonUrban-NonUrban - Urban-Urban distance","p adj")
  Beta_Test$Substrate <- unique_substrate
  metadata_subset <- as.data.frame(sample_data(physeq_retain_subset))
  beta_dist <- phyloseq::distance(physeq_retain_subset, "chao")
  if(length(unique(metadata_subset$Area))>1){
    tmp <- as.data.frame(TukeyHSD(betadisper(beta_dist,group=metadata_subset$Area))$group)
    Beta_Test$`NonUrban-NonUrban - Urban-Urban distance` <- tmp$diff
    Beta_Test$`p adj` <- tmp$`p adj`
  }

  Beta_Test$Number_Urban_Samples <- sum(sample_data(physeq_retain_subset)$Area=="Urban")
  Beta_Test$Number_NonUrban_Samples <- sum(sample_data(physeq_retain_subset)$Area=="Non-Urban")
  
  Beta_Tests[[i]] <- Beta_Test
  i=i+1
  print(unique_substrate)
}
Beta_Tests <- rbindlist(Beta_Tests)

#Generate file containing information on each unique taxon, which sampling substrate it’s found on,
#and whether it’s found in urban areas, non-urban areas, or both.
i=1
Taxa_All_df <- c()
for(unique_substrate in unique_substrates){
  physeq_retain_subset <- subset_samples(physeq_retain,Substrate==unique_substrate)
  
  #Get sample coordinates
  Sample_Points <- data.frame(sample_data(physeq_retain_subset)[,c("Sample.ID","Latitude","Longitude","Substrate")])
  Spatial_Points <- sf::st_as_sf(Sample_Points, coords = c("Longitude","Latitude"),crs=4326)
  
  #Check which samples are found within a MSA.
  MSA_Points <- sf::st_join(Spatial_Points,MSA_Boundaries,suffix=c("_Samples","_Urban"))
  MSA_Points_df <- as.data.frame(MSA_Points[,c("Sample.ID","NAME20","Substrate")])
  MSA_Points_df$geometry <- NULL
  MSA_Points_df$Area <- ifelse(is.na(MSA_Points_df$NAME20),"Non-Urban","Urban")
  Urban_Samples <- MSA_Points_df[!is.na(MSA_Points_df$NAME20),"Sample.ID"]
  NonUrban_Samples <- MSA_Points_df[is.na(MSA_Points_df$NAME20),"Sample.ID"]
  
  #Get the lists of urban, non-urban, and common taxa
  All_Taxa <- as.data.frame(otu_table(physeq_retain_subset))
  Urban_Taxa <- All_Taxa[,colnames(All_Taxa) %in% Urban_Samples]
  Urban_Taxa <- Urban_Taxa[rowSums(Urban_Taxa)>0,]
  Urban_Taxa <- rownames(Urban_Taxa)
  #
  NonUrban_Taxa <- All_Taxa[,colnames(All_Taxa) %in% NonUrban_Samples]
  NonUrban_Taxa <- NonUrban_Taxa[rowSums(NonUrban_Taxa)>0,]
  NonUrban_Taxa <- rownames(NonUrban_Taxa)
  #
  Common_Taxa <- intersect(Urban_Taxa,NonUrban_Taxa)
  Urban_Taxa_tmp <- Urban_Taxa[!Urban_Taxa %in% NonUrban_Taxa]
  NonUrban_Taxa_tmp <- NonUrban_Taxa[!NonUrban_Taxa %in% Urban_Taxa]
  Urban_Taxa <- Urban_Taxa_tmp
  NonUrban_Taxa <- NonUrban_Taxa_tmp
  #
  Urban_Taxa <- as.data.frame(Urban_Taxa)
  colnames(Urban_Taxa) <- c("Taxon")
  Urban_Taxa <- separate(data = Urban_Taxa, col = "Taxon", sep=";",into = TaxonomicRanks)
  Urban_Taxa[Urban_Taxa == 'NA'] <- NA
  if(nrow(Urban_Taxa)>0){
    Urban_Taxa$Area <- "Urban"
    Urban_Taxa$Substrate <- unique_substrate
  } else{
    Urban_Taxa <- data.frame(matrix(nrow=1,ncol=9))
    colnames(Urban_Taxa) <- c(TaxonomicRanks,"Area","Substrate")
  }
  #
  NonUrban_Taxa <- as.data.frame(NonUrban_Taxa)
  colnames(NonUrban_Taxa) <- c("Taxon")
  NonUrban_Taxa <- separate(data = NonUrban_Taxa, col = "Taxon", sep=";",into = TaxonomicRanks)
  NonUrban_Taxa[NonUrban_Taxa == 'NA'] <- NA
  if(nrow(NonUrban_Taxa)>0){
    NonUrban_Taxa$Area <- "NonUrban"
    NonUrban_Taxa$Substrate <- unique_substrate
  } else{
    NonUrban_Taxa <- data.frame(matrix(nrow=1,ncol=9))
    colnames(NonUrban_Taxa) <- c(TaxonomicRanks,"Area","Substrate")
  }
  #
  Common_Taxa <- as.data.frame(Common_Taxa)
  colnames(Common_Taxa) <- c("Taxon")
  Common_Taxa <- separate(data = Common_Taxa, col = "Taxon", sep=";",into = TaxonomicRanks)
  Common_Taxa[Common_Taxa == 'NA'] <- NA
  if(nrow(Common_Taxa)>0){
    Common_Taxa$Area <- "NonUrban"
    Common_Taxa$Substrate <- unique_substrate
  } else{
    Common_Taxa <- data.frame(matrix(nrow=1,ncol=9))
    colnames(Common_Taxa) <- c(TaxonomicRanks,"Area","Substrate")
  }
  #
  Taxa_df <- rbind(Urban_Taxa,NonUrban_Taxa,Common_Taxa)
  Taxa_All_df[[i]] <- Taxa_df
  i=i+1
}
Taxa_All_df <- rbindlist(Taxa_All_df)
write.table(Taxa_All_df,paste(Primer,"_TaxaBySubstrateAndArea.txt",sep=""), row.names=FALSE, sep="\t",quote = FALSE)
Taxa_All_df <- as.data.frame(Taxa_All_df)

#Get unique lists of purely urban and non-urban taxa.
#Get sample coordinates
Sample_Points <- data.frame(sample_data(physeq_retain)[,c("Sample.ID","Latitude","Longitude","Substrate")])
Spatial_Points <- sf::st_as_sf(Sample_Points, coords = c("Longitude","Latitude"),crs=4326)
#Check which samples are found within a MSA.
MSA_Points <- sf::st_join(Spatial_Points,MSA_Boundaries,suffix=c("_Samples","_Urban"))
MSA_Points_df <- as.data.frame(MSA_Points[,c("Sample.ID","NAME20","Substrate")])
MSA_Points_df$geometry <- NULL
MSA_Points_df$Area <- ifelse(is.na(MSA_Points_df$NAME20),"Non-Urban","Urban")
Urban_Samples <- MSA_Points_df[!is.na(MSA_Points_df$NAME20),"Sample.ID"]
NonUrban_Samples <- MSA_Points_df[is.na(MSA_Points_df$NAME20),"Sample.ID"]
#Get the lists of urban, non-urban, and common taxa
All_Taxa <- as.data.frame(otu_table(physeq_retain))
Urban_Taxa <- All_Taxa[,colnames(All_Taxa) %in% Urban_Samples]
Urban_Taxa <- Urban_Taxa[rowSums(Urban_Taxa)>0,]
Urban_Taxa <- rownames(Urban_Taxa)
NonUrban_Taxa <- All_Taxa[,colnames(All_Taxa) %in% NonUrban_Samples]
NonUrban_Taxa <- NonUrban_Taxa[rowSums(NonUrban_Taxa)>0,]
NonUrban_Taxa <- rownames(NonUrban_Taxa)
#Get taxa common to both area categories
Common_Taxa <- intersect(Urban_Taxa,NonUrban_Taxa)
Urban_Taxa_tmp <- Urban_Taxa[!Urban_Taxa %in% NonUrban_Taxa]
NonUrban_Taxa_tmp <- NonUrban_Taxa[!NonUrban_Taxa %in% Urban_Taxa]
Urban_Taxa <- Urban_Taxa_tmp
NonUrban_Taxa <- NonUrban_Taxa_tmp
#Get unique urban taxa
Urban_Taxa <- as.data.frame(Urban_Taxa)
colnames(Urban_Taxa) <- c("Taxon")
Urban_Taxa <- separate(data = Urban_Taxa, col = "Taxon", sep=";",into = TaxonomicRanks)
Urban_Taxa[Urban_Taxa == 'NA'] <- NA
Urban_Taxa <- Urban_Taxa[!duplicated(Urban_Taxa),]
#Get unique non-urban taxa
NonUrban_Taxa <- as.data.frame(NonUrban_Taxa)
colnames(NonUrban_Taxa) <- c("Taxon")
NonUrban_Taxa <- separate(data = NonUrban_Taxa, col = "Taxon", sep=";",into = TaxonomicRanks)
NonUrban_Taxa[NonUrban_Taxa == 'NA'] <- NA
NonUrban_Taxa <- NonUrban_Taxa[!duplicated(NonUrban_Taxa),]

#Find total taxonomic richness per sample location and export it.
#Create presence/absence taxa by sample table.
All_OTU_Table <- as.data.frame(otu_table(physeq_retain))
All_OTU_Table <- All_OTU_Table %>% mutate_all(~ ifelse(. > 0, 1, 0))
#Find total taxonomic richness and merge into sample metadata for export
All_Taxa_Richness <- as.data.frame(colSums(All_OTU_Table))
colnames(All_Taxa_Richness) <- "Richness"
All_Taxa_Richness$Sample.ID <- rownames(All_Taxa_Richness)
All_Taxa_Richness <- dplyr::left_join(metadata_total,All_Taxa_Richness)
write.table(All_Taxa_Richness,paste(Primer,"_TotalRichness.csv",sep=""), row.names=FALSE, sep=",",quote = FALSE)

#Create presence/absence taxa by sample table for urban taxa.
Urban_OTU_Table <- All_OTU_Table[,colnames(All_OTU_Table) %in% Urban_Samples]
Urban_OTU_Table$Taxonomic_Path <- rownames(Urban_OTU_Table)
Urban_OTU_Table <- separate(data = Urban_OTU_Table, col = "Taxonomic_Path", sep=";",into = TaxonomicRanks)
Urban_OTU_Table[Urban_OTU_Table == "NA"] <- NA
Urban_OTU_Table <- dplyr::inner_join(Urban_Taxa,Urban_OTU_Table)
#Find urban taxonomic richness and merge into sample metadata for export
Urban_OTU_Table <- as.data.frame(Urban_OTU_Table)
Urban_Taxa_Richness <- as.data.frame(colSums(Urban_OTU_Table[,!colnames(Urban_OTU_Table) %in% TaxonomicRanks]))
colnames(Urban_Taxa_Richness) <- "Richness"
Urban_Taxa_Richness$Sample.ID <- rownames(Urban_Taxa_Richness)
Urban_Taxa_Richness <- dplyr::left_join(Urban_Taxa_Richness,metadata_total)
write.table(Urban_Taxa_Richness,paste(Primer,"_UrbanRichness.csv",sep=""), row.names=FALSE, sep=",",quote = FALSE)

#Create presence/absence taxa by sample table for nonurban taxa.
NonUrban_OTU_Table <- All_OTU_Table[,!colnames(All_OTU_Table) %in% Urban_Samples]
NonUrban_OTU_Table$Taxonomic_Path <- rownames(NonUrban_OTU_Table)
NonUrban_OTU_Table <- separate(data = NonUrban_OTU_Table, col = "Taxonomic_Path", sep=";",into = TaxonomicRanks)
NonUrban_OTU_Table[NonUrban_OTU_Table == "NA"] <- NA
NonUrban_OTU_Table <- dplyr::inner_join(NonUrban_Taxa,NonUrban_OTU_Table)
#Find urban taxonomic richness and merge into sample metadata for export
NonUrban_OTU_Table <- as.data.frame(NonUrban_OTU_Table)
NonUrban_Taxa_Richness <- as.data.frame(colSums(NonUrban_OTU_Table[,!colnames(NonUrban_OTU_Table) %in% TaxonomicRanks]))
colnames(NonUrban_Taxa_Richness) <- "Richness"
NonUrban_Taxa_Richness$Sample.ID <- rownames(NonUrban_Taxa_Richness)
NonUrban_Taxa_Richness <- dplyr::left_join(NonUrban_Taxa_Richness,metadata_total)
write.table(NonUrban_Taxa_Richness,paste(Primer,"_NonUrbanRichness.csv",sep=""), row.names=FALSE, sep=",",quote = FALSE)
