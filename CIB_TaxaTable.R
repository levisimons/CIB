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

#Select a primer
Primers <- c("12S_MiFish_U","16S_Bacteria","18S_Euk","CO1_Metazoa","ITS1_Fungi","ITS2_Plants","vert12S")
Primer <- "CO1_Metazoa"

#Aggregate tronko-assign results to a particular taxonomic level.
TaxonomicRanks <- c("superkingdom","phylum","class","order","family","genus","species")
TaxonomicRank <- "family"

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
      contamdf.prev_prevalence_filtered = filter_taxa(contamdf.prev_threshold_filtered,function(x) sum(x > 1) > 1, TRUE)
      contamdf.prev_no_blanks = subset_samples(contamdf.prev_prevalence_filtered, Sample_or_Control == "True Sample")
      contamdf.prev_min_threshold_a = prune_taxa(taxa_sums(contamdf.prev_no_blanks) > 19, contamdf.prev_no_blanks)
      contamdf.prev_min_threshold_p = prune_samples(sample_sums(contamdf.prev_min_threshold_a)>=20, contamdf.prev_min_threshold_a)
      ##
      Invertebrate_Phyla <- c("Arthropoda","Mollusca","Annelida","Echinodermata","Platyhelminthes","Nematoda","Onychophora","Rotifera","Brachiopoda","Tardigrades","Dicyemida","Bryozoa")
      contamdf.prev_min_threshold_p = subset_taxa(contamdf.prev_min_threshold_p, phylum %in% Invertebrate_Phyla)
      ##
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
#Aggregate tronko-assign results into a single table.
TronkoTables <- rbind.fill(TronkoTables)
TronkoTables[is.na(TronkoTables)] <- 0
#Aggregate project metadata results into a single table.
MetadataTables <- rbind.fill(MetadataTables)
#Aggregate summary statistics
Project_Samples <- rbind.fill(Project_Samples)
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
TronkoTables <- TronkoTables[,c(TaxonomicRank,SampleIDs)]
TronkoTables <- TronkoTables[!is.na(TronkoTables[,TaxonomicRank]),]
TronkoTables <- TronkoTables[TronkoTables[,TaxonomicRank]!="NA",]
names(TronkoTables)[names(TronkoTables) == TaxonomicRank] <- "Taxon"
TronkoTables <- aggregate(.~Taxon, data = TronkoTables, FUN = sum)
rownames(TronkoTables) <- TronkoTables$Taxon
TronkoTables$Taxon <- NULL

#Merge metadata with tronko-assign data.
TronkoTables <- as.data.frame(t(TronkoTables))
TronkoTables$fastqid <- rownames(TronkoTables)
rownames(TronkoTables) <- NULL
MergedData <- dplyr::left_join(MetadataTables,TronkoTables,by=c("fastqid"))

#Get list of taxa names.
taxa_initial <- colnames(MergedData)[!(colnames(MergedData) %in% c("Sample ID","Sample Date","Sample Type","Latitude","Longitude","Fastq Forward Reads Filename","projectid","fastqid","Sample_or_Control","is.neg"))]
#Get the list of taxa prevalent enough to run species distribution models on.
MergedData[,taxa_initial][MergedData[,taxa_initial] > 0] <- 1
#taxa_retain <- taxa_initial[colSums(MergedData[taxa_initial],na.rm=T) >= 30]
taxa_retain <- taxa_initial

#Merge in extracted environmental metadata with biodiversity data.
Biodiversity_Environment <- dplyr::left_join(MergedData,Environmental_Metadata,by=c("Sample ID"="name","projectid"="projectid","Latitude"="Latitude","Longitude"="Longitude","Sample Date"="Sample_Date"),multiple = "all")
Environmental_Variables <- colnames(Environmental_Metadata)[!(colnames(Environmental_Metadata) %in% c("name","Sample_Date","Latitude","Longitude","Spatial_Uncertainty","projectid"))]

#Export biodiversity spatial data for analysis in ArcGIS Pro
Biodiversity_Environment_Subset <- Biodiversity_Environment[,c(taxa_retain,"Latitude","Longitude",Environmental_Variables)]
Biodiversity_Environment_Subset$Richness <- rowSums(Biodiversity_Environment_Subset[,taxa_retain])
require(sf)
tmp <- Biodiversity_Environment_Subset[,c("Latitude","Longitude","Richness")]
tmp <- tmp[complete.cases(tmp),]
write.table(tmp,"InvertebrateFamilies.csv",quote=FALSE,sep=",",row.names = FALSE)
# Create a simple feature (sf) object
sf <- st_as_sf(tmp, coords = c("Longitude", "Latitude"), crs = 4326)
# Save the sf object as a shapefile
st_write(sf, "InvertebrateFamilies.shp")

#Run random forest models and aggregate results.
accuracy_summary <- list()
importance_summary <- list()
set.seed(1)
k=1
for(j in 1:length(taxa_retain)){
  for(i in 1:10){
    #Construct input for random forest modeling.
    Taxon <- taxa_retain[j]
    RF_Input <- Biodiversity_Environment[,c(Taxon,Environmental_Variables)]
    names(RF_Input)[names(RF_Input) == Taxon] <- "Taxon"
    #Remove largely absent environmental variables.
    RF_Input<- RF_Input[, colMeans(is.na(RF_Input)) <= 0.5]
    RF_Input <- RF_Input[complete.cases(RF_Input),]
    #Set factor variables
    Factor_Variables <- c("Taxon","grtgroup","biome_type","Landform","ECO_NAME","REALM","HYBAS_ID","ENDO","COAST","ORDER")
    RF_Input <- RF_Input %>% mutate_at(vars(one_of(Factor_Variables)), factor)
    #Construct a training and testing sets.
    RF_Input_Present <- RF_Input[RF_Input[,"Taxon"]==1,]
    RF_Input_Absent <- RF_Input[RF_Input[,"Taxon"]==0,]
    if(nrow(RF_Input_Present)>=30 & nrow(RF_Input_Absent)>=30){
      group <- kfold(RF_Input_Present,5)
      Training_Present <- RF_Input_Present[group!=1,]
      Testing_Present <- RF_Input_Present[group==1,]
      group <- kfold(RF_Input_Absent,5)
      Training_Absent <- RF_Input_Absent[group!=1,]
      Testing_Absent <- RF_Input_Absent[group==1,]
      Training_Set <- rbind(Training_Present,Training_Absent)
      Testing_Set <- rbind(Testing_Present,Testing_Absent)
      
      #Run the random forest model using the training data.
      colnames(Training_Set) <- gsub(" ","_",colnames(Training_Set))
      colnames(Testing_Set) <- gsub(" ","_",colnames(Testing_Set))
      rf1 <- randomForest(Taxon ~ . , data=Training_Set)
      # Predict the testing set with the trained model
      predictions <- predict(rf1, Testing_Set, type = "class")
      # Accuracy and other metrics
      RF_metrics <- as.data.frame(confusionMatrix(predictions, Testing_Set$Taxon)[["byClass"]])
      colnames(RF_metrics) <- c("summary")
      test <- as.data.frame(matrix(nrow=1,ncol=2))
      colnames(test) <- c("Taxon","TSS")
      test$TSS <- RF_metrics["Sensitivity","summary"]+RF_metrics["Specificity","summary"]-1
      test$Taxon <- Taxon
      accuracy_summary[[k]] <- test[,c("Taxon","TSS")]
      tmp <- as.data.frame(importance(rf1))
      tmp$Variable <- rownames(tmp)
      tmp$Taxon <- Taxon
      importance_summary[[k]] <- tmp
      print(paste(j,Taxon,i))
      k=k+1
    }
  } 
}
#Calculate the mean, and standard deviation, on the accuracy of the RF model.
accuracy_total <- rbind.fill(accuracy_summary)
accuracy_total <- accuracy_total[complete.cases(accuracy_total),]
accuracy_total <- ddply(accuracy_total, .(Taxon), summarize,  Mean_TSS=mean(TSS), SD_TSS=sd(TSS))
#Calculate the mean, and standard deviation, on the variable importances of the RF model.
importance_total <- rbind.fill(importance_summary)
importance_total <- importance_total[complete.cases(importance_total),]
importance_total <- ddply(importance_total, .(Taxon,Variable), summarize,  Mean_MeanDecreaseGini=mean(MeanDecreaseGini), SD_MeanDecreaseGini=sd(MeanDecreaseGini))

###
require(leaflet)
require(sf)
require(viridis)
require(dplyr)
require(htmlwidgets)
require(mapview)
# Map eDNA locations.
MergedData <- MergedData[!is.na(MergedData$Longitude),]
m <- leaflet(data=MergedData) %>% addTiles() %>% addProviderTiles(providers$CartoDB.Positron) %>%
  addCircleMarkers(~Longitude, ~Latitude, radius=2, color = "green") %>%
  setView(lng = mean(MergedData$Longitude), mean(MergedData$Latitude), zoom = 5)
m
#
n <- 88
taxon <- taxa_retain[n]
tmp <- MergedData[,c("Latitude","Longitude",taxon)]
tmp[,taxon][is.na(tmp[,taxon])] <- 0
tmp[,taxon][tmp[,taxon] > 0] <- 1 
tmp[,taxon] <- as.factor(tmp[,taxon])
pal <- colorFactor(palette = c("red","blue"),levels=tmp[,taxon])
m <- leaflet(data=tmp) %>% addTiles() %>% addProviderTiles(providers$CartoDB.Positron) %>%
  addCircleMarkers(~Longitude, ~Latitude, radius=2^as.numeric(tmp[[taxon]]),opacity=0.15,color = ~pal(na.omit(tmp[[taxon]]))) %>%
  addLegend(position = 'topright',pal=pal,values=~tmp[[taxon]],title = paste(taxon,"detected<br>")) %>%
  setView(lng = mean(tmp$Longitude), mean(tmp$Latitude), zoom = 5)
m
