rm(list=ls())
require(dplyr)
require(plyr)
require(tidyr)
require(data.table)
require(stringr)
require(pbmcapply)
require(taxonbridge)#Make sure taxonkit is installed: conda install -c bioconda taxonkit

wd <- ""

setwd(wd)

#Download GBIF taxonomic backbone
download.file(url="https://hosted-datasets.gbif.org/datasets/backbone/current/backbone.zip",destfile=paste(taxonomy_home,"backbone.zip",sep="/"))
system("unzip -o backbone.zip")
system("rm backbone.zip")
#Download NCBI taxonomic backbone
download.file(url="https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdmp.zip",destfile=paste(taxonomy_home,"taxdmp.zip",sep="/"))
system("unzip -o taxdmp.zip")
system("rm taxdmp.zip")
#Load NCBI taxonomic backbone into All.lineages.tsv.gz
system(paste("taxonkit list --data-dir=",taxonomy_home," --ids 1 | taxonkit lineage --show-lineage-taxids --show-lineage-ranks --show-rank --show-name --data-dir=",taxonomy_home," | taxonkit reformat --taxid-field 1 --data-dir=",taxonomy_home," -o All.lineages.tsv.gz",sep=""))
#Load combined NCBI and GBIF taxonomies.
custom_taxonomy <- load_taxonomies(paste(taxonomy_home,"Taxon.tsv",sep="/"), paste(taxonomy_home,"All.lineages.tsv.gz",sep="/"))
#Export initial GBIF-NCBI Taxonbridge
write.table(custom_taxonomy, paste(taxonomy_home,"ncbi_gbif_backbone_full.tsv",sep="/"), sep = "\t", col.names = TRUE, row.names = FALSE)

#NCBI taxonomic ranks.
TaxonomicRanks <- c("superkingdom","phylum","class","order","family","genus","species")
#GBIF taxonomic ranks.
GBIF_ranks <- c("species", "genus", "family", "order", "class", "phylum", "kingdom")

#Read in NCBI-GBIF Taxonbridge 
Taxon_Bridge <- fread(input="ncbi_gbif_backbone_full.tsv",sep="\t")
Taxon_Bridge <- Taxon_Bridge[Taxon_Bridge$taxonomicStatus!="doubtful"]
#Preferentially using GBIF entries with an accepted taxonomy for the GBIF-NCBI backbone
Taxon_Bridge_Accepted <- Taxon_Bridge %>%
  group_by(canonicalName, taxonRank) %>%
  slice(if ("accepted" %in% taxonomicStatus) which.max(taxonomicStatus == "accepted") else 1) %>%
  ungroup()
Taxon_Bridge_Accepted <- Taxon_Bridge_Accepted[Taxon_Bridge_Accepted$taxonRank %in% GBIF_ranks,]
write.table(Taxon_Bridge_Accepted,"Taxon_Bridge_Accepted.tsv",quote=FALSE,sep="\t",row.names = FALSE)
Taxon_Bridge_Accepted <- fread(input="Taxon_Bridge_Accepted.tsv",sep="\t")

#Read in Open Tree of Life (OTL) backbone https://tree.opentreeoflife.org/about/taxonomy-version/ott3.5
OTL_taxonomy <- fread(input="otl_taxonomy.tsv",sep="\t")
OTL_taxonomy <- OTL_taxonomy[,c("uid","parent_uid","name","rank","sourceinfo")]

#Get NCBI and GBIF taxonomic ranks and names from OTL
OTL_db <- OTL_taxonomy
OTL_db <- OTL_db %>% mutate(rank = ifelse(rank == "domain", "superkingdom", rank))
#Retain entries with at least a NCBI ID
OTL_db <- OTL_db %>% filter(str_detect(sourceinfo, "ncbi:\\d+"))
#Extract the first returns for GBIF and NCBI IDs, along with how many of them are associated with a particular OTL ID.
OTL_standardized <- OTL_db[OTL_db$rank %in% TaxonomicRanks,]
OTL_standardized <- OTL_standardized %>% mutate(ncbi_count = str_count(sourceinfo, "ncbi:\\d+"))
OTL_standardized <- OTL_standardized %>% mutate(ncbi_id = str_extract(sourceinfo, "ncbi:\\d+") %>% str_extract("\\d+"))
OTL_standardized <- OTL_standardized %>% mutate(gbif_count = str_count(sourceinfo, "gbif:\\d+"))
OTL_standardized <- OTL_standardized %>% mutate(gbif_id = str_extract(sourceinfo, "gbif:\\d+") %>% str_extract("\\d+"))
OTL_standardized$gbif_id <- as.integer(OTL_standardized$gbif_id)
OTL_standardized$ncbi_id <- as.integer(OTL_standardized$ncbi_id)
#Add in GBIF IDs and their taxonomic status from Taxon Bridge to check for missing entries downstream.
OTL_standardized <- dplyr::left_join(OTL_standardized,Taxon_Bridge[,c("taxonID","taxonRank","taxonomicStatus")],by=c("gbif_id"="taxonID"))
OTL_standardized <- OTL_standardized[!duplicated(OTL_standardized),]
OTL_standardized$gbif_id <- as.integer(OTL_standardized$gbif_id)
OTL_standardized$ncbi_id <- as.integer(OTL_standardized$ncbi_id)

write.table(OTL_standardized,"OTL_standardized.tsv",quote=FALSE,sep="\t",row.names = FALSE)
OTL_standardized <- fread(input="OTL_standardized.tsv",sep="\t")

#Function to find the most taxonomically resolved unique link between NCBI and GBIF
#taxonomic IDs within the Open Tree of Life database
select_rows <- function(df) {
  result <- c()
  for (i in 1:nrow(df)) {
    current_row <- df[i, ]
    original_name <- current_row$name
    original_rank <- current_row$rank
    original_ncbi_id <- current_row$ncbi_id
    final_rank <- current_row$rank
    # Check if ncbi_count and gbif_count both are equal to 1
    if (current_row$ncbi_count == 1 & current_row$gbif_count == 1) {
      current_row$name <- original_name
      current_row$rank <- original_rank
      current_row$ncbi_id <- original_ncbi_id
      result[[i]] <- current_row
    } else {
      # Keep checking until ncbi_count and gbif_count both are equal to 1
      while (!(current_row$ncbi_count == 1 & current_row$gbif_count == 1)) {
        current_row <- df[df$uid == current_row$parent_uid, ]
        if (nrow(current_row) == 0) break
        current_row$name <- original_name
        current_row$rank <- original_rank
        current_row$ncbi_id <- original_ncbi_id
      }
      result[[i]] <- current_row
      print(paste(i,"of",nrow(df)))
    }
  }
  return(result)
}

#Find the most taxonomically resolved unique link between NCBI and GBIF
#taxonomic IDs within the Open Tree of Life database
OTL_GBIF_updated <- select_rows(OTL_standardized)
OTL_GBIF_updated <- rbind.fill(OTL_GBIF_updated)
#Clear doubtful and depreciated GBIF taxonomic entries.
OTL_GBIF_updated <- OTL_GBIF_updated[!is.na(OTL_GBIF_updated$taxonRank),]
#Rename columns
names(OTL_GBIF_updated)[names(OTL_GBIF_updated) == "name"] <- "ncbi_name"
names(OTL_GBIF_updated)[names(OTL_GBIF_updated) == "rank"] <- "ncbi_rank"
names(OTL_GBIF_updated)[names(OTL_GBIF_updated) == "taxonRank"] <- "gbif_rank"
names(OTL_GBIF_updated)[names(OTL_GBIF_updated) == "name"] <- "ncbi_name"

#Store results
write.table(OTL_GBIF_updated,"OTL_GBIF_updated.tsv",quote=FALSE,sep="\t",row.names = FALSE)
OTL_GBIF_updated <- fread(input="OTL_GBIF_updated.tsv",sep="\t")

#NCBI node and name descriptions: https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump_readme.txt
#Download https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdmp.zip (2024-06-21 version)
#Read in ncbi nodes
ncbi_nodes <- fread("nodes.dmp",sep="\t")
ncbi_nodes <- ncbi_nodes[,c("V1","V3","V5")]
ncbi_nodes <- ncbi_nodes[!duplicated(ncbi_nodes),]
colnames(ncbi_nodes) <- c("ncbi_id","parent_id","ncbi_rank")
#Read in ncbi names
ncbi_names <- fread("names.dmp",sep="\t")
ncbi_names <- ncbi_names[!is.na(ncbi_names$V5) & ncbi_names$V7=="scientific name",]
ncbi_names <- ncbi_names[,c("V1","V3")]
ncbi_names <- ncbi_names[!duplicated(ncbi_names),]
colnames(ncbi_names) <- c("ncbi_id","ncbi_name")
#Merge names and nodes
ncbi_taxa_named_nodes <- dplyr::left_join(ncbi_nodes,ncbi_names)

#Take NCBI database and reformat it so that all of the taxonomic ranks are arranged in
#order columns
ncbi_taxa <- c()
ncbi_taxa_rows <- nrow(ncbi_taxa_named_nodes)
i=1
ncbi_taxa[[i]] <- ncbi_taxa_named_nodes
while(ncbi_taxa_rows>1){
  i=i+1
  ncbi_taxa[[i]] <- ncbi_taxa[[i-1]][(ncbi_id) %in% (parent_id), ]
  ncbi_taxa_rows <- nrow(ncbi_taxa[[i]])
  print(paste(i,ncbi_taxa_rows))
}
i_max <- i
for(i in 1:i_max){
  ncbi_taxa[[i]] <- setnames(ncbi_taxa[[i]], old = names(ncbi_taxa[[i]]), new = c(paste("ncbi_id_",(i),sep=""),paste("ncbi_id_",(i+1),sep=""),paste("ncbi_rank_",(i),sep=""),paste("ncbi_name_",(i),sep="")))
  if(i==1){ncbi_taxa_full <- ncbi_taxa[[i]]}
  if(i>1){
    ncbi_taxa_full <- dplyr::left_join(ncbi_taxa_full,ncbi_taxa[[i]])
  }
}

#Define a function to only retain standard taxonomic ranks.
rank_cols <- colnames(ncbi_taxa_full)[grep("ncbi_rank",colnames(ncbi_taxa_full))]
is_valid_column <- function(column) {
  any(column %in% c("species", "genus", "family", "order", "class", "phylum", "kingdom", "superkingdom"))
}

#Define a function to only retain columns related to the names and IDs of
#NCBI entries with standard taxonomic ranks.
ncbi_fix <- function(row_num) {
  j <- row_num
  valid_ranks <- rank_cols[sapply(ncbi_taxa_full[j, ..rank_cols], is_valid_column)]
  if(length(valid_ranks)>0){
    valid_indices <- as.integer(gsub("[^0-9]", "", valid_ranks))
    valid_names <- paste0("ncbi_name_",valid_indices,sep="")
    valid_ids <- paste0("ncbi_id_",valid_indices,sep="")
    tmp1 <- ncbi_taxa_full[j, ..valid_names]
    colnames(tmp1) <- paste0("ncbi_",unlist(ncbi_taxa_full[j, ..valid_ranks]),sep="")
    tmp2 <- ncbi_taxa_full[j, ..valid_ids]
    colnames(tmp2) <- paste0("ncbi_",unlist(ncbi_taxa_full[j, ..valid_ranks]),"_id",sep="")
    tmp <- cbind(tmp1,tmp2)
    ncbi_taxa_reformatted <- tmp
    return(ncbi_taxa_reformatted)
  }
}

#Convert NCBI database to only have IDs and names for entries with standard taxonomic ranks.
k_min <- 1
k_delta <- 10000
k <- 1
k_max <- k*k_delta
ncbi_taxa_reformatted <- c()
while(k_max < nrow(ncbi_taxa_full)){
  #Get full taxonomic paths for each taxon in NCBI
  tmp <- pbmclapply(k_min:k_max, ncbi_fix,mc.cores=detectCores())
  tmp <- rbind.fill(tmp)
  ncbi_taxa_reformatted[[k]] <- tmp
  print(paste(k,k_min,k_max))
  k <- k+1
  k_min <- k_max+1
  k_max <- k*k_delta
}
k_max <- nrow(ncbi_taxa_full)
print(paste(k,k_min,k_max))
tmp <- pbmclapply(k_min:k_max, ncbi_fix,mc.cores=detectCores())
tmp <- rbind.fill(tmp)
ncbi_taxa_reformatted[[k]] <- tmp
#Get full taxonomic paths for each taxon in NCBI
ncbi_taxa_reformatted <- rbind.fill(ncbi_taxa_reformatted)
#Generate a NCBI ID column
ncbi_taxa_reformatted <- ncbi_taxa_reformatted %>% mutate(ncbi_id = coalesce(ncbi_species_id, ncbi_genus_id, ncbi_family_id, ncbi_order_id, ncbi_class_id, ncbi_phylum_id, ncbi_kingdom_id, ncbi_superkingdom_id))
#Export results
write.table(ncbi_taxa_reformatted,"ncbi_taxa_full.tsv",quote=FALSE,sep="\t",row.names = FALSE)
ncbi_taxa_full <- fread(input="ncbi_taxa_full.tsv",sep="\t")

#Merge Open Tree of Life database entries to the full GBIF taxonomic paths.
GBIF_NCBI <- dplyr::left_join(OTL_GBIF_updated[,c("ncbi_name","ncbi_rank","ncbi_id","gbif_id","gbif_rank")],
                              Taxon_Bridge_Accepted[,c("taxonID","canonicalName","taxonRank","specificEpithet","genericName","family","order","class","phylum","kingdom")],
                              by=c("gbif_id"="taxonID"))

#Mege in NCBI backbone.
GBIF_NCBI <- dplyr::left_join(GBIF_NCBI,ncbi_taxa_full,by=c("ncbi_id"="ncbi_id"))
#Remove extraneous columns
drop_cols <- c("canonicalName","taxonRank","specificEpithet","genericName")
GBIF_NCBI <- GBIF_NCBI %>% dplyr::select(-one_of(drop_cols))
GBIF_NCBI <- GBIF_NCBI[!duplicated(GBIF_NCBI),]
#Clear out entries with NCBI entries with no corresponding GBIF entries
GBIF_NCBI <- GBIF_NCBI[rowSums(is.na(GBIF_NCBI[, ..GBIF_ranks])) != length(GBIF_ranks)]

#Export
write.table(GBIF_NCBI,"GBIF_NCBI_export.tsv",quote=FALSE,sep="\t",row.names = FALSE)
GBIF_NCBI <- fread(input="GBIF_NCBI_export.tsv",sep="\t")

#Pull in full GBIF taxonomic data set and determine the most accepted taxon keys for each
#name and rank.
Taxon_GBIF <- fread(input="Taxon.tsv",sep="\t",select=c("taxonID","taxonRank","taxonomicStatus","parentNameUsageID","canonicalName","species","genus","family","order","class","phylum","kingdom"))
Taxon_GBIF <- Taxon_GBIF[Taxon_GBIF$taxonomicStatus!="doubtful" & Taxon_GBIF$taxonRank %in% c("species","genus","family","order","class","phylum","kingdom"),]
#Preferentially select accepted taxonomic assignments.
Taxon_GBIF <- Taxon_GBIF %>%
  group_by(canonicalName, taxonRank) %>%
  slice(if ("accepted" %in% taxonomicStatus) which.max(taxonomicStatus == "accepted") else 1) %>%
  ungroup()
Taxon_GBIF <- setDT(Taxon_GBIF)
Taxon_GBIF <- Taxon_GBIF[Taxon_GBIF$canonicalName!="",]
Taxon_GBIF <- Taxon_GBIF[, species := ifelse(taxonRank == "species", canonicalName, NA)]

#Pull species entries from GBIF, assign species keys
Taxon_GBIF_species <- Taxon_GBIF[taxonRank == "species"]
Taxon_GBIF_species[, speciesKey := taxonID]
#Pull genus entries from GBIF, assign genus keys
Taxon_GBIF_genus <- Taxon_GBIF[taxonRank == "genus"]
Taxon_GBIF_genus[, genusKey := taxonID]
#Pull family entries from GBIF, assign family keys
Taxon_GBIF_family <- Taxon_GBIF[taxonRank == "family"]
Taxon_GBIF_family[, familyKey := taxonID]
#Pull order entries from GBIF, assign order keys
Taxon_GBIF_order <- Taxon_GBIF[taxonRank == "order"]
Taxon_GBIF_order[, orderKey := taxonID]
#Pull class entries from GBIF, assign class keys
Taxon_GBIF_class <- Taxon_GBIF[taxonRank == "class"]
Taxon_GBIF_class[, classKey := taxonID]
#Pull phylum entries from GBIF, assign phylum keys
Taxon_GBIF_phylum <- Taxon_GBIF[taxonRank == "phylum"]
Taxon_GBIF_phylum[, phylumKey := taxonID]
#Pull kingdom entries from GBIF, assign kingdom keys
Taxon_GBIF_kingdom <- Taxon_GBIF[taxonRank == "kingdom"]
Taxon_GBIF_kingdom[, kingdomKey := taxonID]
#Assign higher order taxon keys
Taxon_GBIF_species <- Taxon_GBIF_species[, genusKey := ifelse(parentNameUsageID %in% na.omit(unique(Taxon_GBIF_genus$taxonID)), parentNameUsageID, NA)]
Taxon_GBIF_genus <- Taxon_GBIF_genus[, familyKey := ifelse(parentNameUsageID %in% na.omit(unique(Taxon_GBIF_family$taxonID)), parentNameUsageID, NA)]
Taxon_GBIF_family <- Taxon_GBIF_family[, orderKey := ifelse(parentNameUsageID %in% na.omit(unique(Taxon_GBIF_order$taxonID)), parentNameUsageID, NA)]
Taxon_GBIF_order <- Taxon_GBIF_order[, classKey := ifelse(parentNameUsageID %in% na.omit(unique(Taxon_GBIF_class$taxonID)), parentNameUsageID, NA)]
Taxon_GBIF_class <- Taxon_GBIF_class[, phylumKey := ifelse(parentNameUsageID %in% na.omit(unique(Taxon_GBIF_phylum$taxonID)), parentNameUsageID, NA)]
Taxon_GBIF_phylum <- Taxon_GBIF_phylum[, kingdomKey := ifelse(parentNameUsageID %in% na.omit(unique(Taxon_GBIF_kingdom$taxonID)), parentNameUsageID, NA)]
#Merge in higher order keys to GBIF entries known down to species.
tmp <- Taxon_GBIF_genus[,c("family","familyKey")]
tmp <- tmp[!duplicated(tmp),]
tmp <- tmp[complete.cases(tmp),]
Taxon_GBIF_species <- dplyr::left_join(Taxon_GBIF_species,tmp)
tmp <- Taxon_GBIF_family[,c("order","orderKey")]
tmp <- tmp[!duplicated(tmp),]
tmp <- tmp[complete.cases(tmp),]
Taxon_GBIF_species <- dplyr::left_join(Taxon_GBIF_species,tmp)
tmp <- Taxon_GBIF_order[,c("class","classKey")]
tmp <- tmp[!duplicated(tmp),]
tmp <- tmp[complete.cases(tmp),]
Taxon_GBIF_species <- dplyr::left_join(Taxon_GBIF_species,tmp)
tmp <- Taxon_GBIF_class[,c("phylum","phylumKey")]
tmp <- tmp[!duplicated(tmp),]
tmp <- tmp[complete.cases(tmp),]
Taxon_GBIF_species <- dplyr::left_join(Taxon_GBIF_species,tmp)
tmp <- Taxon_GBIF_phylum[,c("kingdom","kingdomKey")]
tmp <- tmp[!duplicated(tmp),]
tmp <- tmp[complete.cases(tmp),]
Taxon_GBIF_species <- dplyr::left_join(Taxon_GBIF_species,tmp)
#Merge in higher order keys to GBIF entries known down to genus.
tmp <- Taxon_GBIF_family[,c("order","orderKey")]
tmp <- tmp[!duplicated(tmp),]
tmp <- tmp[complete.cases(tmp),]
Taxon_GBIF_genus <- dplyr::left_join(Taxon_GBIF_genus,tmp)
tmp <- Taxon_GBIF_order[,c("class","classKey")]
tmp <- tmp[!duplicated(tmp),]
tmp <- tmp[complete.cases(tmp),]
Taxon_GBIF_genus <- dplyr::left_join(Taxon_GBIF_genus,tmp)
tmp <- Taxon_GBIF_class[,c("phylum","phylumKey")]
tmp <- tmp[!duplicated(tmp),]
tmp <- tmp[complete.cases(tmp),]
Taxon_GBIF_genus <- dplyr::left_join(Taxon_GBIF_genus,tmp)
tmp <- Taxon_GBIF_phylum[,c("kingdom","kingdomKey")]
tmp <- tmp[!duplicated(tmp),]
tmp <- tmp[complete.cases(tmp),]
Taxon_GBIF_genus <- dplyr::left_join(Taxon_GBIF_genus,tmp)
#Merge in higher order keys to GBIF entries known down to family.
tmp <- Taxon_GBIF_order[,c("class","classKey")]
tmp <- tmp[!duplicated(tmp),]
tmp <- tmp[complete.cases(tmp),]
Taxon_GBIF_family <- dplyr::left_join(Taxon_GBIF_family,tmp)
tmp <- Taxon_GBIF_class[,c("phylum","phylumKey")]
tmp <- tmp[!duplicated(tmp),]
tmp <- tmp[complete.cases(tmp),]
Taxon_GBIF_family <- dplyr::left_join(Taxon_GBIF_family,tmp)
tmp <- Taxon_GBIF_phylum[,c("kingdom","kingdomKey")]
tmp <- tmp[!duplicated(tmp),]
tmp <- tmp[complete.cases(tmp),]
Taxon_GBIF_family <- dplyr::left_join(Taxon_GBIF_family,tmp)
#Merge in higher order keys to GBIF entries known down to order.
tmp <- Taxon_GBIF_class[,c("phylum","phylumKey")]
tmp <- tmp[!duplicated(tmp),]
tmp <- tmp[complete.cases(tmp),]
Taxon_GBIF_order <- dplyr::left_join(Taxon_GBIF_order,tmp)
tmp <- Taxon_GBIF_phylum[,c("kingdom","kingdomKey")]
tmp <- tmp[!duplicated(tmp),]
tmp <- tmp[complete.cases(tmp),]
Taxon_GBIF_order <- dplyr::left_join(Taxon_GBIF_order,tmp)
#Merge in higher order keys to GBIF entries known down to class.
tmp <- Taxon_GBIF_phylum[,c("kingdom","kingdomKey")]
tmp <- tmp[!duplicated(tmp),]
tmp <- tmp[complete.cases(tmp),]
Taxon_GBIF_class <- dplyr::left_join(Taxon_GBIF_class,tmp)
#Merge GBIF tables together to export taxonomic rank by taxonomic key rank table
gbif_taxa_full <- rbind.fill(Taxon_GBIF_species,Taxon_GBIF_genus,Taxon_GBIF_family,Taxon_GBIF_order,Taxon_GBIF_class,Taxon_GBIF_phylum,Taxon_GBIF_kingdom)
gbif_taxa_full <- gbif_taxa_full[,c("species","genus","family","order","class","phylum","kingdom","speciesKey","genusKey","familyKey","orderKey","classKey","phylumKey","kingdomKey")]
gbif_taxa_full <- gbif_taxa_full[!duplicated(gbif_taxa_full),]
write.table(gbif_taxa_full,"gbif_taxa_full.tsv",quote=FALSE,sep="\t",row.names = FALSE)
gbif_taxa_full <- fread(input="gbif_taxa_full.tsv",sep="\t")
gbif_taxa_full <- setDT(gbif_taxa_full)
for (col in names(gbif_taxa_full)[sapply(gbif_taxa_full, is.character)]) {
  gbif_taxa_full[get(col) == "", (col) := NA]
}

#Clean up missing genus key entries
missing_genusKey_rows <- gbif_taxa_full[!is.na(genus) & is.na(genusKey), ]
tmp <- Taxon_GBIF[Taxon_GBIF$taxonRank=="genus",c("genus","taxonID","taxonomicStatus")]
tmp <- tmp[complete.cases(tmp),]
tmp <- tmp[!duplicated(tmp),]
tmp <- tmp %>%
  group_by(genus) %>%
  slice(if ("accepted" %in% taxonomicStatus) which.max(taxonomicStatus == "accepted") else 1) %>%
  ungroup() %>% dplyr::rename(genusKey = taxonID) %>% select(genus,genusKey)
missing_genusKey_rows <- dplyr::left_join(missing_genusKey_rows[, !names(missing_genusKey_rows) %in% c("genusKey"), with = FALSE],tmp)
missing_genusKey_rows <- missing_genusKey_rows[!is.na(missing_genusKey_rows$genus),]
#Clean up missing family key entries
missing_familyKey_rows <- gbif_taxa_full[!is.na(family) & is.na(familyKey), ]
tmp <- Taxon_GBIF[Taxon_GBIF$taxonRank=="family",c("family","taxonID","taxonomicStatus")]
tmp <- tmp[complete.cases(tmp),]
tmp <- tmp[!duplicated(tmp),]
tmp <- tmp %>%
  group_by(family) %>%
  slice(if ("accepted" %in% taxonomicStatus) which.max(taxonomicStatus == "accepted") else 1) %>%
  ungroup() %>% dplyr::rename(familyKey = taxonID) %>% select(family,familyKey)
missing_familyKey_rows <- dplyr::left_join(missing_familyKey_rows[, !names(missing_familyKey_rows) %in% c("familyKey"), with = FALSE],tmp)
missing_familyKey_rows <- missing_familyKey_rows[!is.na(missing_familyKey_rows$familyKey),]
#Clean up missing order key entries
missing_orderKey_rows <- gbif_taxa_full[!is.na(order) & is.na(orderKey), ]
tmp <- Taxon_GBIF[Taxon_GBIF$taxonRank=="order",c("order","taxonID","taxonomicStatus")]
tmp <- tmp[complete.cases(tmp),]
tmp <- tmp[!duplicated(tmp),]
tmp <- tmp %>%
  group_by(order) %>%
  slice(if ("accepted" %in% taxonomicStatus) which.max(taxonomicStatus == "accepted") else 1) %>%
  ungroup() %>% dplyr::rename(orderKey = taxonID) %>% select(order,orderKey)
missing_orderKey_rows <- dplyr::left_join(missing_orderKey_rows[, !names(missing_orderKey_rows) %in% c("orderKey"), with = FALSE],tmp)
missing_orderKey_rows <- missing_orderKey_rows[!is.na(missing_orderKey_rows$orderKey),]
#Clean up missing class key entries
missing_classKey_rows <- gbif_taxa_full[!is.na(class) & is.na(classKey), ]
tmp <- Taxon_GBIF[Taxon_GBIF$taxonRank=="class",c("class","taxonID","taxonomicStatus")]
tmp <- tmp[complete.cases(tmp),]
tmp <- tmp[!duplicated(tmp),]
tmp <- tmp %>%
  group_by(class) %>%
  slice(if ("accepted" %in% taxonomicStatus) which.max(taxonomicStatus == "accepted") else 1) %>%
  ungroup() %>% dplyr::rename(classKey = taxonID) %>% select(class,classKey)
missing_classKey_rows <- dplyr::left_join(missing_classKey_rows[, !names(missing_classKey_rows) %in% c("classKey"), with = FALSE],tmp)
missing_classKey_rows <- missing_classKey_rows[!is.na(missing_classKey_rows$classKey),]
#Clean up missing phylum key entries
missing_phylumKey_rows <- gbif_taxa_full[!is.na(phylum) & is.na(phylumKey), ]
tmp <- Taxon_GBIF[Taxon_GBIF$taxonRank=="phylum",c("phylum","taxonID","taxonomicStatus")]
tmp <- tmp[complete.cases(tmp),]
tmp <- tmp[!duplicated(tmp),]
tmp <- tmp %>%
  group_by(phylum) %>%
  slice(if ("accepted" %in% taxonomicStatus) which.max(taxonomicStatus == "accepted") else 1) %>%
  ungroup() %>% dplyr::rename(phylumKey = taxonID) %>% select(phylum,phylumKey)
missing_phylumKey_rows <- dplyr::left_join(missing_phylumKey_rows[, !names(missing_phylumKey_rows) %in% c("phylumKey"), with = FALSE],tmp)
missing_phylumKey_rows <- missing_phylumKey_rows[!is.na(missing_phylumKey_rows$phylumKey),]
#Clean up missing kingdom key entries
missing_kingdomKey_rows <- gbif_taxa_full[!is.na(kingdom) & is.na(kingdomKey), ]
tmp <- Taxon_GBIF[Taxon_GBIF$taxonRank=="kingdom",c("kingdom","taxonID","taxonomicStatus")]
tmp <- tmp[complete.cases(tmp),]
tmp <- tmp[!duplicated(tmp),]
tmp <- tmp %>%
  group_by(kingdom) %>%
  slice(if ("accepted" %in% taxonomicStatus) which.max(taxonomicStatus == "accepted") else 1) %>%
  ungroup() %>% dplyr::rename(kingdomKey = taxonID) %>% select(kingdom,kingdomKey)
missing_kingdomKey_rows <- dplyr::left_join(missing_kingdomKey_rows[, !names(missing_kingdomKey_rows) %in% c("kingdomKey"), with = FALSE],tmp)
missing_kingdomKey_rows <- missing_kingdomKey_rows[!is.na(missing_kingdomKey_rows$kingdomKey),]
#Merged rows which had missing key entries, but which were fixed.
missing_Key_rows <- rbind.fill(missing_genusKey_rows,missing_familyKey_rows,missing_orderKey_rows,missing_classKey_rows,missing_phylumKey_rows,missing_kingdomKey_rows)
#Merge fixed entries back into the GBIF taxonomic rank by taxonomic key rank table
gbif_taxa_full <- gbif_taxa_full[!missing_Key_rows, on = .(species, genus, family, order, class, phylum, kingdom)]
write.table(gbif_taxa_full,"gbif_taxa_full.tsv",quote=FALSE,sep="\t",row.names = FALSE)
gbif_taxa_full <- fread(input="gbif_taxa_full.tsv",sep="\t")
gbif_taxa_full <- setDT(gbif_taxa_full)
for (col in names(gbif_taxa_full)[sapply(gbif_taxa_full, is.character)]) {
  gbif_taxa_full[get(col) == "", (col) := NA]
}
write.table(gbif_taxa_full,"gbif_taxa_full.tsv",quote=FALSE,sep="\t",row.names = FALSE)
