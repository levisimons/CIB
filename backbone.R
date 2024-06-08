rm(list=ls())
require(dplyr)
require(data.table)

#Generated NCBI-GBIF taxonomy backbone initially here: https://github.com/eDNA-Explorer/eDNAExplorer/blob/main/Backbone_generator.R
custom_taxonomy <- fread(file="ncbi_gbif_backbone_full.tsv",sep="\t")

#Remove doubtful taxonomies
filtered_taxonomy <- custom_taxonomy[custom_taxonomy$taxonomicStatus!="doubtful",]
#Keep only the first unique entry of taxonID.
filtered_taxonomy <- filtered_taxonomy %>% distinct(taxonID, .keep_all = TRUE)
#Standardize GBIF species naming to be genus+species.
taxonomicRanks <- c("species","genus","family","order","class","phylum","kingdom")
names(filtered_taxonomy)[names(filtered_taxonomy) == 'genericName'] <- 'genus'
names(filtered_taxonomy)[names(filtered_taxonomy) == 'specificEpithet'] <- 'species'
tmp <- filtered_taxonomy[,c("taxonID","genus","species")]
tmp <- tmp[complete.cases(tmp),]
tmp <- tmp[!duplicated(tmp),]
tmp$species <- paste(tmp$genus,tmp$species)
filtered_taxonomy$species <- NULL
filtered_taxonomy <- dplyr::left_join(filtered_taxonomy,tmp[,c("taxonID","species")])
filtered_taxonomy <- as.data.frame(filtered_taxonomy)

#Read in taxonomy NCBI - GBIF synonyms from OTL https://tree.opentreeoflife.org/about/taxonomy-version/ott3.5
OTL_taxonomy <- fread(file="otl_taxonomy.tsv",sep="\t")
OTL_bridge <- as.data.frame(OTL_taxonomy$sourceinfo)
colnames(OTL_bridge) <- c("sourceinfo")
#Get rows which have both gbif and ncbi entries.
OTL_bridge <- OTL_bridge %>% 
  filter(grepl("gbif",sourceinfo)) %>%
  filter(grepl("ncbi",sourceinfo))

#Create a data frame linking GBIF and NCBI ids.
combinations <- c()
for(i in 1:nrow(OTL_bridge)){
  OTL_bridge_subset <- OTL_bridge[i,]
  OTL_bridge_subset_list <- as.list(unlist(strsplit(OTL_bridge_subset, ",")))
  gbif_list <- c()
  for(j in 1:length(OTL_bridge_subset_list)){
    if(grepl("gbif",OTL_bridge_subset_list[j])){
      gbif_list <- c(gbif_list,gsub("gbif:","",OTL_bridge_subset_list[j]))
    }
  }
  ncbi_list <- c()
  for(j in 1:length(OTL_bridge_subset_list)){
    if(grepl("ncbi",OTL_bridge_subset_list[j])){
      ncbi_list <- c(ncbi_list,gsub("ncbi:","",OTL_bridge_subset_list[j]))
    }
  }
  combinations[[i]] <- as.data.frame(expand.grid(taxonID = gbif_list, ncbi_id = ncbi_list))
  print(paste(i,nrow(OTL_bridge)))
}
OTL_GBIF_NCBI <- rbindlist(combinations, use.names=TRUE, fill=TRUE)
OTL_GBIF_NCBI <- as.data.frame(OTL_GBIF_NCBI)
OTL_GBIF_NCBI$taxonID <- as.integer(as.character(OTL_GBIF_NCBI$taxonID))
OTL_GBIF_NCBI$ncbi_id <- as.integer(as.character(OTL_GBIF_NCBI$ncbi_id))

#Find the gbif ids for where a ncbi one is missing
gbif_missing_ncbi <- na.omit(filtered_taxonomy[is.na(filtered_taxonomy$ncbi_id),"taxonID"])
#Find additional ncbi ids from the OTL backbone.
additional_ncbi <- as.numeric(as.character(na.omit(OTL_GBIF_NCBI[na.omit(OTL_GBIF_NCBI$taxonID) %in% gbif_missing_ncbi,"ncbi_id"])))
ncbi_to_add <- filtered_taxonomy[na.omit(filtered_taxonomy$ncbi_id),]
ncbi_to_add <- ncbi_to_add[ncbi_to_add$ncbi_id %in% additional_ncbi,colnames(ncbi_to_add)[sapply(colnames(ncbi_to_add), function(x) grepl("ncbi", x, ignore.case = TRUE))]]
#Add in additional gbif ids.
ncbi_to_add <- dplyr::left_join(ncbi_to_add,OTL_GBIF_NCBI,multiple="all",relationship = "many-to-many")
#Combine back with entries with corresponding gbif ids.
ncbi_to_add <- dplyr::left_join(ncbi_to_add,filtered_taxonomy[filtered_taxonomy$taxonID %in% ncbi_to_add$taxonID,colnames(filtered_taxonomy)[sapply(colnames(filtered_taxonomy), function(x) !grepl("ncbi", x, ignore.case = TRUE))]]) 
ncbi_to_add <- ncbi_to_add[,colnames(filtered_taxonomy)]
#Add additional entries into taxonomy table.
filtered_taxonomy_expanded <- rbind(filtered_taxonomy,ncbi_to_add)

#Create taxonomic keys for GBIF taxonomies.
filtered_taxonomy_withKeys <- filtered_taxonomy_expanded[filtered_taxonomy_expanded$taxonRank %in% taxonomicRanks,]
for(taxonomicRank in taxonomicRanks){
  taxon_tmp <- filtered_taxonomy_withKeys[filtered_taxonomy_withKeys$taxonRank==taxonomicRank,c("taxonID","taxonomicStatus",taxonomicRank)]
  taxon_tmp <- taxon_tmp[!duplicated(taxon_tmp),]
  taxon_tmp <- taxon_tmp[complete.cases(taxon_tmp),]
  #Preferentially retain taxa with accepted taxonomies.
  taxon_tmp <- taxon_tmp %>% group_by(!!sym(taxonomicRank)) %>% arrange(desc(taxonomicStatus == 'accepted')) %>% slice(1) %>% ungroup()
  # Group by 'taxonID' and taxon, count occurrences
  grouped <- taxon_tmp %>% group_by(!!sym(taxonomicRank),taxonID) %>% summarize(count = n())
  # Find the most common 'vernacularName' for each 'taxonID'
  most_common <- grouped %>% group_by(taxonID) %>% slice_max(order_by = count, n = 1) %>% ungroup()
  #Keep only the first value of 'vernacularName' for each 'taxonID'
  first_common <- most_common %>% group_by(!!sym(taxonomicRank)) %>% slice(1) %>% select(-count)
  first_common$taxonID <- as.numeric(first_common$taxonID)
  taxon_tmp <- first_common
  colnames(taxon_tmp)[which(names(taxon_tmp) == "taxonID")] <- paste(taxonomicRank,"Key",sep="")
  filtered_taxonomy_withKeys <- dplyr::left_join(filtered_taxonomy_withKeys,taxon_tmp,multiple="all")
}

#Retain certain columns and then remove duplicate rows.
taxonomy_export <- filtered_taxonomy_withKeys[,c("taxonID","canonicalName","taxonRank","kingdom","phylum","class","order","family","genus","species","ncbi_id","ncbi_rank","ncbi_kingdom","ncbi_phylum","ncbi_class","ncbi_order","ncbi_family","ncbi_genus","ncbi_species","speciesKey","genusKey","familyKey","orderKey","classKey","phylumKey","kingdomKey")]
taxonomy_export <- taxonomy_export[!is.na(taxonomy_export$ncbi_id),]
taxonomy_export <- taxonomy_export %>% distinct()
#Retain the most common value of the gbif id for a given ncbi id.
taxonomy_export <- taxonomy_export %>%
  group_by(ncbi_id, taxonID) %>%
  summarise(count = n()) %>%
  arrange(desc(count), .by_group = TRUE) %>%
  filter(rank(desc(count), ties.method = "first") == 1) %>%
  ungroup() %>%
  left_join(taxonomy_export, by = c("ncbi_id", "taxonID")) %>%
  select(-count)

#Export NCBI-GBIF table.
taxonomy_export <- as.data.frame(taxonomy_export)
write.table(taxonomy_export,"GBIF_NCBI.csv",quote=FALSE,sep=",",row.names = FALSE)
