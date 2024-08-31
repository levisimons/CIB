rm(list=ls())
require(plyr)
require(tidyr)
require(data.table)
require(sf)

wd <- ""

setwd(wd)

#Read in California boundaries for spatial filtering.
California <- st_read(paste(wd,"ca-state-boundary/CA_State_TIGER2016.shp",sep="/"))
California <- st_transform(California,crs=4326)

#Source BOLD data queried from https://www.boldsystems.org/index.php/Public_SearchTerms
#Read in BOLD fungal occurrences, filter them to be within California, and generate a unique taxonomic table.
BOLD_Fungi_Input <- fread(input="BOLD_Fungi_Download.txt",sep="\t",fill=T)
BOLD_Fungi <-  BOLD_Fungi_Input[,c("phylum_name","class_name","order_name","family_name","genus_name","species_name","lat","lon")]
BOLD_Fungi <- BOLD_Fungi[!is.na(BOLD_Fungi$lat) & !is.na(BOLD_Fungi$lon),]
BOLD_Fungi_Points <- st_as_sf(BOLD_Fungi, coords = c("lon", "lat"), crs = 4326, agr = "constant")
BOLD_Fungi_Points <- st_filter(BOLD_Fungi_Points,California)
BOLD_Fungi <- st_drop_geometry(BOLD_Fungi_Points)
BOLD_Fungi <- as.data.frame(BOLD_Fungi[,c("phylum_name","class_name","order_name","family_name","genus_name","species_name")])
BOLD_Fungi <- BOLD_Fungi %>% mutate_all(na_if,"")
BOLD_Fungi <- BOLD_Fungi[!duplicated(BOLD_Fungi),]
BOLD_Fungi <- BOLD_Fungi %>% mutate(BOLD_name = coalesce(species_name, genus_name, family_name, order_name,class_name, phylum_name))
BOLD_Fungi <- BOLD_Fungi %>% mutate(BOLD_rank = case_when(!is.na(species_name) ~ "species",!is.na(genus_name) ~ "genus",!is.na(family_name) ~ "family",!is.na(order_name) ~ "order",!is.na(class_name) ~ "class",!is.na(phylum_name) ~ "phylum",TRUE ~ NA_character_))
BOLD_Fungi <- BOLD_Fungi %>% dplyr::rename(species=species_name,genus=genus_name,family=family_name,order=order_name,class=class_name,phylum=phylum_name)
write.table(BOLD_Fungi,"BOLD_Fungi.tsv",quote=FALSE,sep="\t",row.names = FALSE)

#Read in BOLD invertebrate occurrences, filter them to be within California, and generate a unique taxonomic table.
BOLD_Invertebrate_Input <- fread(input="BOLD_Invertebrates_Download.txt",sep="\t",fill=T)
BOLD_Invertebrate <- BOLD_Invertebrate_Input[,c("phylum_name","class_name","order_name","family_name","genus_name","species_name","lat","lon")]
BOLD_Invertebrate <- BOLD_Invertebrate[!is.na(BOLD_Invertebrate$lat) & !is.na(BOLD_Invertebrate$lon),]
BOLD_Invertebrate_Points <- st_as_sf(BOLD_Invertebrate, coords = c("lon", "lat"), crs = 4326, agr = "constant")
BOLD_Invertebrate_Points <- st_filter(BOLD_Invertebrate_Points,California)
BOLD_Invertebrate <- st_drop_geometry(BOLD_Invertebrate_Points)
BOLD_Invertebrate <- as.data.frame(BOLD_Invertebrate[,c("phylum_name","class_name","order_name","family_name","genus_name","species_name")])
BOLD_Invertebrate <- BOLD_Invertebrate %>% mutate_all(na_if,"")
BOLD_Invertebrate <- BOLD_Invertebrate[!duplicated(BOLD_Invertebrate),]
BOLD_Invertebrate <- BOLD_Invertebrate %>% mutate(BOLD_name = coalesce(species_name, genus_name, family_name, order_name,class_name, phylum_name))
BOLD_Invertebrate <- BOLD_Invertebrate %>% mutate(BOLD_rank = case_when(!is.na(species_name) ~ "species",!is.na(genus_name) ~ "genus",!is.na(family_name) ~ "family",!is.na(order_name) ~ "order",!is.na(class_name) ~ "class",!is.na(phylum_name) ~ "phylum",TRUE ~ NA_character_))
BOLD_Invertebrate <- BOLD_Invertebrate %>% dplyr::rename(species=species_name,genus=genus_name,family=family_name,order=order_name,class=class_name,phylum=phylum_name)
write.table(BOLD_Invertebrate,"BOLD_Invertebrates.tsv",quote=FALSE,sep="\t",row.names = FALSE)
