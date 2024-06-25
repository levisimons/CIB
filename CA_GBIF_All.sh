#!/bin/bash

#This script determines the unique taxa from the following GBIF query:
#GBIF.org (25 June 2024) GBIF Occurrence Download https://doi.org/10.15468/dl.u8mz2u
# Define the path to the zipped file and the file inside the zip
#zip_file_path="path/to/your/file.zip"
zip_file_path="0093168-240506114902167.zip"
#file_name="your_file.tsv"
file_name="0093168-240506114902167.csv"

# Use unzip to read the file without extracting it, select the relevant columns, and get unique rows
unzip -p "$zip_file_path" "$file_name" | \
  awk -F'\t' '
  BEGIN {
    OFS = "\t";
  }
  NR == 1 {
    for (i = 1; i <= NF; i++) {
      if ($i == "species") species_col = i;
      if ($i == "genus") genus_col = i;
      if ($i == "family") family_col = i;
      if ($i == "order") order_col = i;
      if ($i == "class") class_col = i;
      if ($i == "phylum") phylum_col = i;
      if ($i == "kingdom") kingdom_col = i;
    }
    print "species", "genus", "family", "order", "class", "phylum", "kingdom";
    next;
  }
  {
    key = $species_col OFS $genus_col OFS $family_col OFS $order_col OFS $class_col OFS $phylum_col OFS $kingdom_col;
    if (!seen[key]++) {
      print $species_col, $genus_col, $family_col, $order_col, $class_col, $phylum_col, $kingdom_col;
    }
  }' > CA_GBIF_All_Taxa.tsv

echo "Unique rows have been saved to CA_GBIF_Taxa_All.tsv"
