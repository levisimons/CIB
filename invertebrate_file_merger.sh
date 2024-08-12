# List of CSV files to combine
# iNaturalist_Arachnida.csv is from https://www.inaturalist.org/observations/export?quality_grade=any&identifications=most_agree&iconic_taxa%5B%5D=Arachnida&place_id=14&lrank=genus&d2=2024-06-25
# iNaturalist_Insecta1.csv is from http://www.inaturalist.org/observations/export?quality_grade=any&identifications=most_agree&iconic_taxa%5B%5D=Insecta&place_id=14&lrank=genus&d2=2020-01-01
# iNaturalist_Insecta2.csv is from http://www.inaturalist.org/observations/export?quality_grade=any&identifications=most_agree&iconic_taxa%5B%5D=Insecta&place_id=14&lrank=genus&d1=2020-01-02&d2=2022-01-01 
# iNaturalist_Insecta3.csv is from http://www.inaturalist.org/observations/export?quality_grade=any&identifications=most_agree&iconic_taxa%5B%5D=Insecta&place_id=14&lrank=genus&d1=2022-01-02&d2=2024-06-25 
# iNaturalist_Mollusca.csv is from https://www.inaturalist.org/observations/export?quality_grade=any&identifications=most_agree&iconic_taxa%5B%5D=Mollusca&place_id=14&lrank=genus&d2=2024-06-25 
file_list="iNaturalist_Arachnida.csv iNaturalist_Insecta1.csv iNaturalist_Insecta2.csv iNaturalist_Insecta3.csv iNaturalist_Mollusca.csv"

# Output file
output_file="iNaturalist_Invertebrates.csv"

# Use the first file as the header
head -n 1 iNaturalist_Arachnida.csv > $output_file

# Combine the rest of the files without headers
for file in $file_list; do
  tail -n +2 $file >> $output_file
done
