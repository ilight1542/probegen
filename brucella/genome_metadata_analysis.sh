## number of unique species
grep "Brucella" brucella_genomes_metadata.tsv | awk -F'\t' 'NR > 1 {print $1}' | awk '{print $1, $2}' | sort | uniq > brucella_species.txt

while read line; do
    count=$(grep -c "$line" brucella_genomes_metadata.tsv)
    echo "$line" $count >> brucella_species_counts.txt
done <brucella_species.txt

egrep "Complete|Chromosome" brucella_genomes_metadata.tsv > brucella_genomes_metadata_chorm_compl.tsv

while read line; do
    count=$(grep -c "$line" brucella_genomes_metadata_chorm_compl.tsv)
    echo "$line" $count >> brucella_species_counts_chrom_compl.txt
done <brucella_species.txt
