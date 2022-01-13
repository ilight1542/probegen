## number of unique species
grep "Brucella" brucella_genomes_metadata.tsv | awk -F'\t' 'NR > 1 {print $1}' | awk '{print $1, $2}' | sort | uniq > brucella_species.txt

while read line; do
    count=$(grep -c "$line" brucella_genomes_metadata.tsv)
    echo "$line" $count >> brucella_species_counts.txt
done <brucella_species.txt

egrep "Complete|Chromosome" brucella_genomes_metadata.tsv > brucella_genomes_metadata_chorm_compl.tsv

touch brucella_genomes_metadata_chorm_compl.tsv.filtered
for i in $(seq 1 $(wc -l < brucella_genomes_metadata_chorm_compl.tsv))
    do
    current=$(awk -v i="$i" -F'\t' 'NR == i {print $2}' brucella_genomes_metadata_chorm_compl.tsv)
    if (( $(grep -c "$current" brucella_genomes_metadata_chorm_compl.tsv) == 1 )) ; then
        grep "$current" brucella_genomes_metadata_chorm_compl.tsv >> brucella_genomes_metadata_chorm_compl.tsv.filtered
    elif (( $(grep -c "$current" brucella_genomes_metadata_chorm_compl.tsv.filtered) == 0 )) ; then
        grep "$current" brucella_genomes_metadata_chorm_compl.tsv | head -n 1 >> brucella_genomes_metadata_chorm_compl.tsv.filtered
    else
    fi
done


while read line; do
    count=$(grep -c "$line" brucella_genomes_metadata_chorm_compl.tsv)
    echo "$line" $count >> brucella_species_counts_chrom_compl.txt
done <brucella_species.txt
