## dataset comes from 
https://www.ncbi.nlm.nih.gov/datasets/genomes/?taxon=234&utm_source=genome&utm_medium=referral
## NCBI genomes database

## number of unique species
grep "Brucella" brucella_genomes_metadata.tsv | awk -F'\t' 'NR > 1 {print $1}' | awk '{print $1, $2}' | sort | uniq > brucella_species.txt

while read line; do
    count=$(grep -c "$line" brucella_genomes_metadata.tsv)
    echo "$line" $count >> brucella_species_counts.txt
done <brucella_species.txt

egrep "Complete|Chromosome" brucella_genomes_metadata.tsv > brucella_genomes_metadata_chorm_compl.tsv

## removing any doubled records (eg from both refseq and genbank)
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
    echo "$line", $count >> brucella_species_counts.txt
done <brucella_species.txt


while read line; do
    count=$(grep -c "$line" brucella_genomes_metadata_chorm_compl.tsv.filtered)
    echo "$line", $count  >> brucella_species_counts_filtered.txt
done <brucella_species.txt