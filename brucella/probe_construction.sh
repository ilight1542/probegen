#### downloading
metadata_path=~/probe_design/brucella/brucella_genomes_metadata_chorm_compl.tsv.filtered

awk -F '\t' '{print $3}' ${metadata_path} > ~/probe_design/brucella/chrom_compl.txt
datasets download genome accession --inputfile ~/probe_design/brucella/chrom_compl.txt --exclude-gff3 --exclude-protein --exclude-rna

#### formatting and running dustmasker

for fasta in */GC*.fna
    do
    ## dustmasker
    dustmasker -in ${fasta} -out ${fasta}.dusted.windows -window 60 -outfmt acclist
    ## removeing all line ends
    awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' ${fasta} > ${fasta}.formatted
done

cat */GC*.dusted.windows >> brucella_genus_chrom_compl.dusted_windows

#### finally running the tiling
path=/ptmp/iclight/brucella_seqs/chrom_compl_no_duplicates/ncbi_dataset/data

## no masking
python3 tiling.py -i ${path}/*/GC*.fna.formatted -l 60 -s 5 --convert_n --randseed 100\
    --out_name brucella_genus_chrom_compl_probes.stepsize5.length60

## only masked windows
python3 tiling.py -i ${path}/*/GC*.fna.formatted -l 60 -s 5 --convert_n --randseed 100\
    --out_name brucella_genus_chrom_compl_probes.stepsize5.length60.masking10.txt \
    --masked_regions ${path}/brucella_genus_chrom_compl.dusted_windows --masked_cutoff 10

## masked windows + reverse_complement
python3 tiling.py \
    -i ${path}/*/GC*.fna.formatted -l 60 -s 5 --convert_n --randseed 100\
    --out_name brucella_genus_chrom_compl_probes.stepsize5.length60.masking10.reverse_complement \
    --masked_regions ${path}/brucella_genus_chrom_compl.dusted_windows \
    --masked_cutoff 10 --reverse_complement --outfmt all

## running cd-hit-est
cd-hit-est -M 8000000 -i brucella_genus_chrom_compl_probes.stepsize5.length60.masking10.reverse_complement.fasta \
    -G 0 -c 1 -o brucella_genus_chrom_compl_probes.stepsize5.length60.masking10.reverse_complement.fasta.clusters \
    -n 11 -aL 0.98 -AL 2 -AS 2 -uL 0 -uS 0 -d 0

cd-hit-est -i test_full.fasta \
-G 0 -c 1 -o test_full_out.clusters \
-n 11 -aL 0.98 -AL 2 -AS 2 -uL 0 -uS 0 -d 0

