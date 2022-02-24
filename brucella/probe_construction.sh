#### downloading
metadata_path=~/probe_design/brucella/brucella_genomes_metadata_chorm_compl.tsv.filtered

awk -F '\t' '{print $3}' ${metadata_path} > ~/probe_design/brucella/chrom_compl.txt
datasets download genome accession --inputfile ~/probe_design/brucella/chrom_compl.txt --exclude-gff3 --exclude-protein --exclude-rna

#### formatting and running dustmasker


## getting number of bases that are ambiguous and not doing these (non- ATCG bases)
mkdir final_set
for i in GC*/; do
    trimmed_i=$(echo $i | sed 's/\///g')
    total=$(egrep "^[^>]" ${i}/GC*_genomic.fna | wc -c)
    value=$(egrep "^[^>]" ${i}/GC*_genomic.fna | egrep "[^ATCG]" -o | wc -l)
    a=$( echo "${total}*5/1000" | bc )
    if (( $value < $a )); then ## put genome into final set of genomes text (for reference) and move genome to subdirectory
        echo ${i} >> ~/probe_design/brucella/final_genomes.txt
        mv ${i}/GC*_genomic.fna final_set
    fi
done


cd final_set
for fasta in GC*.fna
    do
    ## dustmasker
    dustmasker -in ${fasta} -out ${fasta}.dusted.windows -window 52 -outfmt acclist
done

cat GC*.dusted.windows >> brucella_genus_chrom_compl.dusted_windows
## pre formatting/gathering input files done

#### finally running the tiling
cd ~/probe_design/brucella
conda activate py
path=/ptmp/iclight/brucella_seqs/chrom_compl_no_duplicates/ncbi_dataset/data/final_set

## no masking
python3 tiling.py -i ${path}/GC*.fna -l 52 -s 5 --convert_n --randseed 100\
    --out_name brucella_genus_chrom_compl_probes.stepsize5.length52

## only masked windows
python3 tiling.py -i ${path}/GC*.fna -l 52 -s 5 --convert_n --randseed 100 \
    --out_name brucella_genus_chrom_compl_probes.stepsize5.length52.masking10.txt \
    --masked_regions ${path}/brucella_genus_chrom_compl.dusted_windows --masked_cutoff 10

## masked windows + reverse_complement #NOTE: REV COMPLEMENT NOT NEEDED, FUNCTIONALITY IN CD-HIT
python3 tiling.py \
    -i ${path}/*/GC*.fna -l 52 -s 5 --convert_n --randseed 100\
    --out_name brucella_genus_chrom_compl_probes.stepsize5.length52.masking10.reverse_complement \
    --masked_regions ${path}/brucella_genus_chrom_compl.dusted_windows \
    --masked_cutoff 10 --reverse_complement --outfmt fasta

## running cd-hit-est ## TODO: run cd hit on masked output
cd-hit-est -M 8000000 -i brucella_genus_chrom_compl_probes.stepsize5.length52.masking10.reverse_complement.fasta \
    -G 0 -c 1 -o brucella_genus_chrom_compl_probes.stepsize5.length52.masking10.reverse_complement.fasta.clusters \
    -n 11 -aL 0.98 -AL 2 -AS 2 -uL 0 -uS 0 -d 0

## testing cd-hit-est
cd-hit-est -i test_full.fasta \
-G 0 -c 1 -o test_full_out.clusters \
-n 11 -aL 0.98 -AL 2 -AS 2 -uL 0 -uS 0 -d 0


## just probes + rev complement
python3 tiling.py \
    -i ${path}/*/GC*.fna -l 52 -s 5 --convert_n --randseed 100\
    --out_name brucella_genus_chrom_compl_probes.stepsize5.length52.masking10.reverse_complement \
    --masked_regions ${path}/brucella_genus_chrom_compl.dusted_windows \
    --masked_cutoff 10 --reverse_complement --outfmt fasta

## just probes
python3 tiling.py \
    -i ${path}/*/GC*.fna -l 52 -s 5 --convert_n --randseed 100\
    --out_name brucella_genus_chrom_compl_probes.stepsize5.length52.masking10 \
    --masked_regions ${path}/brucella_genus_chrom_compl.dusted_windows \
    --masked_cutoff 100 --outfmt txt

## just probes + masking
python3 tiling.py \
    -i ${path}/*/GC*.fna -l 52 -s 5 --convert_n --randseed 100\
    --out_name brucella_genus_chrom_compl_probes.stepsize5.length52.masking10_update \
    --masked_regions ${path}/brucella_genus_chrom_compl.dusted_windows \
    --masked_cutoff 10 --outfmt fasta



## with prior_probe branch
python3 tiling.py \
    -i ${path}/GC*.fna -l 52 -s 5 --randseed 100\
    --out_name brucella_genus_chrom_compl_probes.stepsize5.length52.masking10.update \
    --masked_regions ${path}/brucella_genus_chrom_compl.dusted_windows \
    --masked_cutoff 10 --outfmt fasta

python3 tiling.py \
    -i ${path}/GC*.fna -l 52 -s 6 --randseed 100\
    --out_name brucella_genus_chrom_compl_probes.stepsize6.length52.masking10 \
    --masked_regions ${path}/brucella_genus_chrom_compl.dusted_windows \
    --masked_cutoff 10 --outfmt fasta

python3 tiling.py \
    -i ${path}/GC*.fna -l 52 -s 7 --randseed 100\
    --out_name brucella_genus_chrom_compl_probes.stepsize7.length52.masking10 \
    --masked_regions ${path}/brucella_genus_chrom_compl.dusted_windows \
    --masked_cutoff 10 --outfmt fasta

python3 tiling.py \
    -i ${path}/GC*.fna -l 52 -s 10 --randseed 100\
    --out_name brucella_genus_chrom_compl_probes.stepsize10.length52.masking10 \
    --masked_regions ${path}/brucella_genus_chrom_compl.dusted_windows \
    --masked_cutoff 10 --outfmt fasta

python3 tiling.py \
    -i ${path}/GC*.fna -l 52 -s 12 --randseed 100\
    --out_name brucella_genus_chrom_compl_probes.stepsize12.length52.masking10 \
    --masked_regions ${path}/brucella_genus_chrom_compl.dusted_windows \
    --masked_cutoff 10 --outfmt fasta

conda activate basic_genomics
clustering_output_name=brucella_genus_chrom_compl_probes.stepsize10.length52.masking10.fasta
cd_hit_output_name=${clustering_output_name}.aL_096_AL2_AS_2_uL04_uS04.clusters
## actual run
cd-hit-est -T 8 -M 32000 -i ${clustering_output_name} \
    -r 1 -G 0 -c 0.96 -o ${cd_hit_output_name} \
    -n 11 -aL 0.96 -AL 2 -AS 2 -uL 0.4 -uS 0.4 -d 0


clustering_output_name=brucella_genus_chrom_compl_probes.stepsize12.length52.masking10.fasta
cd_hit_output_name=${clustering_output_name}.aL_096_AL2_AS_2_uL04_uS04.clusters
## actual run
cd-hit-est -T 8 -M 32000 -i ${clustering_output_name} \
    -r 1 -G 0 -c 0.96 -o ${cd_hit_output_name} \
    -n 11 -aL 0.96 -AL 2 -AS 2 -uL 0.4 -uS 0.4 -d 0

clustering_output_name=brucella_genus_chrom_compl_probes.stepsize12.length52.masking10.fasta
cd_hit_output_name=${clustering_output_name}.c_94_aL_094_AL_3_AS_3_uL04_uS04.clusters
## actual run
cd-hit-est -T 8 -M 32000 -i ${clustering_output_name} \
    -r 1 -G 0 -c 0.94 -o ${cd_hit_output_name} \
    -n 9 -aL 0.94 -AL 3 -AS 3 -uL 0.4 -uS 0.4 -d 0


clustering_output_name=brucella_genus_chrom_compl_probes.stepsize10.length52.masking10.fasta
cd_hit_output_name=${clustering_output_name}.c_94_aL_094_AL_3_AS_3_uL04_uS04.clusters
## actual run
cd-hit-est -T 8 -M 32000 -i ${clustering_output_name} \
    -r 1 -G 0 -c 0.94 -o ${cd_hit_output_name} \
    -n 9 -aL 0.94 -AL 3 -AS 3 -uL 0.4 -uS 0.4 -d 0



## 28.01.2022
## trying looser cd hit algo with reverse complement as well

## allows 1 mismatch within 50 bp run
## test ##tests work OK on 52 bp probes
cd-hit-est -i test_full.fasta \
-G 0 -c 0.96 -o test_full_out_AL2_AS2_uL2_uS2.clusters_rev_comp \
-n 11 -aL 0.98 -AL 2 -AS 2 -uL 0.98 -uS 0.98 -d 0 -U 2


## should be only 1 mismatch in matching region ## BUT TEST THIS WITH 52 bp PROBES
cd-hit-est -M 8000 -i brucella_genus_chrom_compl_probes.stepsize5.length52.masking10.reverse_complement.fasta \
    -G 0 -c 0.98 -o brucella_genus_chrom_compl_probes.stepsize5.length52.masking10.reverse_complement.fasta.clusters \
    -n 11 -AL 2 -AS 2 -uL 0.98 -uS 0.98 -d 0 -U 2

## should be 2 mismatches in matching region ## BUT TEST THIS WITH 52 bp PROBES
cd-hit-est -M 8000 -i brucella_genus_chrom_compl_probes.stepsize5.length52.masking10.reverse_complement.fasta \
    -G 0 -c 0.96 -o brucella_genus_chrom_compl_probes.stepsize5.length52.masking10.reverse_complement.fasta.clusters \
    -n 11 -AL 2 -AS 2 -uL 0.4 -uS 0.4 -d 0

clustering_output_name=brucella_genus_chrom_compl_probes.stepsize5.length52.masking10.fasta
cd_hit_output_name=${clustering_output_name}.aL_096_AL2_AS_2_uL04_uS04.clusters

## actual run

cd-hit-est -T 8 -M 32000 -i ${clustering_output_name} \
    -r 1 -G 0 -c 0.96 -o ${cd_hit_output_name} \
    -n 11 -aL 0.96 -AL 2 -AS 2 -uL 0.4 -uS 0.4 -d 0

cd_hit_output_name=${clustering_output_name}.aL_098_AL2_AS_2_uL04_uS04.clusters

cd-hit-est -T 8 -M 32000 -i ${clustering_output_name} \
    -r 1 -G 0 -c 0.98 -o ${cd_hit_output_name} \
    -n 11 -aL 0.98 -AL 2 -AS 2 -uL 0.2 -uS 0.2 -d 0


## updated
clustering_output_name=brucella_genus_chrom_compl_probes.stepsize5.length52.masking10.update.fasta
cd_hit_output_name=${clustering_output_name}.aL_096_AL2_AS_2_uL04_uS04.clusters
## actual run
cd-hit-est -T 8 -M 32000 -i ${clustering_output_name} \
    -r 1 -G 0 -c 0.96 -o ${cd_hit_output_name} \
    -n 11 -aL 0.96 -AL 2 -AS 2 -uL 0.4 -uS 0.4 -d 0

clustering_output_name=brucella_genus_chrom_compl_probes.stepsize6.length52.masking10.fasta
cd_hit_output_name=${clustering_output_name}.aL_096_AL2_AS_2_uL04_uS04.clusters
## actual run
cd-hit-est -T 8 -M 32000 -i ${clustering_output_name} \
    -r 1 -G 0 -c 0.96 -o ${cd_hit_output_name} \
    -n 11 -aL 0.96 -AL 2 -AS 2 -uL 0.4 -uS 0.4 -d 0


clustering_output_name=brucella_genus_chrom_compl_probes.stepsize7.length52.masking10.fasta
cd_hit_output_name=${clustering_output_name}.aL_096_AL2_AS_2_uL04_uS04.clusters
## actual run
cd-hit-est -T 8 -M 32000 -i ${clustering_output_name} \
    -r 1 -G 0 -c 0.96 -o ${cd_hit_output_name} \
    -n 11 -aL 0.96 -AL 2 -AS 2 -uL 0.4 -uS 0.4 -d 0


clustering_output_name=brucella_genus_chrom_compl_probes.stepsize5.length52.masking10.update.fasta
cd_hit_output_name=${clustering_output_name}.aL_094_AL3_AS2_uL01_uS01.clusters
## actual run
cd-hit-est -T 8 -M 32000 -i ${clustering_output_name} \
    -r 1 -G 0 -c 0.94 -o ${cd_hit_output_name} \
    -n 10 -aL 0.94 -AL 3 -AS 3 -uL 0.1 -uS 0.1 -d 0





#############################################################
####### final used run ##############
############################################
conda activate py
path=/ptmp/iclight/brucella_seqs/chrom_compl_no_duplicates/ncbi_dataset/data/final_set

python3 tiling.py \
    -i ${path}/GC*.fna -l 52 -s 11 --randseed 100\
    --out_name brucella_genus_chrom_compl_probes.stepsize11.length52.masking10 \
    --masked_regions ${path}/brucella_genus_chrom_compl.dusted_windows \
    --masked_cutoff 10 --outfmt fasta


python3 tiling.py \
    -i ${path}/GC*.fna -l 52 -s 11 --randseed 100\
    --out_name brucella_genus_chrom_compl_probes.stepsize11.length52 \
    --outfmt txt

conda activate basic_genomics
clustering_output_name=brucella_genus_chrom_compl_probes.stepsize11.length52.masking10.fasta
cd_hit_output_name=${clustering_output_name}.c_94_aL_094_AL_3_AS_3_uL04_uS04.clusters
## actual run
cd-hit-est -T 8 -M 32000 -i ${clustering_output_name} \
    -r 1 -G 0 -c 0.94 -o ${cd_hit_output_name} \
    -n 9 -aL 0.94 -AL 3 -AS 3 -uL 0.4 -uS 0.4 -d 0


## putting on adapter and formatting
clustering_output_name=brucella_genus_chrom_compl_probes.stepsize11.length52.masking10.fasta
index=1
for i in $(egrep "[ATCG]" ${clustering_output_name}.c_94_aL_094_AL_3_AS_3_uL04_uS04.clusters)
    do
    echo ">${clustering_output_name}.c_94_aL_094_AL_3_AS_3_uL04_uS04.clusters_probe_${index}    ${i}CACTGCGG" >> brucella_genus_chrom_compl_probes.stepsize11.length52.masking10.post_clustering_aL94_AL3_AS3.fasta
    index=$(( index + 1 ))
done