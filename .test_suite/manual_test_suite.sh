#! /usr/bin/env bash

/*
-I [file name]
-O [file name]
--percent-ambiguous-base-threshold. Default=1
--run-dustmasker - Turn on dustmasker of low-complexity regions. Default=true
--length - length of probes to generate (without adapter). Default=52
--step-size - stepsize of probes to generate. Default=0
--randseed - Seed for random pick of replacement of non-ATCG bases in probe. Default=100
--run-clustering - Run clustering of similar probes using cd-hit-est. Default=true
--minimum-percent-identity - Minimum percent identity on aligned region for cd-hit-est to cluster a probe with a representative. Default=95
--max-terminal-mismatches - Maximum number of terminal mismatches (leading or tail) for cd-hit-est to consider an alignment. Default=2
--masked-threshold - Maximum number of reads in a probe that can be masked by dustmasker before the probe is not considered. Default=10
--adapter-seq - Adaptersequence to be appended to all probes. No default.
*/

# FAILING INPUT PATH - OK
../probegen.sh -I test_in_fake_path.txt -O test_out.txt
# FAILING INPUT - OK
../probegen.sh -O test_out.txt
# FAILING OUTPUT - OK
../probegen.sh -I test_in_fake_path.txt

# FAILING PARAMS at parameter checks - OK
step_size=0 length=0 minimumpercentidentity=0 maxterminalmismatches=10 percentambiguousbasethreshold=101 maskedthreshold=10 #SHOULD FAIL ON ALL PARAMS
../probegen.sh -I test_in.txt -O test_out.txt --percent-ambiguous-base-threshold ${percentambiguousbasethreshold} --length ${length} --step-size ${step_size} --minimum-percent-identity ${minimumpercentidentity} --max-terminal-mismatches ${maxterminalmismatches} --masked-threshold ${maskedthreshold}

# RUNS, NO ADAPTER SEQ ADDING - OK
../probegen.sh -I test_in_true_path.txt -O no_adapters.fasta 

# RUNS, NO ADAPTER SEQ ADDING, LENGTH CHANGE - OK
../probegen.sh -I test_in_true_path.txt -O no_adapters.fasta -l 100

# RUNS, NO ADAPTER SEQ ADDING, STEP SIZE CHANGE - OK
../probegen.sh -I test_in_true_path.txt -O no_adapters.fasta -s 100

# RUNS, NO ADAPTER SEQ ADDING, INTERNAL MISMATCHES CHANGE - OK
../probegen.sh -I test_in_true_path.txt -O no_adapters.fasta --max-terminal-mismatches 10

# RUNS, NO ADAPTER SEQ ADDING, NO clustering - OK
../probegen.sh -I test_in_true_path.txt -O no_adapters.fasta -C

# RUNS, NO ADAPTER SEQ ADDING, NO DUSTMASKER - OK
../probegen.sh -I test_in_true_path.txt -O no_adapters.fasta -D

# RUNS, NO ADAPTER SEQ ADDING, NO GENOME QUAL FILTER - OK
../probegen.sh -I test_in_true_path.txt -O no_adapters.fasta -F
