#! /usr/bin/env bash

/*
-I [file name]
-O [file name]
--percentambiguousbasethreshold. Default=1
--rundustmasker - Turn on dustmasker of low-complexity regions. Default=true
--length - length of probes to generate (without adapter). Default=52
--stepsize - stepsize of probes to generate. Default=0
--randseed - Seed for random pick of replacement of non-ATCG bases in probe. Default=100
--runclustering - Run clustering of similar probes using cd-hit-est. Default=true
--minimumpercentidentity - Minimum percent identity on aligned region for cd-hit-est to cluster a probe with a representative. Default=95
--maxterminalmismatches - Maximum number of terminal mismatches (leading or tail) for cd-hit-est to consider an alignment. Default=2
--maskedthreshold - Maximum number of reads in a probe that can be masked by dustmasker before the probe is not considered. Default=10
--adapterseq - Adaptersequence to be appended to all probes. No default.
*/

# FAILING INPUT PATH
../probegen.sh -I test_in_fake_path.txt -O test_out.txt


# FAILING PARAMS at parameter checks
step_size=0 length=0 minimumpercentidentity=0 maxterminalmismatches=10 percentambiguousbasethreshold=101 maskedthreshold=10 #SHOULD FAIL ON ALL PARAMS
../probegen.sh -I test_in.txt -O test_out.txt --percentambiguousbasethreshold ${percentambiguousbasethreshold} --length ${length} --stepsize=${step_size} --minimumpercentidentity ${minimumpercentidentity} --maxterminalmismatches ${maxterminalmismatches} --maskedthreshold ${maskedthreshold}
step_size=0 length=0 minimumpercentidentity=80 maxterminalmismatches=10 #SHOULD FAIL
../probegen.sh -I test_in.txt -O test_out.txt --percentambiguousbasethreshold ${percentambiguousbasethreshold} --length ${length} --stepsize=${step_size} --minimumpercentidentity ${minimumpercentidentity} --maxterminalmismatches ${maxterminalmismatches} --maskedthreshold ${maskedthreshold}
step_size=0 length=0 minimumpercentidentity=100 maxterminalmismatches=10 #SHOULD FAIL
../probegen.sh -I test_in.txt -O test_out.txt --percentambiguousbasethreshold ${percentambiguousbasethreshold} --length ${length} --stepsize=${step_size} --minimumpercentidentity ${minimumpercentidentity} --maxterminalmismatches ${maxterminalmismatches} --maskedthreshold ${maskedthreshold}
step_size=0 length=1 minimumpercentidentity=0 maxterminalmismatches=10 #SHOULD FAIL
../probegen.sh -I test_in.txt -O test_out.txt --percentambiguousbasethreshold ${percentambiguousbasethreshold} --length ${length} --stepsize=${step_size} --minimumpercentidentity ${minimumpercentidentity} --maxterminalmismatches ${maxterminalmismatches} --maskedthreshold ${maskedthreshold}
step_size=0 length=1 minimumpercentidentity=80 maxterminalmismatches=10 #SHOULD FAIL
../probegen.sh -I test_in.txt -O test_out.txt --percentambiguousbasethreshold ${percentambiguousbasethreshold} --length ${length} --stepsize=${step_size} --minimumpercentidentity ${minimumpercentidentity} --maxterminalmismatches ${maxterminalmismatches} --maskedthreshold ${maskedthreshold}
step_size=0 length=1 minimumpercentidentity=100 maxterminalmismatches=10 #SHOULD FAIL
../probegen.sh -I test_in.txt -O test_out.txt --percentambiguousbasethreshold ${percentambiguousbasethreshold} --length ${length} --stepsize=${step_size} --minimumpercentidentity ${minimumpercentidentity} --maxterminalmismatches ${maxterminalmismatches} --maskedthreshold ${maskedthreshold}
step_size=0 length=100 minimumpercentidentity=0 maxterminalmismatches=10 #SHOULD FAIL
../probegen.sh -I test_in.txt -O test_out.txt --percentambiguousbasethreshold ${percentambiguousbasethreshold} --length ${length} --stepsize=${step_size} --minimumpercentidentity ${minimumpercentidentity} --maxterminalmismatches ${maxterminalmismatches} --maskedthreshold ${maskedthreshold}
step_size=0 length=100 minimumpercentidentity=80 maxterminalmismatches=10 #SHOULD FAIL
../probegen.sh -I test_in.txt -O test_out.txt --percentambiguousbasethreshold ${percentambiguousbasethreshold} --length ${length} --stepsize=${step_size} --minimumpercentidentity ${minimumpercentidentity} --maxterminalmismatches ${maxterminalmismatches} --maskedthreshold ${maskedthreshold}
step_size=0 length=100 minimumpercentidentity=100 maxterminalmismatches=10 #SHOULD FAIL
../probegen.sh -I test_in.txt -O test_out.txt --percentambiguousbasethreshold ${percentambiguousbasethreshold} --length ${length} --stepsize=${step_size} --minimumpercentidentity ${minimumpercentidentity} --maxterminalmismatches ${maxterminalmismatches} --maskedthreshold ${maskedthreshold}
