#! /usr/bin/env bash

## defaults, will be overwritten if redefined in options parsing
script_name=$(basename $0)
prog_version=0.1.0

# threshold ambiguous bases for genome to be excluded
percentambiguousbasethreshold="1"

# dustmasking options
runmasking=true

# probe clustering options
runclustering=true
length="52"
stepsize="1"
minimumpercentidentity="95"
maxterminalmismatches="2"
adapterseq=""
randseed="100"
maskedthreshold="10"

## functions

usage() { # Function: Print a help message.
  printf "Usage: $script_name [-O <OUTPUT DIR>] [-I </path/to/folder/with/eager/loops> ] [-n </path/to/malt/flagging/list.txt>, see postprocessing.AMPS ] [-d <damage cutoff on first/last read>, see postprocessing.AMPS] [-e <default read ratio>, see postprocessing.AMPS] [-a <ancient read ratio>, see postprocessing.AMPS] [-s <Sequencing strategy (PE vs SE)>, see postprocessing.AMPS] [-t <NUM_THREADS>] -h(elp) -v(erbose) -V(ersion)" 1>&2 
  # Will make output NEXT to rma6 file supplying -r will only export major ranks
}

exit_abnormal() { # print usage and exit
    usage
    exit 1
}

help() { # print help, explanation for all parameters
    printf "
    $(basename ${script_name^^})
    $script_name - Generate unique probes from reference genomes

    SYNOPSIS
    $script_name [-I </path/to/folder/with/eager/loops> ] [-O <OUTPUT DIR>]  [OPTIONAL ARGUMENTS]...

    DESCRIPTION
        $script_name 
        All requirements filled by custom yaml for unix systems (see https://github.com/ilight1542/eager_postprocessing/tree/main/screening)
    
    OPTIONS
        Mandatory:
        -I [file name]
            --genomepathsfile - File cotaining a line separated list of reference genomes that should be considered to create the probeset. All files must be gzipped!
        -O [file name]
            --output - Output name for completed probe set

        Optional:
        -b
            --percentambiguousbasethreshold. Default=1
        -M
            --rundustmasker - Turn on dustmasker of low-complexity regions. Default=true
        -l
            --length - length of probes to generate (without adapter). Default=52
        -s
            --stepsize - stepsize of probes to generate. Default=0
        -r
            --randseed - Seed for random pick of replacement of non-ATCG bases in probe. Default=100
        -C
            --runclustering - Run clustering of similar probes using cd-hit-est. Default=true
        -i
            --minimumpercentidentity - Minimum percent identity on aligned region for cd-hit-est to cluster a probe with a representative. Default=95
        -t
            --maxterminalmismatches - Maximum number of terminal mismatches (leading or tail) for cd-hit-est to consider an alignment. Default=2
        -m
            --maskedthreshold - Maximum number of reads in a probe that can be masked by dustmasker before the probe is not considered. Default=10
        -a
            --adapterseq - Adaptersequence to be appended to all probes. No default.
        -h      
            print this help message
        -v      
            make execution verbose
        -V      
            print version

    AUTHOR
        Ian Light-Maka (ilight1542@gmail.com)

    VERSION

        ${prog_version}
    \n    
    "
}

## argument parsing ## 
while getopts ":I:ObMlsrCitmavh" options; do         # Loop: Get the next option;
                                          # use silent error checking
  case "${options}" in                    # 
    I)genomepathsfile=${OPTARG};;
    O)output=${OPTARG};;
    b)percentambiguousbasethreshold=${OPTARG};;
    M)runmasking=true;;
    l)length=${OPTARG};;
    s)stepsize=${OPTARG};;
    r)randseed=${OPTARG};;
    C)runclustering=true;;
    i)minimumpercentidentity=${OPTARG};;
    t)maxterminalmismatches=${OPTARG};;
    m)maskedthreshold=${OPTARG};;
    a)adapterseq=${OPTARG};;
    v)version=true;;
    h)help=true;;
    :)                                    # If expected argument omitted:
      echo "Error: -${OPTARG} requires an argument."
      exit_abnormal                       # Exit abnormally.
      ;;
    *)                                    # If unknown (any other) option:
      exit_abnormal                       # Exit abnormally.
      ;;
  esac
done

## validation ##

if [[ $version == true ]]; then
  printf $prog_version
  exit 0
fi

if [[ $OPTIND -eq 1 ]]; then echo "No options were passed!" && exit_abnormal; fi
shift $((OPTIND-1))

if [[ $help == true ]]; then help && exit 0; fi


# BEGIN SCRIPT 
if [[ ! ${genomepathsfile} ]] ; then
    printf "Predownloaded genome paths file must be given! Please supply a text file with a path to each reference genome"
    exit 1
fi

# Check all file paths
date=$(date)
echo "PROBEGEN - ${date}: Running pipeline checker"
if ! python3 pipeline_checker.py -i ${genomepathsfile} ; then
    exit 1
fi

### Removing genomes with rates of ambiguous bases that are too high
date=$(date)
echo "PROBEGEN - ${date}: Removing fastas with ambiguous bases above ${percentambiguousbasethreshold} percent"
cat ${genomepathsfile} | while read fasta_path; do
total_bases=$(zcat ${fasta_path} | grep -E "^[^>]" | tr -d \\n | wc -c)
ambiguous_bases_observed=$(zcat ${fasta_path} | grep -E "^[^>]" | grep -E "[^ATCG]" -o | wc -l)
max_ambiguous_bases=$( echo "${total_bases}*${percentambiguousbasethreshold}/100" | bc )
if (( $ambiguous_bases_observed < $max_ambiguous_bases )); then ## put genome into final set of genomes text (for reference) and move genome to subdirectory
    echo ${fasta_path} >> final_genomes.txt
fi
done

#### formatting and running dustmasker, and tiling the genomes into probes
date=$(date)
echo "PROBEGEN - ${date}: Generating probe set"
if ( ${runmasking} ) ; then
    date=$(date)
    echo "     - ${date}: Running dustmasker"
    cat final_genomes.txt | while read fasta_path
        do
            ## dustmasker
            zcat ${fasta_path} | dustmasker -out temp.dusted.windows -window ${length} -outfmt acclist
            cat temp.dusted.windows >> dusted_genomes.fasta
    done
    date=$(date)
    echo "     - ${date}: Running clustering"
    if ( ${runclustering} ) ; then
        python3 tiling.py -i final_genomes.txt -l ${length} -s ${stepsize} --randseed ${randseed}\
            --masked_regions dusted_genomes.fasta --masked_cutoff ${maskedthreshold}\
            --out_name tiling_out.txt
        file_to_add_apaters=${output}_cdhit.fasta
        python3 execute_clustering.py -f tiling_out.txt -o ${output}_cdhit.fasta -l ${length} -i ${minimumpercentidentity} -t ${maxterminalmismatches}
    else
        # check if reverse complement is already present in growing probe set
        python3 tiling.py -i final_genomes.txt -l ${length} -s ${stepsize} --randseed ${randseed}\
            --masked_regions dusted_genomes.fasta --masked_cutoff ${maskedthreshold}\
            --reverse_complement --out_name tiling_out.txt
        file_to_add_apaters=tiling_out.txt
    fi
else 
    date=$(date)
    echo "     - ${date}: Running clustering"
    if ( ${runclustering} ) ; then
        python3 tiling.py -i final_genomes.txt -l ${length} -s ${stepsize} --randseed ${randseed}\
            --out_name tiling_out.txt
        file_to_add_apaters=${output}_cdhit.fasta
        python3 execute_clustering.py -f tiling_out.txt -o ${file_to_add_apaters} -l ${length} -i ${minimumpercentidentity} -t ${maxterminalmismatches}
    else
            # check if reverse complement is already present in growing probe set
        python3 tiling.py -i final_genomes.txt -l ${length} -s ${stepsize} --randseed ${randseed}\
            --out_name tiling_out.txt
        file_to_add_apaters=tiling_out.txt
    fi
        
fi 

## Add adapter sequence (if provided) and format
date=$(date)
echo "     - ${date}: Adding adapters and formatting"

index=1
for i in $(grep -E "[ATCG]" ${file_to_add_apaters})
    do
    echo ">${output}_probe_${index}    ${i}${adapterseq}" >> ${output}.fasta
    index=$(( index + 1 ))
done

## Done!
date=$(date)
echo "Probe generation completed on ${date}"
