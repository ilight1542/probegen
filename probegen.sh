#! /usr/bin/env bash

## defaults, will be overwritten if redefined in options parsing
script_name=$(basename $0)
prog_version=0.1.0

# threshold ambiguous bases for genome to be excluded
runfiltering=true
percentambiguousbasethreshold="5"

# dustmasking options
rundustmasker=true

# probe clustering options
runclustering=true
length="52"
stepsize="5"
minimumpercentidentity="95"
maxterminalmismatches="2"
adapterseq=""
randseed="100"
maskedthreshold="10"

## functions

# Function to display script usage
usage() {
    echo "Usage: $0 -I <file> -O <file> [options]"
    echo "Mandatory parameters:"
    echo "  -I, --genome-paths-file <file>                   Set the genome paths file (mandatory)"
    echo "  -O, --output <file>                              Set the output file (mandatory)"
    echo "Options:"
    echo "  -F, --disable-genome-qual-filtering              Disable removing genomes with given proportion of ambigious bases (default: filtering run)"
    echo "  -b, --percent-ambiguous-base-threshold <value>   Set the percent ambiguous base threshold for genome inclusion (default: 5)"
    echo "  -D, --disable-dustmasker                            Disable dustmasker (default: dustmasker run)"
    echo "  -l, --length <value>                             Set the length for probes (default: 52)"
    echo "  -s, --step-size <value>                          Set the step size for probes (default: 5)"
    echo "  -m, --masked-threshold <value>                   Set the masked threshold for probe inclusion (default: 10)"
    echo "  -r, --randseed <value>                           Set the random seed for probe generation (default: 100)"
    echo "  -C, --disable-clustering                         Disable clustering (default: clustering run)"
    echo "  -i, --minimum-percent-identity <value>           Set the minimum percent identity for clustering (default: 95)"
    echo "  -t, --max-terminal-mismatches <value>            Set the max terminal mismatches for clustering (default: 2)"
    echo "  -a, --adapter-seq <seq>                          Set the adapter sequence for appending to probes (default: none)"
    echo "  -v, --version                                    Show version"
    echo "  -h, --help                                       Show help"
    exit 1
}
# Parse command-line options
while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        -I|--genome-paths-file)
            genomepathsfile="$2"
            shift
            shift
            ;;
        -O|--output)
            output="$2"
            shift
            shift
            ;;
        -F|--disable-genome-qual-filtering)
            runfiltering=false
            shift
            ;;
        -b|--percent-ambiguous-base-threshold)
            percentambiguousbasethreshold="$2"
            shift
            shift
            ;;
        -D|--disable-dustmasker)
            rundustmasker=false
            shift
            ;;
        -l|--length)
            length="$2"
            shift
            shift
            ;;
        -s|--step-size)
            stepsize="$2"
            shift
            shift
            ;;
        -m|--masked-threshold)
            maskedthreshold="$2"
            shift
            shift
            ;;
        -r|--randseed)
            randseed="$2"
            shift
            shift
            ;;
        -C|--disable-clustering)
            runclustering=false
            shift
            ;;
        -i|--minimum-percent-identity)
            minimumpercentidentity="$2"
            shift
            shift
            ;;
        -t|--max-terminal-mismatches)
            maxterminalmismatches="$2"
            shift
            shift
            ;;
        -a|--adapter-seq)
            adapterseq="$2"
            shift
            shift
            ;;
        -v|--version)
            version=true
            shift
            ;;
        -h|--help)
            help=true
            shift
            ;;
        *)
            echo "Unknown option: $key"
            usage
            ;;
    esac
done

help() { # print help, explanation for all parameters
    printf  "
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
            --genomepathsfile - File cotaining a line separated list of reference genomes that should be considered to create the probeset. All reference genome files must be gzipped!
        -O [file name]
            --output - Output name for completed probe set

        Optional:
        -F
            --disable-genome-qual-filtering - Turn off filtering genomes which have too many ambiguous basecalls (set by --percent-ambiguous-base-threshold) Default: filtering enabled
        -b
            --percent-ambiguous-base-threshold. Default=5
        -D
            --disable-dustmasker - Turn off dustmasker of low-complexity regions. Default=Dustmasker masking enabled
        -l
            --length - Length of probes to generate (without adapter). Default=52
        -s
            --step-size - Step size of probes to generate. Default=5
        -m
            --masked-threshold - Maximum number of reads in a probe that can be masked by dustmasker before the probe is not considered. Default=10
        -r
            --randseed - Seed for random pick of replacement of non-ATCG bases in probe. Default=100
        -C
            --disable-clustering - Turn off clustering of similar probes using cd-hit-est. Default=Clustering enabled
        -i
            --minimum-percent-identity - Minimum percent identity on aligned region for cd-hit-est to cluster a probe with a representative. Default=95
        -t
            --max-terminal-mismatches - Maximum number of terminal mismatches (leading or tail) for cd-hit-est to consider an alignment. Default=2
        -a
            --adapter-seq - Adapter sequence to be appended to all probes.
        -h     
            --help - Print this help message
        -v      
            --version - Print version

    AUTHOR
        Ian Light-Maka (ilight1542@gmail.com)

    VERSION
        ${prog_version}
    \n    
    "
}

## validation ##
if [[ $version == true ]]; then
  printf $prog_version
  exit 0
fi

if [[ $help == true ]]; then help && exit 0; fi

# BEGIN SCRIPT 
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

# Check all file paths
date=$(date)
echo "PROBEGEN - ${date}: Running pipeline checker"
if ! python3 ${SCRIPT_DIR}/pipeline_checker.py -i ${genomepathsfile} --stepsize ${stepsize} --length ${length} --maxterminalmismatches ${maxterminalmismatches} --maskedthreshold ${maskedthreshold} --minimumpercentidentity ${minimumpercentidentity} --percentambiguousbasethreshold ${percentambiguousbasethreshold}; then
    echo "ERROR - Pipeline_checker failed!" && exit 1
fi

### Removing genomes with rates of ambiguous bases that are too high
if ( ${runfiltering} ); then
    date=$(date)
    echo "PROBEGEN - ${date}: Removing fastas with ambiguous bases above ${percentambiguousbasethreshold} percent"
    if [[ -e final_genomes.txt ]] ; then
        rm final_genomes.txt 
    fi
    cat ${genomepathsfile} | while read fasta_path; do
        total_bases=$(zcat ${fasta_path} | grep -E "^[^>]" | tr -d \\n | wc -c)
        ambiguous_bases_observed=$(zcat ${fasta_path} | grep -E "^[^>]" | grep -E "[^ATCG]" -o | wc -l)
        max_ambiguous_bases=$( echo "${total_bases}*${percentambiguousbasethreshold}/100" | bc )
        if (( $ambiguous_bases_observed < $max_ambiguous_bases )); then ## put genome into final set of genomes text (for reference) and move genome to subdirectory
            echo ${fasta_path} >> final_genomes.txt
        fi
    done
    probegen_input=final_genomes.txt
else 
    probegen_input=${genomepathsfile}
fi

#### formatting and running dustmasker, and tiling the genomes into probes
tiling_output=tiling_out.fasta
if ( ${rundustmasker} ) ; then
    date=$(date)
    echo "PROBEGEN - ${date}: Running dustmasker"
    cat ${probegen_input} | while read fasta_path
        do
            ## dustmasker
            zcat ${fasta_path} | dustmasker -out temp.dusted.windows -window ${length} -outfmt acclist
            if [ $? -eq 1 ]; then
                echo "ERROR - Dustmasker failed!" && exit 1
            fi
            cat temp.dusted.windows >> dusted_genomes.fasta
    done
    if ( ${runclustering} ) ; then
        date=$(date)
        echo "PROBEGEN - ${date}: Generating probe set"
        python3 ${SCRIPT_DIR}/tiling.py -i ${probegen_input} -l ${length} -s ${stepsize} --randseed ${randseed}\
            --masked_regions dusted_genomes.fasta --masked_cutoff ${maskedthreshold}\
            --out_name ${tiling_output} --outfmt fasta
        file_to_add_apaters=$(echo ${tiling_output} | sed 's/\.fasta/_cd_hit\.fasta/g')
        date=$(date)
        echo "     - ${date}: Running clustering"
        python3 ${SCRIPT_DIR}/execute_clustering.py -f ${tiling_output} -o ${file_to_add_apaters} -l ${length} -i ${minimumpercentidentity} -t ${maxterminalmismatches}
    else
        # check if reverse complement is already present in growing probe set
        date=$(date)
        echo "PROBEGEN - ${date}: Generating probe set"
        python3 ${SCRIPT_DIR}/tiling.py -i ${probegen_input} -l ${length} -s ${stepsize} --randseed ${randseed}\
            --masked_regions dusted_genomes.fasta --masked_cutoff ${maskedthreshold}\
            --reverse_complement --out_name ${tiling_output} --outfmt fasta
        file_to_add_apaters=${tiling_output}
    fi
else 
    if ( ${runclustering} ) ; then
        date=$(date)
        echo "PROBEGEN - ${date}: Generating probe set"
        python3 ${SCRIPT_DIR}/tiling.py -i ${probegen_input} -l ${length} -s ${stepsize} --randseed ${randseed}\
            --out_name ${tiling_output} --outfmt fasta
        file_to_add_apaters=$(echo ${tiling_output} | sed 's/\.fasta/_cd_hit\.fasta/g')
        date=$(date)
        echo "     - ${date}: Running clustering"
        python3 ${SCRIPT_DIR}/execute_clustering.py -f ${tiling_output} -o ${file_to_add_apaters} -l ${length} -i ${minimumpercentidentity} -t ${maxterminalmismatches}
    else
        # check if reverse complement is already present in growing probe set
        date=$(date)
        echo "PROBEGEN - ${date}: Generating probe set"
        python3 ${SCRIPT_DIR}/tiling.py -i ${probegen_input} -l ${length} -s ${stepsize} --randseed ${randseed}\
            --out_name ${tiling_output} --outfmt fasta
        file_to_add_apaters=${tiling_output}
    fi
        
fi 

## Add adapter sequence (if provided) and format
if [[ ! -z "${adapterseq}" ]]; then
    date=$(date)
    echo "     - ${date}: Adding adapters to probes"

    while IFS= read -r line; do
        if [[ ! $line =~ ^\> ]]; then
            line+="$adapterseq"
        fi
        echo "$line"
    done < "$file_to_add_apaters" > "${output}"
else
    cp ${tiling_output} ${output}
fi

## Done!
date=$(date)
echo "Probe generation completed on ${date}"
