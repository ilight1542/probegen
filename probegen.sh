#! /usr/bin/env bash

## defaults, will be overwritten if redefined in options parsing
script_name=$(basename $0)
prog_version=0.1.0

# downloading genomes
skipdownload=false
taxonid=""
assemblylevel="chromosome,complete_genome"
assemblysource="refseq"
genomepathsfile=""
output="probeset"

# threshold ambiguous bases
ambiguousbasethreshold=""

# dustmasking options
runmasking=true

# probe clustering options
length="52"
stepsize="5"
mismatches=""
adapterseq=""

threads=1


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
    NAME
      
    $script_name - collect parallelized eager runs and postprocess/reanalyse malt(extract) outputs 
    SYNOPSIS
      
    $script_name [-O <OUTPUT DIR>] [-I </path/to/folder/with/eager/loops> ] [OPTIONAL ARGUMENTS]...

    DESCRIPTION

        $script_name is a bash script that collects and reanalyses outputs from parallelized nf-core/eager runs when malt/maltextract is ran.
        Additionally, James Fellows Yates' rma-tabuliser is ran. Optional arguments modify malt-extract postprocessing scripts to output summary pdfs of candidate hits.

        Requires (from rma-tabuliser): MEGAN (>= v6.21.7) to be installed on your system, and the contents of the tools/ directory (in the 
        MEGAN installation path) to be in your \$PATH. (Tip: the bioconda version of MEGAN puts these tools already in 
        your path). rma-tabuliser in your \$PATH variable (from https://github.com/jfy133/rma-tabuliser/blob/main/rma-tabuliser)
        Requires (from custom_postprocessing.AMPS.r): R (tested on version 3.3.1 +), R libraries: parallel, getopt, gridBase, gridExtra Optional: R libraries: jsonlite (if you want generate json-formatted raw heatmap data)

        All requirements filled by custom yaml for unix systems (see https://github.com/ilight1542/eager_postprocessing/tree/main/screening)
    
    OPTIONS

        Mandatory:

        -O [PATH]       
            Output directory for summary files, and malt-extract reanalysis, will create the path if not present already
        -I [PATH]       
            Input path to (parallel) eager processing, eg where to find files to reanalyze/summarize
            essentially at least one level above /results/.. directory from eager output
            parallel eager processes should be \$I/[run 1, 2, 3...]/results
        -n [file PATH]      
            path to node file for malt-extract; List (\\n separated) of nodes to be reported on (aka input species/node list used for MALTextract). 
            custom_postprocessing.AMPS.r option
        
        Optional:

        -d [damage cutoff]      
            Cutoff threshold for end of read damage for outputting plot. 
            Default: 0, no cutoff is used. Range [0,1). custom_postprocessing.AMPS.r option
        -i [read distribution cutoff]
            Distribution of reads required for top reference for output of candidate profile. Default 0, no cutoff used.
            Range [0,1). custom_postprocessing.AMPS.r option
        -e [default read edit ratio]         
            Ratio for default read edit distances. Default: 0.9, strong declining edit distance required. 
            Range [0,1). custom_postprocessing.AMPS.r option
        -a [ancient read edit ratio]        
            Ratio for ancient read edit distances. Default:0.8, fairly strong declining edit distance required. 
            Range [0,1). custom_postprocessing.AMPS.r option
        -s [sequence type]      
            Paired end or single end for calculating damage cutoff (if it allows damage on either end to satisfy the condition).
            Default SE. Options (SE, PE). custom_postprocessing.AMPS.r option
        -t [NUM_THREADS]        
            Number of threads to use, default 1
        -S      
            skip rma6-tabuliser
        -f
            firm cutoff outputs, output maltextract pdf summary only if all thresholds are exceeded (stp3 files only)
            custom_postprocessing.AMPS.r option   
        -h      
            print this help message
        -v      
            make execution verbose
        -V      
            print version

    AUTHOR
        
        Ian Light (ilight1542@gmail.com), rma-tabuliser: James A. Fellows Yates (jfy133@gmail.com)

    VERSION

        $prog_version
    \n    
    "
}

## argument parsing ## 
while getopts ":O:IbMlsmTcarvVh" options; do         # Loop: Get the next option;
                                          # use silent error checking
  case "${options}" in                    # 
    O)output=${OPTARG};;
    I)genomepathsfile=${OPTARG};;
    b)ambiguousbasethreshold=${OPTARG};;
    M)runmasking=true;;
    l)length=${OPTARG};;
    s)stepsize=${OPTARG};;
    m)mismatches=${OPTARG};;
    T)maskedthreshold=${OPTARG};;
    c)convertn=true;;
    a)adapterseq=${OPTARG};;
    r)randseed=${OPTARG};;
    v)version=true;;
    V)verbose="-v";;
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


if (!genomepathsfile) {
    printf "Predownloaded genome paths file must be given! Please supply a text file with a path to each reference genome"
    exit 1
}

# Check all file paths
python3 pipeline_checker.py -i ${genomepathsfile}
if ! python python3 pipeline_checker.py -i ${genomepathsfile} ; then
    exit 1
fi

### Removing genomes with rates of ambiguous bases that are too high
cat ${genomepathsfile} | while read fasta_path; do
    total_bases=$(zcat ${fasta_path} | egrep "^[^>]" ${i} | tr -d \\n | wc -c)
    ambiguous_bases_observed=$(zcat ${fasta_path} egrep "^[^>]" | egrep "[^ATCG]" -o | wc -l)
    max_ambiguous_bases=$( echo "${total_bases}*${ambiguousbasethreshold}" | bc )
    if (( $ambiguous_bases_observed < $max_ambiguous_bases )); then ## put genome into final set of genomes text (for reference) and move genome to subdirectory
        echo ${fasta_path} >> final_genomes.txt
    fi
done

#### formatting and running dustmasker
if (runmasking) {
    cat final_genomes.txt | while read fasta_path
    do
        ## dustmasker
        dustmasker -in ${fasta_path} -out temp.dusted.windows -window ${length} -outfmt acclist
        cat temp.dusted.windows >> dusted_genomes.fasta
    done

    python3 tiling.py -i ${path}/GC*.fna -l ${length} -s ${stepsize} --convert_n --randseed 100\
        masked_regions ${path}/dusted_genomes.fasta --masked_cutoff 10\
        --out_name tiling_out.txt
}
else {
    python3 tiling.py -i ${path}/GC*.fna -l ${length} -s ${stepsize} --convert_n --randseed 100\
    --out_name tiling_out.txt
}
fi 

# Cluster probes based on similarity thresholds
cd-hit-est -T 8 -M 32000 -i tiling_out.txt \
    -r 1 -G 0 -c 0.94 -o ${output}_cdhit.clusters \
    -n 9 -aL 0.94 -AL 3 -AS 3 -uL 0.4 -uS 0.4 -d 0

## Add adapter sequence (if provided) and format
index=1
for i in $(egrep "[ATCG]" ${output}_cdhit.clusters)
    do
    echo ">${output}_cdhit.clusters_${index}    ${i}${adapterseq}" >> ${output}.fasta
    index=$(( index + 1 ))
done

