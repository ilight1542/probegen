# PROBEGEN - Generate unique probes from reference genomes
Developed by Ian Light-Maka - 2023

### DESCRIPTION
PROBEGEN - Generate unique probes from reference genomes by a sliding window approach. Filtering of low-quality genomes, and clustering of similar probes.

### INSTALLATION and REQUIREMENTS
Downloadable from Zenodo or [github](https://github.com/ilight1542/probegen) 
#### Minimum requirements:
- Conda
- Unix operating system (tested)
  
All dependencies by custom conda environment yaml for unix systems. Please run `conda env create --name probegen --file=probegen.yaml`

Then prior to running the pipeline, activate the conda env with `conda activate probegen`

## USAGE
    probegen [-I </path/to/folder/with/eager/loops> ] [-O <OUTPUT DIR>] [OPTIONAL ARGUMENTS]...

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
            --percent-ambiguous-base-threshold. Default=1
        -D
            --disable-dustmasker - Turn off dustmasker of low-complexity regions. Default=Dustmasker masking enabled
        -l
            --length - Length of probes to generate (without adapter). Default=52
        -s
            --step-size - Step size of probes to generate. Default=0
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
        1.0.0

### Tips and tricks
#### Run quality filtering and dustmasking! 
These two steps are optional but will help reduce poor-quality probes produced. One could also consider running [Conterminator](https://github.com/steineggerlab/conterminator), but this is outside the scope of this pipeline.

#### Clustering parameters have large impacts on the final probe amounts!
Decreasing `--minimum-percent-identity` or increasing `--max-terminal-mismatches` will quite rapidly reduce the total number of probes output, since more and more distinct probes will be clustered together. You could try running different clustering thresholds by calling the execute_clustering.py script independent of the whole-pipeline.

#### For downloading genomes for input, consider using [NCBI datasets](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/) tool.

For instance, if you have a list of NCIB accessions you want to be included in the probes in a file called MY_ACCESSIONS.txt, try:
`datasets download genome accession --inputfile MY_ACCESSIONS.txt --exclude-gff3 --exclude-protein --exclude-rna --exclude-genomic-cds --filename genomes_accession_file.zip`

Or, if you have a given taxid you want to use with specific assembly source (eg refseq) and assembly level (eg complete), try: `datasets download genome taxon ${taxonid} --assembly-source {assemblysource} --assembly-level ${assemblylevel} --exclude-gff3 --exclude-protein --exclude-rna --exclude-genomic-cds --filename genomes_taxon_file.zip`