# PROBEGEN - Generate unique probes from reference genomes
Developed by Ian Light-Maka - 2023

Please cite: Ian Light, Alexander Herbig, and Felix M Key. (2023) PROBEGEN (v1.0.0). Zenodo. [10.5281/zenodo.10159353](https://doi.org/10.5281/zenodo.10159353)

### DESCRIPTION
PROBEGEN - Generate unique probes from reference genomes with a sliding window approach. Filtering of low-quality genomes, masking low-complexity regions, and clustering highly similar probes.

### Example usage:
#### Standard usage: 

`./probegen.sh -I paths_to_reference_genomes.txt -O my_probes.fasta`

#### No quality filtering of references (-F), no masking low-complexity regions (-D) and no clustering (-C) of probes based on high identity:

`./probegen.sh -I paths_to_reference_genomes.txt -O my_probes.fasta -F -D -C`

#### Set the percent identity threshold on aligned region to 85% for similar probes to be clustered:

`./probegen.sh -I paths_to_reference_genomes.txt -O my_probes.fasta -i 85`

#### Example outputs: 
A fasta file with unique probes on each line for the some *Y. pestis* reference genomes. Note: the probes are not ordered according to their appearance in the reference genome!
```
>0
GAATTATTTTGGTCAGGGGGCGTTATTACTGAAAAACCCTGAAGCCATCAAA
>1
CGAACTCAACAGGATGCGGAGTCATACGCACGCGATCCATTTCTTCTCTACG
>2
AACAAATACCATCTGAGCGATACGCTCACCGGGTTCGATGGTGAAAGGCTGC
>3
ATATTCTGAGCTGCCTAAACCAACCGCCCCAAAGCGTACTTGGGATAAATCA
>4
ATGATTGCGGTGGATGGTTGGTTACAGCAAGAGCCAGAACCGTTGGTGCGTG
>5
ACAGGTACAAAGGCAGACTCGCAGTCTAATTAACGACGACCTTCAGCAATGG
...
```

Or, with unique probes with an adapter sequence of AAAAAAA appended to the end. Note: the probes are not ordered according to their appearance in the reference genome!
```
>0
GAATTATTTTGGTCAGGGGGCGTTATTACTGAAAAACCCTGAAGCCATCAAAAAAAAAAA
>1
CGAACTCAACAGGATGCGGAGTCATACGCACGCGATCCATTTCTTCTCTACGAAAAAAAA
>2
AACAAATACCATCTGAGCGATACGCTCACCGGGTTCGATGGTGAAAGGCTGCAAAAAAAA
>3
ATATTCTGAGCTGCCTAAACCAACCGCCCCAAAGCGTACTTGGGATAAATCAAAAAAAAA
>4
ATGATTGCGGTGGATGGTTGGTTACAGCAAGAGCCAGAACCGTTGGTGCGTGAAAAAAAA
>5
ACAGGTACAAAGGCAGACTCGCAGTCTAATTAACGACGACCTTCAGCAATGGAAAAAAAA
...
```

### INSTALLATION and REQUIREMENTS
Downloadable from Zenodo or [github](https://github.com/ilight1542/probegen) 
#### Minimum requirements:
- Conda
- Unix operating system (Other OS may work but have not been checked)
  
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
            --randseed - Seed for random pick of replacement of non-ATCG bases in probes with ambiguous bases. Default=100
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
These two steps are optional but they will help reduce poor-quality probes produced. One could also consider running [Conterminator](https://github.com/steineggerlab/conterminator), but this is outside the scope of this pipeline.

#### Clustering parameters have large impacts on the final probe amounts!
Decreasing `--minimum-percent-identity` or increasing `--max-terminal-mismatches` will quite rapidly reduce the total number of probes output, since more and more distinct probes will be clustered together. You could try running different clustering thresholds by calling the execute_clustering.py script independent of the whole-pipeline.

#### For downloading genomes for input, consider using [NCBI datasets](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/) tool.

For instance, if you have a list of NCIB accessions you want to be included in the probes in a file called MY_ACCESSIONS.txt, try:
`datasets download genome accession --inputfile MY_ACCESSIONS.txt --exclude-gff3 --exclude-protein --exclude-rna --exclude-genomic-cds --filename genomes_accession_file.zip`

Or, if you have a given taxid you want to use with specific assembly source (eg refseq) and assembly level (eg complete), try: `datasets download genome taxon ${taxonid} --assembly-source {assemblysource} --assembly-level ${assemblylevel} --exclude-gff3 --exclude-protein --exclude-rna --exclude-genomic-cds --filename genomes_taxon_file.zip`

### References
Sayers EW, Bolton EE, Brister JR, et al. Database resources of the national center for biotechnology information. Nucleic Acids Res. 2022;50(D1):D20-D26. https://doi.org/10.1093/nar/gkab1112

Steinegger, M., Salzberg, S.L. Terminating contamination: large-scale search identifies more than 2,000,000 contaminated entries in GenBank. Genome Biol 21, 115 (2020). https://doi.org/10.1186/s13059-020-02023-1

Walt, Stéfan van der, S. Chris Colbert, and Gaël Varoquaux. 2011. “The NumPy Array: A Structure for Efficient Numerical Computation.” Computing in Science & Engineering 13 (2): 22–30. https://doi.org/10.1109/MCSE.2011.37.

McKinney, Wes, and Others. 2010. “Data Structures for Statistical Computing in Python.” In Proceedings of the 9th Python in Science Conference, 445:51–56. Austin, TX. http://conference.scipy.org/proceedings/scipy2010/pdfs/mckinney.pdf
