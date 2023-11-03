# Probe design scripts

Scripts and code for tiling of genomes for capture probe design.

Required input:
-I input genomes paths file -- path to each genome


For downloading genomes for input, consider using NCBI datasets tool

`datasets download genome accession --inputfile ${accessionpathfile} --exclude-gff3 --exclude-protein --exclude-rna --exclude-genomic-cds --filename genomes_accession_file.zip`

or `datasets download genome taxon ${taxonid} --assembly-source {assemblysource} --assembly-level ${assemblylevel} --exclude-gff3 --exclude-protein --exclude-rna --exclude-genomic-cds --filename genomes_taxon_file.zip`