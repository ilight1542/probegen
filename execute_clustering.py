import argparse
import subprocess

parser = argparse.ArgumentParser(description='Helper script for checking existance of all reference genomes. Fields, -i')
parser.add_argument('-f', '--fasta', metavar='fasta', \
    required=True, help="fasta with unclustered probes")
parser.add_argument('-o', '--output', metavar='fasta', \
    required=True, help="output fasta with clustered probes")
parser.add_argument('-l', '--length', metavar='length of probe', \
    required=True, help="Lenght of probes to create")
parser.add_argument('-i', '--percent_identity_aligned_portion', help='Minimum percent identity for aligned portion of reads (excluding terminal mismatches)', \
    required=True)
parser.add_argument('-t', '--max_terminal_mismatches', \
    required=True, help='Max terminal mismatches \
    \\n Important parameter, since a 1bp insertion at the start of two otherwise identical fastas would cause 2x the probes if step size is even. \
    However, setting a maximum termianl of 1bp merges all the probes. Thus setting minimum overlap at least 1 or 2 is recommended. \
    \\n Larger values will cause larger indels or start/end SNPs to be ignored and merged if the aligned portion is above the percent_identity_aligned_portion threshold.', \
    default=2)

args=parser.parse_args()
length=int(args.length)
percent_identity_aligned_portion=float(args.percent_identity_aligned_portion)
max_terminal_mismatches=int(args.max_terminal_mismatches)


def main(input_fasta, output_fasta, length, percent_identity_aligned_portion, max_terminal_mismatches):
    min_overlap=length-max_terminal_mismatches
    minimumpercentidentity=percent_identity_aligned_portion/100
    subprocess.run([ f"cd-hit-est -i {input_fasta} -r 1 -G 0 -o {output_fasta} -n 9 -A {min_overlap} -d 0 -c {minimumpercentidentity} -p 1"  ],shell=True)

if __name__ == '__main__':
    main(args.fasta, args.output, length, percent_identity_aligned_portion,max_terminal_mismatches)

