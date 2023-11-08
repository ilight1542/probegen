import argparse
import sys
import os

parser = argparse.ArgumentParser(description='Helper script for checking existance of all reference genomes. Fields, -i')
parser.add_argument('-i', '--input', metavar='input genomes text file', \
    required=True, help="Text file of input genome paths to be made into probe")

args=parser.parse_args()

def main(input_fasta_paths_file):
    not_found_fasta=[]
    with open(input_fasta_paths_file) as f:
        for path in f:
            path=path.strip()
            if not os.path.isfile(f'{path}'):
                not_found_fasta.append(f'{path}')
    if len(not_found_fasta) > 0:
        for e in not_found_fasta:
            print('ERROR: The following fasta file paths do not exist. Please correct them within the file paths input.')
            print(e)
        sys.exit(1)
    else: sys.exit(0)

if __name__ == '__main__':
    main(args.input)

