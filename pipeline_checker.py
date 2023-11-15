import argparse
import sys
import os

parser = argparse.ArgumentParser(description='Helper script for checking existance of all reference genomes. Fields, -i')
parser.add_argument('-i', '--input', metavar='input genomes text file', \
    required=True, help="Text file of input genome paths to be made into probe")
parser.add_argument('--stepsize', metavar='stepsize', \
    required=True, help="stepsize for probe")
parser.add_argument('--length',required=True)
parser.add_argument('--maxterminalmismatches',required=True)
parser.add_argument('--maskedthreshold',required=True)
parser.add_argument('--minimumpercentidentity',required=True)
parser.add_argument('--percentambiguousbasethreshold', required=True)

args=parser.parse_args()

def check_fasta_paths(input_fasta_paths_file):
    not_found_fasta=[]
    with open(input_fasta_paths_file) as f:
        for path in f:
            path=path.strip()
            if not os.path.isfile(f'{path}'):
                not_found_fasta.append(f'{path}')
    if len(not_found_fasta) > 0:
        print(f'INPUT ERROR: The following fasta file paths declared in {input_fasta_paths_file} do not exist. Please correct them within the file paths input.')
        for e in not_found_fasta:
            print(f'    {e}')
        return False
    else: 
        return True

def check_length_bounds(length):
    if int(length) < 2:
        raise ValueError('PARAMETER ERROR: Length of probe must be >=2!')
    
def check_maxterminalmismatches_bounds(maxterminalmismatches,length):
    if int(maxterminalmismatches) >= int(length):
        raise ValueError('PARAMETER ERROR: Max terminal mismatches cannot be equal to or greater than the total lenght of the probe!')

def check_maskedthreshold_bounds(maskedthreshold,length):
    if int(maskedthreshold) >= int(length):
        raise ValueError('PARAMETER ERROR: Max masked bases allowed in probe cannot be equal to or greater than the total lenght of the probe!')

def check_step_size_bounds(stepsize):
    if int(stepsize) < 1:
        raise ValueError('PARAMETER ERROR: Step size must a value >=1!')

def check_minpercentidentity_bounds(minimumpercentidentity):
    if int(minimumpercentidentity) < 80 or int(minimumpercentidentity) >= 100:
        raise ValueError('PARAMETER ERROR: Minimum percent identity for probe clustering must be >=80% and <100%')

def check_percentambiguousbasethreshold_bounds(percentambiguousbasethreshold):
    if int(percentambiguousbasethreshold) < 0 or int(percentambiguousbasethreshold) > 100:
        raise ValueError('PARAMETER ERROR: Percent ambiguous bases allowed in fasta must be >=0% and <=100%')

def main(input_fasta_paths_file,stepsize,length,maxterminalmismatches,maskedthreshold,minimumpercentidentity,percentambiguousbasethreshold):
    error_list=[]
    paths_ok=check_fasta_paths(input_fasta_paths_file)
    try: check_step_size_bounds(stepsize)
    except ValueError as e:
        print('error encountered')
        error_list.append(e)
    try: check_length_bounds(length)
    except ValueError as e:
        error_list.append(e)
    try: check_maxterminalmismatches_bounds(maxterminalmismatches,length)
    except ValueError as e:
        error_list.append(e)
    try: check_maskedthreshold_bounds(maskedthreshold,length)
    except ValueError as e:
        error_list.append(e)
    try: check_minpercentidentity_bounds(minimumpercentidentity)
    except ValueError as e:
        error_list.append(e)
    try: check_percentambiguousbasethreshold_bounds(percentambiguousbasethreshold)
    except ValueError as e:
        error_list.append(e)
    if len(error_list) > 0:
        for e in error_list:
            print(e)
        sys.exit(1)
    if not paths_ok:
        sys.exit(1)
    sys.exit(0)
    

if __name__ == '__main__':
    main(args.input,args.stepsize,args.length,args.maxterminalmismatches,args.maskedthreshold,args.minimumpercentidentity,args.percentambiguousbasethreshold)
