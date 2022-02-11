import argparse
import math
import random
from Bio import SeqIO

parser = argparse.ArgumentParser(description='tile a set of genomes for capture probe. mandatory fields, -i, -l, -s')
parser.add_argument('-i', '--input', metavar='input genomes', \
    required=True, nargs='+', help="input genomes to be made into probe")
parser.add_argument('-l', '--length', metavar="length of probe",\
    required=True, nargs=1, help="length of probe (before any addition of adapters")
parser.add_argument('-s', '--step_size', metavar="step size",\
    required=True, nargs=1, help="number of bases between probes")
parser.add_argument('--out_name', nargs=1, help="output file name")
parser.add_argument('--masked_regions', metavar="dustmasker output acclist",\
    required = False, nargs=1, help="fasta header --> low complexity regions, single file")
parser.add_argument('--masked_cutoff', metavar="masked probe percentage",\
    required = False, nargs=1, help="percentage of probe that can be masked before discarding (eg 10 == up to 10 percent of read can be masked)")
parser.add_argument('--reverse_complement', action='store_true', help="filter cross probe reverse complements (default = no filtering)")
parser.add_argument('--outfmt', metavar="txt, fasta, all", help="output format, .txt or .fasta")
parser.add_argument('--convert_n', action='store_true', help="convert N to random base for probe set")
parser.add_argument('--randseed', nargs=1, help="random seed for consistent probe construction across multiple runs")
args=parser.parse_args()

## required argumets
if args.outfmt:
    outfmt=args.outfmt
else: outfmt="txt"
length=int(args.length[0])
step_size=int(args.step_size[0])

random.seed=(args.randseed[0])

## sliding window approach with set corresponding to masked regions
## eg set(10,11,12,..,16) for window 10..70 has 6 masked reads
## then just set.remove(10) if 10 in set and add 71 if 71 in masked set

## functions for parsing optional arguments
def parse_masked(masked_input):
    masked_dict_sets={}
    with open(masked_input, newline="\n") as f:
        for line in f:
            (key,start,end)=line.split('\t')
            key=key.split(' ')[0][1:] ## format to match seqIO record.id
            start=int(start)
            end=int(end)
            if key in masked_dict_sets:
                for i in range(start, end+1):
                    masked_dict_sets[key].add(i)
            else:
                masked_dict_sets[key] = set()
                for i in range(start, end+1):
                    masked_dict_sets[key].add(i)
    return masked_dict_sets

def reverse_complement(tile, convert_n=False):
    """generates reverse complement of tile, and will also mask any N (or other ambiguous base calls) as a random base (of possible bases)"""
    rev_comp_tile=""
    comp={"A":"T","T":"A","C":"G","G":"C","N":"N"}
    bases=["A","T","C","G"]
    non_atcg={"N": bases, "R":["A","G"], "Y":["C","T"], "K":["G","T"], "M":["A","C"],"S":["C","G"],"W":["A","T"],"B":["T","C","G"],"D":["A","T","G"],"H":["A","T","C"],"V":["A","C","G"]}
    for i in range(len(tile),0,-1):
        idx=i-1
        if convert_n:
            to_add=tile[idx]
            if to_add not in bases: 
                to_add=random.choice(non_atcg[to_add])
                tile=tile[:idx]+to_add+tile[idx+1:]
            rev_comp_tile+=comp[to_add]
        else:
            rev_comp_tile+=comp[tile[idx]]
    return [rev_comp_tile,tile]

def remove_from_set(set, start, step_size, len, masked_set):
    """updates set counting how many bases in current tile are marked as low complexity"""
    to_remove=range(start,start+step_size)
    to_add=range(start+len,start+len+step_size)
    for i in to_remove:
        if i in set:
            set.remove(i)
    for i in to_add:
        if i in masked_set:
            set.add(i)

def tiling_masked(input_fastas, length, step_size, masked_cutoff, masked_dict_sets, check_reverse_complement=False, convert_n=False):
    """tiling for genomes when dustmasker should be taken into account, slower than tiling() but both are linear in time complexity"""
    tiling_set=set()
    tiling_set_reverse_complements={}
    for current_input in input_fastas:
        with open(current_input) as f:
            for record in SeqIO.parse(f, "fasta"):
            ## get section of fasta ready for tiling
                current_header=record.id
                line=str(record.seq)
                iter_needed=math.ceil((len(line)-length)/int(step_size))

            ## initialize set of masked bases for this chromosome
                start=0  ## used only if masking regions
                masked_set=set()
                for i in range(start,length): 
                    if i in masked_dict_sets[current_header]: 
                        masked_set.add(i) ## add index to masked set if in masked_dict_sets (preprocessing)

            ## begin processing tiles
                for iterations in range(iter_needed):
                    tile=line[:length] ## get tile of given length
                    if convert_n or check_reverse_complement:
                        [rev_comp_tile,rc_conv_n]=reverse_complement(tile, convert_n)
                        if convert_n:
                            tile=rc_conv_n
                    if check_reverse_complement and rev_comp_tile in tiling_set:
                        pass
                    elif tile in tiling_set: ## check if already in set
                        pass
                    else: ## check if the len of the set is above cutoff, if not add tile, if caring about reverse complement add rev to dict with key(tile) --> value(rev_compl)
                        if len(masked_set) < ((masked_cutoff/100)*length):
                            tiling_set.add(tile)
                            if check_reverse_complement:
                                tiling_set_reverse_complements[tile]=rev_comp_tile
                        remove_from_set(masked_set, start, step_size, length, masked_dict_sets[current_header])
                        start+=step_size
                    line=line[step_size:]
    return [tiling_set,tiling_set_reverse_complements]


## main function to create probe set
def tiling(input_fastas, length, step_size, check_reverse_complement=False, convert_n=True):
    """tiling, just base genome, no masking checking"""
    tiling_set=set()
    for current_input in input_fastas:
        # iterate over each input file (fasta without linebreaks)
         with open(current_input) as f:
            for record in SeqIO.parse(f, "fasta"):
            ## get section of fasta ready for tiling
                current_header=record.id
                line=str(record.seq)
                iter_needed=math.ceil((len(line)-length)/int(step_size))
            ## begin processing tiles
                for iterations in range(iter_needed):
                    tile=line[:length] ## get tile of given length
                    if convert_n or check_reverse_complement:
                        rc_conv_n=reverse_complement(tile, convert_n)
                        if convert_n:
                            tile=rc_conv_n[1]
                    if check_reverse_complement and rc_conv_n[0] in tiling_set:
                        pass
                    elif tile in tiling_set: ## check if already in set
                        pass
                    else:
                        if len(tile) == length: ## 
                            tiling_set.add(tile)
                    line=line[step_size:]
    return tiling_set

def write_fasta(tiling_set,filename="probe_set.fasta", reverse_complement=False, tuple=False):
    """
    output in fasta format with accession headers being a growing index, if reverse_complement==True, will output forward probes first, then the corresponding reverse complements with the same index+_rev_complement
    useful for cd-hit clustering to remove any probe with reverse complement getting into clusters
    """
    if '.fasta' not in filename:
        filename+='.fasta'
    index=0
    if reverse_complement:
        with open(filename, "w") as file_out:
            probes=list(tiling_set[1].keys())
            for line in probes:
                file_out.write(">"+str(index)+"\n")
                file_out.write(line+"\n")
                index+=1
            index=0
            for line in probes:
                file_out.write(">"+str(index)+"_rev_complement"+"\n")
                file_out.write(tiling_set[1][line]+"\n")
                index+=1
    else:
        if tuple:
            with open(filename, "w") as file_out:
                for line in list(tiling_set[0]):
                    file_out.write(">"+str(index)+"\n")
                    file_out.write(line+"\n")
                    index+=1
        else:
            probes=tiling_set
            with open(filename, "w") as file_out:
                for line in probes:
                    file_out.write(">"+str(index)+"\n")
                    file_out.write(line+"\n")
                    index+=1


def write_text(tiling_set,filename="probe_set.txt",reverse_complement=False, tuple=False):
    """output tiles as a single probe per line file"""
    if '.txt' not in filename:
            filename+='.txt'
    if tuple:
        if reverse_complement:
            probes=list(tiling_set[1].keys)+list(tiling_set[1].values)
        else: probes=list(tiling_set[0])
    else: probes=tiling_set
    with open(filename, "w") as file_out:
        file_out.write('\n'.join(probes))

def write_output(tiling_set,filename="probe_set.txt", outfmt="txt",reverse_complement=False, tuple=False):
    """sends tiling set (either just base tiling set or set with reverse complement (dict)) to output writing functions"""
    if outfmt=="txt":
        write_text(tiling_set,filename,reverse_complement,tuple)
    elif outfmt=="fasta":
        write_fasta(tiling_set,filename,reverse_complement,tuple)
    elif outfmt=="all":
        write_text(tiling_set,filename,reverse_complement,tuple)
        write_fasta(tiling_set,filename,reverse_complement,tuple)
        

if __name__ == '__main__':
    if args.masked_regions:
        print("tiling masked")
        masked_regions=parse_masked(args.masked_regions[0])
        write_output(tiling_masked(args.input, length, step_size, int(args.masked_cutoff[0]), masked_regions, args.reverse_complement, args.convert_n), args.out_name[0], outfmt,args.reverse_complement, tuple=True)
    else:
        write_output(tiling(args.input, length, step_size, args.reverse_complement, args.convert_n),args.out_name[0], outfmt, tuple=False)