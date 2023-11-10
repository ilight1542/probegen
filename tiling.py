import argparse
import math
import random
from Bio import SeqIO
import gzip

parser = argparse.ArgumentParser(description='tile a set of genomes for capture probe. mandatory fields, -i, -l, -s')
parser.add_argument('-i', '--input', metavar='input genomes text file (line sep)', \
    required=True, help="input genomes text file to be made into probe")
parser.add_argument('-l', '--length', metavar="length of probe",\
    required=True, help="length of probe (before any addition of adapters")
parser.add_argument('-s', '--step_size', metavar="step size",\
    required=True, help="number of bases between probes")
parser.add_argument('--out_name', help="output file name")
parser.add_argument('--masked_regions', metavar="dustmasker output acclist",\
    required = False, help="fasta header --> low complexity regions, single file")
parser.add_argument('--masked_cutoff', metavar="masked probe percentage",\
    required = False, help="percentage of probe that can be masked before discarding (eg 10 == up to 10 percent of read can be masked), default =10")
parser.add_argument('--reverse_complement', action='store_true', help="filter cross probe reverse complements (default = no filtering done, functionality contained within cd-hit)")
parser.add_argument('--output_reverse_complement', action='store_true',help="output both probe and a probe's reverse complement to fasta or text file (default = output only probe)")
parser.add_argument('--outfmt', metavar="txt, fasta, all", help="output format, .txt or .fasta, default =txt")
parser.add_argument('--no_convert_n', action='store_true', help="do not convert N to random base for probe set")
parser.add_argument('--randseed', help="random seed for consistent probe construction across multiple runs")
parser.add_argument('--add_to_existing', help="probe set which we want to add new probes to from genomes marked at input, \
    shoudld be in either fasta format with a header with probe name followed by the probe, or txt with a probe per line")
args=parser.parse_args()

## required argumets
## sliding window approach with set corresponding to masked regions
## eg set(10,11,12,..,16) for window 10..70 has 7 masked basepairs
## then just set.remove(10) if 10 in set and add 71 if 71 in masked set

## functions for parsing optional arguments
def parse_masked(masked_input):
    masked_dict_sets={}
    with open(masked_input, newline="\n") as f:
        for line in f:
            (key,start,end)=line.split('\t')
            key=key.split(' ')[0][1:] ## format to match seqIO record.id, just accession without >
            start=int(start)
            end=int(end)
            if key in masked_dict_sets:
                to_union=set(range(start, end+1))
                masked_dict_sets[key].union(to_union)
            else:
                masked_dict_sets[key] = set(range(start, end+1))
    return masked_dict_sets

def parse_probes(prior_probes):
    """parse prior probe set to add new probes to"""
    set_to_output=set()
    if prior_probes.split(".")[-1] in ["fasta", "fna", "fa"]:
        with open(prior_probes) as f:
            for probe in f:
                set_to_output.add(probe.seq)
    else:
        with open(prior_probes, newline="\n") as f:
            for line in f:
                set_to_output.add(line.strip())
    return set_to_output

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

def tiling_masked(input_fastas, length, step_size, masked_cutoff, masked_dict_sets, check_reverse_complement=False, output_reverse_complement=False, convert_n=False,add_probes=None):
    """tiling for genomes when dustmasker should be taken into account, slower than tiling() but both are linear in time complexity"""
    if add_probes is None:
        tiling_set=set()
    else:
        tiling_set=add_probes
    #reverse complement checking not necessary if using cd-hit, reverse complement removal done by cd-hit-est
    tiling_set_reverse_complements={}
    with open(input_fastas) as fasta_f:
        input_fasta_paths = [line.rstrip() for line in fasta_f]
    for current_input in input_fasta_paths:
        with gzip.open(current_input,'rt') as f:
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
                            if output_reverse_complement:
                                tiling_set_reverse_complements[tile]=rev_comp_tile
                        remove_from_set(masked_set, start, step_size, length, masked_dict_sets[current_header])
                        start+=step_size
                    line=line[step_size:]
    if output_reverse_complement:                
        return [tiling_set,tiling_set_reverse_complements]
    else:
        return tiling_set

def tiling(input_fastas, length, step_size, check_reverse_complement=False, output_reverse_complement=False, convert_n=True, add_probes=None):
    """tiling, just base genome, no masking checking"""
    if add_probes is None:
        tiling_set=set()
    else:
        tiling_set=add_probes
    tiling_set_reverse_complements={}
    with open(input_fastas) as fasta_f:
        input_fasta_paths = [line.rstrip() for line in fasta_f]
    for current_input in input_fasta_paths:
        # iterate over each input file (fasta without linebreaks)
        with gzip.open(current_input,'rt') as f:
            for record in SeqIO.parse(f, "fasta"):
            ## get section of fasta ready for tiling
                current_header=record.id
                line=str(record.seq)
                iter_needed=math.ceil((len(line)-length)/int(step_size))
            ## begin processing tiles
                for iterations in range(iter_needed):
                    tile=line[:length] ## get tile of given length
                    if convert_n or check_reverse_complement or output_reverse_complement:
                        [rev_comp_tile,rc_conv_n]=reverse_complement(tile, convert_n)
                        if convert_n:
                            tile=rc_conv_n
                    if check_reverse_complement and rev_comp_tile in tiling_set:
                        pass
                    else:
                        if len(tile) == length: ## check if already in set
                            tiling_set.add(tile)
                        if output_reverse_complement:
                            tiling_set_reverse_complements[tile] = rev_comp_tile
                    line=line[step_size:]
    if output_reverse_complement:                
        return [tiling_set,tiling_set_reverse_complements]
    else:
        return tiling_set

def write_fasta(tiling_set,filename="probe_set.fasta", output_reverse_complement=False):
    """
    output in fasta format with accession headers being a growing index, if reverse_complement==True, will output forward probes first, then the corresponding reverse complements with the same index+_rev_complement
    useful for cd-hit clustering to remove any probe with reverse complement getting into clusters
    """
    if '.fasta' not in filename:
        filename+='.fasta'
    index=0
    if output_reverse_complement:
        with open(filename, "w") as file_out:
            probes=list(tiling_set[1].keys()) ## parse probes
            for line in probes:
                file_out.write(">"+str(index)+"\n")
                file_out.write(line+"\n")
                index+=1
            index=0
            for line in probes: ## parse probe's complement
                file_out.write(">"+str(index)+"_rev_complement"+"\n")
                file_out.write(tiling_set[1][line]+"\n")
                index+=1
    else:
        probes=tiling_set
        with open(filename, "w") as file_out:
            for line in probes:
                file_out.write(">"+str(index)+"\n")
                file_out.write(line+"\n")
                index+=1

def write_text(tiling_set,filename="probe_set.txt", output_reverse_complement=False):
    """output tiles as a single probe per line file"""
    if '.txt' not in filename:
            filename+='.txt'
    if output_reverse_complement:
        probes=list(tiling_set[1].keys)+list(tiling_set[1].values)
    else: probes=tiling_set
    with open(filename, "w") as file_out:
        file_out.write('\n'.join(probes))

def write_output(tiling_set,filename="probe_set.txt", outfmt="txt", output_reverse_complement=False):
    """sends tiling set (either just base tiling set or set with reverse complement (dict)) to output writing functions"""
    if outfmt=="txt":
        write_text(tiling_set,filename,output_reverse_complement)
    elif outfmt=="fasta":
        write_fasta(tiling_set,filename,output_reverse_complement)
    elif outfmt=="all":
        write_text(tiling_set,filename,output_reverse_complement)
        write_fasta(tiling_set,filename,output_reverse_complement)

def main(args):
    """executing function, will call other helper functions to create probe set and output"""
    ## parse arguments for input to varisou functions 
    ## mandatory arg parsing
    length=int(args.length)
    step_size=int(args.step_size)

    ## optional arg parsing that have defaults and are not automatically parsed 
    ## NOTE: this may be implementable in argparser at top, but I am lazy
    ## convert n
    if args.no_convert_n:
        convert_n=False
    else: convert_n=True
    ## output format
    if args.outfmt:
        outfmt_parsed=args.outfmt
    else: 
        print("No option passed to output format (--outfmt), using default output of txt")
        outfmt_parsed="txt"
    ## rand seed
    if args.randseed:
        random.seed(args.randseed)
    else: 
        print("No option passed to random seed (--randseed), using default value of 100")
        random.seed(100)
    
    if args.add_to_existing is None:
        add_probes=None
    else:
        add_probes=parse_probes(args.add_to_existing)

    #######################
    ## execute functions ##
    #######################
    if args.masked_regions is None:
        write_output(
            tiling(args.input, length, step_size,\
            args.reverse_complement, args.output_reverse_complement, convert_n, add_probes),\
        args.out_name, outfmt_parsed, args.output_reverse_complement)
    else:
        masked_regions_parsed=parse_masked(args.masked_regions)
        ## check masked_cutoff set, if not use default of 10%
        if args.masked_cutoff is None:
            masked_cutoff_parsed = 10
            print("No option passed to masked cutoff (--masked_cutoff), using default value of 10%")
        else: masked_cutoff_parsed = int(args.masked_cutoff)
        write_output(
            tiling_masked(args.input, length, step_size,\
            masked_cutoff_parsed, masked_regions_parsed,\
            args.reverse_complement, args.output_reverse_complement, convert_n, add_probes),\
        args.out_name, outfmt_parsed, args.output_reverse_complement)

if __name__ == '__main__':
    main(args)
    
