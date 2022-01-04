import argparse
import math

parser = argparse.ArgumentParser(description='tile a set of genomes for probe')
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
args=parser.parse_args()

## required argumets
length=int(args.length[0])
step_size=int(args.step_size[0])

## functions for parsing optional arguments
def parse_masked(masked_input):
    masked_dict={}
    with open(masked_input, newline="\n") as f:
        for line in f:
            (key,start,end)=line.split('\t')
            start=int(start)
            end=int(end)
            if key in masked_dict:
                masked_dict[key].append((start,end))
            else:
                masked_dict[key] = []
                masked_dict[key].append((start,end))
    return masked_dict

def is_overlapping(start,length,segments):
    ## return number of reads overlapping with masked regions for each probe
    ## essentially just Binary search since our lists of intervals are sorted
    start_index=0
    end_index=len(segments)
    found=False
    end = start + length
    while not found and start_index != end_index:
        search_index=(start_index+end_index)//2
        if start >= segments[search_index][1]:
            start_index=search_index+1
        elif end <= segments[search_index][0]:
            end_index=search_index
        else:
            overlap = min(end,segments[search_index][1]) - max(start,segments[search_index][0])
            found=True
            return overlap
    return 0

def write_output(tiling_set,filename="probe_set.txt"):
    if '.txt' not in filename:
        filename+='.txt'
    with open(filename, "w") as file_out:
        file_out.write('\n'.join(list(tiling_set)))

## main function to create probe set
def tiling(input_fastas, length,step_size, masked_cutoff=10, overlap=None):
    print(overlap)
    print((masked_cutoff/100)*length)
    tiling_set=set()
    new_uniq_tiles_per_ref={}
    for current_input in input_fastas:
        # iterate over each input file (fasta without linebreaks)
        with open(current_input, newline="\n") as f:
            for line in f:
                if line.startswith(">"): 
                    ## start of new accession in fasta
                    current_header=str(line)[:-1]
                    new_uniq_tiles_per_ref[current_header] = 0
                else:
                    ## bases
                    iter_needed=math.ceil(len(line)/int(step_size))
                    start=0  ## used only if masking regions
                    for iterations in range(iter_needed):
                        tile=line[:length] ## get tile of given length
                        if tile in tiling_set: ## check if already in set
                            pass
                        else:
                            if tile == "" or tile =="\n":
                                pass
                            else:
                                if overlap != None: ## overlap dictionary provided
                                    current_overlap=is_overlapping(start, length, overlap[current_header])
                                    if current_overlap > 30:
                                        print(tile)
                                    if len(tile) == length and current_overlap <= ((masked_cutoff/100)*length):
                                        tiling_set.add(tile)
                                        new_uniq_tiles_per_ref[current_header]+=1
                                    start+=step_size
                                else:
                                    if len(tile) == length: ## 
                                        tiling_set.add(tile)
                                        new_uniq_tiles_per_ref[current_header]+=1
                        line=line[step_size:]
    return (tiling_set,new_uniq_tiles_per_ref)

if __name__ == '__main__':
    if args.masked_regions:
        masked_regions=parse_masked(args.masked_regions[0])
    write_output(tiling(args.input, length, step_size, int(args.masked_cutoff[0]), masked_regions)[0], args.out_name[0])
