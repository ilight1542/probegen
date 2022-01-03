import argparse
import math
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description='tile a set of genomes for probe')
parser.add_argument('-i', '--input', metavar="input fasta files (remove all but first linebreaks)",\
     required=True, nargs='+', help="input genomes to be made into probe")
parser.add_argument('-l', '--length', metavar="length of probe",\
     required=True, nargs=1, help="length of probe (before any addition of adapters")
parser.add_argument('-s', '--step_size', metavar="step size of tiles",\
     required=True, nargs=1, help="number of bases between probes")
args=parser.parse_args()

length=int(args.length[0])
step_size=int(args.step_size[0])

tiling_set=set()
new_uniq_tiles_per_ref={}
for current_input in args.input:
    with open(current_input, newline="\n") as f:
        i=0
        for line in f:
            print("line")
            if line.startswith(">"):
                new_uniq_tiles_per_ref[str(line)] = 0
                current_header=str(line)
            else:
                tiles_needed=math.ceil(len(line)//int(step_size))
                for iterations in range(tiles_needed):
                    tile=line[:length]
                    if tile in tiling_set:
                        pass
                    else:
                        if tile == "" or tile =="\n":
                            pass
                        else:
                            tiling_set.add(tile)
                            new_uniq_tiles_per_ref[current_header]+=1
                    line=line[length:]

print(len(tiling_set))
print(new_uniq_tiles_per_ref)

with open("probe_set.txt", "w") as file_out:
    file_out.write('\n'.join(list(tiling_set)))