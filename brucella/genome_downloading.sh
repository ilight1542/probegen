## probe design



## once downloaded

for fasta in */*.fna
    do
    awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' ${fasta} > ${fasta}.formatted
done