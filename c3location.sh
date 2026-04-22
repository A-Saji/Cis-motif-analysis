#!/bin/bash

output="fimo_c3_positions.tsv"

echo -e "Gene\tMotif_ID\tSequence\tStart\tEnd\tStrand\tp_value\tq_value" > "$output"

for dir in GRMZM*
do
    gene=$(basename "$dir")

    fimo_file="fimo_results/$gene/fimo.tsv"

    if [ ! -f "$fimo_file" ]; then
        continue
    fi

    awk -v g="$gene" 'BEGIN{FS=OFS="\t"}
    NR==1 {next}
    {
        print g, $2, $3, $4, $5, $6, $8, $9
    }' "$fimo_file" >> "$output"

done
