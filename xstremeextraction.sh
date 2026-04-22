#!/bin/bash

output="meme_positions.tsv"

echo -e "Gene\tMotif_ID\tSequence\tStrand\tStart\tp_value" > "$output"

for dir in GRMZM*
do
    gene=$(basename "$dir")
    meme_file="$dir/${gene}_OG_xst_out/meme_out/meme.txt"

    if [ ! -f "$meme_file" ]; then
        echo "Skipping $gene"
        continue
    fi

    echo "Processing $gene"

    awk -v g="$gene" '
    BEGIN {FS="[[:space:]]+"}

    # Step 1: capture motif ID
    /Motif .* MEME-[0-9]+/ {
        for (i=1; i<=NF; i++) {
            if ($i ~ /^MEME-/) {
                motif = $i
            }
        }
        next
    }

    # Step 2: detect header separator (start of table)
    /^[-]+[[:space:]]+[-]+/ {
        read_data = 1
        next
    }

    # Step 3: detect full dashed line (end of table)
    /^-+$/ && read_data == 1 {
        read_data = 0
        next
    }

    # Step 4: extract data rows
    read_data == 1 && NF >= 4 {
        seq = $1
        strand = $2
        start = $3
        pval = $4

        # skip header row accidentally captured
        if (seq != "Sequence") {
            print g, motif, seq, strand, start, pval
        }
    }

    ' "$meme_file" >> "$output"

done
