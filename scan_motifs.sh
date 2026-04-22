#!/bin/bash

mkdir -p fimo_results

output="combined_motif_table.tsv"

echo -e "Gene\tMotif_ID\tMotif_Name\tTotal_C3_promoters\tC3_promoters_with_motif\tSequences_with_motif\tBest_pvalue\tBest_qvalue\tOEL_present\tChala_present" > "$output"

for gene_dir in GRMZM*
do
    gene=$(basename "$gene_dir")

    motif_file="$gene_dir/${gene}_OG_xst_out/xstreme.txt"
    C3_fasta="$gene_dir/${gene}_OG_C3.fasta"

    if [ ! -f "$motif_file" ]; then
        echo "Skipping $gene (no motifs)"
        continue
    fi

    echo "Processing $gene"

    clean_fasta="$gene_dir/${gene}_OG_C3_clean.fasta"

    # --- Fix FASTA headers so FIMO keeps only the first token ---
    awk '/^>/ {print $1; next} {print}' "$C3_fasta" > "$clean_fasta"

    total_promoters=$(grep -c "^>" "$clean_fasta")

    # detect ortholog presence
    grep -q "OEL" "$clean_fasta" && oel_status="YES" || oel_status="NO"
    grep -q "Chala" "$clean_fasta" && chala_status="YES" || chala_status="NO"

    # run FIMO
    fimo --thresh 1e-6 --oc "fimo_results/$gene" "$motif_file" "$clean_fasta"

    grep "^MOTIF" "$motif_file" | grep "MEME-" > tmp_meme_motifs.txt

    while read line
    do
        motif_name=$(echo "$line" | awk '{print $2}')
        motif_id=$(echo "$line" | awk '{print $3}')

        if [ -f "fimo_results/$gene/fimo.tsv" ]; then

            seqs=$(awk -v m="$motif_id" 'NR>1 && $2==m {print $3}' fimo_results/$gene/fimo.tsv | sort | uniq)

            count=$(echo "$seqs" | grep -c .)

            seq_list=$(echo "$seqs" | paste -sd ";" -)

            best_p=$(awk -v m="$motif_id" 'NR>1 && $2==m {print $8}' fimo_results/$gene/fimo.tsv | sort -g | head -1)

            best_q=$(awk -v m="$motif_id" 'NR>1 && $2==m {print $9}' fimo_results/$gene/fimo.tsv | sort -g | head -1)

        else
            count=0
            seq_list="NA"
            best_p="NA"
            best_q="NA"
        fi

        echo -e "$gene\t$motif_id\t$motif_name\t$total_promoters\t$count\t$seq_list\t$best_p\t$best_q\t$oel_status\t$chala_status" >> "$output"

    done < tmp_meme_motifs.txt

done

rm tmp_meme_motifs.txt
