import os
import pandas as pd
import re
# ----------------------------
# INPUT
# ----------------------------
table = pd.read_csv("allmotifs.csv", sep="\t")
print(table.columns)
output_file = "final_selected.meme"

# ----------------------------
# Write MEME header
# ----------------------------
with open(output_file, "w") as out:
    out.write("""MEME version 4

ALPHABET= ACGT

strands: + -

Background letter frequencies:
A 0.25 C 0.25 G 0.25 T 0.25

""")

    # ----------------------------
    # Process each gene
    # ----------------------------
    import re

    for gene in table['Gene'].unique():
    
        meme_path = f"{gene}/{gene}_OG_xst_out/combined.meme"
    
        if not os.path.exists(meme_path):
            print(f"Missing: {meme_path}")
            continue
    
        # clean motif list
        motifs_to_extract = (
            table[table['Gene'] == gene]['Motif_ID']
            .astype(str)
            .str.strip()
            .tolist()
        )
    
        print(f"Checking {gene}: {motifs_to_extract}")
        
        with open(meme_path) as f:
            lines = f.readlines()
    
        i = 0
        while i < len(lines):
            line = lines[i]
        
            if line.startswith("MOTIF"):
            
            # extract MEME-X using regex
                match = re.search(r"MEME-\d+", line)
                if not match:
                    i += 1
                    continue
            
                motif_name = match.group().strip()
            
                print(f"{gene}: Found {motif_name}")
            
                if motif_name in motifs_to_extract:
                
                    new_name = f"{gene}_{motif_name}"
                
                # write new header
                    out.write(f"MOTIF {new_name}\n")
                
                    i += 1
                
                # copy everything until next MOTIF
                    while i < len(lines) and not lines[i].startswith("MOTIF"):
                        out.write(lines[i])
                        i += 1
                
                    continue  # skip increment
            
            i += 1

print("Done! Output: final_selected.meme")
