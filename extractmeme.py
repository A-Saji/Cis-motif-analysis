import sys

meme_file = "final_selected.meme"
motif_list_file = "test_C4_shifted_motif_uid.txt"
output_file = "test_C4_shifted_motif_pwms.meme"

# --- load motif IDs ---
with open(motif_list_file) as f:
    motifs_to_keep = set(line.strip() for line in f)

# --- parse MEME file ---
with open(meme_file) as f:
    lines = f.readlines()

out = []
keep = False
current_motif = None

for line in lines:
    
    if line.startswith("MOTIF"):
        parts = line.strip().split()
        current_motif = parts[1]  # e.g. GRMZM2G088242_MEME-9
        
        if current_motif in motifs_to_keep:
            keep = True
        else:
            keep = False
    
    if keep:
        out.append(line)

# --- write output ---
with open(output_file, "w") as f:
    f.writelines(out)

print(f"Extracted {len(motifs_to_keep)} motifs (if present) to {output_file}")
