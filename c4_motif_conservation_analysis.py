#!/usr/bin/env python3
"""
Extract motifs conserved across C4 plant species from XSTREME/MEME output.

C4 Plant Species Analyzed:
- GRMZM: Zea mays (Maize)
- Pahal: Panicum hallii (Hall's panicgrass)
- Sevir: Setaria viridis (Green foxtail)
- Urofu: Urochloa fusca
- ELECO: Eleusine coracana (Finger millet)
- Pavag: Paspalum vaginatum (Seashore paspalum)

Usage:
    python c4_motif_conservation_analysis.py <gene_name>
    
Example:
    python c4_motif_conservation_analysis.py GRMZM2G000171
    
The script will look for files in:
    <BASE_DIR>/upstream_seqs/<gene_name>/<gene_name>_xst_out/
"""

import pandas as pd
import re
import sys
from pathlib import Path


# C4 plant species dictionary
C4_PLANTS_DICT = {
    'Zea mays': 'GRMZM',
    'Panicum halli': 'Pahal',
    'Setaria viridis': 'Sevir',
    'Urochloa fusca': 'Urofu',
    'Eleusine coracana': 'ELECO',
    'Paspalum vaginatum': 'Pavag'
}

# Species prefixes to check for conservation
C4_SPECIES = list(C4_PLANTS_DICT.values())  # ['GRMZM', 'Pahal', 'Sevir', 'Urofu', 'ELECO', 'Pavag']


def parse_xstreme_motifs(xstreme_path):
    """
    Extract motif IDs, consensus sequences, and E-values from XSTREME output.
    
    Args:
        xstreme_path: Path to xstreme.txt file
        
    Returns:
        dict: Dictionary with motif_id as key and dict with consensus and e_value as values
    """
    with open(xstreme_path, 'r') as f:
        lines = f.readlines()
    
    motif_info = {}
    
    for i, line in enumerate(lines):
        # Look for lines with "MOTIF" and either "MEME-" or "STREME-"
        if "MOTIF" in line and ("MEME-" in line or "STREME-" in line):
            parts = line.split()
            
            # Extract consensus sequence (first element after MOTIF)
            if len(parts) >= 3 and parts[0] == "MOTIF":
                consensus = parts[1]
                motif_id = parts[2]
                
                # Look for E-value in the next few lines
                e_value = "N/A"
                for j in range(i+1, min(i+5, len(lines))):
                    if "E=" in lines[j]:
                        # Extract E-value
                        e_match = re.search(r'E=\s*([\d.e+-]+)', lines[j])
                        if e_match:
                            e_value = e_match.group(1)
                        break
                
                motif_info[motif_id] = {
                    'consensus': consensus,
                    'e_value': e_value
                }
    
    return motif_info


def extract_motif_blocks(meme_path):
    """
    Extract BL (BLOCKS format) sections from MEME output file.
    Each block contains sequences that have the motif.
    
    Args:
        meme_path: Path to meme.txt file
        
    Returns:
        dict: Dictionary with motif_id as key and list of sequence IDs as values
    """
    with open(meme_path, 'r') as f:
        lines = f.readlines()
    
    motif_blocks = {}
    current_motif = None
    in_block = False
    
    for i, line in enumerate(lines):
        # Look for the section header that contains "MEME-X in BLOCKS format"
        if "in BLOCKS format" in line:
            match = re.search(r'(MEME-\d+)', line)
            if match:
                current_motif = match.group(1)
                # Note: don't start collecting yet, wait for BL line
        
        # Detect start of actual BLOCKS data
        elif line.startswith("BL   MOTIF") and current_motif:
            motif_blocks[current_motif] = []
            in_block = True
        
        # End of block marked by "//"
        elif line.strip() == "//" and in_block:
            in_block = False
            current_motif = None
        
        # Extract sequence IDs within the block
        elif in_block and current_motif:
            # Skip the BL header line and separator lines
            if not line.startswith("BL") and not line.startswith("-"):
                # Parse sequence ID (first column before parenthesis)
                parts = line.split()
                if len(parts) > 0 and "(" in line:
                    seq_id = parts[0]
                    motif_blocks[current_motif].append(seq_id)
    
    return motif_blocks


def check_species_presence(sequence_ids, required_species):
    """
    Check if all required C4 species are represented in the sequence list.
    
    Args:
        sequence_ids: List of sequence identifiers
        required_species: List of species prefixes to check for
        
    Returns:
        tuple: (bool: all species present, dict: species presence)
    """
    species_present = {species: False for species in required_species}
    
    for seq_id in sequence_ids:
        for species in required_species:
            if seq_id.startswith(species):
                species_present[species] = True
    
    all_present = all(species_present.values())
    
    return all_present, species_present


def analyze_motif_conservation(xstreme_path, meme_path):
    """
    Main analysis function to identify conserved motifs across C4 species.
    
    Args:
        xstreme_path: Path to xstreme.txt file
        meme_path: Path to meme.txt file
        
    Returns:
        pd.DataFrame: Results showing motif conservation across species
    """
    print(f"Analyzing {len(C4_SPECIES)} C4 plant species: {', '.join(C4_SPECIES)}")
    print(f"Reading XSTREME output: {xstreme_path}")
    print(f"Reading MEME output: {meme_path}\n")
    
    # Extract motif information from XSTREME
    motif_info = parse_xstreme_motifs(xstreme_path)
    print(f"Found {len(motif_info)} motifs in XSTREME output:")
    for motif_id, info in motif_info.items():
        print(f"  - {motif_id}: {info['consensus']} (E={info['e_value']})")
    print()
    
    # Extract sequence information for each motif from MEME
    motif_blocks = extract_motif_blocks(meme_path)
    print(f"Extracted sequence information for {len(motif_blocks)} motifs from MEME output\n")
    
    # Analyze each motif for species conservation
    results = []
    
    for motif_id in motif_info.keys():
        if motif_id in motif_blocks:
            sequences = motif_blocks[motif_id]
            all_present, species_dict = check_species_presence(sequences, C4_SPECIES)
            
            # Create result row
            result = {
                'Motif_ID': motif_id,
                'Consensus_Sequence': motif_info[motif_id]['consensus'],
                'E_value': motif_info[motif_id]['e_value'],
                'Width': len(motif_info[motif_id]['consensus']),
                'All_C4_Species_Present': all_present,
                'Total_Sequences': len(sequences),
                'Species_Count': sum(species_dict.values())
            }
            
            # Add individual species presence
            for species in C4_SPECIES:
                result[f'{species}_present'] = species_dict[species]
            
            # Add list of sequences
            result['Sequence_IDs'] = '; '.join(sequences)
            
            results.append(result)
    
    # Create DataFrame
    if len(results) == 0:
        print("WARNING: No motifs found in MEME blocks!")
        return pd.DataFrame()
    
    df = pd.DataFrame(results)
    
    # Reorder columns for better readability
    base_cols = ['Motif_ID', 'Consensus_Sequence', 'E_value', 'Width', 
                 'All_C4_Species_Present', 'Species_Count', 'Total_Sequences']
    species_cols = [f'{species}_present' for species in C4_SPECIES]
    other_cols = ['Sequence_IDs']
    
    df = df[base_cols + species_cols + other_cols]
    
    # Sort by Species_Count (descending) and E_value (ascending)
    df['E_value_numeric'] = pd.to_numeric(df['E_value'], errors='coerce')
    df = df.sort_values(['Species_Count', 'E_value_numeric'], ascending=[False, True])
    df = df.drop('E_value_numeric', axis=1)
    
    return df


def main():
    """Main function to run the analysis with gene name as argument."""
    
    # Check command line arguments
    if len(sys.argv) < 2:
        print("ERROR: Gene name not provided!")
        print("\nUsage: python c4_motif_conservation_analysis.py <gene_name>")
        print("\nExample:")
        print("  python c4_motif_conservation_analysis.py GRMZM2G000171")
        sys.exit(1)
    
    # Get gene name from command line
    gene_name = sys.argv[1]
    
    # Construct file paths
    # CONFIGURE YOUR BASE PATH HERE:
    # Option 1: Standard upstream sequences
    # base_path = f"<BASE_DIR>/upstream_seqs/{gene_name}/{gene_name}_xst_out/"
    
    # Option 2: Control upstream sequences
    # base_path = f"<BASE_DIR>/upstream_seqs_control/{gene_name}/{gene_name}_OG_xst_out/"
    
    # Option 3: Without UTR upstream sequences (currently active)
    base_path = f"<BASE_DIR>/without_UTR/upstream_seqs/{gene_name}/{gene_name}_OG_xst_out/"
    
    xstreme_path = f"{base_path}/xstreme.txt"
    meme_path = f"{base_path}/meme_out/meme.txt"
    
    print("="*80)
    print(f"C4 PLANT MOTIF CONSERVATION ANALYSIS")
    print(f"Gene: {gene_name}")
    print("="*80)
    print()
    
    # Check if files exist
    if not Path(xstreme_path).exists():
        print(f"ERROR: XSTREME file not found: {xstreme_path}")
        sys.exit(1)
    
    if not Path(meme_path).exists():
        print(f"ERROR: MEME file not found: {meme_path}")
        sys.exit(1)
    
    # Run analysis
    results_df = analyze_motif_conservation(xstreme_path, meme_path)
    
    if len(results_df) == 0:
        print("No results to display. Check your input files.")
        sys.exit(1)
    
    # Display results
    print("="*80)
    print("MOTIF CONSERVATION ANALYSIS RESULTS")
    print("="*80)
    
    # Display simplified table (without sequence IDs for readability)
    display_cols = ['Motif_ID', 'Consensus_Sequence', 'E_value', 'Width', 
                    'All_C4_Species_Present', 'Species_Count', 'Total_Sequences'] + \
                   [f'{species}_present' for species in C4_SPECIES]
    print(results_df[display_cols].to_string(index=False))
    print()
    
    # Summary statistics
    conserved_motifs = results_df[results_df['All_C4_Species_Present'] == True]
    print("="*80)
    print("SUMMARY")
    print("="*80)
    print(f"Total motifs analyzed: {len(results_df)}")
    print(f"Motifs conserved in ALL {len(C4_SPECIES)} C4 species: {len(conserved_motifs)}")
    if len(results_df) > 0:
        print(f"Conservation rate: {len(conserved_motifs)/len(results_df)*100:.1f}%")
    
    # Show distribution by species count
    print("\nMotif distribution by number of species:")
    species_count_dist = results_df['Species_Count'].value_counts().sort_index(ascending=False)
    for count, num_motifs in species_count_dist.items():
        print(f"  {count}/{len(C4_SPECIES)} species: {num_motifs} motif(s)")
    print()
    
    if len(conserved_motifs) > 0:
        print("Conserved motif IDs (present in ALL species):")
        for _, row in conserved_motifs.iterrows():
            print(f"  ✓ {row['Motif_ID']:10s} {row['Consensus_Sequence']:20s} (E={row['E_value']})")
    else:
        print("No motifs found that are conserved across all C4 species.")
    print()
    
    # Save to CSV in the same directory as input files
    output_dir = base_path
    output_file = f"{output_dir}/{gene_name}_motif_conservation_results.csv"
    results_df.to_csv(output_file, index=False)
    print(f"Full results saved to: {output_file}")
    
    # Save conserved motifs only
    if len(conserved_motifs) > 0:
        conserved_file = f"{output_dir}/{gene_name}_conserved_motifs_only.csv"
        conserved_motifs.to_csv(conserved_file, index=False)
        print(f"Conserved motifs saved to: {conserved_file}")
    
    print()
    print("="*80)
    print("Analysis complete!")
    print("="*80)


if __name__ == "__main__":
    main()
