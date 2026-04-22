#!/usr/bin/env python3
"""
Split FASTA files containing C3 and C4 plant sequences into separate files.
Organizes output into folders by gene ID.
"""

import os
import glob
from pathlib import Path

# Dictionary mapping gene prefixes to C3/C4 categories
plant_categories = {
    'Bradi': 'C3',      # Brachypodium distachyon
    'Chala': 'C3',      # Chlamydomonas (algae)
    'ELECO': 'C4',      # Eleusine coracana (finger millet)
    'LOC_Os': 'C3',     # Oryza sativa (rice)
    'Pahal': 'C4',      # Paspalum haumanii
    'Pavag': 'C4',      # Panicum virgatum (switchgrass)
    'Sevir': 'C4',      # Setaria viridis
    'Urofu': 'C4',      # Urochloa
    'GRMZM': 'C4',      # Zea mays (maize)
    'OEL': 'C3'         # Dichanthelium oligosanthes
}

def get_plant_category(header):
    """Determine if a sequence is from a C3 or C4 plant based on the header."""
    for prefix, category in plant_categories.items():
        if header.startswith(prefix):
            return category
    return None

def parse_fasta(filename):
    """Parse a FASTA file and return sequences grouped by C3/C4 category."""
    sequences = {'C3': [], 'C4': []}

    with open(filename, 'r') as f:
        current_header = None
        current_seq = []

        for line in f:
            line = line.strip()
            if not line:
                continue

            if line.startswith('>'):
                # Save previous sequence if exists
                if current_header is not None:
                    category = get_plant_category(current_header)
                    if category:
                        sequences[category].append({
                            'header': current_header,
                            'sequence': ''.join(current_seq)
                        })
                    else:
                        print(f"Warning: Unknown plant prefix in header: {current_header}")

                # Start new sequence
                current_header = line[1:]  # Remove '>'
                current_seq = []
            else:
                current_seq.append(line)

        # Don't forget the last sequence
        if current_header is not None:
            category = get_plant_category(current_header)
            if category:
                sequences[category].append({
                    'header': current_header,
                    'sequence': ''.join(current_seq)
                })
            else:
                print(f"Warning: Unknown plant prefix in header: {current_header}")

    return sequences

def write_fasta(sequences, filename):
    """Write sequences to a FASTA file."""
    with open(filename, 'w') as f:
        for seq in sequences:
            f.write(f">{seq['header']}\n")
            f.write(f"{seq['sequence']}\n")

def process_file(fasta_file):
    """Process a single FASTA file and split it into C3/C4 files."""
    # Extract gene ID from filename (e.g., GRMZM2G000743 from GRMZM2G000743_OG.fasta)
    filename = os.path.basename(fasta_file)
    gene_id = filename.replace('_OG.fasta', '')

    print(f"Processing {gene_id}...")

    # Parse the FASTA file
    sequences = parse_fasta(fasta_file)

    # Create directory for this gene
    gene_dir = Path(gene_id)
    gene_dir.mkdir(exist_ok=True)

    # Write C3 and C4 sequences to separate files
    c3_count = len(sequences['C3'])
    c4_count = len(sequences['C4'])

    if c3_count > 0:
        write_fasta(sequences['C3'], gene_dir / 'C3.fasta')
        print(f"  - Wrote {c3_count} C3 sequences to {gene_id}/C3.fasta")

    if c4_count > 0:
        write_fasta(sequences['C4'], gene_dir / 'C4.fasta')
        print(f"  - Wrote {c4_count} C4 sequences to {gene_id}/C4.fasta")

    return c3_count, c4_count

def main():
    """Main function to process all FASTA files."""
    # Find all FASTA files
    fasta_files = glob.glob('*_OG.fasta')

    if not fasta_files:
        print("No FASTA files found matching pattern '*_OG.fasta'")
        return

    print(f"Found {len(fasta_files)} FASTA files to process\n")

    total_c3 = 0
    total_c4 = 0

    for fasta_file in sorted(fasta_files):
        c3_count, c4_count = process_file(fasta_file)
        total_c3 += c3_count
        total_c4 += c4_count

    print(f"\n{'='*60}")
    print(f"Processing complete!")
    print(f"Total C3 sequences: {total_c3}")
    print(f"Total C4 sequences: {total_c4}")
    print(f"Created {len(fasta_files)} gene directories")
    print(f"{'='*60}")

if __name__ == '__main__':
    main()
