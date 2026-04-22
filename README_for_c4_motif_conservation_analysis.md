# C4 Motif Conservation Analysis

This script extracts and analyzes motifs conserved across C4 plant species from XSTREME/MEME output files.

## C4 Plant Species Analyzed

- **GRMZM**: Zea mays (Maize)
- **Pahal**: Panicum hallii (Hall's panicgrass)
- **Sevir**: Setaria viridis (Green foxtail)
- **Urofu**: Urochloa fusca
- **ELECO**: Eleusine coracana (Finger millet)
- **Pavag**: Paspalum vaginatum (Seashore paspalum)

## Configuration

Before running the script, you need to configure the base directory path in the `main()` function.

**Line ~260-270** in `c4_motif_conservation_analysis.py`:

Replace `<BASE_DIR>` with your actual data directory path. Three options are provided:

```python
# Option 1: Standard upstream sequences
# base_path = f"<BASE_DIR>/upstream_seqs/{gene_name}/{gene_name}_xst_out/"

# Option 2: Control upstream sequences
# base_path = f"<BASE_DIR>/upstream_seqs_control/{gene_name}/{gene_name}_OG_xst_out/"

# Option 3: Without UTR upstream sequences (currently active)
base_path = f"<BASE_DIR>/without_UTR/upstream_seqs/{gene_name}/{gene_name}_OG_xst_out/"
```

**Example configuration:**
```python
base_path = f"/home/user/data/c4_motif/without_UTR/upstream_seqs/{gene_name}/{gene_name}_OG_xst_out/"
```

## Expected Directory Structure

The script expects the following directory structure:

```
<BASE_DIR>/
├── upstream_seqs/               # Option 1
│   └── <gene_name>/
│       └── <gene_name>_xst_out/
│           ├── xstreme.txt
│           └── meme_out/
│               └── meme.txt
├── upstream_seqs_control/       # Option 2
│   └── <gene_name>/
│       └── <gene_name>_OG_xst_out/
│           ├── xstreme.txt
│           └── meme_out/
│               └── meme.txt
└── without_UTR/                 # Option 3
    └── upstream_seqs/
        └── <gene_name>/
            └── <gene_name>_OG_xst_out/
                ├── xstreme.txt
                └── meme_out/
                    └── meme.txt
```

## Usage

```bash
python c4_motif_conservation_analysis.py <gene_name>
```

### Example

```bash
python c4_motif_conservation_analysis.py GRMZM2G000171
```

## Input Files

The script requires two input files from XSTREME/MEME analysis:

1. **xstreme.txt**: XSTREME output containing motif consensus sequences and E-values
2. **meme.txt**: MEME output in BLOCKS format containing sequence information

## Output Files

The script generates two CSV files in the same directory as the input files:

1. **`<gene_name>_motif_conservation_results.csv`**: Complete results for all motifs
2. **`<gene_name>_conserved_motifs_only.csv`**: Results filtered to motifs conserved across ALL C4 species

## Output Columns

- **Motif_ID**: Identifier for the motif (e.g., MEME-1, STREME-1)
- **Consensus_Sequence**: Consensus sequence of the motif
- **E_value**: Statistical significance of the motif
- **Width**: Length of the consensus sequence
- **All_C4_Species_Present**: Boolean indicating if motif is present in all 6 C4 species
- **Species_Count**: Number of C4 species containing the motif
- **Total_Sequences**: Total number of sequences containing the motif
- **<Species>_present**: Boolean columns for each C4 species (GRMZM, Pahal, Sevir, Urofu, ELECO, Pavag)
- **Sequence_IDs**: Semicolon-separated list of all sequence IDs containing the motif

## Requirements

- Python 3.6+
- pandas

Install requirements:
```bash
pip install pandas
```

## License

[Add your license here]

## Citation

[Add citation information if applicable]
