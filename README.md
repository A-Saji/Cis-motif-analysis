# Cis-Motif Analysis Pipeline (C4 vs C3 Grasses)

This repository contains scripts used to identify, extract, classify, and analyze cis-regulatory motifs associated with C4 evolution across grass species.

The workflow integrates orthogroup-based upstream extraction, MEME/XSTREME motif discovery, FIMO scanning, and comparative motif conservation + positional analysis.


**Workflow overview**

Orthogroups
   ↓
Upstream sequence extraction
   ↓
FASTA splitting (C3 vs C4)
   ↓
Motif discovery (MEME / XSTREME)
   ↓
Motif extraction & filtering
   ↓
FIMO scanning (C3 background)
   ↓
Motif classification (presence + synteny)
   ↓
Visualization & TF annotation

**Dependencies**
Python:
Python 3 (some scripts use Python 2 ⚠️)
pandas
pathlib
re

R:
tidyverse
ggplot2
dplyr
ggpattern
tidygraph
ggraph

External Tools:
MEME Suite (meme, xstreme, fimo)
Unix tools (awk, grep, bash)

## Pipeline Steps

### 1. Extract Upstream Sequences
Extract upstream regions using genome FASTA + GFF files.

python extract_upstream_of_orthologs.py
Read the log files to identify any sequences for which upstream could not be found and add them in manually

## 2. Split FASTA into C3 and C4

python split_fasta.py
Outputs:

C3.fasta
C4.fasta

## 3. Motif Discovery (MEME/XSTREME)
c4_motif_conservation_analysis.py
Run using MEME Suite.

## 4. Extract Motif Positions
bash xstremeextraction.sh

## 5. Scan C3 Promoters (FIMO)
bash scan_motifs.sh

## 6. Extract FIMO Positions
bash c3location.sh

## 7. Filter and Merge Selected Motifs
python finalmeme.py
python extractmeme.py

Add in the header from a meme file for background frequencies in extractmeme.py output

## 8. Motif Classification (Presence + Synteny)

Run R scripts:

source("motifsynternyreal.R")

## 9. Visualization & TF Annotation
source("motifplots.R")
source("motifssyntenyclusteringannotation.R")

the annotation files comes from running tomtom using meme file against jaspar database.
