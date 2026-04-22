#!/usr/bin/env python2
# -*- coding: utf-8 -*-

# --- Dictionaries mapping species keys to genome and GFF files ---
key2genome = {
    'GRM': 'maize/Zmays_284_AGPv3.fa',
    'AC1': 'maize/Zmays_284_AGPv3.fa',
    'AC2': 'maize/Zmays_284_AGPv3.fa',
    'Bra': 'Bdistachyon/Bdistachyon_556_v3.0.fa',
    'Cha': 'Claxumv1.1/Claxum_690_v1.0.fa',
    'OEL': 'Doligosanthusv2/GCA_001633215.2_ASM163321v2_genomic.fna',
    'ELE': 'Ecoracanav1.1/Ecoracana_560_v1.0.fa',
    'LOC': 'rice/Osativa_204_v7.0.fa',
    'Pah': 'Phalliiv3.2/Phallii_590_v3.0.fa',
    'Pav': 'Pvaginatumv3.1/Pvaginatum_672_v3.0.fa',
    'Sev': 'Sviridisv4.1/Sviridis_726_v4.0.fa',
    'Uro': 'Ufuscav1.1/Ufusca_669_v1.0.fa'
}

key2gff = {
    'GRM': 'maize/Zmays_284_Ensembl-18_2010-01-MaizeSequence.gene.gff3',
    'AC1': 'maize/Zmays_284_Ensembl-18_2010-01-MaizeSequence.gene.gff3',
    'AC2': 'maize/Zmays_284_Ensembl-18_2010-01-MaizeSequence.gene.gff3',
    'Bra': 'Bdistachyon/Bdistachyon_556_v3.2.gene.gff3',
    'Cha': 'Claxumv1.1/Claxum_690_v1.1.gene.gff3',
    'OEL': 'Doligosanthusv2/genomic.gff',
    'ELE': 'Ecoracanav1.1/Ecoracana_560_v1.1.gene.gff3',
    'LOC': 'rice/Osativa_204_v7.0.gene.gff3',
    'Pah': 'Phalliiv3.2/Phallii_590_v3.2.gene.gff3',
    'Pav': 'Pvaginatumv3.1/Pvaginatum_672_v3.1.gene.gff3',
    'Sev': 'Sviridisv4.1/Sviridis_726_v4.1.gene.gff3',
    'Uro': 'Ufuscav1.1/Ufusca_669_v1.1.gene.gff3'
}

# --- Complement base dictionary (handles ambiguous bases too) ---
compBase = {
    'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N',
    'R': 'Y', 'Y': 'R', 'S': 'S', 'W': 'W',
    'K': 'M', 'M': 'K', 'B': 'V', 'D': 'H', 'H': 'D', 'V': 'B'
}

# --- Open log file for missing/failed IDs ---
missing_log = open("missing_IDs.log", "w")

# --- Process orthogroups ---
with open("ortholog_cleaned_positive.csv") as fh:
    for othgline in fh:
        if 'GRM' not in othgline:
            continue

        print "\n\n**Orthogroup**\n" + othgline
        othgline = othgline.replace('"', '')
        ids = othgline.strip().split(',')
        id_count = 0
        fout = None

        for i in ids:
            if i == "":
                continue

            # First ID: create new output file
            if id_count == 0:
                id_count += 1
                fout = open(i + '_OG.fasta', 'w')
                continue

            id_count += 1

            # --- Automatically detect species key ---
            key = None
            for k in key2gff:
                if i.startswith(k):
                    key = k
                    break

            if key is None:
                msg = "**Warning: No species key found for ID {}. Skipping.**".format(i)
                print msg
                missing_log.write(msg + "\n")
                continue

            print "**Processing ID:** {} → Species key: {}".format(i, key)

            # --- Search for gene in GFF ---
            gffh = open(key2gff[key])
            chrm = ""
            strand = ""
            start = 0
            end = 0

            for gffline in gffh:
                if (i in gffline) and ('CDS' in gffline):
                    ftrarr = gffline.split()
                    chrm = ftrarr[0]
                    start = int(ftrarr[3])
                    end = int(ftrarr[4])
                    strand = ftrarr[6]
                    break
            gffh.close()

            if (start == 0) or (end == 0):
                msg = "**Ortholog ID not found in GFF:** {}".format(i)
                print msg
                missing_log.write(msg + "\n")
                output = ">{}\tError! NO details found for this ID in GFF file\n".format(i)
                fout.write(output)
                continue

            # --- Extract sequence from FASTA ---
            fastah = open(key2genome[key])
            flag = 0
            loc = 0
            chrm_1 = chrm
            chrm = '>' + chrm

            if strand == '+':
                upstr_s = start - 1500
                upstr_e = start
            else:
                upstr_s = end
                upstr_e = end + 1500

            upstr = ""
            revcomp = ""

            for fasline in fastah:
                if flag == 1:
                    if (loc + len(fasline.strip())) >= upstr_e:
                        upstr += fasline.strip()
                        if strand == '+':
                            output = ">{}\t{}:{}-{}\t{}\n{}\n".format(i, chrm_1, start, end, strand, upstr)
                            fout.write(output)
                        else:
                            revcomp = ""
                            upstr = upstr.upper()
                            for idx in range(len(upstr) - 1, -1, -1):
                                base = upstr[idx].upper()
                                if base not in compBase:
                                    print "Warning: unexpected base '{}' in {}, replaced with 'N'".format(base, i)
                                revcomp += compBase.get(base, 'N')
                            output = ">{}\t{}:{}-{}\t{}\n{}\n".format(i, chrm_1, start, end, strand, revcomp)
                            fout.write(output)
                        flag = 0
                        loc = 0
                        upstr = ""
                        revcomp = ""
                        break

                    elif (loc + len(fasline.strip())) >= upstr_s:
                        upstr += fasline.strip()
                        loc += len(fasline.strip())

                    elif (loc + len(fasline.strip())) < upstr_s:
                        loc += len(fasline.strip())

                elif (flag == 0) and (chrm in fasline):
                    flag = 1

            fastah.close()

        if fout:
            fout.close()

# --- Close missing log file ---
missing_log.close()

print "\nDone. Any missing or skipped IDs have been logged in 'missing_IDs.log'."

