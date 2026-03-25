#!/usr/bin/env python3
"""
CRISPR Cas9 gRNA Design Tool  (Python 3)
=========================================
Identifies gRNA candidates from hg38 exonic FASTA sequences and reports
off-target sites scored with the MIT/Hsu position-weighted score
(Hsu et al., Nature Biotechnology 2013).

Usage:
    python Guide_Design_tool.py [--input FILE] [--output FILE] [--max-mm N]
"""

import re
import argparse
from cfd_scorer import calc_cfd_score, find_off_targets

# --- Configuration -----------------------------------------------------------

PAM_PATTERNS = {"AAG", "GGG", "CGG", "TGG", "GAG", "TAG", "CAG"}
GC_THRESHOLD = 60   # minimum % G or C in seed region to keep candidate
SEED_START   = 10   # 0-based start of seed window within the 20nt gRNA
SEED_END     = 17   # 0-based end (exclusive)

# --- FASTA parser ------------------------------------------------------------

def parse_fasta(filepath):
    """
    Load a FASTA file into a list of (header_str, gene_name, sequence) tuples.
    Handles multi-line sequences.
    """
    records = []
    current_header = current_gene = None
    current_seq = []

    with open(filepath, "r") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_header is not None:
                    records.append((current_header, current_gene, "".join(current_seq)))
                raw = line.lstrip(">")
                parts = raw.split("|")
                current_gene = parts[1] if len(parts) > 1 else parts[0]
                current_header = raw
                current_seq = []
            else:
                current_seq.append(line.upper())

    if current_header is not None:
        records.append((current_header, current_gene, "".join(current_seq)))

    return records

# --- gRNA candidate design ---------------------------------------------------

def find_grna_candidates(records):
    """
    Slide a 20nt window across each sequence.
    Keep windows whose last 3nt match a PAM pattern AND whose seed region
    has >= GC_THRESHOLD % G or >= GC_THRESHOLD % C.
    """
    candidates = []

    for header, gene_name, seq in records:
        for pos in range(len(seq) - 20):
            window = seq[pos:pos + 20]
            if window[-3:] not in PAM_PATTERNS:
                continue

            seed    = window[SEED_START:SEED_END]
            seq_len = len(seed)
            pct_g   = (seed.count("G") / seq_len) * 100
            pct_c   = (seed.count("C") / seq_len) * 100

            if int(pct_g) >= GC_THRESHOLD or int(pct_c) >= GC_THRESHOLD:
                candidates.append({
                    "header":   header,
                    "gene":     gene_name,
                    "guide":    window,
                    "position": pos,
                    "pct_c":    int(pct_c),
                    "pct_g":    int(pct_g),
                })

    return candidates

# --- Main --------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="CRISPR Cas9 gRNA Design Tool")
    parser.add_argument("--input",  default="hg38_new.fasta",
                        help="Input FASTA file (default: hg38_new.fasta)")
    parser.add_argument("--output", default="hg38_matched_guideRNA.txt",
                        help="Output TSV file (default: hg38_matched_guideRNA.txt)")
    parser.add_argument("--max-mm", type=int, default=4,
                        help="Max mismatches for off-target search (default: 4)")
    args = parser.parse_args()

    print(f"[1/3] Loading genome from {args.input} ...")
    records = parse_fasta(args.input)
    print(f"      Loaded {len(records)} sequence(s)")

    print("[2/3] Finding gRNA candidates ...")
    candidates = find_grna_candidates(records)
    print(f"      Found {len(candidates)} candidate(s)")

    genome_seqs = [seq for _, _, seq in records]

    print(f"[3/3] Scoring off-targets (max mismatches={args.max_mm}) ...")
    header_line = (
        "Header\tGene\tGuide_Sequence\tPct_C\tPct_G\t"
        "Num_OffTargets\tMin_CFD_Score\n"
    )

    with open(args.output, "w") as fout:
        fout.write(header_line)

        for cand in candidates:
            off_targets = find_off_targets(
                cand["guide"], genome_seqs, max_mismatches=args.max_mm
            )
            cfd_scores = [
                calc_cfd_score(cand["guide"], ot["sequence"], ot["pam"])
                for ot in off_targets
            ]
            min_cfd       = round(min(cfd_scores), 6) if cfd_scores else "NA"
            n_off_targets = len(off_targets)

            row = (
                f"{cand['header']}\t{cand['gene']}\t{cand['guide']}\t"
                f"{cand['pct_c']}\t{cand['pct_g']}\t"
                f"{n_off_targets}\t{min_cfd}\n"
            )
            fout.write(row)
            print(row.strip())

    print(f"\nResults written to {args.output}")


if __name__ == "__main__":
    main()
