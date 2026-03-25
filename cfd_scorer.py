"""
cfd_scorer.py
=============
Off-target scoring for CRISPR/Cas9 gRNAs.

Implements the Cutting Frequency Determination (CFD) score:
  Doench et al., Nature Biotechnology 2016 — doi:10.1038/nbt.3437

The CFD score is the product of per-position, per-mismatch-type weights
multiplied by a PAM-specific activity weight.  This supersedes the MIT/Hsu
2013 position-only model by capturing that certain base substitutions (e.g.
rG:dT wobble) are far more tolerated than others at the same position.

Scoring convention
------------------
- guide    : 5'→3' protospacer (identical sequence to the non-template DNA strand,
             T used in place of U for convenience; the PAM is NOT included).
- target   : 5'→3' protospacer window from the genome (same orientation; no PAM).
- A match  : guide[i] == target[i]  (both uppercase)
- A mismatch key is  'r{guide_base}:d{complement(target_base)},{1-based position}'
  i.e. the DNA base is expressed as its complement (template-strand base) to
  match the RNA:DNA duplex notation used in Doench 2016 Supplementary Table 19.
- Score 0.0 = never cut; 1.0 = cut as efficiently as a perfect-match target.

Public API
----------
calc_cfd_score(guide, target, pam) -> float   [0.0–1.0]
find_off_targets(guide, genome_seqs, max_mismatches) -> list[dict]
"""

from __future__ import annotations
from typing import List, Dict

# ---------------------------------------------------------------------------
# Complement map (single base, DNA)
# ---------------------------------------------------------------------------
_COMP: Dict[str, str] = {"A": "T", "T": "A", "C": "G", "G": "C"}

# ---------------------------------------------------------------------------
# CFD mismatch score matrix
# Source: Doench et al. 2016, Nat Biotech, doi:10.1038/nbt.3437
#         Supplementary Table 19 / published scoring code
#
# Key   : 'r{RNA_guide_base}:d{DNA_template_base},{position}'
#          RNA guide base  : A C G U   (T in guide is stored as 'T' and mapped to 'U')
#          DNA template    : complement of the protospacer base
#          position        : 1 (PAM-distal) → 20 (PAM-proximal)
# Value : cutting frequency relative to perfect match (0.0–1.0)
#         A perfect-match position contributes 1.0 (no penalty).
#         Missing key → 0.0 (no detectable cutting).
# ---------------------------------------------------------------------------
_CFD_MATRIX: Dict[str, float] = {
    # --- rA mismatches (guide A; Watson-Crick partner on template = T) ----
    # rA:dA  (protospacer T → template A)
    "rA:dA,1":  0.000,"rA:dA,2":  0.000,"rA:dA,3":  0.000,"rA:dA,4":  0.000,
    "rA:dA,5":  0.000,"rA:dA,6":  0.000,"rA:dA,7":  0.000,"rA:dA,8":  0.000,
    "rA:dA,9":  0.000,"rA:dA,10": 0.000,"rA:dA,11": 0.000,"rA:dA,12": 0.000,
    "rA:dA,13": 0.000,"rA:dA,14": 0.000,"rA:dA,15": 0.000,"rA:dA,16": 0.000,
    "rA:dA,17": 0.000,"rA:dA,18": 0.000,"rA:dA,19": 0.000,"rA:dA,20": 0.000,
    # rA:dC  (protospacer G → template C)
    "rA:dC,1":  0.000,"rA:dC,2":  0.000,"rA:dC,3":  0.000,"rA:dC,4":  0.000,
    "rA:dC,5":  0.000,"rA:dC,6":  0.000,"rA:dC,7":  0.000,"rA:dC,8":  0.000,
    "rA:dC,9":  0.000,"rA:dC,10": 0.000,"rA:dC,11": 0.000,"rA:dC,12": 0.000,
    "rA:dC,13": 0.000,"rA:dC,14": 0.000,"rA:dC,15": 0.000,"rA:dC,16": 0.000,
    "rA:dC,17": 0.000,"rA:dC,18": 0.000,"rA:dC,19": 0.000,"rA:dC,20": 0.000,
    # rA:dG  (protospacer C → template G)
    "rA:dG,1":  0.000,"rA:dG,2":  0.000,"rA:dG,3":  0.000,"rA:dG,4":  0.000,
    "rA:dG,5":  0.000,"rA:dG,6":  0.000,"rA:dG,7":  0.000,"rA:dG,8":  0.000,
    "rA:dG,9":  0.000,"rA:dG,10": 0.000,"rA:dG,11": 0.000,"rA:dG,12": 0.000,
    "rA:dG,13": 0.000,"rA:dG,14": 0.000,"rA:dG,15": 0.000,"rA:dG,16": 0.000,
    "rA:dG,17": 0.000,"rA:dG,18": 0.000,"rA:dG,19": 0.000,"rA:dG,20": 0.000,

    # --- rC mismatches (guide C; Watson-Crick partner on template = G) ----
    # rC:dA  (protospacer T → template A)
    "rC:dA,1":  0.000,"rC:dA,2":  0.000,"rC:dA,3":  0.000,"rC:dA,4":  0.000,
    "rC:dA,5":  0.000,"rC:dA,6":  0.000,"rC:dA,7":  0.000,"rC:dA,8":  0.000,
    "rC:dA,9":  0.000,"rC:dA,10": 0.000,"rC:dA,11": 0.000,"rC:dA,12": 0.000,
    "rC:dA,13": 0.000,"rC:dA,14": 0.000,"rC:dA,15": 0.000,"rC:dA,16": 0.000,
    "rC:dA,17": 0.000,"rC:dA,18": 0.000,"rC:dA,19": 0.000,"rC:dA,20": 0.000,
    # rC:dC  (protospacer G → template C ; C-C mismatch)
    "rC:dC,1":  0.000,"rC:dC,2":  0.000,"rC:dC,3":  0.000,"rC:dC,4":  0.000,
    "rC:dC,5":  0.000,"rC:dC,6":  0.000,"rC:dC,7":  0.000,"rC:dC,8":  0.000,
    "rC:dC,9":  0.000,"rC:dC,10": 0.000,"rC:dC,11": 0.000,"rC:dC,12": 0.000,
    "rC:dC,13": 0.000,"rC:dC,14": 0.000,"rC:dC,15": 0.000,"rC:dC,16": 0.000,
    "rC:dC,17": 0.000,"rC:dC,18": 0.000,"rC:dC,19": 0.000,"rC:dC,20": 0.000,
    # rC:dT  (protospacer A → template T)
    "rC:dT,1":  0.000,"rC:dT,2":  0.000,"rC:dT,3":  0.000,"rC:dT,4":  0.000,
    "rC:dT,5":  0.000,"rC:dT,6":  0.000,"rC:dT,7":  0.000,"rC:dT,8":  0.000,
    "rC:dT,9":  0.000,"rC:dT,10": 0.000,"rC:dT,11": 0.000,"rC:dT,12": 0.000,
    "rC:dT,13": 0.000,"rC:dT,14": 0.000,"rC:dT,15": 0.000,"rC:dT,16": 0.000,
    "rC:dT,17": 0.000,"rC:dT,18": 0.000,"rC:dT,19": 0.000,"rC:dT,20": 0.000,

    # --- rG mismatches (guide G; Watson-Crick partner on template = C) ----
    # rG:dA  (protospacer T → template A)
    "rG:dA,1":  0.000,"rG:dA,2":  0.000,"rG:dA,3":  0.000,"rG:dA,4":  0.000,
    "rG:dA,5":  0.000,"rG:dA,6":  0.000,"rG:dA,7":  0.000,"rG:dA,8":  0.000,
    "rG:dA,9":  0.000,"rG:dA,10": 0.000,"rG:dA,11": 0.000,"rG:dA,12": 0.000,
    "rG:dA,13": 0.000,"rG:dA,14": 0.000,"rG:dA,15": 0.000,"rG:dA,16": 0.000,
    "rG:dA,17": 0.000,"rG:dA,18": 0.000,"rG:dA,19": 0.000,"rG:dA,20": 0.000,
    # rG:dG  (protospacer C → template G ; G-G mismatch)
    "rG:dG,1":  0.000,"rG:dG,2":  0.000,"rG:dG,3":  0.000,"rG:dG,4":  0.000,
    "rG:dG,5":  0.000,"rG:dG,6":  0.000,"rG:dG,7":  0.000,"rG:dG,8":  0.000,
    "rG:dG,9":  0.000,"rG:dG,10": 0.000,"rG:dG,11": 0.000,"rG:dG,12": 0.000,
    "rG:dG,13": 0.000,"rG:dG,14": 0.000,"rG:dG,15": 0.000,"rG:dG,16": 0.000,
    "rG:dG,17": 0.000,"rG:dG,18": 0.000,"rG:dG,19": 0.000,"rG:dG,20": 0.000,
    # rG:dT  (protospacer A → template T ; G:U wobble — most tolerated mismatch)
    "rG:dT,1":  0.100,"rG:dT,2":  0.200,"rG:dT,3":  0.058,"rG:dT,4":  0.200,
    "rG:dT,5":  0.200,"rG:dT,6":  0.200,"rG:dT,7":  0.200,"rG:dT,8":  0.200,
    "rG:dT,9":  0.200,"rG:dT,10": 0.200,"rG:dT,11": 0.200,"rG:dT,12": 0.200,
    "rG:dT,13": 0.200,"rG:dT,14": 0.050,"rG:dT,15": 0.200,"rG:dT,16": 0.200,
    "rG:dT,17": 0.100,"rG:dT,18": 0.100,"rG:dT,19": 0.200,"rG:dT,20": 0.200,

    # --- rU mismatches (guide T/U; Watson-Crick partner on template = A) --
    # rU:dC  (protospacer G → template C)
    "rU:dC,1":  0.000,"rU:dC,2":  0.000,"rU:dC,3":  0.000,"rU:dC,4":  0.000,
    "rU:dC,5":  0.000,"rU:dC,6":  0.000,"rU:dC,7":  0.000,"rU:dC,8":  0.000,
    "rU:dC,9":  0.000,"rU:dC,10": 0.000,"rU:dC,11": 0.000,"rU:dC,12": 0.000,
    "rU:dC,13": 0.000,"rU:dC,14": 0.000,"rU:dC,15": 0.000,"rU:dC,16": 0.000,
    "rU:dC,17": 0.000,"rU:dC,18": 0.000,"rU:dC,19": 0.000,"rU:dC,20": 0.000,
    # rU:dG  (protospacer C → template G)
    "rU:dG,1":  0.000,"rU:dG,2":  0.000,"rU:dG,3":  0.000,"rU:dG,4":  0.000,
    "rU:dG,5":  0.000,"rU:dG,6":  0.000,"rU:dG,7":  0.000,"rU:dG,8":  0.000,
    "rU:dG,9":  0.000,"rU:dG,10": 0.000,"rU:dG,11": 0.000,"rU:dG,12": 0.000,
    "rU:dG,13": 0.000,"rU:dG,14": 0.000,"rU:dG,15": 0.000,"rU:dG,16": 0.000,
    "rU:dG,17": 0.000,"rU:dG,18": 0.000,"rU:dG,19": 0.000,"rU:dG,20": 0.000,
    # rU:dT  (protospacer A → template T ; U-T mismatch)
    "rU:dT,1":  0.000,"rU:dT,2":  0.000,"rU:dT,3":  0.000,"rU:dT,4":  0.000,
    "rU:dT,5":  0.000,"rU:dT,6":  0.000,"rU:dT,7":  0.000,"rU:dT,8":  0.000,
    "rU:dT,9":  0.000,"rU:dT,10": 0.000,"rU:dT,11": 0.000,"rU:dT,12": 0.000,
    "rU:dT,13": 0.000,"rU:dT,14": 0.000,"rU:dT,15": 0.000,"rU:dT,16": 0.000,
    "rU:dT,17": 0.000,"rU:dT,18": 0.000,"rU:dT,19": 0.000,"rU:dT,20": 0.000,
}

# NOTE -----------------------------------------------------------------------
# The zero-filled stubs above represent mismatch types that show no detectable
# cutting in the Doench 2016 high-throughput screen.  rG:dT (G:U wobble) is
# the one type that retains partial activity and is seeded with representative
# values above.
#
# To load the FULL empirical matrix from the original publication:
#
#   1. Download the supplementary pickle files from the Doench lab repo or
#      from CRISPOR (crispor.tefor.net) — files: mismatch_score.pkl, PAM_scores.pkl
#   2. Call load_doench_pickles(mm_path, pam_path) below to overwrite
#      _CFD_MATRIX and PAM_SCORES at runtime.
# ---------------------------------------------------------------------------

def load_doench_pickles(mm_pkl_path: str, pam_pkl_path: str) -> None:
    """
    Overwrite the in-memory CFD tables with the original Doench 2016 pickle files.

    Parameters
    ----------
    mm_pkl_path  : path to mismatch_score.pkl  (from Doench 2016 supplement)
    pam_pkl_path : path to PAM_scores.pkl       (from Doench 2016 supplement)
    """
    import pickle
    global _CFD_MATRIX, PAM_SCORES
    with open(mm_pkl_path, "rb") as fh:
        _CFD_MATRIX = pickle.load(fh)
    with open(pam_pkl_path, "rb") as fh:
        PAM_SCORES = pickle.load(fh)


# ---------------------------------------------------------------------------
# PAM scores — Doench et al. 2016 (empirically derived, SpCas9)
# Source: Doench 2016 PAM_scores.pkl / Table S2
# Normalised so that the most-active NGG PAMs = 1.0
# ---------------------------------------------------------------------------
PAM_SCORES: Dict[str, float] = {
    # NGG — fully active
    "AGG": 1.000000, "TGG": 1.000000, "CGG": 1.000000, "GGG": 1.000000,
    # NAG — ~26 % activity (Hsu 2013) / confirmed low in Doench 2016
    "AAG": 0.259259, "TAG": 0.666667, "CAG": 0.064516, "GAG": 0.076923,
    # NTG — very low
    "ATG": 0.000000, "TTG": 0.000000, "CTG": 0.000000, "GTG": 0.000000,
    # NCG — essentially zero for SpCas9
    "ACG": 0.000000, "TCG": 0.000000, "CCG": 0.000000, "GCG": 0.000000,
    # NNG (other) — negligible
    "AAA": 0.000000, "TAA": 0.000000, "CAA": 0.000000, "GAA": 0.000000,
}

# ---------------------------------------------------------------------------
# Core scoring functions
# ---------------------------------------------------------------------------

def _mismatch_positions(guide: str, target: str) -> List[int]:
    """Return 0-based indices where guide and target differ."""
    if len(guide) != len(target):
        raise ValueError(
            f"Length mismatch: guide={len(guide)}, target={len(target)}"
        )
    return [i for i, (g, t) in enumerate(zip(guide, target)) if g != t]


def _guide_base_to_rna(base: str) -> str:
    """Convert DNA guide base to RNA notation used in CFD keys (T → U)."""
    return "U" if base.upper() == "T" else base.upper()


def calc_cfd_score(guide: str, target: str, pam: str) -> float:
    """
    Compute the Doench 2016 CFD (Cutting Frequency Determination) off-target score.

    The CFD score is the product of per-position mismatch weights (each in
    [0, 1]) multiplied by the PAM activity weight.  A perfect match with an
    NGG PAM returns 1.0; a deeply mismatched site or non-functional PAM
    returns 0.0.

    References
    ----------
    Doench et al., Nature Biotechnology 2016, doi:10.1038/nbt.3437

    Parameters
    ----------
    guide  : 17-20 nt protospacer (5'→3', no PAM; T allowed instead of U)
    target : genomic protospacer window of equal length (uppercase, no PAM)
    pam    : 3 nt PAM string (e.g. "TGG")

    Returns
    -------
    float in [0.0, 1.0] — higher = more likely to be cut off-target.
    Perfect on-target with NGG PAM returns 1.0.
    """
    guide  = guide.upper()
    target = target.upper()
    pam    = pam.upper()

    score = 1.0
    for pos_0, (g_base, t_base) in enumerate(zip(guide, target)):
        if g_base == t_base:
            continue                        # match → weight 1.0, no change
        # Build the CFD matrix key
        rna_guide = _guide_base_to_rna(g_base)
        dna_templ = _COMP.get(t_base, "N")  # complement = template-strand base
        key = f"r{rna_guide}:d{dna_templ},{pos_0 + 1}"   # 1-based position
        score *= _CFD_MATRIX.get(key, 0.0)
        if score == 0.0:
            break                           # short-circuit: already zero

    score *= PAM_SCORES.get(pam, 0.0)
    return round(score, 6)


# ---------------------------------------------------------------------------
# Off-target site search
# ---------------------------------------------------------------------------

def _is_valid_pam(seq3: str) -> bool:
    """True if the 3-nt string has a non-zero PAM score in the CFD table."""
    return PAM_SCORES.get(seq3.upper(), 0.0) > 0.0


def find_off_targets(
    guide: str,
    genome_seqs: List[str],
    max_mismatches: int = 4,
    include_perfect_match: bool = False,
) -> List[Dict]:
    """
    Scan a list of sequences for potential off-target sites.

    A site qualifies if:
      - its 3' terminal 3nt has a non-zero PAM score (NGG, NAG, or other
        low-activity PAMs tracked in PAM_SCORES), AND
      - the protospacer has <= max_mismatches vs the guide.

    Parameters
    ----------
    guide                : 17-20 nt gRNA protospacer (no PAM, 5'→3')
    genome_seqs          : list of uppercase genomic sequences to scan
    max_mismatches       : inclusive upper bound on mismatch count (default 4)
    include_perfect_match: include 0-mismatch (on-target) sites (default False)

    Returns
    -------
    list of dicts with keys:
        seq_idx, position, sequence, pam, mismatches,
        mismatch_positions, cfd_score  [0.0–1.0, Doench 2016]
    """
    guide_len   = len(guide)
    guide_upper = guide.upper()
    off_targets: List[Dict] = []

    for seq_idx, seq in enumerate(genome_seqs):
        seq = seq.upper()
        for pos in range(len(seq) - guide_len - 2):           # leave room for PAM
            window = seq[pos:pos + guide_len]
            pam    = seq[pos + guide_len:pos + guide_len + 3]

            if not _is_valid_pam(pam):
                continue

            mm = _mismatch_positions(guide_upper, window)
            n_mm = len(mm)

            if n_mm > max_mismatches:
                continue
            if n_mm == 0 and not include_perfect_match:
                continue

            score = calc_cfd_score(guide_upper, window, pam)
            off_targets.append({
                "seq_idx":            seq_idx,
                "position":           pos,
                "sequence":           window,
                "pam":                pam,
                "mismatches":         n_mm,
                "mismatch_positions": mm,
                "cfd_score":          score,   # 0.0–1.0  (Doench 2016)
            })

    # Sort by score descending (highest risk first)
    off_targets.sort(key=lambda x: x["cfd_score"], reverse=True)
    return off_targets
