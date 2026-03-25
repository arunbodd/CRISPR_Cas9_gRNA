# CRISPR Cas9 gRNA Design & Off-Target Scoring Pipeline

A Python 3 pipeline for designing SpCas9 guide RNAs (gRNAs) from human genome build 38 (hg38) exonic sequences and scoring predicted off-target sites using the Cutting Frequency Determination (CFD) model.

---

## Overview

This tool identifies candidate 20 nt gRNA sequences from hg38 exonic FASTA input by applying PAM compatibility and seed-region GC content filters tuned for SpCas9 specificity. Each candidate is then evaluated for off-target activity across the input sequences using the empirically derived CFD scoring model (Doench et al. 2016), which accounts for both mismatch type and position within the protospacer.

---

## Features

- **gRNA Candidate Design** — sliding 20 nt window scan with NGG/NAG PAM filtering and ≥60% GC content in the seed region (positions 10–17)
- **Off-Target Site Search** — genome-wide scan for protospacer windows with ≤4 mismatches and a valid PAM
- **CFD Off-Target Scoring** — Doench 2016 position- and mismatch-type-weighted scoring; output range 0.0–1.0 (1.0 = on-target efficiency)
- **Expanded PAM Table** — NGG (full activity), NAG, NTG, and NCG PAM weights from Doench 2016 Table S2
- **TSV Output** — scored results with gene, guide sequence, GC%, off-target count, and minimum CFD score
- **CLI Interface** — `--input`, `--output`, and `--max-mm` flags

---

## Requirements

- Python ≥ 3.8
- No third-party dependencies (standard library only)

---

## Usage

```bash
python Guide_Design_tool.py \
    --input  hg38_exons.fasta \
    --output results.txt \
    --max-mm 4
```

| Argument    | Default                    | Description                                      |
|-------------|----------------------------|--------------------------------------------------|
| `--input`   | `hg38_new.fasta`           | Input multi-sequence FASTA (hg38 exonic regions) |
| `--output`  | `hg38_matched_guideRNA.txt`| Output TSV file                                  |
| `--max-mm`  | `4`                        | Maximum mismatches for off-target search         |

### Output columns

| Column          | Description                                          |
|-----------------|------------------------------------------------------|
| Header          | FASTA header (Ensembl gene/exon IDs, coordinates)    |
| Gene            | Gene symbol                                          |
| Guide_Sequence  | 20 nt gRNA protospacer                               |
| Pct_C           | % cytosine in seed region (pos 10–17)                |
| Pct_G           | % guanine in seed region (pos 10–17)                 |
| Num_OffTargets  | Number of off-target sites found (≤ max-mm)          |
| Min_CFD_Score   | Lowest CFD score among off-target sites (0.0–1.0)    |

---

## Off-Target Scoring Model

Off-target scores are computed using the **CFD (Cutting Frequency Determination)** model:

$$\text{CFD} = \prod_{i=1}^{20} W_{r,d,i} \times W_{\text{PAM}}$$

where $W_{r,d,i}$ is the empirical weight for RNA base $r$ mismatched against DNA template base $d$ at position $i$, and $W_{\text{PAM}}$ is the PAM activity weight. A score of **1.0** indicates the site is predicted to be cleaved as efficiently as the on-target; **0.0** indicates no predicted cleavage.

**Key properties of the model:**
- The rG:dT (G:U wobble) mismatch is the most tolerated substitution and retains partial activity at all positions
- PAM-proximal mismatches (positions 17–20) incur heavier penalties than PAM-distal mismatches
- NGG PAM weight = 1.0; NAG PAM weights vary by trinucleotide context (0.06–0.67)

To replace the built-in matrix with the original Broad Institute empirical data:
```python
from cfd_scorer import load_doench_pickles
load_doench_pickles("mismatch_score.pkl", "PAM_scores.pkl")
```

---

## Pipeline Architecture

The pipeline architecture is documented in `architecture.d2`. Render with:

```bash
brew install d2
d2 architecture.d2 architecture.svg
```

The diagram shows the data flow from FASTA input → gRNA candidate design → off-target scoring → scored TSV output, with a dashed future node for an SML-based predictor (DNABERT-2 fine-tuned on GUIDE-seq/CIRCLE-seq data).

---

## References

1. **Doench, J.G. et al.** (2016). Optimized sgRNA design to maximize activity and minimize off-target effects of CRISPR-Cas9. *Nature Biotechnology*, 34, 184–191. https://doi.org/10.1038/nbt.3437

2. **Hsu, P.D. et al.** (2013). DNA targeting specificity of RNA-guided Cas9 nucleases. *Nature Biotechnology*, 31, 827–832. https://doi.org/10.1038/nbt.2647

3. **Ran, F.A. et al.** (2013). Genome engineering using the CRISPR-Cas9 system. *Nature Protocols*, 8, 2281–2308. https://doi.org/10.1038/nprot.2013.143

4. **Ensembl** (hg38 / GRCh38). https://www.ensembl.org

---

## License

MIT License

Copyright (c) 2026 Arun Boddapati

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
