# Tissue-resolved single-cell expression of a query gene in NSCLC T cells

Reproducible analysis testing whether a **query gene** is upregulated in
**tumor-infiltrating** T cells versus **peripheral blood** and **adjacent normal**
tissue, using a single-cell T-cell atlas of non-small-cell lung cancer.

**Dataset:** [GSE99254](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE99254)
— Guo *et al.*, *Nature Medicine* (2018). FACS-sorted T cells profiled by SMART-seq2
from tumor / adjacent-normal / peripheral blood of 14 treatment-naive NSCLC patients
(12,346 cells; 12,210 with a parseable tissue+subtype label).

The processed TPM matrix is downloaded directly from GEO on first run and cached in `data/`.
Tissue and T-cell subtype are decoded from the cell-name prefix (e.g. `TTC5-0616A` =
**T**umor, **C**D8): `T/N/P` = Tumor / Normal / Peripheral; `C/H/R/Y/S` =
CD8 / CD4-helper / Treg / CD4-other / CD8-other.

## Query gene

The default query gene is **`PDCD1`** (PD-1) — a canonical T-cell activation/exhaustion
marker, which serves as a positive control for the pipeline. Change the `GENE`
parameter in the notebook config cell (or `--gene` on the CLI script) to analyze any gene.

## Contents

| File | Description |
|------|-------------|
| `PDCD1_GSE99254_analysis.ipynb` | End-to-end notebook (executed, with outputs) |
| `run_analysis.py` | Command-line equivalent: `python run_analysis.py --gene PDCD1` |
| `requirements.txt` | Python dependencies |
| `results/` | Figure (PNG) + summary and p-value CSVs |

## Reproduce

```bash
pip install -r requirements.txt

# Notebook:
jupyter notebook PDCD1_GSE99254_analysis.ipynb

# or command line (any gene):
python run_analysis.py --gene PDCD1
python run_analysis.py --gene UBASH3A
```

## Method

- Expression normalized as `log2(TPM/10 + 1)` (the convention used in the source study).
- Cell-level **Mann-Whitney U** (two-sided): Tumor vs Blood and Tumor vs Normal,
  computed overall and **within each T-cell subtype** (to rule out that a whole-tissue
  difference is driven by shifted subtype composition).
- Memory-efficient: only the single query-gene row is streamed from the ~330 MB matrix,
  so the analysis runs on a laptop.

**Caveat.** The numeric suffix in cell names is a within-sort-group serial, **not** a
stable patient ID recoverable from this file, so tests are cell-level rather than
patient-paired pseudobulk.

## Example result (PDCD1 / PD-1)

Mean `log2(TPM/10+1)` rises Blood **0.43** -> Normal **1.09** -> Tumor **1.96**, with
tumor enrichment significant in every subtype (CD8, CD4-helper, Treg, CD4-other) — the
expected exhaustion-marker biology.
