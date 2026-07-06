#!/usr/bin/env python
"""
Tissue-resolved single-cell expression of a query gene in NSCLC T cells (GSE99254).

Tests whether a query gene is upregulated in tumor-infiltrating T cells vs
peripheral blood and adjacent normal tissue. Default gene: PDCD1 (PD-1).

Usage:
    python run_analysis.py --gene PDCD1
    python run_analysis.py --gene UBASH3A --outdir results
"""
import argparse, os, re, gzip, urllib.request
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu

GEO_ACCESSION = "GSE99254"
TPM_URL = ("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE99nnn/"
           "GSE99254/suppl/GSE99254_NSCLC.TCell.S12346.TPM.txt.gz")
TISSUE = {"N": "Normal", "P": "Peripheral", "T": "Tumor"}
SUBTYPE = {"C": "CD8", "H": "CD4_Th", "R": "Treg", "Y": "CD4_other", "S": "CD8_other"}
ORDER = ["Peripheral", "Normal", "Tumor"]
LABEL = {"Peripheral": "Blood", "Normal": "Normal", "Tumor": "Tumor"}
SUBS = ["CD8", "CD4_Th", "Treg", "CD4_other"]
SUBLAB = {"CD8": "CD8", "CD4_Th": "CD4 helper", "Treg": "Treg", "CD4_other": "CD4 other"}


def download(data_dir):
    os.makedirs(data_dir, exist_ok=True)
    path = os.path.join(data_dir, "GSE99254_TPM.txt.gz")
    if not os.path.exists(path):
        print("Downloading TPM matrix (~330 MB) ...")
        urllib.request.urlretrieve(TPM_URL, path)
    return path


def load_gene(path, gene):
    """Stream only the header and the single query-gene row (low memory)."""
    with gzip.open(path, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            a = line.find("\t"); b = line.find("\t", a + 1)
            if line[a + 1:b] == gene:
                row = line.rstrip("\n").split("\t")
                cells = header[2:]
                return cells, pd.Series(np.array(row[2:], dtype=float), index=cells), row[0]
    raise SystemExit(f"Gene {gene!r} not found in matrix")


def parse(cell):
    m = re.match(r"^([NPT])T([CHRYS])", cell)
    return (TISSUE[m.group(1)], SUBTYPE[m.group(2)]) if m else (None, None)


def mw(a, b):
    return mannwhitneyu(a, b, alternative="two-sided")[1] if len(a) >= 3 and len(b) >= 3 else np.nan


def stars(p):
    if p is None or np.isnan(p): return "n/a"
    return "***" if p < 1e-3 else "**" if p < 1e-2 else "*" if p < 5e-2 else "ns"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gene", default="PDCD1")
    ap.add_argument("--data-dir", default="data")
    ap.add_argument("--outdir", default="results")
    args = ap.parse_args()
    gene = args.gene
    os.makedirs(args.outdir, exist_ok=True)

    path = download(args.data_dir)
    cells, expr, gid = load_gene(path, gene)
    print(f"{gene} found (geneID={gid}) across {len(cells)} cells")

    meta = pd.DataFrame({"cell": cells})
    meta[["tissue", "subtype"]] = pd.DataFrame([parse(c) for c in cells], index=meta.index)
    meta["TPM"] = expr.values
    meta["log2"] = np.log2(meta["TPM"] / 10 + 1)
    m = meta.dropna(subset=["tissue"]).copy()
    print(f"labeled T cells: {len(m)} / {len(meta)}")

    # summary + p-value tables
    srows = []
    for comp, sel in [("All T cells", m)] + [(SUBLAB[s], m[m.subtype == s]) for s in SUBS]:
        for t in ORDER:
            v = sel[sel.tissue == t]
            srows.append([comp, LABEL[t], len(v), round(v.log2.mean(), 3),
                          round(v.log2.median(), 3), round((v.TPM > 0).mean() * 100, 1)])
    summary = pd.DataFrame(srows, columns=["Compartment", "Tissue", "n_cells",
                                           "mean_log2", "median_log2", "pct_expressing"])
    prows = []
    for comp, sel in [("All T cells", m)] + [(SUBLAB[s], m[m.subtype == s]) for s in SUBS]:
        T = sel[sel.tissue == "Tumor"].log2; P = sel[sel.tissue == "Peripheral"].log2; N = sel[sel.tissue == "Normal"].log2
        prows.append([comp, mw(T, P), mw(T, N)])
    pvals = pd.DataFrame(prows, columns=["Compartment", "p_Tumor_vs_Blood", "p_Tumor_vs_Normal"])
    summary.to_csv(f"{args.outdir}/{gene}_{GEO_ACCESSION}_summary.csv", index=False)
    pvals.to_csv(f"{args.outdir}/{gene}_{GEO_ACCESSION}_pvalues.csv", index=False)
    print(summary.to_string(index=False))
    print(pvals.to_string(index=False))

    # figure
    mpl.rcParams.update({"savefig.dpi": 300, "savefig.bbox": "tight", "font.size": 8,
                         "axes.titlesize": 8.5, "axes.spines.top": False,
                         "axes.spines.right": False, "axes.titlelocation": "left"})
    colors = {"Peripheral": "#7fb3d5", "Normal": "#95a5a6", "Tumor": "#c0392b"}
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.0), gridspec_kw={"width_ratios": [1, 1.55]})
    ax = axes[0]
    data = [m[m.tissue == t].log2.values for t in ORDER]
    parts = ax.violinplot(data, positions=range(3), showextrema=False, widths=0.85)
    for i, pc in enumerate(parts["bodies"]):
        pc.set_facecolor(colors[ORDER[i]]); pc.set_alpha(0.55); pc.set_edgecolor("none")
    ymax = max(v.max() for v in data)
    for i, t in enumerate(ORDER):
        v = m[m.tissue == t].log2; q1, med, q3 = v.quantile([.25, .5, .75])
        ax.plot([i - .18, i + .18], [med, med], color="k", lw=1.8, zorder=5)
        ax.plot([i, i], [q1, q3], color="k", lw=1.0, zorder=4)
        ax.text(i, ymax * 0.94, f"{v.mean():.2f}", ha="center", va="bottom", fontsize=6, color=colors[t])
    ax.set_xticks(range(3)); ax.set_xticklabels([LABEL[t] for t in ORDER])
    ax.set_ylabel(f"{gene}  log$_2$(TPM/10+1)"); ax.set_ylim(-0.3, ymax * 1.05)
    ax.set_title(f"{gene} across tissue compartments")
    pAP = mw(m[m.tissue=="Tumor"].log2, m[m.tissue=="Peripheral"].log2)
    pAN = mw(m[m.tissue=="Tumor"].log2, m[m.tissue=="Normal"].log2)
    def brack(x1, x2, y, txt):
        ax.plot([x1, x1, x2, x2], [y, y+ymax*0.02, y+ymax*0.02, y], color="k", lw=0.8)
        ax.text((x1+x2)/2, y+ymax*0.025, txt, ha="center", va="bottom", fontsize=6)
    brack(1, 2, ymax*0.78, f"{stars(pAN)} T vs N"); brack(0, 2, ymax*0.87, f"{stars(pAP)} T vs P")
    ax = axes[1]
    x = np.arange(len(SUBS)); w = 0.26
    for j, t in enumerate(ORDER):
        means = [m[(m.subtype==s)&(m.tissue==t)].log2.mean() for s in SUBS]
        cis = [1.96*m[(m.subtype==s)&(m.tissue==t)].log2.std()/np.sqrt(len(m[(m.subtype==s)&(m.tissue==t)])) for s in SUBS]
        ax.bar(x + (j-1)*w, means, w, yerr=cis, capsize=2, color=colors[t], label=LABEL[t], error_kw=dict(lw=0.8))
    ax.set_xticks(x); ax.set_xticklabels([SUBLAB[s] for s in SUBS])
    ax.set_ylabel(f"{gene}  log$_2$(TPM/10+1)")
    ax.set_title(f"{gene}: tumor vs blood/normal within T-cell subtypes")
    ax.legend(frameon=False, loc="upper left", ncol=3, fontsize=6, columnspacing=1.0, handlelength=1.2)
    ymaxB = 0
    for i, s in enumerate(SUBS):
        tops = [m[(m.subtype==s)&(m.tissue==t)].log2.mean()+1.96*m[(m.subtype==s)&(m.tissue==t)].log2.std()/np.sqrt(len(m[(m.subtype==s)&(m.tissue==t)])) for t in ORDER]
        top = max(tops); ymaxB = max(ymaxB, top)
        pTP = mw(m[(m.subtype==s)&(m.tissue=="Tumor")].log2, m[(m.subtype==s)&(m.tissue=="Peripheral")].log2)
        ax.text(i + w, top + 0.03, stars(pTP), ha="center", fontsize=7)
    ax.set_ylim(0, ymaxB * 1.20)
    for a, l in zip(axes, "ab"):
        a.text(-0.08, 1.02, l, transform=a.transAxes, fontweight="bold", fontsize=11, va="bottom")
    fig.suptitle(f"{gene} in tumor-infiltrating T cells — NSCLC ({GEO_ACCESSION}, {len(m):,} T cells)",
                 x=0.01, ha="left", fontsize=8.5, weight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    outfig = f"{args.outdir}/{gene}_{GEO_ACCESSION}.png"
    fig.savefig(outfig)
    print("saved", outfig)


if __name__ == "__main__":
    main()
