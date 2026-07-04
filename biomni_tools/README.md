# Lung T-cell scRNA-seq tools for Biomni

Packages the T-cell subset reprocessing pipeline (QC → normalize → HVG/PCA →
Harmony batch correction → UMAP → Leiden clustering → subset annotation →
marker heatmap/dot plot → donor-level stats → reference-atlas comparison),
originally built for the CELLxGENE Census LuCA lung cancer atlas, as
[Biomni](https://github.com/snap-stanford/Biomni) tools.

## Files

- `lung_tcell_scrna_tools.py` — 9 plain, type-annotated Python functions, one
  pipeline step each. Every function reads/writes files on disk and returns a
  short human-readable log string describing what it did (Biomni's
  convention: the return string is the LLM-visible "observation", not the
  data itself — that stays in the h5ad/csv/png files).
- `tool_description.py` — the declarative schema (`name`, `description`,
  `required_parameters`, `optional_parameters`) Biomni's tool retriever and
  registry use to decide when to surface each tool and how to validate a call.

## Tools

| Function | Purpose |
|---|---|
| `fetch_cellxgene_census_cells` | Pull a balanced cell subset from a CELLxGENE Census dataset by `dataset_id` |
| `load_h5ad_from_url` | Download an AnnData directly from a URL (e.g. a CELLxGENE Explorer download link) and optionally filter/downsample it |
| `qc_normalize_cluster_umap` | QC filter → normalize → HVG/PCA → Harmony → neighbors/UMAP → Leiden |
| `annotate_clusters_by_label_map` | Map cluster IDs to curated biological subset labels |
| `donor_level_composition_stats` | Per-donor subset fractions + cross-condition test (avoids pseudoreplication) |
| `plot_umap_by_group` | UMAP scatter colored by any categorical obs column |
| `plot_marker_gene_heatmap` | Mean-expression (z-scored) heatmap of marker genes across groups |
| `plot_marker_gene_dotplot` | Dot plot: size = % expressing, color = mean expression |
| `analyze_gene_expression_by_group` | Single-gene expression summary across groups + Kruskal-Wallis + donor-level check |
| `compare_clusters_to_reference_labels` | ARI / NMI / majority-label purity against a reference/published atlas |

Each function generalizes the concrete LuCA pipeline: any Census
`dataset_id`, any cell types/tissue/disease filter, any marker panel, any
manually curated cluster→label map.

## Integration into a Biomni checkout

Biomni supports two ways to add tools (see `CONTRIBUTION.md` in the Biomni
repo for the authoritative reference):

**A. Static registration (contributing to Biomni itself)**
1. Copy `lung_tcell_scrna_tools.py` to `biomni/tool/lung_tcell_scrna.py`.
2. Merge the schema list in `tool_description.py` into
   `biomni/tool/tool_description/lung_tcell_scrna.py` (or add it as a new
   module and register it in `biomni/tool/tool_description/__init__.py`
   alongside the other domain modules).
3. Add the new description module to the aggregation dict Biomni's registry
   builds from (`database_config.yaml` / the tool-description `__init__.py`
   import list, depending on Biomni version) so `ToolRegistry` picks it up
   at agent construction time.

**B. Dynamic registration on a running agent (no repo edit needed)**
```python
from biomni.agent import A1
import lung_tcell_scrna_tools as t
from tool_description import lung_tcell_scrna_tools as schemas

agent = A1(path="./biomni_data")
schema_by_name = {s["name"]: s for s in schemas}
for fn_name in [
    "fetch_cellxgene_census_cells", "load_h5ad_from_url",
    "qc_normalize_cluster_umap",
    "annotate_clusters_by_label_map", "donor_level_composition_stats",
    "plot_umap_by_group", "plot_marker_gene_heatmap",
    "plot_marker_gene_dotplot", "analyze_gene_expression_by_group",
    "compare_clusters_to_reference_labels",
]:
    fn = getattr(t, fn_name)
    agent.add_tool(fn, tool_description=schema_by_name[fn_name])
```
Use path B for quick experimentation in a notebook/session; use path A when
contributing the tools back to a shared Biomni install other users share.

## Example: prompt-driven data loading from a direct h5ad URL

When the user hands the agent a specific file location instead of a Census
`dataset_id` — e.g. pasting a CELLxGENE Explorer "Download Dataset" link —
prompt the agent along these lines:

> "Load the T cells from this dataset: https://datasets.cellxgene.cziscience.com/<uuid>.h5ad
> — keep only CD4-positive, CD8-positive, and regulatory T cells, downsample
> to 10,000 cells, then run QC, clustering, and UMAP."

Biomni's retriever should match this to `load_h5ad_from_url` (not
`fetch_cellxgene_census_cells`, which expects a Census `dataset_id` rather
than a URL) followed by `qc_normalize_cluster_umap`. Called directly:

```python
from lung_tcell_scrna_tools import load_h5ad_from_url

print(load_h5ad_from_url(
    url="https://datasets.cellxgene.cziscience.com/<uuid>.h5ad",
    output_h5ad="luca_tcells_raw.h5ad",
    cell_types=["CD4-positive, alpha-beta T cell", "CD8-positive, alpha-beta T cell",
                "regulatory T cell"],
    max_cells=10000))
```
Note: `max_cells` here downsamples uniformly at random; use
`fetch_cellxgene_census_cells` instead if you need the per-disease-group
balanced downsampling used in the original LuCA analysis.

## Example: reproduce the LuCA T-cell analysis end to end

```python
from lung_tcell_scrna_tools import (
    fetch_cellxgene_census_cells, qc_normalize_cluster_umap,
    annotate_clusters_by_label_map, plot_umap_by_group,
    plot_marker_gene_heatmap, plot_marker_gene_dotplot,
    donor_level_composition_stats, analyze_gene_expression_by_group,
)

print(fetch_cellxgene_census_cells(
    dataset_id="232f6a5a-a04c-4758-a6e8-88ab2e3a6e69",
    cell_types=["CD4-positive, alpha-beta T cell", "CD8-positive, alpha-beta T cell",
                "regulatory T cell", "T cell"],
    tissue="lung", n_per_group=3500, group_by="disease",
    output_h5ad="luca_tcells_raw.h5ad"))

print(qc_normalize_cluster_umap(
    input_h5ad="luca_tcells_raw.h5ad", output_h5ad="luca_tcells_processed.h5ad",
    batch_key="assay", leiden_resolution=1.0))

label_map = {
    "0": "CD4 T resting", "1": "CD4 Naive/Tcm", "3": "Regulatory T (Treg)",
    "4": "CD8 Effector Memory (GZMK+)", "6": "CD8 Terminal Effector (GZMH+GNLY+)",
    "8": "CD8 Exhausted/Tumor-reactive (CXCL13+GZMB+)", "9": "Proliferating T (MKI67+)",
    # ... remaining clusters
}
print(annotate_clusters_by_label_map(
    input_h5ad="luca_tcells_processed.h5ad", output_h5ad="luca_tcells_annotated.h5ad",
    cluster_key="leiden", label_map=label_map,
    low_quality_labels=["2", "5", "13", "14", "16"]))

markers = ["CCR7", "SELL", "TCF7", "IL7R", "FOXP3", "CTLA4", "IL2RA",
           "GZMK", "GZMA", "GZMH", "GZMB", "GNLY", "NKG7", "PRF1",
           "CXCL13", "HAVCR2", "PDCD1", "LAG3", "TIGIT", "MKI67", "KLRB1", "KLRC1"]

print(plot_umap_by_group(input_h5ad="luca_tcells_annotated.h5ad", color_col="subset",
                          output_png="umap_by_cluster.png"))
print(plot_marker_gene_heatmap(input_h5ad="luca_tcells_annotated.h5ad", marker_genes=markers,
                                groupby="subset", output_png="marker_heatmap.png"))
print(plot_marker_gene_dotplot(input_h5ad="luca_tcells_annotated.h5ad", marker_genes=markers,
                                groupby="subset", output_png="marker_dotplot.png"))
print(donor_level_composition_stats(
    input_h5ad="luca_tcells_annotated.h5ad", donor_col="donor_id", condition_col="disease",
    subset_col="subset", subsets_of_interest=["CD8 Exhausted/Tumor-reactive (CXCL13+GZMB+)"],
    reference_condition="normal", output_csv="donor_composition.csv"))
print(analyze_gene_expression_by_group(
    input_h5ad="luca_tcells_annotated.h5ad", gene="UBASH3A", groupby="subset",
    donor_col="donor_id", output_csv="ubash3a_by_subset.csv", output_png="ubash3a_by_subset.png"))
```

## Requirements

`scanpy`, `harmonypy`, `cellxgene_census`, `requests`, `pandas`, `numpy`,
`scipy`, `scikit-learn`, `matplotlib` — the same stack used in the original
`repro_scripts/` pipeline this module was refactored from.
