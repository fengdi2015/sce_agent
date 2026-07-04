"""
Biomni tool module: T-cell subset scRNA-seq analysis.

Follows the Biomni tool convention (see snap-stanford/Biomni CONTRIBUTION.md):
plain, type-annotated Python functions, one analysis step per function, each
returning a human-readable log string describing what was done and where
outputs were written (Biomni agents inspect this return string as the tool's
"observation"). Side effects (files written) are the actual deliverable;
the return string is documentation of those side effects, not the data
itself, so it stays short regardless of dataset size.

Drop this file in as biomni/tool/lung_tcell_scrna.py in a Biomni checkout,
or register the functions directly on a running agent with
`agent.add_tool(<function>)` (Biomni's dynamic tool-addition path, no
tool_registry.py edit required). See tool_description.py in this directory
for the declarative schema Biomni's retriever uses to decide when to surface
each tool, and README.md for both integration paths.

Generalizes the concrete lung-cancer LuCA-atlas pipeline this module was
built from into reusable steps: any CELLxGENE Census dataset_id, any set of
diseases/tissues/cell types, any marker panel, any subset annotation.
"""

from __future__ import annotations

import os
from typing import Optional


# ---------------------------------------------------------------------------
# 1. Data acquisition
# ---------------------------------------------------------------------------

def fetch_cellxgene_census_cells(
    dataset_id: str,
    cell_types: list[str],
    tissue: Optional[str] = None,
    diseases: Optional[list[str]] = None,
    n_per_group: int = 3500,
    group_by: str = "disease",
    output_h5ad: str = "census_cells_raw.h5ad",
    random_seed: int = 0,
) -> str:
    """Pull a balanced cell subset from a CELLxGENE Census dataset into an AnnData h5ad.

    Filters a named Census dataset to the requested cell type(s) (and,
    optionally, tissue and disease values), then downsamples to at most
    `n_per_group` cells per distinct value of `group_by` (e.g. balanced
    across disease groups) so downstream clustering isn't dominated by
    whichever group happens to have the most cells.

    Args:
        dataset_id: CELLxGENE Census dataset_id (UUID) to pull from, e.g. the
            LuCA lung cancer atlas ("232f6a5a-a04c-4758-a6e8-88ab2e3a6e69").
        cell_types: Cell Ontology `cell_type` labels to keep, e.g.
            ["CD4-positive, alpha-beta T cell", "CD8-positive, alpha-beta T cell",
             "regulatory T cell"].
        tissue: If given, restrict to this `tissue` value (e.g. "lung").
        diseases: If given, restrict to these `disease` values (e.g.
            ["normal", "lung adenocarcinoma"]). If None, all disease values
            present for the matched cells are kept.
        n_per_group: Max cells sampled per distinct value of `group_by`.
        group_by: obs column used to balance the downsample (default "disease").
        output_h5ad: Path to write the resulting raw (unprocessed) AnnData.
        random_seed: Seed for the per-group sampling.

    Returns:
        Log string with cell counts before/after balancing and the output path.
    """
    import numpy as np
    import pandas as pd
    import cellxgene_census

    census = cellxgene_census.open_soma(census_version="stable")

    value_filter = f"dataset_id == '{dataset_id}'"
    obs_df = cellxgene_census.get_obs(
        census, "homo_sapiens",
        value_filter=value_filter,
        column_names=["soma_joinid", "cell_type", "disease", "tissue", "assay", "sex"],
    )

    mask = obs_df.cell_type.isin(cell_types)
    if tissue is not None:
        mask &= obs_df.tissue.astype(str) == tissue
    if diseases is not None:
        mask &= obs_df.disease.astype(str).isin(diseases)
    sub = obs_df[mask].copy()
    n_before = len(sub)

    np.random.seed(random_seed)
    parts = []
    for _, g in sub.groupby(sub[group_by].astype(str)):
        n = min(n_per_group, len(g))
        parts.append(g.sample(n=n, random_state=random_seed))
    sel = pd.concat(parts)

    adata = cellxgene_census.get_anndata(
        census, organism="homo_sapiens",
        obs_coords=sel.soma_joinid.tolist(),
        column_names={"obs": ["soma_joinid", "cell_type", "disease", "tissue", "assay",
                                "sex", "dataset_id", "self_reported_ethnicity", "donor_id"]},
    )
    adata.write_h5ad(output_h5ad)

    return (f"Fetched {n_before} matching cells from dataset {dataset_id}; "
            f"balanced-downsampled (max {n_per_group}/{group_by}) to {adata.n_obs} cells "
            f"x {adata.n_vars} genes. Wrote raw AnnData to {output_h5ad}.")


def load_h5ad_from_url(
    url: str,
    output_h5ad: str = "loaded_dataset.h5ad",
    cell_types: Optional[list[str]] = None,
    cell_type_col: str = "cell_type",
    tissue: Optional[str] = None,
    tissue_col: str = "tissue",
    max_cells: Optional[int] = None,
    random_seed: int = 0,
    timeout: int = 600,
) -> str:
    """Download an AnnData object from a direct URL and (optionally) subset it.

    Use this instead of `fetch_cellxgene_census_cells` when the user hands
    Biomni a specific dataset location — e.g. a CELLxGENE Explorer "Download"
    link, a GEO/figshare/Zenodo direct file URL, or any other hosted
    `.h5ad` — rather than a Census `dataset_id` to query. Streams the file
    to disk (no full in-memory buffering of the download), then applies the
    same optional cell-type/tissue filter and downsample used elsewhere in
    this module so it slots into the same pipeline.

    Args:
        url: Direct HTTP(S) URL to an `.h5ad` file (must resolve to raw
            bytes, not an HTML landing page — e.g. CELLxGENE Explorer's
            "Download Dataset" link, not the dataset's browse-page URL).
        output_h5ad: Path to save the downloaded (and possibly subset)
            AnnData.
        cell_types: If given, keep only cells whose `cell_type_col` value is
            in this list, e.g. ["CD4-positive, alpha-beta T cell",
            "CD8-positive, alpha-beta T cell", "regulatory T cell"].
        cell_type_col: obs column holding cell-type labels, checked only if
            `cell_types` is given.
        tissue: If given, keep only cells whose `tissue_col` value equals
            this string.
        tissue_col: obs column holding tissue labels, checked only if
            `tissue` is given.
        max_cells: If given and fewer than the filtered cell count, randomly
            downsample to this many cells (uniform, not stratified — use
            `fetch_cellxgene_census_cells`'s balanced downsampling if you
            need per-group balance).
        random_seed: Seed for the downsample.
        timeout: Download timeout in seconds.

    Returns:
        Log string with the downloaded file size, cell/gene counts before
        and after filtering, and the output path.
    """
    import os
    import numpy as np
    import requests
    import scanpy as sc

    raw_path = output_h5ad + ".download"
    with requests.get(url, stream=True, timeout=timeout) as r:
        r.raise_for_status()
        content_type = r.headers.get("content-type", "")
        if "html" in content_type.lower():
            raise ValueError(
                f"URL {url!r} returned content-type {content_type!r}, which looks like an "
                "HTML page rather than an h5ad file. Use a direct file-download link "
                "(e.g. CELLxGENE Explorer's 'Download Dataset' URL), not a browse-page URL."
            )
        with open(raw_path, "wb") as f:
            for chunk in r.iter_content(chunk_size=8 * 1024 * 1024):
                f.write(chunk)
    size_mb = os.path.getsize(raw_path) / 1e6

    adata = sc.read_h5ad(raw_path)
    os.remove(raw_path)
    n_before = adata.n_obs

    if cell_types is not None:
        if cell_type_col not in adata.obs.columns:
            raise ValueError(f"'{cell_type_col}' not found in obs columns: {list(adata.obs.columns)}")
        adata = adata[adata.obs[cell_type_col].astype(str).isin(cell_types)].copy()
    if tissue is not None:
        if tissue_col not in adata.obs.columns:
            raise ValueError(f"'{tissue_col}' not found in obs columns: {list(adata.obs.columns)}")
        adata = adata[adata.obs[tissue_col].astype(str) == tissue].copy()
    n_after_filter = adata.n_obs

    if max_cells is not None and adata.n_obs > max_cells:
        rng = np.random.default_rng(random_seed)
        idx = rng.choice(adata.n_obs, size=max_cells, replace=False)
        adata = adata[np.sort(idx)].copy()

    adata.write_h5ad(output_h5ad)

    return (f"Downloaded {size_mb:.1f} MB from {url} ({n_before} cells x {adata.n_vars} genes). "
            f"After cell_type/tissue filter: {n_after_filter} cells. "
            f"After downsample (max_cells={max_cells}): {adata.n_obs} cells. "
            f"Wrote AnnData to {output_h5ad}.")


# ---------------------------------------------------------------------------
# 2. QC, normalization, batch-corrected clustering + UMAP
# ---------------------------------------------------------------------------

def qc_normalize_cluster_umap(
    input_h5ad: str,
    output_h5ad: str = "processed.h5ad",
    batch_key: Optional[str] = "assay",
    min_genes: int = 200,
    min_cells: int = 3,
    max_total_counts_quantile: float = 0.995,
    drop_assay_values: Optional[list[str]] = ("Smart-seq2",),
    n_top_genes: int = 2000,
    n_pcs: int = 50,
    n_neighbors: int = 15,
    leiden_resolution: float = 1.0,
) -> str:
    """QC filter, normalize, batch-correct, cluster, and UMAP-embed an AnnData.

    Pipeline: filter cells (min genes) and genes (min cells) -> drop cells
    above the `max_total_counts_quantile` total-count percentile (doublet-ish
    outliers) -> optionally drop specified assay values (e.g. non-UMI
    Smart-seq2 counts, which are not magnitude-comparable to 10x UMI counts)
    -> log1p-normalize -> HVG selection -> PCA -> Harmony batch correction on
    `batch_key` (if given) -> neighbor graph -> UMAP -> Leiden clustering.

    Args:
        input_h5ad: Path to raw AnnData (expects a `feature_name` var column
            if var_names are not already gene symbols, and raw counts in .X).
        output_h5ad: Path to write the processed AnnData. Adds/uses:
            layers["counts"], layers["lognorm"], obsm["X_pca"],
            obsm["X_pca_harmony"] (if batch_key set), obsm["X_umap"],
            obs["leiden"]. Sets .raw to the full log-normalized matrix.
        batch_key: obs column for Harmony batch correction; set to None to
            skip batch correction and cluster directly on PCA space.
        min_genes: Minimum genes/cell (sc.pp.filter_cells).
        min_cells: Minimum cells/gene (sc.pp.filter_genes).
        max_total_counts_quantile: Upper total-count quantile cutoff.
        drop_assay_values: obs["assay"] values to exclude entirely (needs an
            "assay" column); pass None to keep all assays.
        n_top_genes: Number of highly variable genes for PCA.
        n_pcs: Number of principal components.
        n_neighbors: Neighbors for the kNN graph feeding UMAP/Leiden.
        leiden_resolution: Leiden clustering resolution.

    Returns:
        Log string with cell/gene counts after each filtering step, cluster
        sizes, and the output path.
    """
    import numpy as np
    import scanpy as sc

    adata = sc.read_h5ad(input_h5ad)
    if "feature_name" in adata.var.columns:
        adata.var_names = adata.var["feature_name"].astype(str)
        adata.var_names_make_unique()
        adata.var.index.name = None
    adata.layers["counts"] = adata.X.copy()

    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True, percent_top=None)

    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_genes(adata, min_cells=min_cells)
    upper = adata.obs.total_counts.quantile(max_total_counts_quantile)
    adata = adata[adata.obs.total_counts <= upper].copy()
    if drop_assay_values and "assay" in adata.obs.columns:
        adata = adata[~adata.obs.assay.astype(str).isin(list(drop_assay_values))].copy()
    n_after_qc = adata.n_obs

    adata.X = adata.layers["counts"].copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.layers["lognorm"] = adata.X.copy()

    hvg_kwargs = {"n_top_genes": n_top_genes}
    if batch_key is not None and batch_key in adata.obs.columns:
        hvg_kwargs["batch_key"] = batch_key
    sc.pp.highly_variable_genes(adata, **hvg_kwargs)
    adata.raw = adata
    adata_hvg = adata[:, adata.var.highly_variable].copy()
    sc.pp.scale(adata_hvg, max_value=10)
    sc.tl.pca(adata_hvg, n_comps=n_pcs, svd_solver="arpack")
    adata.obsm["X_pca"] = adata_hvg.obsm["X_pca"]

    rep = "X_pca"
    if batch_key is not None and batch_key in adata.obs.columns:
        import harmonypy
        ho = harmonypy.run_harmony(adata.obsm["X_pca"], adata.obs, [batch_key], max_iter_harmony=20)
        Zc = np.asarray(ho.Z_corr)
        adata.obsm["X_pca_harmony"] = Zc.T if Zc.shape[0] != adata.n_obs else Zc
        rep = "X_pca_harmony"

    sc.pp.neighbors(adata, use_rep=rep, n_neighbors=n_neighbors)
    sc.tl.umap(adata)
    sc.tl.leiden(adata, resolution=leiden_resolution, key_added="leiden", flavor="igraph", n_iterations=2)

    adata.write_h5ad(output_h5ad)
    cluster_counts = adata.obs.leiden.value_counts().sort_index()

    return (f"QC: {n_after_qc} cells retained after filtering (min_genes={min_genes}, "
            f"total_counts<={max_total_counts_quantile:.3f} quantile"
            f"{', assay filter applied' if drop_assay_values else ''}). "
            f"Batch correction: {'Harmony on ' + batch_key if rep == 'X_pca_harmony' else 'none'}. "
            f"Leiden (res={leiden_resolution}) found {cluster_counts.shape[0]} clusters, "
            f"sizes {dict(cluster_counts)}. Wrote processed AnnData with new UMAP embedding "
            f"(obsm['X_umap']) to {output_h5ad}.")


# ---------------------------------------------------------------------------
# 3. Cluster -> biological subset annotation
# ---------------------------------------------------------------------------

def annotate_clusters_by_label_map(
    input_h5ad: str,
    output_h5ad: str,
    cluster_key: str,
    label_map: dict[str, str],
    low_quality_labels: Optional[list[str]] = None,
    new_col: str = "subset",
) -> str:
    """Map cluster IDs (e.g. Leiden) to manually curated biological subset labels.

    Any cluster present in the data but missing from `label_map` is labeled
    "Cluster <id> (unannotated)" rather than silently dropped, so re-running
    clustering at a different resolution never loses cells.

    Args:
        input_h5ad: Path to a clustered AnnData (must have `cluster_key` in .obs).
        output_h5ad: Path to write the annotated AnnData.
        cluster_key: obs column holding cluster IDs, e.g. "leiden".
        label_map: Mapping from cluster ID (str) to subset label, e.g.
            {"3": "Regulatory T (Treg)", "10": "CD8 Exhausted (CXCL13+)"}.
        low_quality_labels: Subset of `label_map` VALUES (or cluster-id keys)
            to flag as low quality in a new obs["is_low_quality"] boolean
            column, e.g. clusters driven by ambient RNA, stress response, or
            doublets. If None, obs["is_low_quality"] is all False.
        new_col: Name of the obs column to write the mapped label into.

    Returns:
        Log string with the resulting label distribution and output path.
    """
    import scanpy as sc

    adata = sc.read_h5ad(input_h5ad)
    present = adata.obs[cluster_key].astype(str).unique()
    full_map = dict(label_map)
    for c in present:
        full_map.setdefault(c, f"Cluster {c} (unannotated)")
    adata.obs[new_col] = adata.obs[cluster_key].astype(str).map(full_map).astype("category")

    if low_quality_labels:
        lq = set(low_quality_labels)
        adata.obs["is_low_quality"] = (
            adata.obs[cluster_key].astype(str).isin(lq) | adata.obs[new_col].astype(str).isin(lq)
        )
    else:
        adata.obs["is_low_quality"] = False

    adata.write_h5ad(output_h5ad)
    counts = adata.obs[new_col].value_counts()

    return (f"Annotated {adata.n_obs} cells into {counts.shape[0]} labels via '{new_col}': "
            f"{dict(counts)}. Flagged {int(adata.obs['is_low_quality'].sum())} cells as "
            f"low quality. Wrote {output_h5ad}.")


# ---------------------------------------------------------------------------
# 4. Donor-level composition + group comparison
# ---------------------------------------------------------------------------

def donor_level_composition_stats(
    input_h5ad: str,
    donor_col: str,
    condition_col: str,
    subset_col: str,
    subsets_of_interest: list[str],
    exclude_low_quality: bool = True,
    min_cells_per_donor: int = 10,
    reference_condition: Optional[str] = None,
    output_csv: str = "donor_level_composition.csv",
) -> str:
    """Compute per-donor subset fractions and test them across condition groups.

    Aggregating to one row per donor (rather than testing at the cell level)
    avoids pseudoreplication: cells from the same donor are correlated, so a
    cell-level test overstates significance whenever a condition group is
    dominated by a few donors with many captured cells.

    Args:
        input_h5ad: Path to an annotated AnnData (needs `donor_col`,
            `condition_col`, `subset_col` in .obs).
        donor_col: obs column identifying the donor/sample of origin.
        condition_col: obs column identifying the comparison groups (e.g.
            "disease").
        subset_col: obs column with the cell-subset / cluster label.
        subsets_of_interest: Subset labels to compute per-donor % for, e.g.
            ["CD8 Exhausted/Tumor-reactive (CXCL13+GZMB+)", "Regulatory T (Treg)"].
        exclude_low_quality: If True and obs["is_low_quality"] exists, drop
            those cells first.
        min_cells_per_donor: Minimum cells for a donor to be included (donors
            with fewer cells give unstable fraction estimates).
        reference_condition: If given, this condition's donors are compared
            against all other conditions pooled via Mann-Whitney U
            (two-sided) for each subset in `subsets_of_interest`. If None,
            no statistical test is run — only the composition table is written.
        output_csv: Path to write the donor-level composition table.

    Returns:
        Log string with donor/cell counts and (if `reference_condition` set)
        the Mann-Whitney U statistic and p-value per subset.
    """
    import scanpy as sc
    import pandas as pd
    from scipy.stats import mannwhitneyu

    adata = sc.read_h5ad(input_h5ad)
    if exclude_low_quality and "is_low_quality" in adata.obs.columns:
        adata = adata[~adata.obs["is_low_quality"]].copy()

    def _row(g):
        out = {"n_cells": len(g)}
        for s in subsets_of_interest:
            out[f"pct_{s}"] = 100 * (g[subset_col] == s).mean()
        return pd.Series(out)

    donor_ct = adata.obs.groupby([donor_col, condition_col], observed=True).apply(
        _row, include_groups=False
    ).reset_index()
    donor_ct = donor_ct[donor_ct.n_cells >= min_cells_per_donor]
    donor_ct.to_csv(output_csv, index=False)

    log = (f"{adata.n_obs} cells across {donor_ct[donor_col].nunique()} donors "
           f"(>= {min_cells_per_donor} cells/donor) written to {output_csv}.")

    if reference_condition is not None:
        is_ref = donor_ct[condition_col].astype(str) == reference_condition
        ref_grp, other_grp = donor_ct[is_ref], donor_ct[~is_ref]
        test_lines = []
        for s in subsets_of_interest:
            col = f"pct_{s}"
            if len(ref_grp[col]) < 2 or len(other_grp[col]) < 2:
                continue
            stat, p = mannwhitneyu(other_grp[col], ref_grp[col], alternative="two-sided")
            test_lines.append(
                f"  {s}: {reference_condition} median={ref_grp[col].median():.2f}%, "
                f"other median={other_grp[col].median():.2f}%, U={stat:.1f}, p={p:.2e}"
            )
        log += "\nMann-Whitney U (other conditions vs. " + reference_condition + "):\n" + "\n".join(test_lines)

    return log


# ---------------------------------------------------------------------------
# 5. Visualization: UMAP
# ---------------------------------------------------------------------------

def plot_umap_by_group(
    input_h5ad: str,
    color_col: str,
    output_png: str = "umap_by_group.png",
    exclude_low_quality: bool = True,
    point_size: float = 3.0,
    title: Optional[str] = None,
) -> str:
    """Scatter-plot an AnnData's UMAP embedding colored by a categorical obs column.

    Args:
        input_h5ad: Path to an AnnData with obsm["X_umap"] already computed
            (see `qc_normalize_cluster_umap`).
        color_col: obs column to color points by (e.g. "subset", "leiden",
            "disease").
        output_png: Path to save the figure.
        exclude_low_quality: If True and obs["is_low_quality"] exists, drop
            those cells first.
        point_size: Marker size passed to matplotlib scatter.
        title: Plot title; defaults to "UMAP colored by <color_col>".

    Returns:
        Log string with the number of cells/groups plotted and output path.
    """
    import numpy as np
    import scanpy as sc
    import matplotlib.pyplot as plt

    adata = sc.read_h5ad(input_h5ad)
    if exclude_low_quality and "is_low_quality" in adata.obs.columns:
        adata = adata[~adata.obs["is_low_quality"]].copy()

    groups = adata.obs[color_col].astype(str).value_counts().index.tolist()
    palette = plt.get_cmap("tab20")(np.linspace(0, 1, len(groups)))
    color_map = dict(zip(groups, palette))
    coords = adata.obsm["X_umap"]

    fig, ax = plt.subplots(figsize=(8, 6.5))
    for g in groups:
        m = adata.obs[color_col].astype(str) == g
        ax.scatter(coords[m, 0], coords[m, 1], s=point_size, alpha=0.7, color=color_map[g], label=g, linewidths=0)
    ax.set_xlabel("UMAP1"); ax.set_ylabel("UMAP2")
    ax.set_title(title or f"UMAP colored by {color_col}", loc="left")
    ax.legend(fontsize=6, markerscale=3, loc="upper left", bbox_to_anchor=(1.0, 1.0), frameon=False)
    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)
    fig.tight_layout()
    fig.savefig(output_png, dpi=200, bbox_inches="tight")

    return f"Plotted {adata.n_obs} cells across {len(groups)} '{color_col}' groups to {output_png}."


# ---------------------------------------------------------------------------
# 6. Visualization: marker heatmap + dot plot
# ---------------------------------------------------------------------------

def _mean_and_pct_expr(adata, marker_genes, groupby):
    import pandas as pd
    marker_genes = [g for g in marker_genes if g in adata.raw.var_names]
    if not marker_genes:
        raise ValueError("None of the requested marker_genes were found in adata.raw.var_names.")
    X = adata.raw[:, marker_genes].X
    expr = pd.DataFrame(X.toarray() if hasattr(X, "toarray") else X,
                         columns=marker_genes, index=adata.obs_names)
    expr[groupby] = adata.obs[groupby].values
    mean_expr = expr.groupby(groupby, observed=True).mean()
    pct_expr = (expr.drop(columns=groupby) > 0).groupby(expr[groupby], observed=True).mean() * 100
    return marker_genes, mean_expr, pct_expr.loc[mean_expr.index]


def plot_marker_gene_heatmap(
    input_h5ad: str,
    marker_genes: list[str],
    groupby: str,
    output_png: str = "marker_heatmap.png",
    exclude_low_quality: bool = True,
    zscore: bool = True,
) -> str:
    """Plot a mean-marker-expression heatmap across groups (e.g. cell subsets).

    Args:
        input_h5ad: Path to an AnnData with .raw set to log-normalized values
            (see `qc_normalize_cluster_umap`, which sets .raw automatically).
        marker_genes: Gene symbols to plot as heatmap columns; any not found
            in adata.raw.var_names are silently skipped.
        groupby: obs column defining heatmap rows (e.g. "subset").
        output_png: Path to save the figure.
        exclude_low_quality: If True and obs["is_low_quality"] exists, drop
            those cells first.
        zscore: If True, z-score each gene's mean expression across groups
            (highlights relative enrichment); if False, plot raw mean
            log-normalized expression.

    Returns:
        Log string with genes actually plotted (vs. requested) and output path.
    """
    import scanpy as sc
    import matplotlib.pyplot as plt

    adata = sc.read_h5ad(input_h5ad)
    if exclude_low_quality and "is_low_quality" in adata.obs.columns:
        adata = adata[~adata.obs["is_low_quality"]].copy()

    genes_used, mean_expr, _ = _mean_and_pct_expr(adata, marker_genes, groupby)
    order = adata.obs[groupby].value_counts().index
    mean_expr = mean_expr.loc[[g for g in order if g in mean_expr.index]]

    plot_vals = mean_expr
    label = "mean log-norm expr"
    if zscore:
        plot_vals = (mean_expr - mean_expr.mean(axis=0)) / (mean_expr.std(axis=0) + 1e-9)
        label = "z-score"

    fig, ax = plt.subplots(figsize=(max(6, 0.4 * len(genes_used)), max(4, 0.3 * len(mean_expr))))
    vmax = 2 if zscore else float(plot_vals.values.max())
    vmin = -2 if zscore else 0
    im = ax.imshow(plot_vals.values, cmap="RdBu_r" if zscore else "viridis", vmin=vmin, vmax=vmax, aspect="auto")
    ax.set_xticks(range(len(genes_used))); ax.set_xticklabels(genes_used, rotation=90, fontsize=7)
    ax.set_yticks(range(len(plot_vals))); ax.set_yticklabels(plot_vals.index, fontsize=7)
    ax.set_title(f"Marker expression by {groupby}", loc="left")
    cbar = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
    cbar.set_label(label, fontsize=7)
    fig.tight_layout()
    fig.savefig(output_png, dpi=200, bbox_inches="tight")

    skipped = set(marker_genes) - set(genes_used)
    note = f" (skipped, not found: {sorted(skipped)})" if skipped else ""
    return f"Plotted {len(genes_used)}/{len(marker_genes)} marker genes x {len(plot_vals)} '{groupby}' groups{note} to {output_png}."


def plot_marker_gene_dotplot(
    input_h5ad: str,
    marker_genes: list[str],
    groupby: str,
    output_png: str = "marker_dotplot.png",
    exclude_low_quality: bool = True,
) -> str:
    """Dot plot of marker genes across groups: dot size = % expressing, color = mean expression.

    Args:
        input_h5ad: Path to an AnnData with .raw set to log-normalized values.
        marker_genes: Gene symbols to plot; any not found are silently skipped.
        groupby: obs column defining dot-plot rows (e.g. "subset").
        output_png: Path to save the figure.
        exclude_low_quality: If True and obs["is_low_quality"] exists, drop
            those cells first.

    Returns:
        Log string with genes actually plotted (vs. requested) and output path.
    """
    import scanpy as sc
    import matplotlib.pyplot as plt

    adata = sc.read_h5ad(input_h5ad)
    if exclude_low_quality and "is_low_quality" in adata.obs.columns:
        adata = adata[~adata.obs["is_low_quality"]].copy()

    genes_used, mean_expr, pct_expr = _mean_and_pct_expr(adata, marker_genes, groupby)
    order = adata.obs[groupby].value_counts().index
    mean_expr = mean_expr.loc[[g for g in order if g in mean_expr.index]]
    pct_expr = pct_expr.loc[mean_expr.index]

    fig, ax = plt.subplots(figsize=(max(6, 0.4 * len(genes_used)), max(4, 0.3 * len(mean_expr))))
    xs, ys, sizes, colors = [], [], [], []
    for i, s in enumerate(mean_expr.index):
        for j, g in enumerate(genes_used):
            xs.append(j); ys.append(i)
            sizes.append(pct_expr.loc[s, g] * 6)
            colors.append(mean_expr.loc[s, g])
    sc_plot = ax.scatter(xs, ys, s=sizes, c=colors, cmap="viridis", edgecolors="black", linewidths=0.3)
    ax.set_xticks(range(len(genes_used))); ax.set_xticklabels(genes_used, rotation=90, fontsize=7)
    ax.set_yticks(range(len(mean_expr))); ax.set_yticklabels(mean_expr.index, fontsize=7)
    ax.invert_yaxis()
    ax.set_title(f"Marker expression dot plot by {groupby}", loc="left")
    for pct_val in [10, 30, 60]:
        ax.scatter([], [], s=pct_val * 6, c="grey", edgecolors="black", linewidths=0.3, label=f"{pct_val}%")
    ax.legend(title="% expressing", loc="upper left", bbox_to_anchor=(1.12, 1.0), fontsize=7, title_fontsize=7.5, frameon=False)
    cbar = fig.colorbar(sc_plot, ax=ax, fraction=0.03, pad=0.09)
    cbar.set_label("mean expr (all cells)", fontsize=7)
    fig.tight_layout()
    fig.savefig(output_png, dpi=200, bbox_inches="tight")

    skipped = set(marker_genes) - set(genes_used)
    note = f" (skipped, not found: {sorted(skipped)})" if skipped else ""
    return f"Plotted {len(genes_used)}/{len(marker_genes)} marker genes x {len(mean_expr)} '{groupby}' groups{note} to {output_png}."


# ---------------------------------------------------------------------------
# 7. Single-gene expression analysis across subsets
# ---------------------------------------------------------------------------

def analyze_gene_expression_by_group(
    input_h5ad: str,
    gene: str,
    groupby: str,
    donor_col: Optional[str] = None,
    min_cells_per_donor: int = 5,
    exclude_low_quality: bool = True,
    output_csv: str = "gene_expression_by_group.csv",
    output_png: Optional[str] = "gene_expression_by_group.png",
) -> str:
    """Summarize one gene's expression across groups, with an optional donor-level check.

    Computes, per group: mean/median log-normalized expression, % of cells
    expressing (count > 0), and mean expression among only the expressing
    cells. Runs a Kruskal-Wallis test for a difference in expression across
    groups. If `donor_col` is given, also aggregates to donor-level means
    (>= `min_cells_per_donor` cells) as a check that the cell-level signal
    isn't driven by one or two donors, and plots both a dot-plot summary and
    a donor-level strip plot.

    Args:
        input_h5ad: Path to an AnnData with .raw / a "lognorm" layer set.
        gene: Gene symbol to analyze.
        groupby: obs column defining the groups (e.g. "subset").
        donor_col: obs column identifying donors, for the donor-level check;
            skip donor-level validation if None.
        min_cells_per_donor: Minimum cells per donor to include in the
            donor-level check.
        exclude_low_quality: If True and obs["is_low_quality"] exists, drop
            those cells first.
        output_csv: Path to write the per-group summary table.
        output_png: Path to save the figure; pass None to skip plotting.

    Returns:
        Log string with the per-group summary and Kruskal-Wallis result.
    """
    import numpy as np
    import scanpy as sc
    import pandas as pd
    from scipy.stats import kruskal

    adata = sc.read_h5ad(input_h5ad)
    if exclude_low_quality and "is_low_quality" in adata.obs.columns:
        adata = adata[~adata.obs["is_low_quality"]].copy()

    if gene not in adata.var_names:
        raise ValueError(f"Gene {gene!r} not found in adata.var_names.")
    gene_idx = adata.var.index.get_loc(gene)
    layer = "lognorm" if "lognorm" in adata.layers else None
    mat = adata.layers[layer] if layer else adata.X
    expr = mat[:, gene_idx]
    expr = expr.toarray().ravel() if hasattr(expr, "toarray") else np.asarray(expr).ravel()
    adata.obs[f"{gene}_expr"] = expr

    cols = [groupby, f"{gene}_expr"] + ([donor_col] if donor_col else [])
    df = adata.obs[cols].copy()

    summary = df.groupby(groupby, observed=True)[f"{gene}_expr"].agg(
        mean_expr="mean", median_expr="median",
        pct_expressing=lambda x: (x > 0).mean() * 100, n_cells="size",
    ).sort_values("mean_expr", ascending=False)
    mean_expressing = df[df[f"{gene}_expr"] > 0].groupby(groupby, observed=True)[f"{gene}_expr"].mean()
    summary["mean_expr_among_expressers"] = mean_expressing.reindex(summary.index)
    summary.to_csv(output_csv)

    groups_for_test = [g[f"{gene}_expr"].values for _, g in df.groupby(groupby, observed=True)]
    stat, p = kruskal(*groups_for_test)

    log = (f"{gene} expression across {summary.shape[0]} '{groupby}' groups written to {output_csv}. "
           f"Kruskal-Wallis: H={stat:.2f}, p={p:.3e}.")

    donor_summary = None
    if donor_col is not None:
        donor_summary = df.groupby([groupby, donor_col], observed=True).agg(
            mean_expr=(f"{gene}_expr", "mean"), n=(f"{gene}_expr", "size")
        ).reset_index()
        donor_summary = donor_summary[donor_summary.n >= min_cells_per_donor]
        log += (f" Donor-level check (>= {min_cells_per_donor} cells/donor): "
                f"{donor_summary[donor_col].nunique()} donors retained.")

    if output_png:
        import matplotlib.pyplot as plt
        order = summary.index.tolist()
        ncols = 2 if donor_summary is not None else 1
        fig, axes = plt.subplots(1, ncols, figsize=(6 * ncols, 5.2), squeeze=False)
        ax_dot = axes[0][0]

        y_pos = np.arange(len(order))
        sizes = summary.loc[order, "pct_expressing"].values * 12
        sc_ = ax_dot.scatter(mean_expressing.reindex(order).values, y_pos, s=sizes,
                              c=mean_expressing.reindex(order).values, cmap="viridis",
                              edgecolors="black", linewidths=0.4, zorder=3)
        ax_dot.set_yticks(y_pos); ax_dot.set_yticklabels(order, fontsize=7.5)
        ax_dot.invert_yaxis()
        ax_dot.set_xlabel(f"Mean log-norm {gene} expr (among {gene}+ cells)")
        ax_dot.set_title(f"{gene} expression by {groupby}", loc="left")
        ax_dot.grid(axis="x", alpha=0.3, zorder=0)
        fig.colorbar(sc_, ax=ax_dot, fraction=0.04, pad=0.02).set_label("mean expr (expressers)", fontsize=7)

        if donor_summary is not None:
            ax_donor = axes[0][1]
            rng = np.random.default_rng(0)
            for i, s in enumerate(order):
                vals = donor_summary.loc[donor_summary[groupby] == s, "mean_expr"].values
                jitter = rng.uniform(-0.15, 0.15, size=len(vals))
                ax_donor.scatter(np.full(len(vals), i) + jitter, vals, s=10, alpha=0.6, color="#4C72B0", linewidths=0)
                if len(vals):
                    med = np.median(vals)
                    ax_donor.plot([i - 0.22, i + 0.22], [med, med], color="black", lw=1.8, zorder=5)
            ax_donor.set_xticks(range(len(order)))
            ax_donor.set_xticklabels(order, rotation=45, ha="right", fontsize=7.5)
            ax_donor.set_ylabel(f"Donor-mean log-norm {gene} expr")
            ax_donor.set_title(f"Donor-level check (n>={min_cells_per_donor} cells/donor)", loc="left")

        fig.tight_layout()
        fig.savefig(output_png, dpi=200, bbox_inches="tight")
        log += f" Plotted to {output_png}."

    return log


# ---------------------------------------------------------------------------
# 8. Concordance against a reference / published atlas
# ---------------------------------------------------------------------------

def compare_clusters_to_reference_labels(
    own_h5ad: str,
    reference_h5ad: str,
    own_label_col: str,
    reference_label_col: str,
    join_key: str,
    output_crosstab_csv: str = "label_crosstab.csv",
    output_mapping_csv: str = "label_mapping.csv",
) -> str:
    """Quantify concordance between this analysis' cluster labels and a reference atlas's labels.

    Matches cells between the two AnnData objects on a shared identifier
    column (e.g. a Census `observation_joinid`, or `obs_names` if both were
    subset from the same source), then computes Adjusted Rand Index (ARI),
    Normalized Mutual Information (NMI), and a weighted majority-label
    purity: the fraction of cells whose own-analysis cluster's single most
    common reference label matches that cell's actual reference label.

    Args:
        own_h5ad: Path to this analysis' annotated AnnData.
        reference_h5ad: Path to the reference/published AnnData (must share
            `join_key` values with `own_h5ad`, typically a subset of a much
            larger reference — only overlapping cells are compared).
        own_label_col: obs column in `own_h5ad` with this analysis' cluster
            or subset labels.
        reference_label_col: obs column in `reference_h5ad` with the
            reference author labels.
        join_key: obs column present in both, used to match cells 1:1.
        output_crosstab_csv: Path to write the full label x label cross-tab.
        output_mapping_csv: Path to write, per own-label, the majority
            reference label and purity.

    Returns:
        Log string with n matched cells, ARI, NMI, and overall purity.
    """
    import scanpy as sc
    import pandas as pd
    from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score

    own = sc.read_h5ad(own_h5ad)
    ref = sc.read_h5ad(reference_h5ad)

    own_df = own.obs[[join_key, own_label_col]].copy()
    ref_df = ref.obs[[join_key, reference_label_col]].copy()
    merged = own_df.merge(ref_df, on=join_key, how="inner", validate="one_to_one")
    n_matched = len(merged)
    if n_matched == 0:
        raise ValueError(f"No cells matched between own_h5ad and reference_h5ad on '{join_key}'.")

    ari = adjusted_rand_score(merged[reference_label_col], merged[own_label_col])
    nmi = normalized_mutual_info_score(merged[reference_label_col], merged[own_label_col])

    ct = pd.crosstab(merged[own_label_col], merged[reference_label_col])
    row_frac = ct.div(ct.sum(axis=1), axis=0)
    mapping_df = pd.DataFrame({
        "majority_reference_label": row_frac.idxmax(axis=1),
        "purity": row_frac.max(axis=1),
        "n_cells": ct.sum(axis=1),
    }).sort_values("purity", ascending=False)
    overall_purity = ct.max(axis=1).sum() / ct.values.sum()

    ct.to_csv(output_crosstab_csv)
    mapping_df.to_csv(output_mapping_csv)

    return (f"Matched {n_matched} cells between own and reference labels via '{join_key}'. "
            f"ARI={ari:.3f}, NMI={nmi:.3f}, weighted majority-label purity={overall_purity:.3f}. "
            f"Crosstab -> {output_crosstab_csv}; per-label mapping -> {output_mapping_csv}.")
