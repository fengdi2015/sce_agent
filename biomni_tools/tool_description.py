"""
Declarative tool schema for lung_tcell_scrna_tools.py, in Biomni's
tool_description format (list of dicts with name/description/
required_parameters/optional_parameters). Biomni's tool registry and
retriever read this structure to decide which tool to surface for a
given task and to validate/construct the call; it is the schema half
of the tool, the function body in lung_tcell_scrna_tools.py is the
implementation half.
"""

lung_tcell_scrna_tools = [
    {
        "name": "fetch_cellxgene_census_cells",
        "description": "Pull a balanced cell subset from a CELLxGENE Census dataset into an AnnData h5ad.",
        "required_parameters": [
            {
                "name": "dataset_id",
                "type": "str",
                "description": "CELLxGENE Census dataset_id (UUID) to pull from, e.g. the LuCA lung cancer atlas (\"232f6a5a-a04c-4758-a6e8-88ab2e3a6e69\").",
            },
            {
                "name": "cell_types",
                "type": "List[str]",
                "description": "Cell Ontology `cell_type` labels to keep, e.g. [\"CD4-positive, alpha-beta T cell\", \"CD8-positive, alpha-beta T cell\", \"regulatory T cell\"].",
            },
        ],
        "optional_parameters": [
            {
                "name": "tissue",
                "type": "str",
                "default": "None",
                "description": "If given, restrict to this `tissue` value (e.g. \"lung\").",
            },
            {
                "name": "diseases",
                "type": "List[str]",
                "default": "None",
                "description": "If given, restrict to these `disease` values (e.g. [\"normal\", \"lung adenocarcinoma\"]). If None, all disease values present for the matched cells are kept.",
            },
            {
                "name": "n_per_group",
                "type": "int",
                "default": "3500",
                "description": "Max cells sampled per distinct value of `group_by`.",
            },
            {
                "name": "group_by",
                "type": "str",
                "default": "'disease'",
                "description": "obs column used to balance the downsample (default \"disease\").",
            },
            {
                "name": "output_h5ad",
                "type": "str",
                "default": "'census_cells_raw.h5ad'",
                "description": "Path to write the resulting raw (unprocessed) AnnData.",
            },
            {
                "name": "random_seed",
                "type": "int",
                "default": "0",
                "description": "Seed for the per-group sampling.",
            },
        ],
    },
    {
        "name": "load_h5ad_from_url",
        "description": "Download an AnnData object from a direct URL and (optionally) subset it.",
        "required_parameters": [
            {
                "name": "url",
                "type": "str",
                "description": "Direct HTTP(S) URL to an `.h5ad` file (must resolve to raw bytes, not an HTML landing page \u2014 e.g. CELLxGENE Explorer's \"Download Dataset\" link, not the dataset's browse-page URL).",
            },
        ],
        "optional_parameters": [
            {
                "name": "output_h5ad",
                "type": "str",
                "default": "'loaded_dataset.h5ad'",
                "description": "Path to save the downloaded (and possibly subset) AnnData.",
            },
            {
                "name": "cell_types",
                "type": "List[str]",
                "default": "None",
                "description": "If given, keep only cells whose `cell_type_col` value is in this list, e.g. [\"CD4-positive, alpha-beta T cell\", \"CD8-positive, alpha-beta T cell\", \"regulatory T cell\"].",
            },
            {
                "name": "cell_type_col",
                "type": "str",
                "default": "'cell_type'",
                "description": "obs column holding cell-type labels, checked only if `cell_types` is given.",
            },
            {
                "name": "tissue",
                "type": "str",
                "default": "None",
                "description": "If given, keep only cells whose `tissue_col` value equals this string.",
            },
            {
                "name": "tissue_col",
                "type": "str",
                "default": "'tissue'",
                "description": "obs column holding tissue labels, checked only if `tissue` is given.",
            },
            {
                "name": "max_cells",
                "type": "int",
                "default": "None",
                "description": "If given and fewer than the filtered cell count, randomly downsample to this many cells (uniform, not stratified \u2014 use `fetch_cellxgene_census_cells`'s balanced downsampling if you need per-group balance).",
            },
            {
                "name": "random_seed",
                "type": "int",
                "default": "0",
                "description": "Seed for the downsample.",
            },
            {
                "name": "timeout",
                "type": "int",
                "default": "600",
                "description": "Download timeout in seconds.",
            },
        ],
    },
    {
        "name": "qc_normalize_cluster_umap",
        "description": "QC filter, normalize, batch-correct, cluster, and UMAP-embed an AnnData.",
        "required_parameters": [
            {
                "name": "input_h5ad",
                "type": "str",
                "description": "Path to raw AnnData (expects a `feature_name` var column if var_names are not already gene symbols, and raw counts in .X).",
            },
        ],
        "optional_parameters": [
            {
                "name": "output_h5ad",
                "type": "str",
                "default": "'processed.h5ad'",
                "description": "Path to write the processed AnnData. Adds/uses: layers[\"counts\"], layers[\"lognorm\"], obsm[\"X_pca\"], obsm[\"X_pca_harmony\"] (if batch_key set), obsm[\"X_umap\"], obs[\"leiden\"]. Sets .raw to the full log-normalized matrix.",
            },
            {
                "name": "batch_key",
                "type": "str",
                "default": "'assay'",
                "description": "obs column for Harmony batch correction; set to None to skip batch correction and cluster directly on PCA space.",
            },
            {
                "name": "min_genes",
                "type": "int",
                "default": "200",
                "description": "Minimum genes/cell (sc.pp.filter_cells).",
            },
            {
                "name": "min_cells",
                "type": "int",
                "default": "3",
                "description": "Minimum cells/gene (sc.pp.filter_genes).",
            },
            {
                "name": "max_total_counts_quantile",
                "type": "float",
                "default": "0.995",
                "description": "Upper total-count quantile cutoff.",
            },
            {
                "name": "drop_assay_values",
                "type": "List[str]",
                "default": "('Smart-seq2',)",
                "description": "obs[\"assay\"] values to exclude entirely (needs an \"assay\" column); pass None to keep all assays.",
            },
            {
                "name": "n_top_genes",
                "type": "int",
                "default": "2000",
                "description": "Number of highly variable genes for PCA.",
            },
            {
                "name": "n_pcs",
                "type": "int",
                "default": "50",
                "description": "Number of principal components.",
            },
            {
                "name": "n_neighbors",
                "type": "int",
                "default": "15",
                "description": "Neighbors for the kNN graph feeding UMAP/Leiden.",
            },
            {
                "name": "leiden_resolution",
                "type": "float",
                "default": "1.0",
                "description": "Leiden clustering resolution.",
            },
        ],
    },
    {
        "name": "annotate_clusters_by_label_map",
        "description": "Map cluster IDs (e.g. Leiden) to manually curated biological subset labels.",
        "required_parameters": [
            {
                "name": "input_h5ad",
                "type": "str",
                "description": "Path to a clustered AnnData (must have `cluster_key` in .obs).",
            },
            {
                "name": "output_h5ad",
                "type": "str",
                "description": "Path to write the annotated AnnData.",
            },
            {
                "name": "cluster_key",
                "type": "str",
                "description": "obs column holding cluster IDs, e.g. \"leiden\".",
            },
            {
                "name": "label_map",
                "type": "Dict[str, str]",
                "description": "Mapping from cluster ID (str) to subset label, e.g. {\"3\": \"Regulatory T (Treg)\", \"10\": \"CD8 Exhausted (CXCL13+)\"}.",
            },
        ],
        "optional_parameters": [
            {
                "name": "low_quality_labels",
                "type": "List[str]",
                "default": "None",
                "description": "Subset of `label_map` VALUES (or cluster-id keys) to flag as low quality in a new obs[\"is_low_quality\"] boolean column, e.g. clusters driven by ambient RNA, stress response, or doublets. If None, obs[\"is_low_quality\"] is all False.",
            },
            {
                "name": "new_col",
                "type": "str",
                "default": "'subset'",
                "description": "Name of the obs column to write the mapped label into.",
            },
        ],
    },
    {
        "name": "donor_level_composition_stats",
        "description": "Compute per-donor subset fractions and test them across condition groups.",
        "required_parameters": [
            {
                "name": "input_h5ad",
                "type": "str",
                "description": "Path to an annotated AnnData (needs `donor_col`, `condition_col`, `subset_col` in .obs).",
            },
            {
                "name": "donor_col",
                "type": "str",
                "description": "obs column identifying the donor/sample of origin.",
            },
            {
                "name": "condition_col",
                "type": "str",
                "description": "obs column identifying the comparison groups (e.g. \"disease\").",
            },
            {
                "name": "subset_col",
                "type": "str",
                "description": "obs column with the cell-subset / cluster label.",
            },
            {
                "name": "subsets_of_interest",
                "type": "List[str]",
                "description": "Subset labels to compute per-donor % for, e.g. [\"CD8 Exhausted/Tumor-reactive (CXCL13+GZMB+)\", \"Regulatory T (Treg)\"].",
            },
        ],
        "optional_parameters": [
            {
                "name": "exclude_low_quality",
                "type": "bool",
                "default": "True",
                "description": "If True and obs[\"is_low_quality\"] exists, drop those cells first.",
            },
            {
                "name": "min_cells_per_donor",
                "type": "int",
                "default": "10",
                "description": "Minimum cells for a donor to be included (donors with fewer cells give unstable fraction estimates).",
            },
            {
                "name": "reference_condition",
                "type": "str",
                "default": "None",
                "description": "If given, this condition's donors are compared against all other conditions pooled via Mann-Whitney U (two-sided) for each subset in `subsets_of_interest`. If None, no statistical test is run \u2014 only the composition table is written.",
            },
            {
                "name": "output_csv",
                "type": "str",
                "default": "'donor_level_composition.csv'",
                "description": "Path to write the donor-level composition table.",
            },
        ],
    },
    {
        "name": "plot_umap_by_group",
        "description": "Scatter-plot an AnnData's UMAP embedding colored by a categorical obs column.",
        "required_parameters": [
            {
                "name": "input_h5ad",
                "type": "str",
                "description": "Path to an AnnData with obsm[\"X_umap\"] already computed (see `qc_normalize_cluster_umap`).",
            },
            {
                "name": "color_col",
                "type": "str",
                "description": "obs column to color points by (e.g. \"subset\", \"leiden\", \"disease\").",
            },
        ],
        "optional_parameters": [
            {
                "name": "output_png",
                "type": "str",
                "default": "'umap_by_group.png'",
                "description": "Path to save the figure.",
            },
            {
                "name": "exclude_low_quality",
                "type": "bool",
                "default": "True",
                "description": "If True and obs[\"is_low_quality\"] exists, drop those cells first.",
            },
            {
                "name": "point_size",
                "type": "float",
                "default": "3.0",
                "description": "Marker size passed to matplotlib scatter.",
            },
            {
                "name": "title",
                "type": "str",
                "default": "None",
                "description": "Plot title; defaults to \"UMAP colored by <color_col>\".",
            },
        ],
    },
    {
        "name": "plot_marker_gene_heatmap",
        "description": "Plot a mean-marker-expression heatmap across groups (e.g. cell subsets).",
        "required_parameters": [
            {
                "name": "input_h5ad",
                "type": "str",
                "description": "Path to an AnnData with .raw set to log-normalized values (see `qc_normalize_cluster_umap`, which sets .raw automatically).",
            },
            {
                "name": "marker_genes",
                "type": "List[str]",
                "description": "Gene symbols to plot as heatmap columns; any not found in adata.raw.var_names are silently skipped.",
            },
            {
                "name": "groupby",
                "type": "str",
                "description": "obs column defining heatmap rows (e.g. \"subset\").",
            },
        ],
        "optional_parameters": [
            {
                "name": "output_png",
                "type": "str",
                "default": "'marker_heatmap.png'",
                "description": "Path to save the figure.",
            },
            {
                "name": "exclude_low_quality",
                "type": "bool",
                "default": "True",
                "description": "If True and obs[\"is_low_quality\"] exists, drop those cells first.",
            },
            {
                "name": "zscore",
                "type": "bool",
                "default": "True",
                "description": "If True, z-score each gene's mean expression across groups (highlights relative enrichment); if False, plot raw mean log-normalized expression.",
            },
        ],
    },
    {
        "name": "plot_marker_gene_dotplot",
        "description": "Dot plot of marker genes across groups: dot size = % expressing, color = mean expression.",
        "required_parameters": [
            {
                "name": "input_h5ad",
                "type": "str",
                "description": "Path to an AnnData with .raw set to log-normalized values.",
            },
            {
                "name": "marker_genes",
                "type": "List[str]",
                "description": "Gene symbols to plot; any not found are silently skipped.",
            },
            {
                "name": "groupby",
                "type": "str",
                "description": "obs column defining dot-plot rows (e.g. \"subset\").",
            },
        ],
        "optional_parameters": [
            {
                "name": "output_png",
                "type": "str",
                "default": "'marker_dotplot.png'",
                "description": "Path to save the figure.",
            },
            {
                "name": "exclude_low_quality",
                "type": "bool",
                "default": "True",
                "description": "If True and obs[\"is_low_quality\"] exists, drop those cells first.",
            },
        ],
    },
    {
        "name": "analyze_gene_expression_by_group",
        "description": "Summarize one gene's expression across groups, with an optional donor-level check.",
        "required_parameters": [
            {
                "name": "input_h5ad",
                "type": "str",
                "description": "Path to an AnnData with .raw / a \"lognorm\" layer set.",
            },
            {
                "name": "gene",
                "type": "str",
                "description": "Gene symbol to analyze.",
            },
            {
                "name": "groupby",
                "type": "str",
                "description": "obs column defining the groups (e.g. \"subset\").",
            },
        ],
        "optional_parameters": [
            {
                "name": "donor_col",
                "type": "str",
                "default": "None",
                "description": "obs column identifying donors, for the donor-level check; skip donor-level validation if None.",
            },
            {
                "name": "min_cells_per_donor",
                "type": "int",
                "default": "5",
                "description": "Minimum cells per donor to include in the donor-level check.",
            },
            {
                "name": "exclude_low_quality",
                "type": "bool",
                "default": "True",
                "description": "If True and obs[\"is_low_quality\"] exists, drop those cells first.",
            },
            {
                "name": "output_csv",
                "type": "str",
                "default": "'gene_expression_by_group.csv'",
                "description": "Path to write the per-group summary table.",
            },
            {
                "name": "output_png",
                "type": "str",
                "default": "'gene_expression_by_group.png'",
                "description": "Path to save the figure; pass None to skip plotting.",
            },
        ],
    },
    {
        "name": "compare_clusters_to_reference_labels",
        "description": "Quantify concordance between this analysis' cluster labels and a reference atlas's labels.",
        "required_parameters": [
            {
                "name": "own_h5ad",
                "type": "str",
                "description": "Path to this analysis' annotated AnnData.",
            },
            {
                "name": "reference_h5ad",
                "type": "str",
                "description": "Path to the reference/published AnnData (must share `join_key` values with `own_h5ad`, typically a subset of a much larger reference \u2014 only overlapping cells are compared).",
            },
            {
                "name": "own_label_col",
                "type": "str",
                "description": "obs column in `own_h5ad` with this analysis' cluster or subset labels.",
            },
            {
                "name": "reference_label_col",
                "type": "str",
                "description": "obs column in `reference_h5ad` with the reference author labels.",
            },
            {
                "name": "join_key",
                "type": "str",
                "description": "obs column present in both, used to match cells 1:1.",
            },
        ],
        "optional_parameters": [
            {
                "name": "output_crosstab_csv",
                "type": "str",
                "default": "'label_crosstab.csv'",
                "description": "Path to write the full label x label cross-tab.",
            },
            {
                "name": "output_mapping_csv",
                "type": "str",
                "default": "'label_mapping.csv'",
                "description": "Path to write, per own-label, the majority reference label and purity.",
            },
        ],
    },
]