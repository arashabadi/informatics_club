# How Rich et al. Aligned Seurat and Scanpy

This note summarizes where Rich et al. explain how they made Seurat and Scanpy as similar as possible, and what exact changes they used in the paper and Zenodo code.

## Short Answer

They aligned Seurat and Scanpy in two different ways:

1. **Matched function arguments**: one package was run with non-default settings chosen to imitate the other package.
2. **Controlled identical input after each step**: after filtering, HVG selection, PCA, graph construction, or clustering, they copied the intermediate output from one package into the other so the next step started from the same object.

The strongest "most similar" comparisons are therefore not Figure S2 alone. Figure S2 only aligns arguments. Figures S4 and S5 are the strongest controls because they combine matched arguments with identical upstream input at each step.

## Where They Explain This

Main paper:

- Results around Figure 1: describes why the default workflows differ and how matching arguments changes the result.
- Table 1: classifies steps as equivalent by default, equivalent with matched arguments, or incompatible.
- STAR Methods / Quantification and Statistical Analysis: describes the controlled-input strategy and the similarity metrics.
- Supplemental Figures 2-5 captions: define the specific alignment experiments.
- Supplemental Tables 1-2: list the Seurat and Scanpy function arguments used to match the other package.

Code:

- `paper/zenodo-codes_pachterlab/pachterlab-RMEJLBASBMP_2024-1de630a/analysis/rmd/Seurat_v_Scanpy.Rmd`
- `paper/zenodo-codes_pachterlab/pachterlab-RMEJLBASBMP_2024-1de630a/analysis/yaml/Supp_Fig2.yaml`
- `paper/zenodo-codes_pachterlab/pachterlab-RMEJLBASBMP_2024-1de630a/analysis/yaml/Supp_Fig4.yaml`
- `paper/zenodo-codes_pachterlab/pachterlab-RMEJLBASBMP_2024-1de630a/analysis/yaml/Supp_Fig5.yaml`
- `paper/zenodo-codes_pachterlab/pachterlab-RMEJLBASBMP_2024-1de630a/analysis/yaml/Fig4_Supp_Fig3_Supp_Fig15_Supp_Fig16.yaml`

## The Four Key Alignment Experiments

### Figure S2: Seurat-like Scanpy Arguments Only

YAML:

```yaml
analysis_methods: "seurat_like"
data_input: "default"
```

Meaning:

- Scanpy is run with Seurat-like settings.
- The workflow is still sequential: each package carries forward its own results.
- This tests how much of the difference is due only to different function arguments.

### Figure S3: Default Arguments With Controlled Seurat Input

YAML:

```yaml
analysis_methods: "default"
data_input: "seurat"
```

Meaning:

- Seurat and Scanpy keep default-like settings.
- But after each step, the data created by Scanpy is replaced by the data created by Seurat.
- This isolates the difference introduced by each individual step when the input to that step is identical.

### Figure S4: Seurat-like Arguments Plus Seurat-Controlled Input

YAML:

```yaml
analysis_methods: "seurat_like"
data_input: "seurat"
```

Meaning:

- Scanpy is made Seurat-like.
- Scanpy also receives Seurat's upstream objects after each relevant step.
- This is the strongest Seurat-direction alignment.

### Figure S5: Scanpy-like Arguments Plus Scanpy-Controlled Input

YAML:

```yaml
analysis_methods: "scanpy_like"
data_input: "scanpy"
```

Meaning:

- Seurat is made Scanpy-like.
- Seurat receives Scanpy's upstream objects after each relevant step.
- This is the strongest Scanpy-direction alignment.

## Exact Settings By Step

### 1. Filtering

They used the same general filtering logic in both packages:

- Same input count matrix when comparing Seurat and Scanpy directly.
- Same UMI knee/inflection cutoff.
- `min_cells = 3`
- `min_features = 200`
- Mitochondrial filtering: `pct_mt < 20`
- Scanpy gene-count upper filter was relaxed to `max_n_genes_by_counts_scanpy = 12000` instead of the usual smaller Scanpy tutorial threshold.

In controlled-input runs:

- `data_input: "seurat"` means Scanpy is subset to Seurat's cells and genes.
- `data_input: "scanpy"` means Seurat is subset to Scanpy's cells and genes.

### 2. Normalization

They used log normalization in both packages.

Important point:

- Seurat default scale factor is `10000`.
- Scanpy's function default can differ, but the Scanpy tutorial they followed used `target_sum = 1e4`.
- Therefore, in most analyses, Scanpy was run with `target_sum = 1e4`.

With the same input matrix and the same scale factor, Seurat and Scanpy log normalization became identical.

### 3. HVG Selection

This is one of the central differences.

Seurat default:

```r
selection.method = "vst"
nfeatures = 2000
mean.cutoff = c(0.1, 8)
dispersion.cutoff = c(1, Inf)
```

Scanpy default-like:

```python
flavor = "seurat"
n_top_genes = None
min_mean = 0.0125
max_mean = 3
min_disp = 0.5
```

To make Scanpy Seurat-like:

```python
flavor = "seurat_v3"
n_top_genes = 2000
```

To make Seurat Scanpy-like:

```r
selection.method = "mean.var.plot"
mean.cutoff = c(0.0125, 3)
dispersion.cutoff = c(0.5, Inf)
nfeatures = Inf
```

In controlled-input runs:

- If `data_input == "seurat"`, the Seurat HVG list is written into Scanpy.
- If `data_input == "scanpy"`, the Scanpy HVG list is written into Seurat.

### 4. Regression And Scaling

Seurat default-like:

```r
vars.to.regress = NULL
scale.max = 10
```

Scanpy default-like:

```python
sc.pp.regress_out(adata, ["total_counts", "pct_mt"])
sc.pp.scale(adata, max_value=None)
```

To make Scanpy Seurat-like:

```python
# omit sc.pp.regress_out
max_value = 10
```

To make Seurat Scanpy-like:

```r
vars.to.regress = c("nCount_RNA", "pct_mt")
scale.max = Inf
```

This matters because regression and clipping change the matrix that PCA sees.

### 5. PCA

Both workflows use PCA after HVG selection and scaling.

Common settings in the YAML runs:

```yaml
seu_num_pcs: 50
scan_num_pcs: 50
pca_seed_seu: 42
pca_seed_scan: 0
```

To make Scanpy Seurat-like:

```python
zero_center = False
```

To make Seurat Scanpy-like:

- Use the Scanpy-like upstream HVGs, regression, and scaling.
- In controlled-input mode, copy Scanpy PCA embeddings into the Seurat object.

In controlled-input runs:

- Seurat PCA embeddings can be copied to `adata.obsm["X_pca"]`.
- Scanpy PCA embeddings can be copied into Seurat as a Seurat PCA reduction.

### 6. KNN/SNN Graph

Seurat default-like:

```r
k.param = 20
```

Scanpy default-like:

```python
n_neighbors = 15
```

To make Scanpy Seurat-like:

```python
n_neighbors = 20
use_rep = "X_pca"
```

To make Seurat Scanpy-like:

```r
k.param = 15
```

In controlled-input runs:

- If `data_input == "seurat"`, the Seurat SNN graph is copied into Scanpy as `adata.obsp["connectivities"]`.
- If `data_input == "scanpy"`, the Scanpy SNN graph is copied into Seurat as `seu@graphs$RNA_snn`.

Important result:

- Matching graph arguments helped only modestly.
- The SNN graph remained one of the hardest steps to reconcile, because Seurat and Scanpy build nearest-neighbor/SNN graphs differently.

### 7. Clustering

Seurat default-like:

```r
algorithm = 1   # Louvain
resolution = 0.8
```

Scanpy default-like:

```python
sc.tl.leiden(...)
resolution = 1
n_iterations = -1
```

To make Scanpy Seurat-like:

```python
sc.tl.louvain(...)
resolution = 0.8
```

To make Seurat Scanpy-like:

```r
algorithm = 4   # Leiden
resolution = 1
```

Important result:

- Leiden could be made identical when the graph/input was matched.
- Louvain still produced differences between Seurat and Scanpy even when inputs and arguments were aligned.
- In the Scanpy-like mode, the code often assigns Scanpy's Leiden clusters directly to Seurat to avoid Seurat's memory-intensive Leiden run.

### 8. UMAP

Seurat default-like:

```r
umap.method = "uwot"
min.dist = 0.3
metric = "cosine"
```

Scanpy default-like:

```python
sc.tl.umap(...)
min_dist = 0.5
```

To make Scanpy Seurat-like:

```python
min_dist = 0.3
```

To make Seurat Scanpy-like:

```r
umap.method = "umap-learn"
min.dist = 0.5
metric = "correlation"
```

Important result:

- UMAP remained partially irreconcilable.
- Even with matched upstream PCA/SNN/cluster inputs, `uwot` and `umap-learn` can preserve different local geometry.

### 9. Differential Expression

This is the most biologically important downstream step in the paper.

Seurat default-like DE:

```r
FindAllMarkers(seu)
```

Seurat uses:

- Wilcoxon test.
- Tie correction.
- Bonferroni adjusted p values.
- Marker filtering before/around testing.
- Filters include logFC threshold, fraction of cells expressing the gene, and p value threshold.

Scanpy default-like DE:

```python
sc.tl.rank_genes_groups(
    adata,
    groupby,
    method = "wilcoxon",
    corr_method = "benjamini-hochberg",
    tie_correct = False,
    pts = True
)
```

To make Scanpy Seurat-like:

```python
method = "wilcoxon"
corr_method = "bonferroni"
tie_correct = True
pts = True
```

Then they manually filtered Scanpy markers to mimic Seurat:

```r
filter(!(pts < 0.01 & pts_rest < 0.01))
filter(!(abs(log_fc) < 0.1))
filter(p_value < 0.01)
```

To make Seurat Scanpy-like:

```r
FindAllMarkers(
    seu,
    logfc.threshold = 0,
    min.pct = 0,
    return.thresh = 1.0001,
    min.cells.group = 1
)
```

Then they used Benjamini-Hochberg adjustment:

```r
p.adjust(..., method = "BH")
```

In controlled-input runs:

- They also copied cluster assignments so Seurat and Scanpy performed DE on the same clusters.

Important results:

- Matching Scanpy to Seurat DE settings improved significant marker overlap to a Jaccard index of about `0.73`.
- Providing the same cluster assignments improved significant marker overlap to about `0.99`.
- The remaining difference was mostly attributed to logFC calculation differences.
- Matching Seurat to Scanpy DE settings was worse, about `0.38`, because Seurat/presto could not turn off Wilcoxon tie correction.

## What Could Be Fully Reconciled

These parts could be made identical or nearly identical:

- Cell and gene filtering, when thresholds and inputs were matched.
- Log normalization, when the same matrix and same scale factor were used.
- HVG selection, when the corresponding algorithm was selected or the HVG list was copied.
- Leiden clustering, when the graph/input was matched.
- Seurat-like DE p values, when Scanpy used tie correction and Bonferroni correction.

## What Could Not Be Fully Reconciled

These remained difficult or impossible to perfectly align:

- SNN graph construction, because the graph-building implementations differ.
- Louvain clustering, even with matched graph/input.
- UMAP geometry, because Seurat and Scanpy use different UMAP implementations and settings.
- Marker filtering behavior, especially when trying to make Seurat behave exactly like Scanpy.
- logFC calculation, because Seurat and Scanpy compute logFC differently.
- Scanpy-like DE p values in Seurat, because Seurat/presto cannot disable Wilcoxon tie correction.

## Practical Interpretation

The paper's message is not simply "Seurat and Scanpy are different." It is more precise:

- Some differences are just defaults and can be reduced by matching arguments.
- Some differences are upstream propagation: early choices in HVGs, scaling, PCA, or graphs change everything downstream.
- Some differences are genuine implementation differences and remain even after careful matching.
- The choice of package, version, and function arguments is part of the scientific method, not just a software detail.

For journal club, the most important explanation is:

> Figure S2 shows that matching arguments makes the workflows more similar, but Figures S4 and S5 show the stricter test: even with matched arguments and identical step-by-step inputs, some parts of the pipeline remain incompatible. Therefore, the main conclusion is that defaults explain part of the problem, but not all of it.

## If I Only Have Scanpy: How To Get The Most Seurat-Like Results

Question:

> If I do not have intermediate outputs from Seurat, what can I do to get the most similar results from Scanpy? What should I do at every step for a new dataset that is only being run in Scanpy?

Short answer:

You cannot reproduce the strongest controlled-input analyses from Figures S4 and S5 without Seurat intermediate outputs. Those figures rely on copying Seurat's cells, genes, HVGs, PCA embeddings, SNN graph, and cluster labels into Scanpy after each step. If you are only running Scanpy, the closest practical strategy is to run a **Seurat-like Scanpy workflow**: use the same count matrix, Seurat-like QC thresholds, Seurat-like HVG selection, no Scanpy regression, Seurat-like scaling/clipping, Seurat-like graph and clustering settings, and Seurat-like differential expression settings.

This is closest to the Figure S2 / Figure S4 direction: **make Scanpy imitate Seurat**.

## Practical Scanpy-Only Protocol For A New Dataset

### 0. Freeze The Analysis Environment

Before running anything, record package versions and random seeds. This matters because the paper shows that package versions alone can change DE results.

```python
import scanpy as sc
import numpy as np
import pandas as pd
import scipy

print("scanpy", sc.__version__)
print("numpy", np.__version__)
print("scipy", scipy.__version__)

RANDOM_STATE = 0
np.random.seed(RANDOM_STATE)
```

Best practice:

- Save the raw count matrix path.
- Save the exact Scanpy version.
- Save the exact filtering thresholds.
- Save the final `.h5ad`.
- Save intermediate `.h5ad` files after QC, HVG selection, PCA, neighbors, clustering, UMAP, and DE.

### 1. Start From The Same Kind Of Raw Count Matrix

Use a raw cell-by-gene count matrix, not already normalized data. The paper's Seurat-vs-Scanpy comparison used the same count matrix for both tools when possible.

Example:

```python
adata = sc.read_10x_mtx(
    "path/to/filtered_or_raw_feature_bc_matrix",
    var_names="gene_ids",
    cache=True,
)

adata.var_names_make_unique()
adata.layers["counts"] = adata.X.copy()
```

If you use `kb`, Cell Ranger, STARsolo, or another matrix generator, record that too. Matrix generation can change downstream results.

### 2. Use Seurat-Like Cell And Gene Filtering

Use the same filtering logic as the paper where possible:

- Minimum genes per cell: `200`
- Minimum cells per gene: `3`
- Mitochondrial percent cutoff: `20`
- Avoid Scanpy tutorial's stricter `n_genes_by_counts < 2500` unless you intentionally want that.

```python
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

# Human mitochondrial genes by symbol. If using Ensembl IDs, map or flag MT genes carefully.
adata.var["mt"] = adata.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(
    adata,
    qc_vars=["mt"],
    percent_top=None,
    log1p=False,
    inplace=True,
)

adata = adata[adata.obs["pct_counts_mt"] < 20, :].copy()

# Optional only if you want to mimic the paper's relaxed Scanpy upper gene filter.
# This is much looser than the usual tutorial threshold.
if "n_genes_by_counts" in adata.obs:
    adata = adata[adata.obs["n_genes_by_counts"] < 12000, :].copy()

adata.write("01_after_qc.h5ad")
```

Important:

- If your dataset is not human, mitochondrial gene names may not start with `MT-`.
- For mouse, often use `mt-`.
- If using Ensembl IDs, create the mitochondrial gene list explicitly.

### 3. Select HVGs With Seurat's Default-Like Method

To make Scanpy most Seurat-like, use Scanpy's implementation of Seurat v3-style HVG selection:

```python
sc.pp.highly_variable_genes(
    adata,
    flavor="seurat_v3",
    n_top_genes=2000,
    layer="counts",
)

adata.write("02_after_hvg_selection.h5ad")
```

This step is one of the most important alignments. In the paper:

- Seurat default HVG method is `vst`.
- Scanpy's closest equivalent is `flavor="seurat_v3"`.
- Use `n_top_genes=2000` to match Seurat's usual default.

Note:

- `flavor="seurat_v3"` may require `scikit-misc`.
- This method expects raw counts, so use `layer="counts"` if you stored raw counts.

### 4. Normalize And Log Transform With Seurat-Like Scale

Use total-count normalization to `1e4`, then log1p.

```python
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

adata.raw = adata
adata.write("03_after_normalization_log1p.h5ad")
```

Why:

- Seurat's `NormalizeData` default uses a scale factor of `10000`.
- The paper states that with the same input matrix and this scale factor, Seurat and Scanpy log normalization are effectively identical.

### 5. Keep Only HVGs For PCA

Subset to the Seurat-like HVGs.

```python
adata = adata[:, adata.var["highly_variable"]].copy()
adata.write("04_hvg_subset.h5ad")
```

This matches the usual workflow where PCA is run on the selected variable genes.

### 6. Do Not Regress Out Total Counts Or Mito Percent

For the most Seurat-like Scanpy workflow, **skip** this Scanpy tutorial step:

```python
# Do not do this for Seurat-like Scanpy:
# sc.pp.regress_out(adata, ["total_counts", "pct_counts_mt"])
```

Why:

- Seurat's default in the paper used no regression.
- Scanpy tutorial/default-style workflows often regress total counts and mitochondrial percentage.
- This difference affects scaling, PCA, neighbors, clustering, UMAP, and DE.

### 7. Scale With Seurat-Like Clipping

Use `max_value=10`, matching Seurat's default-like scaling behavior.

```python
sc.pp.scale(adata, max_value=10)
adata.write("05_after_scaling.h5ad")
```

Why:

- Seurat clips scaled values at `10`.
- Scanpy default-like behavior can leave values unclipped.

### 8. Run PCA With A Fixed PC Plan

Use a fixed number of PCs and record it. In Rich et al.'s code, many YAML runs use `50` PCs. If you want to mimic a small Seurat tutorial-style analysis, `10` PCs may be closer to some Seurat defaults, but for this paper's workflow use `50`.

```python
N_PCS = 50

sc.tl.pca(
    adata,
    n_comps=N_PCS,
    svd_solver="arpack",
    zero_center=False,
    random_state=RANDOM_STATE,
)

adata.write("06_after_pca.h5ad")
```

Why:

- The paper's Seurat-like Scanpy setting used `zero_center=False`.
- The selected number of PCs should be fixed and reported.
- Changing the number of PCs can change graph construction and clustering.

### 9. Build A Seurat-Like Neighbor Graph

Use Seurat-like neighbor count:

```python
sc.pp.neighbors(
    adata,
    n_neighbors=20,
    n_pcs=N_PCS,
    use_rep="X_pca",
    random_state=RANDOM_STATE,
)

adata.write("07_after_neighbors.h5ad")
```

Why:

- Seurat default-like `k.param` is `20`.
- Scanpy default-like `n_neighbors` is often `15`.
- Matching this parameter helps, although the paper shows SNN graph construction is still not fully reconcilable.

Important limitation:

- This will not make Scanpy's SNN graph identical to Seurat's.
- Without a Seurat SNN graph to copy, this remains an approximation.

### 10. Cluster With Seurat-Like Louvain Settings

To mimic Seurat's default-like clustering:

```python
sc.tl.louvain(
    adata,
    resolution=0.8,
    random_state=RANDOM_STATE,
    key_added="louvain_seurat_like",
)

adata.obs["seurat_like_clusters"] = adata.obs["louvain_seurat_like"]
adata.write("08_after_clustering.h5ad")
```

Why:

- Seurat default-like clustering in the paper used Louvain with resolution `0.8`.
- Scanpy default-like workflows often use Leiden with resolution `1`.

Important limitation:

- The paper found Louvain remains partly incompatible between Seurat and Scanpy.
- If your goal is biological stability rather than Seurat mimicry, Leiden may be preferable.
- If your goal is maximum Seurat-like behavior, use Louvain `resolution=0.8`.

Optional robust alternative:

```python
sc.tl.leiden(
    adata,
    resolution=0.8,
    random_state=RANDOM_STATE,
    key_added="leiden_resolution_0_8",
)
```

### 11. Run UMAP With Seurat-Like Minimum Distance

Use a Seurat-like `min_dist=0.3`.

```python
sc.tl.umap(
    adata,
    min_dist=0.3,
    random_state=RANDOM_STATE,
)

adata.write("09_after_umap.h5ad")
```

Why:

- Seurat default-like UMAP uses `min.dist=0.3`.
- Scanpy default-like UMAP uses `min_dist=0.5`.

Important limitation:

- Scanpy uses `umap-learn`.
- Seurat often uses `uwot`.
- So even with matched arguments, UMAP geometry may not be identical.

### 12. Run Seurat-Like Differential Expression In Scanpy

Use Wilcoxon, Bonferroni correction, tie correction, and return expression fractions.

```python
sc.tl.rank_genes_groups(
    adata,
    groupby="seurat_like_clusters",
    method="wilcoxon",
    corr_method="bonferroni",
    tie_correct=True,
    pts=True,
    use_raw=True,
)

adata.write("10_after_de.h5ad")
```

Then convert the results to a table and apply Seurat-like marker filtering.

```python
de = sc.get.rank_genes_groups_df(adata, group=None)

# Scanpy column names can vary slightly by version. These are common names:
# group, names, scores, logfoldchanges, pvals, pvals_adj, pct_nz_group, pct_nz_reference

de = de.rename(columns={
    "names": "gene",
    "logfoldchanges": "log_fc",
    "pvals": "p_value",
    "pvals_adj": "p_value_adj",
    "pct_nz_group": "pts",
    "pct_nz_reference": "pts_rest",
})

seurat_like_markers = de[
    ~((de["pts"] < 0.01) & (de["pts_rest"] < 0.01))
].copy()

seurat_like_markers = seurat_like_markers[
    seurat_like_markers["log_fc"].abs() >= 0.1
].copy()

seurat_like_markers = seurat_like_markers[
    seurat_like_markers["p_value"] < 0.01
].copy()

seurat_like_markers.to_csv("seurat_like_scanpy_markers.csv", index=False)
```

Why:

- Seurat default-like DE uses Wilcoxon with tie correction.
- Seurat default-like multiple testing uses Bonferroni correction.
- Seurat filters marker candidates by expression fraction, logFC, and p value.
- Scanpy does not naturally reproduce all of Seurat's marker filtering unless you add it manually.

Important limitation:

- logFC values may still differ from Seurat.
- The paper emphasizes that Seurat and Scanpy compute logFC differently.
- Therefore, use marker overlap and biological interpretation carefully.

### 13. Save A Reproducibility Manifest

Create a small text or YAML file with all key settings.

Example:

```yaml
workflow: scanpy_seurat_like
input_matrix: path/to/matrix
scanpy_version: record_from_runtime
filtering:
  min_genes_per_cell: 200
  min_cells_per_gene: 3
  pct_mito_max: 20
  max_n_genes_by_counts: 12000
hvg:
  flavor: seurat_v3
  n_top_genes: 2000
normalization:
  target_sum: 10000
  log1p: true
regression:
  regress_total_counts: false
  regress_mito_percent: false
scaling:
  max_value: 10
pca:
  n_pcs: 50
  zero_center: false
  svd_solver: arpack
neighbors:
  n_neighbors: 20
  n_pcs: 50
clustering:
  method: louvain
  resolution: 0.8
umap:
  min_dist: 0.3
de:
  method: wilcoxon
  corr_method: bonferroni
  tie_correct: true
  pts: true
  marker_filtering:
    abs_log_fc_min: 0.1
    p_value_max: 0.01
    min_fraction_rule: remove_if_pts_and_pts_rest_both_below_0.01
```

## Minimal Scanpy Code Skeleton

```python
import scanpy as sc
import numpy as np

RANDOM_STATE = 0
N_PCS = 50

adata = sc.read_10x_mtx("path/to/matrix", var_names="gene_ids", cache=True)
adata.var_names_make_unique()
adata.layers["counts"] = adata.X.copy()

# QC
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
adata.var["mt"] = adata.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)
adata = adata[adata.obs["pct_counts_mt"] < 20, :].copy()
adata = adata[adata.obs["n_genes_by_counts"] < 12000, :].copy()

# Seurat-like HVGs
sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=2000, layer="counts")

# Normalize/log
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata

# PCA input
adata = adata[:, adata.var["highly_variable"]].copy()

# No regression for Seurat-like Scanpy
sc.pp.scale(adata, max_value=10)

# PCA
sc.tl.pca(adata, n_comps=N_PCS, svd_solver="arpack", zero_center=False, random_state=RANDOM_STATE)

# Graph
sc.pp.neighbors(adata, n_neighbors=20, n_pcs=N_PCS, use_rep="X_pca", random_state=RANDOM_STATE)

# Seurat-like clustering
sc.tl.louvain(adata, resolution=0.8, random_state=RANDOM_STATE, key_added="seurat_like_clusters")

# UMAP
sc.tl.umap(adata, min_dist=0.3, random_state=RANDOM_STATE)

# DE
sc.tl.rank_genes_groups(
    adata,
    groupby="seurat_like_clusters",
    method="wilcoxon",
    corr_method="bonferroni",
    tie_correct=True,
    pts=True,
    use_raw=True,
)

adata.write("scanpy_seurat_like_final.h5ad")
```

## What To Say In A Presentation

If someone asks, "Can we get Seurat results from Scanpy if we do not have Seurat intermediate objects?", the careful answer is:

> Not exactly. Without Seurat intermediate outputs, we cannot perform the controlled-input alignment used in Figures S4 and S5. But we can make Scanpy as Seurat-like as possible by matching the major Seurat defaults: Seurat v3-style HVG selection, no regression, clipping scaled values at 10, 20 neighbors, Louvain at resolution 0.8, UMAP min_dist 0.3, and Wilcoxon DE with tie correction, Bonferroni correction, and Seurat-like marker filtering. This reduces many default-driven differences, but it cannot eliminate implementation-level differences in graph construction, Louvain clustering, UMAP, and logFC calculation.
