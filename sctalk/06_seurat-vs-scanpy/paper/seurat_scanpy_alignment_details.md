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

