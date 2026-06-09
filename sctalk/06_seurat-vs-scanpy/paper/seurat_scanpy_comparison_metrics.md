# Rich et al. Seurat-vs-Scanpy comparison metrics

This note explains the metrics Rich et al. used to compare Seurat and Scanpy outputs, including the Fig. 2-style downsampling summaries, the formulas behind each metric, and where each metric is computed in the released Zenodo code.

## Code path aliases

The code paths are long, so this note uses these aliases:

| Alias | File |
|---|---|
| `SVS` | `zenodo-codes_pachterlab/pachterlab-RMEJLBASBMP_2024-1de630a/analysis/rmd/Seurat_v_Scanpy.Rmd` |
| `HELPER` | `zenodo-codes_pachterlab/pachterlab-RMEJLBASBMP_2024-1de630a/analysis/scripts/data_analysis_helper.R` |
| `PLOTSTATS` | `zenodo-codes_pachterlab/pachterlab-RMEJLBASBMP_2024-1de630a/analysis/scripts/plotting_and_stats.R` |
| `AGG` | `zenodo-codes_pachterlab/pachterlab-RMEJLBASBMP_2024-1de630a/analysis/rmd/aggregate_plots.Rmd` |
| `BOOT` | `zenodo-codes_pachterlab/pachterlab-RMEJLBASBMP_2024-1de630a/analysis/rmd/bootstrapped_plots.Rmd` |

Line numbers are from the local Zenodo code copy in this folder.

## Big picture

The paper compares the two workflows step by step. At each step, it asks a practical question: if Seurat and Scanpy start from comparable inputs, do they keep the same cells, genes, HVGs, PCs, graphs, clusters, UMAP neighborhoods, marker genes, logFC values, and adjusted p-values?

Most metrics fall into three families:

| Metric family | Used for | Meaning |
|---|---|---|
| Set overlap | cells, genes, HVGs, marker genes, marker gene-cluster pairs, volcano gene sets, graph neighborhoods | "Are the same items selected?" |
| Continuous-value agreement | PCA loading direction, logFC, adjusted p-values | "Are the numeric outputs close?" |
| Partition or label agreement | clustering, UMAP-derived clustering, CellTypist labels | "Do cells receive the same group labels?" |

## Fig. 2 / downsampling metric set

For the downsampling analysis, the authors reran each workflow after read or cell downsampling and compared the downsampled result with the full-size result within the same package. They then asked whether downsampling-induced differences were as large as the default Seurat-vs-Scanpy differences.

The metric list is defined in `AGG:106-230`.

| Step | Summary metric | Ideal value | Default Seurat-vs-Scanpy baseline in code | Code extraction phrase |
|---|---:|---:|---:|---|
| Cell filtering | Cells Jaccard | 1 | 1 | `Cells Jaccard:` |
| Gene filtering | Genes Jaccard | 1 | 1 | `Genes Jaccard:` |
| HVG selection | HVGs Jaccard | 1 | 0.22249151720795 | `HVGs Jaccard:` |
| PCA | Mean PC1-3 loading difference | 0 | 0.417599195296815 | `Mean loading difference of PC1-3:` |
| KNN/SNN graph | Median magnitude log2 SNN degree ratio | 0 | 2.05889368905357 | `Median magnitude of log degree ratio of SNN:` |
| Clustering | Adjusted Rand index | 1 | 0.706349552155871 | `Adjusted Rand index between clusters:` |
| UMAP | Median Jaccard of UMAP-derived KNN neighborhoods | 1 | 0.0638297872340425 | `Median jaccard of UMAP KNN:` |
| Significant marker genes | Marker Genes Jaccard | 1 | 0.615919763464549 | `Marker Genes Jaccard:` |
| All markers, same-input DE | Markers Jaccard | 1 | same-input baseline 0.222991819976808 | `Markers Jaccard:` |
| logFC, same-input DE | Concordance correlation coefficient | 1 | same-input baseline 0.980466856603784 | `logFC CCC:` |
| Adjusted p-value, same-input DE | Fraction crossing p = 0.05 differently | 0 | same-input baseline 0.2 | `Adjusted p value, fraction that flipped across 0.05 threshold:` |

For Fig. 2-style fraction summaries, the code averages each metric over seeds, draws the Seurat-vs-Scanpy baseline as a black dashed line, and finds the downsampling fraction where the within-package curve reaches the baseline within a margin. This is implemented in `AGG:280-318`, `AGG:334-340`, and `AGG:431-459`.

The margin logic is:

```text
if ideal_value == 0:
    target = baseline * (1 + margin)
else:
    target = baseline * (1 - margin)
```

The crossing point is interpolated by `find_intersection()` in `AGG:89-103`.

## Metric formulas and code locations

### 1. Jaccard index

Used for cells, genes, HVGs, significant marker genes, all marker gene-cluster pairs, volcano significant gene sets, SNN neighborhoods, and UMAP-derived KNN neighborhoods.

Formula:

```text
J(A, B) = |A intersect B| / |A union B|
```

Interpretation:

| Value | Meaning |
|---:|---|
| 1 | identical sets |
| 0 | no overlap |
| intermediate | partial overlap |

Code locations:

| Use | Code |
|---|---|
| General vector/list Jaccard helper | `HELPER:224-245` |
| Per-cell KNN/UMAP neighborhood Jaccard | `HELPER:248-258` |
| Euler/UpSet set Jaccard for cells, genes, HVGs, marker genes, markers | `PLOTSTATS:323-420` |
| Cell/gene Jaccard calls after QC | `SVS:552-563` |
| HVG Jaccard calls | `SVS:711-714` |
| Significant marker gene Jaccard calls | `SVS:1518-1540` |
| All marker gene-cluster pair Jaccard calls | `SVS:1711-1726` |
| Volcano gene-set Jaccards | `PLOTSTATS:2367-2424`, called at `SVS:1844` |

Important nuance: for significant marker genes, the set is unique genes across all clusters with adjusted p-value below 0.05. For all markers, the set is not just gene names; it is gene-cluster pairs, so the same gene in two clusters counts as two different marker entries.

### 2. Normalization output equality

Used to check whether the normalized count matrices match after the normalization step.

Formula:

```text
matrix_equal = all.equal(nonzero_values_matrix_1, nonzero_values_matrix_2)
```

This is not one of the main Fig. 2 summary metrics, but it is part of the stepwise workflow comparison.

Code locations:

| Use | Code |
|---|---|
| Seurat-vs-Scanpy normalization matrix equality | `SVS:655-686` |
| Related cell and gene order checks before normalization | `SVS:650-651` |

### 3. PCA loading direction difference

Used for PCA comparison. The paper shows PCA overlays and scree plots, but the Fig. 2-style numeric metric is the mean difference in corresponding PC loadings for PCs 1-3.

For a Seurat loading vector `u` and Scanpy loading vector `v`, the code computes sign-invariant cosine similarity:

```text
cos_sim = abs(sum(u * v) / sqrt(sum(u^2) * sum(v^2)))
```

Then it converts that to a sine-distance-like difference:

```text
pc_diff = sqrt(1 - cos_sim^2)
```

The Fig. 2 summary metric is:

```text
mean_pc1_3_diff = mean(pc_diff_PC1, pc_diff_PC2, pc_diff_PC3)
```

Interpretation:

| Value | Meaning |
|---:|---|
| 0 | same loading direction, allowing sign flips |
| closer to 1 | increasingly different loading directions |

Code locations:

| Use | Code |
|---|---|
| PC loading sine-difference helper | `HELPER:65-80` |
| Intersect genes and build PC-difference table | `HELPER:84-103` |
| PCA variance/scree calculations | `SVS:800-808` |
| Mean PC1-3 loading difference | `SVS:851-872` |
| Fig. 2 extraction phrase and baseline | `AGG:106-230` |

### 4. SNN graph neighborhood Jaccard

Used for graph-level comparison after KNN/SNN construction. For each cell, the packages have a set of graph neighbors. The paper compares overlap cell by cell.

Formula for cell `i`:

```text
J_i = |N_Seurat(i) intersect N_Scanpy(i)| / |N_Seurat(i) union N_Scanpy(i)|
```

Summary values include the median neighborhood Jaccard across cells.

Code locations:

| Use | Code |
|---|---|
| Convert graph matrix to neighbor list | `HELPER:211-220` |
| Per-cell Jaccards for graph neighbor lists | `HELPER:224-245` |
| SNN Jaccard analysis | `SVS:982-1034` |
| Graph Jaccard plotting | `PLOTSTATS:1086-1218` |

### 5. SNN degree ratio

Used with the SNN graph comparison. The degree of a cell is the number of graph neighbors/edges retained for that cell. The paper asks whether one workflow systematically gives each cell more graph connections.

Formula:

```text
D_i = degree_Seurat(i) / degree_Scanpy(i)
log_degree_ratio_i = log2(D_i)
summary = median(abs(log_degree_ratio_i))
```

Interpretation:

| Value | Meaning |
|---:|---|
| 0 | same graph degree |
| 1 | typical two-fold degree difference |
| 2 | typical four-fold degree difference |

Code locations:

| Use | Code |
|---|---|
| Degree ratio and log2 degree ratio | `SVS:1017-1025` |
| Printed summary metric | `SVS:1033-1034` |
| Fig. 2 extraction phrase and baseline | `AGG:106-230` |

The STAR Methods also note a useful upper bound: if two neighborhoods have degrees `d1 < d2`, the maximum possible Jaccard index is `d1 / d2`, because the largest possible intersection is `d1` and the smallest possible union is `d2`.

### 6. Adjusted Rand index

Used for clustering agreement and also for some downstream label comparisons. ARI compares two partitions of the same cells while correcting for chance agreement.

Conceptual formula:

```text
ARI = (RI - expected_RI) / (max_RI - expected_RI)
```

where `RI` is the Rand index, the fraction of cell pairs that are consistently placed together or apart in both clusterings.

Interpretation:

| Value | Meaning |
|---:|---|
| 1 | identical clustering |
| 0 | approximately chance-level agreement |
| below 0 | worse than chance under the ARI model |

Code locations:

| Use | Code |
|---|---|
| Main cluster ARI | `SVS:1144-1170` |
| UMAP-derived KNN Leiden cluster ARI | `SVS:1356-1366` |
| CellTypist label ARI | `SVS:1689` |
| Fig. 2 extraction phrase and baseline | `AGG:106-230` |

### 7. Alluvial cluster mapping

Used for visualizing how clusters correspond across packages. The alluvial plot itself is visual, while ARI is the numeric metric. The paper also describes a cluster-ordering algorithm that uses maximum overlap/Jaccard-style agreement so the plot is easier to read.

Core idea:

```text
For each cluster in group 2:
    find the group 1 cluster with the largest overlap
    reorder group 2 clusters near their best-matching group 1 cluster
```

Code locations:

| Use | Code |
|---|---|
| Build alluvial input table from paired labels | `HELPER:286-292` |
| Sort clusters by agreement | `HELPER:384-464` |
| Main alluvial calls | `SVS:1193-1207` |
| UMAP-derived alluvial calls | `SVS:1377-1389` |
| Plot implementation | `PLOTSTATS:1351-1423` |

### 8. UMAP-derived KNN Jaccard

Used to make UMAP comparison less subjective. Instead of only visually comparing UMAPs, they build a KNN graph in the two-dimensional UMAP coordinate space and compare each cell's UMAP-space neighbors.

Formula for cell `i`:

```text
J_i_UMAP = |KNN_UMAP_Seurat(i) intersect KNN_UMAP_Scanpy(i)| /
           |KNN_UMAP_Seurat(i) union KNN_UMAP_Scanpy(i)|
summary = median_i(J_i_UMAP)
```

Code locations:

| Use | Code |
|---|---|
| Extract UMAP embeddings | `SVS:1284-1285` |
| Build exact KNN in UMAP space | `SVS:1312-1313` |
| Calculate per-cell UMAP KNN Jaccards | `SVS:1318-1326` and `HELPER:248-258` |
| Leiden clustering on UMAP-derived KNN graph | `SVS:1348-1350` |
| UMAP-derived cluster ARI | `SVS:1360-1366` |
| Fig. 2 extraction phrase and baseline | `AGG:106-230` |

### 9. Significant marker gene Jaccard

Used for differential expression/marker selection. This metric compares the set of genes that are significant in at least one cluster.

Formula:

```text
S_Seurat = unique genes with adjusted p < 0.05 in Seurat
S_Scanpy = unique genes with adjusted p < 0.05 in Scanpy
J_marker_genes = |S_Seurat intersect S_Scanpy| / |S_Seurat union S_Scanpy|
```

Interpretation:

| Value | Meaning |
|---:|---|
| 1 | same significant marker genes overall |
| 0 | no shared significant marker genes |

Code locations:

| Use | Code |
|---|---|
| Significant marker filtering and vectorization | `SVS:1518-1540` |
| Euler/UpSet/Jaccard helper | `PLOTSTATS:323-420` |
| Fig. 2 extraction phrase and baseline | `AGG:106-230` |

### 10. All marker gene-cluster pair Jaccard

Used when clusters are matched and the analysis can compare markers cluster by cluster. This is stricter than the significant marker gene metric because it keeps cluster identity attached to the gene.

Formula:

```text
M_Seurat = set of "gene-cluster" marker pairs from Seurat
M_Scanpy = set of "gene-cluster" marker pairs from Scanpy
J_markers = |M_Seurat intersect M_Scanpy| / |M_Seurat union M_Scanpy|
```

Interpretation:

| Value | Meaning |
|---:|---|
| 1 | same genes assigned as markers to the same matched clusters |
| lower values | different marker assignment and/or filtering |

Code locations:

| Use | Code |
|---|---|
| Build gene-cluster marker pairs | `SVS:1711-1726` |
| Euler/UpSet/Jaccard helper | `PLOTSTATS:323-420` |
| Fig. 2 same-input extraction phrase and baseline | `AGG:106-230` |

### 11. logFC concordance correlation coefficient

Used for comparing analogous marker-gene log fold-change values.

Formula from the paper:

```text
CCC = (2 * rho * sigma_x * sigma_y) /
      (sigma_x^2 + sigma_y^2 + (mu_x - mu_y)^2)
```

Equivalent implementation in the code:

```text
CCC = 2 * cov(x, y) / (var(x) + var(y) + (mean(y) - mean(x))^2)
```

where `x` and `y` are matched logFC values from the two workflows.

Interpretation:

| Value | Meaning |
|---:|---|
| 1 | perfect agreement along the identity line |
| near 0 | weak concordance |
| negative | discordant values |

Code locations:

| Use | Code |
|---|---|
| Join analogous Seurat and Scanpy DE rows | `SVS:1730-1753` |
| Filter extreme logFC outliers before CCC | `PLOTSTATS:1591` |
| Compute and print logFC CCC | `PLOTSTATS:1593-1595` |
| CCC function | `PLOTSTATS:1955-1995` and implementation at `PLOTSTATS:1983` |
| Fig. 2 same-input extraction phrase and baseline | `AGG:106-230` |

Additional logFC summaries:

```text
absolute_difference = abs(logFC_Seurat - logFC_Scanpy)
signed_difference = logFC_Seurat - logFC_Scanpy
```

These are summarized by mean, median, and variance in `PLOTSTATS:1560-1650`.

### 12. Adjusted p-value flip fraction

Used for comparing whether a marker gene is called significant in one package but not the other.

Formula:

```text
flip_fraction =
    count((p_adj_Seurat <= 0.05 and p_adj_Scanpy > 0.05) or
          (p_adj_Seurat > 0.05 and p_adj_Scanpy <= 0.05)) /
    number_of_matched_marker_rows
```

Interpretation:

| Value | Meaning |
|---:|---|
| 0 | no marker crosses the 0.05 threshold differently |
| higher values | more significance-threshold disagreement |

Code locations:

| Use | Code |
|---|---|
| Replace exact zero p-values with machine minimum for plotting | `SVS:1769-1775` |
| p-value group counts and flip fraction | `PLOTSTATS:1620-1635` |
| Spearman rank correlation for adjusted p-values | `PLOTSTATS:1637-1648` |
| Fig. 2 same-input extraction phrase and baseline | `AGG:106-230` |

The code also reports the fraction significant in both, significant only in group 1, significant only in group 2, and non-significant in both.

### 13. Adjusted p-value rank correlation

Used as an additional DE p-value agreement check. It is not the main Fig. 2 summary metric, but it is computed in the DE statistics helper.

Formula:

```text
rho_s = SpearmanCorr(rank(p_adj_Seurat), rank(p_adj_Scanpy))
```

Code location:

| Use | Code |
|---|---|
| Spearman correlation on adjusted p-values | `PLOTSTATS:1637-1648` |

### 14. Volcano significant gene sets

Used for downstream biological interpretation in Fig. 4 and supplement-style analyses. The volcano comparison asks whether the two pipelines highlight the same strongly changing genes in each cluster.

Threshold formula:

```text
significant_up = logFC >= 1 and adjusted_p < 0.05
significant_down = logFC <= -1 and adjusted_p < 0.05
significant_any = abs(logFC) >= 1 and adjusted_p < 0.05
```

Per-cluster Jaccard:

```text
J_cluster_c = |V_Seurat_c intersect V_Scanpy_c| /
              |V_Seurat_c union V_Scanpy_c|
```

Global Jaccard:

```text
J_global = |union_c V_Seurat_c intersect union_c V_Scanpy_c| /
           |union_c V_Seurat_c union union_c V_Scanpy_c|
```

Code locations:

| Use | Code |
|---|---|
| Volcano thresholds and plot colors | `PLOTSTATS:2028-2118` |
| Per-cluster and global volcano Jaccards | `PLOTSTATS:2367-2424` |
| Volcano Jaccard call | `SVS:1844` |

### 15. GO/pathway enrichment comparison

Used to ask whether marker-gene differences change downstream pathway interpretation.

Core logic:

```text
genes_for_cluster = genes with abs(logFC) >= 1 and adjusted_p < 0.05
enrichment = Enrichr(genes_for_cluster, database)
pathway_significant = adjusted_pathway_p < 0.05
```

The comparison then counts and plots pathways that are significant in Seurat only, Scanpy only, both, or neither, using the adjusted pathway p-values returned by Enrichr/gget.

Code locations:

| Use | Code |
|---|---|
| Run gget/Enrichr | `PLOTSTATS:2167-2220` |
| Build volcano-to-GO comparisons | `PLOTSTATS:2223-2310` |
| Main call | `SVS:1888` |

### 16. CellTypist downstream label agreement

Used for downstream biological annotation comparison when CellTypist is run on the two workflow outputs.

Metrics:

```text
ARI_celltype = adjustedRandIndex(labels_Seurat, labels_Scanpy)
fraction_agree = mean(labels_Seurat == labels_Scanpy)
```

Code locations:

| Use | Code |
|---|---|
| CellTypist on Scanpy output | `SVS:1543-1581` |
| CellTypist on Seurat output converted to AnnData | `SVS:1621-1652` |
| Align cells and compare labels | `SVS:1663-1692` |

## logFC formulas discussed in the paper

These formulas explain why Seurat and Scanpy can have different logFC values even when p-value testing is aligned. The Zenodo workflow does not reimplement the internal package formulas from scratch; it calls Seurat and Scanpy, extracts their DE tables, and compares the resulting values.

Let `Y_ig` be the log-normalized expression of gene `g` in cell `i`, and let `G1` and `G2` be the two cell groups with sizes `n1` and `n2`.

### Seurat v5 logFC

The paper gives Seurat v5 as:

```text
R_g =
log2((sum_{i in G1}(exp(Y_ig) - 1) + 1) / n1)
-
log2((sum_{i in G2}(exp(Y_ig) - 1) + 1) / n2)
```

Equivalently:

```text
R_g =
log2((1/n1) * sum_{i in G1}(exp(Y_ig) - 1) + 1/n1)
-
log2((1/n2) * sum_{i in G2}(exp(Y_ig) - 1) + 1/n2)
```

Interpretation: Seurat reverses the log transform cell by cell, averages on the unlogged normalized scale, and uses a group-size-dependent pseudocount after division (`1/n_group` in the equivalent form).

Where this appears:

| Use | Location |
|---|---|
| Paper discussion formulas | extracted main text around lines `670-715` |
| Seurat DE call/output extraction in workflow | `SVS:1518-1753` |

### Scanpy logFC

The paper gives Scanpy as:

```text
P_g =
log2(exp((1/n1) * sum_{i in G1}Y_ig) - 1 + epsilon)
-
log2(exp((1/n2) * sum_{i in G2}Y_ig) - 1 + epsilon)
```

with:

```text
epsilon = 1e-9
```

Interpretation: Scanpy averages the logged expression values first, then reverses the transform. That behaves like a geometric-mean-style calculation and can create very large logFC values when one group is close to zero.

Where this appears:

| Use | Location |
|---|---|
| Paper discussion formulas | extracted main text around lines `699-833` |
| Scanpy DE extraction helper | `HELPER:295-336` |
| Joined Seurat-vs-Scanpy DE comparison | `SVS:1730-1753` |

### Seurat v4 logFC

The paper contrasts Seurat v4 with Seurat v5:

```text
R_g =
log2((1/n1) * sum_{i in G1}(exp(Y_ig) - 1) + 1)
-
log2((1/n2) * sum_{i in G2}(exp(Y_ig) - 1) + 1)
```

Interpretation: Seurat v4 adds a larger pseudocount of 1 after averaging. This pushes fold-change ratios closer to 1 and logFC values closer to 0 compared with Seurat v5.

Where this appears:

| Use | Location |
|---|---|
| Paper discussion formulas | extracted main text around lines `753-779` |
| Version comparison scripts | `zenodo-codes_pachterlab/pachterlab-RMEJLBASBMP_2024-1de630a/analysis/rmd/Seurat_version_comparison.Rmd` |

## Bootstrap/null comparison formula

The paper estimates whether a Seurat-vs-Scanpy difference is bigger than expected from pipeline noise. It does this with bootstrapped datasets and combines package-specific deviations from the ideal value.

For a metric with ideal value `I`, bootstrap value from Seurat `X`, and bootstrap value from Scanpy `Y`:

```text
epsilon_X = X - I
epsilon_Y = Y - I
W = I + epsilon_X + epsilon_Y
```

So:

```text
W = I + (X - I) + (Y - I)
```

`W` is the null distribution for the combined amount of package noise expected if both workflows were only varying because of bootstrap sampling.

Code locations:

| Use | Code |
|---|---|
| Bootstrap metric list and baselines | `BOOT:92-206` |
| Choose observed Seurat-vs-Scanpy value | `BOOT:257-264` |
| Compute summed-error null value `W` | `BOOT:281-285` |
| Extract null values | `BOOT:290-293` |
| t-test against observed Seurat-vs-Scanpy value | `BOOT:295-300` |

Note: the STAR Methods text says p-values were calculated with one-tailed t-tests, but the local Zenodo code uses `alternative = "two.sided"` at `BOOT:295-300`.

## Which metrics are higher-is-better vs lower-is-better?

| Metric | Better direction | Ideal |
|---|---|---:|
| Jaccard indices | higher | 1 |
| ARI | higher | 1 |
| CCC | higher | 1 |
| PCA loading difference | lower | 0 |
| Median abs log2 SNN degree ratio | lower | 0 |
| Adjusted p-value flip fraction | lower | 0 |
| Absolute/signed logFC differences | lower absolute difference | 0 |

## One-sentence talk-ready summary

Rich et al. mostly quantify agreement as set overlap, label agreement, or numeric concordance: Jaccard asks whether the same objects were selected, ARI asks whether the same cells were grouped together, CCC asks whether logFC values agree on the identity line, and the downsampling/bootstrap analyses ask whether these Seurat-vs-Scanpy differences are larger than the variability expected from losing reads, losing cells, or resampling the dataset.
