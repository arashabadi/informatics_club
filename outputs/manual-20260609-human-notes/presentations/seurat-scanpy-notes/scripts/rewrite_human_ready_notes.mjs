import path from "node:path";

const ROOT = "/Users/ashi/github/informatics_club";
const ARTIFACT_TOOL =
  "/Users/ashi/.cache/codex-runtimes/codex-primary-runtime/dependencies/node/node_modules/@oai/artifact-tool/dist/artifact_tool.mjs";
const SOURCE = path.join(
  ROOT,
  "sctalk/06_seurat-vs-scanpy/06_sctalk0906_codebase_notes_smooth.pptx",
);
const OUT = path.join(
  ROOT,
  "sctalk/06_seurat-vs-scanpy/06_sctalk0906_ready_to_read_notes.pptx",
);

const { PresentationFile, FileBlob } = await import(ARTIFACT_TOOL);

const TALK_TRACKS = {
  7: "Figure 1 begins with the basic objects that Seurat and Scanpy carry forward. In panel A, the cell and gene filters are mostly a shared setup after the common quality criteria, but the variable-gene choice is where the pipelines split. The small zero in Aii should not be read as Scanpy detecting no genes; it means there are no Scanpy-only genes outside the Seurat gene set at that checkpoint, because the Scanpy workflow carries forward its 3,044 selected genes within the larger Seurat set. Panel Aiii is the main warning: the default HVG methods emphasize different features, and those features become the input for PCA, graphs, and clustering.",
  8: "After feature selection, the workflows enter PCA and graph construction. Panel B shows that the overall PCA structure is still recognizable, but the details change: PC1 explains a different amount of variance, and the eigenvector comparison shows that the coordinate systems are not identical. Panel C matters because clustering will be built from this neighborhood structure. The Jaccard index and degree-ratio plots show that many cells do not keep the same shared-nearest-neighbor relationships between packages. By this stage, Seurat and Scanpy are no longer just using different labels for the same graph; they are passing different graph objects into the next steps.",
  9: "By panels D and E, the upstream differences are visible in the outputs most people would actually discuss. Panel D compares the UMAPs and the alluvial mapping between clusters. The broad structure is related, but cluster assignments and local relationships are not identical. Panel E then moves to marker genes across all clusters, where the default workflows produce different sets of significant markers. That is the practical warning from Figure 1: defaults do not only change an intermediate object. They can change the map, the cluster labels, and the genes that become part of the biological story.",
  10: "Figure 2 puts the Seurat-versus-Scanpy difference onto a calibration scale. The dashed comparison is package choice on the same full dataset, and the curves ask how much downsampling inside one package is needed to reach that level of difference. In panel A, many read-downsampling comparisons reach the package-choice scale only after the data are reduced to less than five percent of reads. Panel B repeats the idea for cell downsampling, where the threshold is usually higher but still often below about a quarter of the cells. The result makes package choice feel concrete: it can resemble a large loss of sequencing or sampling information.",
  11: "Figure 3 shifts the reproducibility problem from package names to version numbers. The top row compares Seurat v5 with v4, where changes in log fold-change behavior affect marker selection and adjusted p values. The middle row does the same for Scanpy v1.9 versus v1.4, and the bottom row adds Cell Ranger v7 versus v6 before the data even enter Seurat or Scanpy. The visual pattern is that differential-expression results can change substantially across versions. For a methods section, the lesson is straightforward: reporting only the package name is not enough when the version can reshape the marker-gene list.",
  12: "Figure 4 connects the computational differences to biological interpretation. Panel A uses CellTypist assignments and shows that a noticeable fraction of cells change labels between the default Seurat and Scanpy pipelines, especially among monocytes and T-cell subtypes. Panels B and C focus on cluster 15, where the Seurat and Scanpy volcano plots call different sets of upregulated and downregulated genes. Panel D makes the limited overlap explicit, and panel E carries the comparison into pathway analysis. The figure turns the earlier workflow differences into a biological consequence: the same dataset can support different cell-type, gene, and pathway narratives depending on defaults.",
  13: "The code archive is organized in two layers. The matrix-generation layer starts from FASTQ files, references, and configuration files, then produces count matrices with kb or Cell Ranger, including downsampling and bootstrap-related steps. The analysis layer sits above that: figure-specific YAML files feed R Markdown notebooks, and helper scripts handle plotting, statistics, differential expression, and version-specific Scanpy behavior. For reproducing a figure, the YAML file is the best entry point because it tells the analysis notebooks which inputs and settings to use. The lower matrix-generation code becomes important when the goal is to regenerate the count matrices from raw sequencing data.",
  17: "Supplemental Figure 1 is the roadmap for the entire study. The top row lays out the major sources of variation the authors test: sequencing depth, sample size, matrix generation, analysis package, package version, and stochasticity. The lower row follows the analysis pipeline from filtering through normalization, HVG selection, PCA, graph construction, clustering, UMAP, and differential expression. The color coding separates sources of variability from the metrics used to compare outputs. This overview is useful because the paper is not making one isolated comparison; it is systematically asking where a single-cell workflow can drift.",
  18: "Supplemental Figure 2 asks whether Scanpy becomes more Seurat-like when the authors give it Seurat-like arguments and then run the whole pipeline sequentially. Panel A begins with filtering and HVG overlap, and the later panels carry those inherited differences into PCA, neighbor graphs, clustering, UMAP, and differential expression. Because each step receives the output of the previous step, the figure includes both local function differences and accumulated history. The comparison shows why argument matching can help without fully solving reproducibility: once early objects differ, downstream steps inherit those differences even when later settings are aligned.",
  19: "Supplemental Figure 3 removes much of the accumulation problem by controlling the input before each step. After each step, the Scanpy object is replaced with the Seurat-created object, so the next function receives the same starting point. That design makes the panels easier to interpret as step-level diagnostics. Differences that remain in filtering, PCA, graph construction, clustering, UMAP, or differential expression are more likely to come from the function behavior or default implementation at that step. Compared with Figure 1, this supplement separates propagation through the workflow from the direct effect of a specific analysis choice.",
  20: "Supplemental Figure 4 combines the controlled-input design with Seurat-like arguments in Scanpy. The controlled input keeps each step from being dominated by earlier differences, while the aligned arguments ask how close the functions can get when the obvious settings are matched. The panels follow the familiar sequence from filtering and HVGs to PCA, SNN graphs, clustering, UMAP, and differential expression. Where agreement improves, parameter choice is a major driver. Where differences remain, the result points toward implementation details, defaults without direct counterparts, or algorithmic behavior that cannot be harmonized by a simple argument change.",
  21: "Supplemental Figure 5 runs the alignment in the other direction by making Seurat as Scanpy-like as possible while controlling the input before each step. That symmetry matters because it keeps the comparison from treating either package as the natural ground truth. The panels again test filtering, PCA, graph construction, clustering, UMAP, and differential expression under controlled conditions. Some differences shrink, but others remain because the packages do not always expose matching settings or use the same internal conventions. The broader message is that reproducibility across tools is not just a matter of copying a parameter list.",
  22: "Supplemental Figure 6 focuses on UMAP-derived graph structure rather than treating UMAP as only a visualization. Each row represents a different comparison mode: default or aligned, and sequential or controlled. Within each row, the density plot summarizes neighborhood overlap, the alluvial plot compares Leiden clusters, and the UMAP panels show the resulting cluster structure. The repeated layout lets the audience see when UMAP neighborhoods remain similar and when they diverge. The important idea is that UMAP can participate in structural comparisons here, not just provide a final picture at the end of the workflow.",
  23: "Supplemental Figure 7 calibrates Seurat against itself after reducing the data to four percent of the original reads. Because the package is fixed, the differences across panels are driven by read depth rather than software choice. Panel A compares retained cells, genes, and HVGs; panels B and C carry the comparison into PCA and neighbor graphs; panel D shows clustering and UMAP; and panel E checks marker-gene behavior. The figure gives a concrete reference scale for the main paper: it shows what a severe read-depth perturbation looks like inside one workflow.",
  24: "Supplemental Figure 8 repeats the four-percent read-downsampling calibration inside Scanpy. The layout mirrors the Seurat version, so the comparison can stay focused on robustness rather than plot interpretation. The panels move from filtering and HVG overlap into PCA, graph neighborhoods, clustering, UMAP, and differential expression. Since both sides of the comparison use Scanpy, any changes reflect the effect of losing reads while holding the package and workflow fixed. Together with the previous supplement, this figure helps separate software-driven variation from the variation expected when sequencing depth is strongly reduced.",
  25: "Supplemental Figure 9 changes the perturbation from reads to cells. In Seurat, the authors keep sixteen percent of cells and compare that reduced dataset with the full dataset. Panel A shows how the retained cells, genes, and HVGs change; panels B and C test whether PCA and graph neighborhoods remain stable; panel D moves into clustering and UMAP; and panel E evaluates differential expression. This comparison asks a different robustness question from read downsampling: whether losing sampled cells alters the structure and markers of the analysis even when reads per remaining cell are not the main variable.",
  26: "Supplemental Figure 10 is the Scanpy counterpart to the sixteen-percent cell-downsampling experiment. The package and workflow stay fixed, while the dataset contains far fewer cells. The repeated panel structure makes the sensitivity easy to follow from filtering through HVG selection, PCA, graph construction, clustering, UMAP, and marker detection. Changes here represent the effect of sample-size loss within Scanpy rather than a Seurat-versus-Scanpy difference. Paired with the Seurat version, the figure helps show whether the two workflows respond similarly when the number of cells, rather than read depth, is sharply reduced.",
  27: "Supplemental Figure 11 summarizes the read-downsampling calibration across many workflow metrics. Each small panel tracks one comparison, including retained cells and genes, HVG overlap, PCA behavior, graph degree, clustering agreement, UMAP-neighborhood overlap, marker-gene overlap, log fold change, and adjusted p-value threshold crossing. Orange represents Seurat, blue represents Scanpy, and the dashed horizontal line marks the full-data Seurat-versus-Scanpy difference. The curves can be read as a scale bar: when a downsampling curve crosses the dashed line, read loss inside one package has reached the same magnitude as switching between packages.",
  28: "Supplemental Figure 12 gives the same threshold-style summary for cell downsampling. The x-axis is the fraction of cells retained, and each panel asks when the reduced dataset becomes as different from the full dataset as Seurat is from Scanpy on the original data. Compared with read downsampling, cell downsampling usually requires keeping a larger fraction before the same level of agreement is preserved. The main-text interpretation is that cell loss is less forgiving than read loss, but many metrics still reach the package-choice scale below roughly a quarter of cells.",
  29: "Supplemental Figure 13 brings the upstream count-matrix software into the reproducibility story by comparing Cell Ranger v7 with v6. The panels follow the usual workflow checkpoints, from filtering and HVG selection into PCA, graph construction, clustering, UMAP, and differential expression. The main paper highlights one specific reason this matters: Cell Ranger v7 includes intron counts by default, while v6 does not. That upstream change can shift the matrix before Seurat or Scanpy begins. The figure broadens the message of the paper: reproducibility depends on the full computational chain, not only the final analysis package.",
  30: "Supplemental Figure 14 gives a compact quality-control view of the input perturbations. Panel A compares UMI and feature distributions between Cell Ranger versions, and panel B compares the full dataset with the four-percent read-downsampled dataset. These distributions sit upstream of filtering, feature selection, and downstream graph construction. When the UMI or feature landscape shifts, later analysis steps can change even if the same high-level workflow is used. The figure is small, but it explains why the paper treats matrix-generation software and read depth as real sources of variation rather than background technical details.",
  31: "Supplemental Figure 15 extends the downstream biological comparison across clusters 0 through 9. Each row is one cluster, and the columns compare Seurat volcano results, Scanpy volcano results, overlap of biologically significant genes, and Enrichr pathway agreement. The volcano colors keep the reading simple: red genes are significantly upregulated, blue genes are significantly downregulated, and gray genes are not significant under the displayed thresholds. The repeated row structure shows that the cluster-15 example in Figure 4 is part of a broader pattern, with package defaults changing gene and pathway calls across multiple clusters.",
  32: "Supplemental Figure 16 continues the same cluster-by-cluster audit for clusters 10 through 18. The layout is identical to the previous supplement, which makes the differences easy to scan: Seurat and Scanpy can produce different volcano patterns, different significant-gene overlaps, and different pathway relationships. The value of this second half is coverage. It supports the claim that the biological consequences are not confined to one selected cluster; the effect of workflow defaults appears across much of the cluster set.",
  33: "Supplemental Figure 17 asks how much disagreement comes from random seed variation. The top half tests Seurat-related randomness in neighbor search, Louvain clustering, uwot UMAP, and UMAP-derived graphs. The bottom half tests analogous Scanpy components, including NNDescent, Leiden, and UMAP-learn. Each metric is compared with the package-level difference so the seed effect has a reference scale. The figure makes the interpretation more balanced: stochasticity can move the analysis, but several Seurat-versus-Scanpy differences are larger than the variation produced by changing seeds.",
  34: "Supplemental Figure 18 gives the graph intuition behind some of the package differences. Panel A shows how approximate KNN construction can add or miss edges, especially around local neighborhoods. Panel B then shows how Seurat-style shared-nearest-neighbor pruning changes which edges survive, with different behavior around hub and peripheral nodes. Scanpy's graph construction is closer to an undirected KNN graph. This schematic helps explain why similar PCA coordinates do not guarantee similar graph neighborhoods. Once the graph differs, clustering and UMAP can diverge even before differential expression begins.",
  35: "Supplemental Figure 19 lays out the graph-construction routes more explicitly. The upper path follows Seurat as it creates KNN and SNN graphs, prunes and stores them, and uses them for clustering or UMAP. The lower path follows the Scanpy route through its corresponding graph objects and downstream functions. Bold labels mark functions or methods, while italic labels mark external packages. The diagram makes the graph step feel architectural rather than cosmetic: the packages do not simply expose two names for the same internal object.",
  36: "Supplemental Figure 20 repeats the default Seurat-versus-Scanpy comparison in a mouse brain dataset. The panels deliberately mirror Figure 1, so the same checkpoints can be compared: filtering and HVGs, PCA, graph neighborhoods, clustering and UMAP, and differential expression. Because the dataset is different from the PBMC example, the result tests generality rather than one dataset-specific artifact. The figure supports the paper's broader claim that workflow divergence can appear in another biological system with a different cell composition.",
  37: "Supplemental Figure 21 repeats the same external-dataset check in non-small cell lung cancer. Again, the panel structure follows the Figure 1 template from early filtering through downstream marker-gene calls. The disease context and cell composition are different, but default Seurat and Scanpy settings still produce measurable differences across the pipeline. This figure strengthens the reproducibility argument by showing that the concern is not limited to PBMCs. The package-default effect remains visible when the biological dataset changes.",
  38: "Supplemental Figure 22 summarizes bootstrap variability across many workflow metrics. Each panel is one consistency measure across bootstrapped datasets generated with different random seeds, with orange for Seurat and blue for Scanpy. The dashed horizontal line again marks the full-data Seurat-versus-Scanpy difference. This gives bootstrap variation the same reference scale used for downsampling. The figure asks whether resampling noise inside a workflow is large enough to explain the package-level comparison. The answer helps place software choice in context alongside data sampling variability.",
  39: "Supplemental Figure 23 turns the bootstrap summary into a concrete example by comparing one Seurat bootstrap run with the original dataset. The panels follow the usual sequence through filtering, HVG selection, PCA, graph neighborhoods, clustering, UMAP, and differential expression. The value of this figure is that it makes the line-plot summary less abstract. A single resampled analysis can perturb the workflow, but those perturbations are being judged against the larger Seurat-versus-Scanpy benchmark used throughout the paper.",
  40: "Table S1 translates the paper's comparison into Seurat settings. The rows follow the workflow from filtering and normalization through HVG selection, scaling, PCA, SNN construction, clustering, UMAP, and differential expression. The color labels indicate whether Seurat matches Scanpy by default, can be made closer with arguments, remains partly incompatible, or has no clean match. The table is a practical methods checklist. For reproducible Seurat analysis, the relevant object is not just the package name; it is the exact function, version, default, and argument choice at each step.",
  41: "Table S2 gives the matching checklist from the Scanpy side. It lists Scanpy functions and defaults across the same pipeline stages and marks where they agree with Seurat, where an argument change helps, and where the behaviors do not have a clean equivalent. The table makes the paper operational. A reproducible method section needs more than a statement that Scanpy was used; it needs the version, function calls, defaults, and any deliberate argument changes that affect filtering, graphs, clustering, UMAP, and differential expression.",
};

const BAD_PATTERNS = [
  /\bstart by\b/i,
  /\borient (?:the audience|the room)\b/i,
  /\bwalk through\b/i,
  /\bpoint out\b/i,
  /\btell (?:them|the audience)\b/i,
  /\bi would\b/i,
  /\buse this slide\b/i,
  /\bthis slide (?:asks|is useful|is important)\b/i,
];

function noteText(slide) {
  return String(slide.speakerNotes?.text || slide.speakerNotes?.getText?.() || "");
}

function splitExistingNote(note) {
  const normalized = String(note || "").trim();
  const legendMatch = normalized.match(/(?:^|\n)Legend:\n([\s\S]*?)(?:\n\nSource:|\nSource:|$)/i);
  const sourceMatch = normalized.match(/(?:^|\n)(Source:[\s\S]*)$/i);
  if (legendMatch) {
    return {
      legend: legendMatch[1].trim(),
      source: sourceMatch ? sourceMatch[1].trim() : "Source: Rich et al., Cell Systems 2026.",
    };
  }

  const blocks = normalized.split(/\n\s*\n/).map((block) => block.trim()).filter(Boolean);
  const sourceIndex = blocks.findIndex((block) => /^Source:/i.test(block));
  const source = sourceIndex >= 0
    ? blocks.slice(sourceIndex).join("\n\n")
    : "Source: Rich et al., Cell Systems 2026.";
  let legendBlocks = sourceIndex >= 0 ? blocks.slice(0, sourceIndex) : blocks;
  if (legendBlocks.length > 1 && /^(Fig|Figure|Table)/i.test(legendBlocks[0])) {
    legendBlocks = legendBlocks.slice(1);
  }
  return {
    legend: legendBlocks.join("\n\n").trim() || normalized,
    source,
  };
}

function setJournalClubNotes(slide, talkTrack) {
  const { legend, source } = splitExistingNote(noteText(slide));
  slide.speakerNotes.setText(
    [
      "Ready-to-read note:",
      talkTrack,
      "",
      "Legend:",
      legend,
      "",
      source,
    ].join("\n"),
  );
}

function extractReadyNote(text) {
  const match = String(text || "").match(/Ready-to-read note:\s*([\s\S]*?)(?:\n\s*Legend:|$)/i);
  return match ? match[1].trim() : "";
}

function auditPresentation(presentation) {
  const missing = [];
  const promptLike = [];
  const tooShort = [];
  for (const slideNumber of Object.keys(TALK_TRACKS).map(Number).sort((a, b) => a - b)) {
    const slide = presentation.slides.items[slideNumber - 1];
    const ready = extractReadyNote(noteText(slide));
    if (!ready) {
      missing.push(slideNumber);
      continue;
    }
    const words = ready.split(/\s+/).filter(Boolean).length;
    if (words < 55) {
      tooShort.push({ slide: slideNumber, words });
    }
    for (const pattern of BAD_PATTERNS) {
      if (pattern.test(ready)) {
        promptLike.push({ slide: slideNumber, pattern: pattern.toString(), text: ready.slice(0, 180) });
      }
    }
  }
  return { missing, promptLike, tooShort };
}

const presentation = await PresentationFile.importPptx(await FileBlob.load(SOURCE));

for (const [slideNumber, talkTrack] of Object.entries(TALK_TRACKS)) {
  const slide = presentation.slides.items[Number(slideNumber) - 1];
  if (!slide) {
    throw new Error(`Missing slide ${slideNumber}`);
  }
  setJournalClubNotes(slide, talkTrack);
}

const preExportAudit = auditPresentation(presentation);
if (preExportAudit.missing.length || preExportAudit.promptLike.length || preExportAudit.tooShort.length) {
  throw new Error(`Pre-export note audit failed: ${JSON.stringify(preExportAudit, null, 2)}`);
}

const pptx = await PresentationFile.exportPptx(presentation);
await pptx.save(OUT);

const reloaded = await PresentationFile.importPptx(await FileBlob.load(OUT));
const postExportAudit = auditPresentation(reloaded);
if (postExportAudit.missing.length || postExportAudit.promptLike.length || postExportAudit.tooShort.length) {
  throw new Error(`Post-export note audit failed: ${JSON.stringify(postExportAudit, null, 2)}`);
}

console.log(
  JSON.stringify(
    {
      source: SOURCE,
      output: OUT,
      slideCount: reloaded.slides.count,
      updatedNotes: Object.keys(TALK_TRACKS).length,
      auditedSlides: Object.keys(TALK_TRACKS).map(Number).sort((a, b) => a - b),
      postExportAudit,
    },
    null,
    2,
  ),
);
