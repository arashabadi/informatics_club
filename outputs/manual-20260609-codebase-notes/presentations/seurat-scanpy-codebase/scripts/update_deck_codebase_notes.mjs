import fs from "node:fs/promises";
import path from "node:path";

const ROOT = "/Users/ashi/github/informatics_club";
const ARTIFACT_TOOL =
  "/Users/ashi/.cache/codex-runtimes/codex-primary-runtime/dependencies/node/node_modules/@oai/artifact-tool/dist/artifact_tool.mjs";
const SOURCE = path.join(ROOT, "sctalk/06_seurat-vs-scanpy/06_sctalk0906.pptx");
const OUT = path.join(ROOT, "sctalk/06_seurat-vs-scanpy/06_sctalk0906_codebase_notes.pptx");
const CODEBASE_IMAGE = path.join(
  ROOT,
  "outputs/manual-20260609-codebase-notes/presentations/seurat-scanpy-codebase/assets/codebase_structure.png",
);

const { PresentationFile, FileBlob } = await import(ARTIFACT_TOOL);

const SLIDE = { width: 1280, height: 720 };
const CODEBASE_FRAME = { left: 72, top: 116, width: 1138, height: 540 };

const TALK_TRACKS = {
  7: "Start this first Figure 1 slide by emphasizing that both workflows begin with the same 10,000 PBMC dataset, but their default filtering and highly variable gene choices already separate the analysis paths. Panel A uses UpSet plots to compare retained cells, detected genes, and selected HVGs. The important point is not just that some sets differ, but that these early object definitions become the input for every later step. This is the first place to tell the audience that software defaults are acting like an experimental variable.",
  8: "This slide continues Figure 1 into the latent-space and graph-construction steps. Panel B compares PCA projections and eigenvector behavior, so it shows whether the two workflows represent the same cells in similar low-dimensional coordinates. Panel C then moves from PCA into the shared-nearest-neighbor graph, using neighborhood Jaccard similarity and degree ratios to show graph-level differences cell by cell. The message is that once HVG and scaling choices differ, PCA and graph topology can diverge even before clustering or biological interpretation begins.",
  9: "Use this last Figure 1 slice to show how upstream computational differences become visible biological outputs. Panel D compares Seurat and Scanpy UMAPs and links clusters through the alluvial plot, which makes disagreements in cluster assignment easy to see. Panel E then asks whether the same marker genes are recovered across clusters. This is the payoff of the figure: defaults do not merely change plots cosmetically; they propagate into clustering, visualization, and the marker genes a reader might use for interpretation.",
  10: "Figure 2 gives a useful effect-size comparison. Panel A asks how much read downsampling is needed before Seurat or Scanpy becomes as different from itself as Seurat is from Scanpy on the full dataset. Panel B asks the same question for cell downsampling. Orange bars are Seurat and blue bars are Scanpy, with metrics spanning filtering, PCA, graph structure, clustering, UMAP, and marker genes. The key framing is that package choice can introduce variation comparable to substantial loss of input data.",
  11: "Figure 3 isolates versioning as another source of variability, focusing on differential expression. The top row compares Seurat v5 and v4, the middle row compares Scanpy v1.9 and v1.4, and the bottom row compares Cell Ranger v7 and v6. Within each row, the panels move from marker-set overlap to log fold-change behavior and adjusted p-value behavior. The practical takeaway is that even when the software name is unchanged, version choices can shift the genes declared significant and therefore alter downstream biological conclusions.",
  12: "Figure 4 translates the computational differences into a biologically motivated downstream analysis. Panel A maps cell-type assignments between the default Seurat and Scanpy pipelines. Panels B and C show volcano plots for cluster 15 in each workflow, using the same log fold-change and p-value thresholds. Panel D compares the biologically significant genes selected from those volcano plots, and panel E asks whether pathway enrichment agrees. This slide is where the paper shows that workflow defaults can change the story told about cell identity, differential expression, and pathways.",
  17: "Use this supplemental overview as the roadmap for the whole paper. The top row separates the sources of variability the authors test: sequencing depth, sample size, count-matrix generation, and analysis package choice. The lower workflow follows the single-cell pipeline from filtering through normalization, HVG selection, scaling, PCA, graph construction, clustering, UMAP, and differential expression. The blue annotations mark the metrics used to compare steps. This slide is useful because it makes clear that the study is not only comparing final UMAPs, but measuring agreement throughout the pipeline.",
  18: "This slide asks what happens when Scanpy is made more Seurat-like while the full pipeline is run sequentially. Read it left to right: panel A compares filtering and HVGs, panel B checks PCA, panel C checks graph similarity, panel D checks clustering and UMAP, and panel E checks differential expression. Because each step receives the output of the previous step, this view captures accumulated workflow differences. The key use of the slide is to show which disagreements persist even after aligning function arguments.",
  19: "This controlled-input supplement is designed to isolate each step from the accumulated differences upstream. The authors run default Seurat and Scanpy, but after each stage they replace Scanpy's intermediate object with the Seurat object so the next step sees matched input. Panels A through E use the same filtering, PCA, graph, clustering, UMAP, and differential-expression readouts as Figure 1. This slide helps distinguish differences caused by a specific function from differences inherited from earlier preprocessing.",
  20: "This slide repeats the controlled-input analysis after changing Scanpy arguments to be more Seurat-like. The structure mirrors the previous figure: filtering and HVGs in panel A, PCA in panel B, SNN graph similarity in panel C, clustering and UMAP in panel D, and differential-expression agreement in panel E. The reason to show it is that argument alignment can reduce some differences, but the controlled design also exposes where package implementations or algorithmic choices still do not behave identically.",
  21: "This figure flips the alignment direction: Seurat is pushed toward Scanpy-like arguments under controlled inputs. Use it as the counterpart to the previous Seurat-like Scanpy comparison. Panels A through E again walk through filtering, PCA, graph construction, clustering or UMAP, and differential expression. The audience should hear that matching arguments is not symmetric or automatically sufficient; each package has defaults, available algorithms, and internal choices that can leave residual differences even when the user tries to harmonize settings.",
  22: "This supplemental UMAP analysis asks whether differences remain when the graph is built from UMAP space and then reclustered. Each row corresponds to one analysis mode: default sequential, aligned sequential, default controlled, and aligned controlled. Within each row, the density plots summarize neighborhood similarity, the alluvial plots compare Leiden clusters, and the UMAP panels show the corresponding cluster assignments. The slide is useful because it checks whether apparent UMAP differences are just visual layout effects or reflect graph and clustering changes as well.",
  23: "Here the authors test Seurat's robustness after reducing reads to 4 percent of the original dataset. The panels follow the same pipeline template: filtering and HVGs, PCA, graph similarity, clustering and UMAP, then differential expression. The comparison is full data versus read-downsampled data within Seurat, so it is a self-consistency benchmark rather than a Seurat-versus-Scanpy comparison. This slide helps calibrate how much technical data reduction is needed before the workflow begins to look meaningfully different from itself.",
  24: "This slide repeats the 4 percent read-downsampling robustness test for Scanpy. Again, the panels track the effect from filtering and HVG selection through PCA, graph construction, clustering, UMAP, and differential expression. Because it mirrors the Seurat version, it lets you compare whether the two ecosystems respond similarly to loss of sequencing depth. The main way to present it is as an internal stability check: how far can the input reads be reduced before the pipeline's own outputs change substantially?",
  25: "This figure switches from read downsampling to cell downsampling within Seurat, retaining only 16 percent of cells. Panel A shows how the retained cell, gene, and HVG sets compare with the full dataset. Panels B and C ask how PCA and graph neighborhoods change, while panel D compares clustering and UMAP. Panel E checks marker-gene overlap. The slide helps separate the effect of fewer molecules per cell from the effect of fewer cells in the population.",
  26: "This is the Scanpy counterpart to the Seurat cell-downsampling experiment. The authors retain 16 percent of cells and then compare the downsampled workflow to the full Scanpy workflow across the same sequence of readouts. The important distinction is that the package is held constant, so any differences reflect sample-size reduction rather than software choice. Use this slide to compare how cell-number loss propagates through Scanpy's filtering, PCA, graph, clustering, UMAP, and marker-gene steps.",
  27: "This summary figure turns the read-downsampling experiments into threshold curves. Each panel is one consistency metric, starting with cells, genes, and HVGs, then moving through PCA loadings, graph degree, clustering, UMAP-neighborhood similarity, marker-gene overlap, fold-change agreement, and adjusted p-value threshold crossing. Orange is Seurat, blue is Scanpy, lighter lines are replicate random seeds, and the dashed reference is the full-data Seurat-versus-Scanpy difference. The slide tells you where downsampling becomes comparable to package-choice variability.",
  28: "This figure is the cell-downsampling version of the threshold summary. The x-axis varies the fraction of cells retained, and each panel tracks one pipeline-consistency metric from filtering through differential expression. Orange and blue separate Seurat and Scanpy, and the dashed reference line keeps the main software-comparison effect in view. The point is to show which parts of the pipeline are more sensitive to losing cells and whether cell downsampling reaches the same magnitude of variation as switching packages.",
  29: "This figure focuses on count-matrix generation by comparing Cell Ranger v7 and v6 before the Seurat/Scanpy analysis. Panels A through E follow the same workflow structure: filtering and HVG overlap, PCA, SNN graph similarity, clustering and UMAP, and marker-gene overlap. Because only the Cell Ranger version changes, this slide shows that upstream matrix-generation software can propagate into downstream single-cell interpretation. It is a reminder that reproducibility depends not only on Seurat or Scanpy, but also on the preprocessing pipeline.",
  30: "This supplemental figure is a compact QC view for two upstream perturbations. Panel A compares Cell Ranger v7 and v6, and panel B compares the full dataset with the 4 percent read-downsampled dataset. The distributions show UMI counts and feature counts before filtering, so this is the place to connect upstream changes to the raw cellular summaries that feed filtering decisions. The slide supports the broader argument that technical preprocessing changes can shift the starting point of the analysis.",
  31: "This slide extends the downstream biology comparison across clusters 0 through 9. Each row is one cluster, and the columns compare Seurat volcano plots, Scanpy volcano plots, UpSet overlap of biologically significant genes, and Enrichr pathway agreement. Red and blue points mark significantly up- and downregulated genes, while gray points are not significant. The presenter takeaway is that the disagreement seen in one highlighted cluster is not a single-cluster anecdote; it recurs across many cluster-level DE analyses.",
  32: "This slide continues the same downstream analysis for clusters 10 through 18. Read each row as a repeated experiment: Seurat DE, Scanpy DE, overlap of significant genes, and Enrichr pathway comparison. Because the layout is identical to the previous slide, the audience can scan for whether disagreement persists across later clusters. The point is that differences in software defaults can spread across the full cluster set, changing not only individual gene calls but also the pathways selected for interpretation.",
  33: "This figure asks how much randomness alone contributes to graph, clustering, and UMAP differences. The top half shows Seurat random seeds across Annoy, Louvain, uwot, and UMAP-derived graph analyses; the bottom half shows the analogous Scanpy seed experiments using NNDescent, Leiden, and UMAP-learn. Each panel tests a place where stochastic algorithms can enter the pipeline. Use this slide to separate random-seed sensitivity from package-default sensitivity: randomness matters, but it is one component of a larger reproducibility problem.",
  34: "This schematic explains why KNN and SNN graph construction can differ between implementations. Panel A shows a small KNN graph example where approximate nearest-neighbor search can add or omit edges. Panel B shows how Seurat's SNN pruning changes the graph, including a hub node and a peripheral node. The contrast with Scanpy is important because Scanpy's SNN behavior is closer to the undirected KNN graph. This slide gives intuition for why graph-based clustering can diverge even with the same cells and PCs.",
  35: "This protocol schematic is the clearest map of how Seurat and Scanpy generate and reuse neighbor graphs. The top path is Seurat and the bottom path is Scanpy. Follow where KNN and SNN graphs are built, pruned, passed into clustering, or reused for UMAP. Bold labels mark package functions or methods, while italic labels mark external packages. The key point for journal club is that the two ecosystems do not just use different defaults; they organize graph objects and downstream calls differently.",
  36: "This supplemental figure tests whether the main Seurat-versus-Scanpy divergence generalizes beyond PBMCs by using a mouse brain dataset. The panels again track filtering, HVGs, PCA, SNN graph similarity, clustering, UMAP, and differential expression. Because the structure mirrors Figure 1, the audience can compare the same metrics in a new biological context. The takeaway is that the default-workflow differences are not limited to one human immune dataset; they also appear in a mouse brain setting.",
  37: "This slide repeats the cross-package default comparison in a non-small cell lung cancer dataset. Panels A through E use the same pipeline diagnostics as the PBMC and mouse brain analyses, moving from filtering and HVGs to PCA, graph structure, clustering, UMAP, and marker genes. Present it as an external dataset check: even when the cell types and biological system change, the computational choices still produce measurable differences across the workflow.",
  38: "This figure evaluates bootstrap variability by changing the random seed used to generate bootstrapped datasets. Each panel summarizes one consistency metric, from detected cells and genes through PCA, graph structure, clustering, UMAP neighborhoods, marker-gene overlap, fold-change concordance, and adjusted p-value threshold crossing. Orange is Seurat, blue is Scanpy, lighter lines are replicate seeds, and the dashed line is the full-data package-comparison reference. The slide asks whether resampling noise reaches the same scale as package-choice differences.",
  39: "This representative bootstrap slide shows what one bootstrapped Seurat dataset looks like compared with the original dataset. The panels use the now-familiar progression: filtering and HVGs, PCA, graph similarity, clustering and UMAP, then differential-expression agreement. Use it to make the abstract bootstrap summary more concrete. The main point is to show how a resampled dataset can alter many downstream outputs even when the analysis package and broad workflow are held constant.",
  40: "Table S1 is the practical checklist for making Seurat behave more like Scanpy. Read the columns as pipeline stages, from filtering through normalization, HVG selection, scaling, PCA, graph construction, clustering, UMAP, and differential expression. The colors tell you whether parameters already agree, can be matched with arguments, are partly incompatible, or are incompatible. This table is useful for the discussion because it turns the paper's comparisons into concrete knobs a user would have to set for reproducibility.",
  41: "Table S2 is the reciprocal checklist for making Scanpy behave more like Seurat. The columns follow the same pipeline stages as Table S1, but now the rows list Scanpy functions and defaults. The color coding shows which settings already align, which require explicit changes, and which do not have clean equivalents. The key message is that harmonization is not just remembering a few parameters; some differences are structural or algorithmic and need to be documented in methods sections."
};

function textOf(shape) {
  return String(shape?.text || "").trim();
}

function updateTitle(slide, oldText, newText) {
  const shape = slide.shapes.items.find((candidate) => textOf(candidate) === oldText);
  if (shape) {
    shape.text = newText;
  }
}

function splitExistingNote(note) {
  const normalized = String(note || "").trim();
  const blocks = normalized.split(/\n\s*\n/).map((block) => block.trim()).filter(Boolean);
  const sourceIndex = blocks.findIndex((block) => /^Source:/i.test(block));
  const source = sourceIndex >= 0 ? blocks.slice(sourceIndex).join("\n\n") : "Source: Rich et al., Cell Systems 2026.";
  let legendBlocks = sourceIndex >= 0 ? blocks.slice(0, sourceIndex) : blocks;
  if (legendBlocks.length > 1 && /^(Fig|Figure|Table)/i.test(legendBlocks[0])) {
    legendBlocks = legendBlocks.slice(1);
  }
  const legend = legendBlocks.join("\n\n").trim() || normalized;
  return { legend, source };
}

function setJournalClubNotes(slide, talkTrack) {
  const note = String(slide.speakerNotes?.text || slide.speakerNotes?.getText?.() || "");
  const { legend, source } = splitExistingNote(note);
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

function fitContain(width, height, box) {
  const scale = Math.min(box.width / width, box.height / height);
  const fittedWidth = width * scale;
  const fittedHeight = height * scale;
  return {
    left: box.left + (box.width - fittedWidth) / 2,
    top: box.top + (box.height - fittedHeight) / 2,
    width: fittedWidth,
    height: fittedHeight,
  };
}

function pngDimensions(buffer) {
  return { width: buffer.readUInt32BE(16), height: buffer.readUInt32BE(20) };
}

const presentation = await PresentationFile.importPptx(await FileBlob.load(SOURCE));

const codeSlide = presentation.slides.items[12];
updateTitle(codeSlide, "Codes", "Codebase");
for (const image of [...codeSlide.images.items]) {
  image.delete();
}
const imageData = await fs.readFile(CODEBASE_IMAGE);
codeSlide.images.add({
  data: imageData,
  contentType: "image/png",
  frame: fitContain(pngDimensions(imageData).width, pngDimensions(imageData).height, CODEBASE_FRAME),
  alt: "Diagram summarizing the Rich et al. reproducibility codebase structure.",
  fit: "contain",
});
codeSlide.speakerNotes.setText(
  [
    "Ready-to-read note:",
    "This slide is a map for reading the reproducibility archive. The repository has two main branches. On the left, matrix_generation is the lower-level workflow that starts from FASTQ files and references, reads config.yaml, optionally downloads inputs, downsamples reads, and produces count matrices with kb or Cell Ranger. On the right, analysis is the figure-producing layer: figure-specific YAML files feed R Markdown notebooks, and helper scripts handle plotting, statistics, differential expression, and version-specific Scanpy behavior. The practical way to navigate the repo is to choose the figure YAML, run the matching Rmd, and only drop into matrix_generation when regenerating matrices from raw FASTQs.",
    "",
    "Legend:",
    "Generated overview of the local Zenodo code archive for Rich et al. The archive separates FASTQ-to-matrix generation from downstream Seurat/Scanpy analysis, with figure-specific YAML files controlling the R Markdown analysis notebooks.",
    "",
    "Source:",
    "Local code archive: sctalk/06_seurat-vs-scanpy/paper/zenodo-codes_pachterlab/pachterlab-RMEJLBASBMP_2024-1de630a.",
  ].join("\n"),
);

for (const [slideNumber, talkTrack] of Object.entries(TALK_TRACKS)) {
  const slide = presentation.slides.items[Number(slideNumber) - 1];
  if (slide) {
    setJournalClubNotes(slide, talkTrack);
  }
}

for (const slide of presentation.slides.items) {
  slide.setViewportSize(SLIDE.width, SLIDE.height);
}

const pptx = await PresentationFile.exportPptx(presentation);
await pptx.save(OUT);

console.log(
  JSON.stringify(
    {
      source: SOURCE,
      output: OUT,
      updatedFigureOrTableNotes: Object.keys(TALK_TRACKS).length,
      codebaseSlide: 13,
      slideCount: presentation.slides.count,
    },
    null,
    2,
  ),
);
