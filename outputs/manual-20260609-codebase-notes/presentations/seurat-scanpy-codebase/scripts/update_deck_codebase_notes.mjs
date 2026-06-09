import fs from "node:fs/promises";
import path from "node:path";

const ROOT = "/Users/ashi/github/informatics_club";
const ARTIFACT_TOOL =
  "/Users/ashi/.cache/codex-runtimes/codex-primary-runtime/dependencies/node/node_modules/@oai/artifact-tool/dist/artifact_tool.mjs";
const SOURCE = path.join(ROOT, "sctalk/06_seurat-vs-scanpy/06_sctalk0906.pptx");
const OUT = path.join(ROOT, "sctalk/06_seurat-vs-scanpy/06_sctalk0906_codebase_notes_smooth.pptx");
const CODEBASE_IMAGE = path.join(
  ROOT,
  "outputs/manual-20260609-codebase-notes/presentations/seurat-scanpy-codebase/assets/codebase_structure.png",
);

const { PresentationFile, FileBlob } = await import(ARTIFACT_TOOL);

const SLIDE = { width: 1280, height: 720 };
const CODEBASE_FRAME = { left: 72, top: 116, width: 1138, height: 540 };

const TALK_TRACKS = {
  7: "I would open Figure 1 by slowing down on panel A, because this is where the pipeline begins to fork before PCA or clustering. The paper says the basic cell and gene filtering criteria themselves did not differ after the shared quality filters. The confusing part is Aii: Scanpy is not detecting zero genes. In this workflow, Scanpy carries forward only the highly variable gene subset for later steps, so the 3,044 Scanpy genes are all inside the Seurat gene set, and the zero is the Scanpy-only intersection. The real package split is clearest in Aiii, where the default HVG algorithms choose very different variable-gene sets.",
  8: "Here the same dataset has moved from feature selection into PCA and graph construction. Panel B shows that the broad PCA shape is recognizable, but the details are not the same: the variance explained by PC1 changes, and the eigenvectors are measurably different. That matters because panel C builds the shared-nearest-neighbor graph from this space. The low neighborhood overlap and degree-ratio shift tell us that the packages are no longer just labeling the same structure differently; they are handing clustering a different graph.",
  9: "This final part of Figure 1 is where the earlier differences become visible to a biologist reading the analysis. In panel D, the UMAPs retain some broad structure, but the alluvial plot shows cluster assignments are not identical, and the paper points out examples where local cluster relationships change. Panel E then moves to marker genes, where Seurat reports substantially more significant markers than Scanpy under default settings. So the main lesson is practical: defaults can change the clusters, the map, and the genes you would discuss.",
  10: "Figure 2 gives the paper a useful sense of scale. The authors ask: if we ignore the Seurat-versus-Scanpy difference, how much data loss would produce a similar amount of variation within one workflow? Panel A answers that most read-downsampling comparisons reach the same scale with less than five percent of reads retained. Panel B shows cell downsampling is more demanding, but often still under about a quarter of cells. The uncomfortable point is that package choice can be as large as throwing away a lot of data.",
  11: "Figure 3 shifts from package choice to version choice, and the warning is very concrete: the differential-expression layer is fragile. The top row compares Seurat v5 with v4, where changes in log fold-change behavior affect marker selection. The middle row compares Scanpy v1.9 with v1.4, where marker filtering changed across releases. The bottom row compares Cell Ranger v7 with v6, where upstream matrix generation also changes DE outputs. I would present this as a version-control slide, not just a DE slide.",
  12: "Figure 4 asks whether these computational differences matter biologically, and the answer is yes. Panel A shows that about eight percent of cells change CellTypist assignment between default Seurat and Scanpy pipelines, especially among monocytes and T-cell subclasses. Panels B and C focus on cluster 15, where Scanpy calls more downregulated genes but fewer upregulated genes than Seurat. Panel D shows the overlap is limited, and panel E shows pathway calls diverge too, even though many pathways remain shared. This is the biological consequence slide.",
  17: "I would use this first supplemental figure as the map of the whole study. The top row names the sources of variation the authors are testing: sequencing depth, sample size, matrix generation, and the analysis package itself. The lower row is the pipeline we keep returning to, from filtering through normalization, HVGs, PCA, graphs, clustering, UMAP, and differential expression. The blue text marks the comparison metrics. This slide is useful because it makes the paper feel systematic rather than a collection of separate plots.",
  18: "This slide asks a simple question: if Scanpy is given Seurat-like arguments, do the two workflows become the same? Because this is a sequential run, each step inherits whatever differences were produced before it. So panel A starts with filtering and HVGs, then the differences flow into PCA, SNN graphs, clustering, UMAP, and DE. I would use this slide to show that matching arguments helps in some places, but it does not magically erase the full pipeline history.",
  19: "This controlled-input figure is a cleaner diagnostic. Instead of letting differences accumulate, the authors replace the Scanpy object with the Seurat object before each step, so each function is tested on the same input. That makes the slide a step-by-step stress test of default behavior. If a difference appears here, it is more likely coming from that specific function or default setting rather than from earlier preprocessing. It is a helpful companion to Figure 1 because it separates propagation from local implementation differences.",
  20: "This slide repeats the controlled-input idea, but now Scanpy is also pushed toward Seurat-like arguments. I would describe it as asking whether the remaining differences are still there after we have been as fair as possible to the two packages. The panels follow the familiar path from filtering and HVGs through PCA, SNN graphs, clustering, UMAP, and DE. When agreement improves, it tells us argument choice matters; when it does not, it points to deeper implementation or algorithmic differences.",
  21: "This is the reverse alignment experiment: Seurat is made as Scanpy-like as the authors can make it, while the input before each step is controlled. That symmetry is important for the journal club discussion, because it prevents us from treating one package as the default truth. The slide shows that harmonization depends on which direction you try to harmonize in, and some settings simply do not have perfect counterparts. In other words, reproducibility is not just copying a parameter list from one package into the other.",
  22: "This UMAP supplement checks whether UMAP differences are only visual or whether they also affect graph-based structure. Each row represents a different analysis mode: default or aligned, sequential or controlled. Within each row, the density plots summarize UMAP-derived neighborhood overlap, the alluvial plots compare Leiden clusters, and the UMAP panels show the clusters on the embedding. I would use this to say that UMAP is not just decoration here; its neighborhood structure can feed back into clustering comparisons.",
  23: "Here the authors ask how stable Seurat is when reads are reduced to four percent of the original dataset. The comparison is Seurat against itself, full data versus read-downsampled data, so this is a calibration slide. The panels walk through the same pipeline checkpoints as Figure 1: filtering, HVGs, PCA, SNN graph, clustering, UMAP, and DE. The reason it matters is that we need a reference scale: how much technical data loss produces differences comparable to switching software?",
  24: "This is the same four-percent read-downsampling experiment, now inside Scanpy. Because the layout matches the Seurat slide, the audience can compare package robustness without learning a new plot grammar. I would present it as a self-consistency test: the package is fixed, the workflow is fixed, and only sequencing depth is reduced. Any changes across filtering, PCA, graphs, clustering, UMAP, or DE are therefore a measure of how much information was lost by downsampling reads.",
  25: "This slide changes the perturbation from reads to cells. In Seurat, the authors keep only sixteen percent of cells and then compare the result with the full dataset. Panel A tells us how cell, gene, and HVG sets change; panels B and C move into PCA and graph neighborhoods; panel D shows clusters and UMAP; and panel E checks marker genes. I would frame this as asking whether fewer sampled cells disturb the analysis in the same way as fewer reads per cell.",
  26: "This is the Scanpy version of the sixteen-percent cell-downsampling test. Again, the package is compared with itself, so we are not yet asking which package is better. We are asking how much the Scanpy workflow moves when the dataset contains far fewer cells. The repeated panel structure makes it easy to see where the pipeline is stable and where it is sensitive. This slide pairs with the previous one to compare how Seurat and Scanpy respond to sample-size loss.",
  27: "This summary figure is where the read-downsampling story becomes quantitative. Each panel is one pipeline metric, from retained cells and genes through HVGs, PCA, graph degree, clustering, UMAP-neighborhood overlap, marker genes, log fold change, and adjusted p-value threshold crossing. Orange is Seurat, blue is Scanpy, and the dashed line is the full-data Seurat-versus-Scanpy difference. I would tell the audience to read the dashed line as the benchmark: when a curve crosses it, downsampling has reached package-choice scale.",
  28: "This is the same threshold summary for cell downsampling. The x-axis is the fraction of cells retained, and each panel asks when that reduced dataset becomes as different from the full dataset as Seurat is from Scanpy. The paper's main text says cell downsampling is still fairly robust, but less so than read downsampling, with many steps comparable below about twenty-five percent of cells. This slide helps separate sequencing-depth robustness from sample-size robustness.",
  29: "This figure brings Cell Ranger into the story. The authors compare v7 and v6 count matrices and then carry those differences through the analysis workflow. Panels A through E show that the matrix-generation version can affect filtering, HVGs, PCA, graph structure, clustering, UMAP, and marker genes. The main paper highlights one key change: Cell Ranger v7 includes intron counts by default, while v6 does not. So this slide reminds us that reproducibility starts before Seurat or Scanpy ever sees the data.",
  30: "This compact QC figure gives a more intuitive view of the upstream perturbations. Panel A compares UMI and feature distributions between Cell Ranger versions, while panel B compares the full dataset with four-percent read downsampling. These are pre-filtering summaries, so they sit upstream of many later decisions. I would use this slide to say: if the input distributions move, filtering and downstream feature selection can move too. It is a small slide, but it explains why upstream software and data depth matter.",
  31: "This slide extends the biological downstream analysis across clusters 0 through 9. Each row is a cluster, and the columns show Seurat volcano results, Scanpy volcano results, overlap of biologically significant genes, and Enrichr pathway agreement. The color coding is straightforward: red and blue are significant up- and downregulated genes, gray is not significant. The point is not to memorize each row, but to see that the Figure 4 cluster-15 story is not an isolated example.",
  32: "This slide continues the same analysis for clusters 10 through 18. I would present it as the second half of the cluster-by-cluster audit. The repeated layout lets the audience scan for a pattern: Seurat and Scanpy often produce different volcano calls, different significant-gene overlaps, and different pathway relationships. Together with the previous slide, it supports the paper's claim that package defaults can change downstream biological interpretation across the cluster set, not just in one hand-picked cluster.",
  33: "This figure asks how much of the disagreement could simply be random seed noise. The top half tests Seurat randomness in neighbor search, Louvain clustering, uwot UMAP, and UMAP-derived graphs. The bottom half tests the analogous Scanpy pieces with NNDescent, Leiden, and UMAP-learn. The main text says random seeds do introduce variation, but several package-level differences are larger than seed effects. So I would use this slide to be fair: randomness matters, but it does not explain everything.",
  34: "This schematic is a small but important graph-intuition slide. Panel A shows how approximate KNN construction can add or miss edges. Panel B then shows Seurat-style SNN pruning, where shared-neighbor structure changes which edges survive, especially around hub and peripheral nodes. The contrast is that Scanpy's method is closer to an undirected KNN graph. This helps explain why two packages can start from similar PCA coordinates yet produce graph neighborhoods that are not the same.",
  35: "This protocol diagram is the clearest way to explain why the graph step is hard to harmonize. The top path shows how Seurat builds and uses KNN and SNN graphs; the bottom path shows the Scanpy route. Follow where graph objects are generated, pruned, stored, and reused for clustering or UMAP. Bold labels are functions or methods, and italic labels are external packages. I would use this slide to say that the disagreement is architectural, not just a single parameter mismatch.",
  36: "This supplement asks whether the main result holds outside the PBMC dataset by repeating the default Seurat-versus-Scanpy comparison in mouse brain. The panels are intentionally familiar: filtering and HVGs, PCA, SNN graph, clustering and UMAP, and DE. Because the figure mirrors Figure 1, the audience can focus on generalization rather than learning new methods. The message is that the paper's concern is not one dataset-specific artifact; similar workflow divergence appears in another biological system.",
  37: "This is the same external-dataset check in non-small cell lung cancer. Again, the authors use the Figure 1 panel structure so we can track the whole pipeline from filtering to DE. I would present this as a robustness check for the paper's claim. Even with a different disease context and cell composition, default Seurat and Scanpy choices still produce measurable differences across intermediate and downstream outputs. That makes the reproducibility issue broader than the PBMC example.",
  38: "This figure looks at bootstrap variability rather than downsampling fraction. Each panel is one consistency metric across bootstrapped datasets generated with different random seeds, and the dashed line again marks the full-data Seurat-versus-Scanpy difference. Orange is Seurat and blue is Scanpy. I would explain it as a second calibration experiment: after asking how much data loss matters, the authors ask how much resampling noise matters. The answer helps place package-choice effects in a more interpretable range.",
  39: "This slide makes the bootstrap summary concrete by showing one representative Seurat bootstrap against the original dataset. The panels move through the same pipeline checkpoints: filtering, HVGs, PCA, graph neighborhoods, clustering, UMAP, and DE. It is useful because line plots can feel abstract; here the audience can see what one resampled analysis actually looks like. The main takeaway is that resampling can perturb the workflow, but it is being measured against the larger package-comparison benchmark.",
  40: "Table S1 is the practical translation of the paper into Seurat settings. Read the columns as the analysis pipeline: filtering, normalization, HVG selection, scaling, PCA, SNN, clustering, UMAP, and DE. The colors tell us whether Seurat agrees with Scanpy by default, can be matched with arguments, is partly incompatible, or is incompatible. I would use this table during discussion: if someone wants reproducible Seurat code, these are the defaults and arguments they need to document.",
  41: "Table S2 is the mirror image for Scanpy. It lists Scanpy functions and defaults across the same pipeline stages, then marks where they already match Seurat, where an argument change can help, and where there is no clean equivalent. This table makes the paper feel very practical. The message is that harmonizing workflows is not just saying 'we used Seurat' or 'we used Scanpy'; the exact defaults, versions, and unmatched behaviors need to be part of the methods."
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
  const existingLegendMatch = normalized.match(/(?:^|\n)Legend:\n([\s\S]*?)(?:\n\nSource:|\nSource:|$)/i);
  const existingSourceMatch = normalized.match(/(?:^|\n)(Source:[\s\S]*)$/i);
  if (existingLegendMatch) {
    return {
      legend: existingLegendMatch[1].trim(),
      source: existingSourceMatch ? existingSourceMatch[1].trim() : "Source: Rich et al., Cell Systems 2026.",
    };
  }
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
for (const shape of [...codeSlide.shapes.items]) {
  const text = textOf(shape);
  if (/github\.com\/pachterlab|zenodo\.org\/records/i.test(text)) {
    shape.delete();
  }
}
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
