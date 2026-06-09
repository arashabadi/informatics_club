import fs from "node:fs/promises";
import path from "node:path";

const ROOT = "/Users/ashi/github/informatics_club";
const ARTIFACT_TOOL =
  "/Users/ashi/.cache/codex-runtimes/codex-primary-runtime/dependencies/node/node_modules/@oai/artifact-tool/dist/artifact_tool.mjs";
const SOURCE = path.join(ROOT, "sctalk/06_seurat-vs-scanpy/06_sctalk0906.pptx");
const BACKUP = path.join(
  ROOT,
  "outputs/manual-20260609-formula-explainers/presentations/seurat-scanpy-formulas/backups/06_sctalk0906.before-formula-explainers.pptx",
);
const ASSET_DIR = path.join(
  ROOT,
  "outputs/manual-20260609-formula-explainers/presentations/seurat-scanpy-formulas/assets",
);

const { PresentationFile, FileBlob } = await import(ARTIFACT_TOOL);

const SLIDE = { width: 1280, height: 720 };
const INSERTIONS = [
  {
    after: 12,
    image: "formula_06_volcano_pathways.png",
    title: "Formula explainer 6: volcano thresholds and pathway input",
    notes:
      "This explainer connects Figure 4 to the simple threshold logic behind volcano plots. The x-axis is log fold change, the y-axis is minus log10 of the adjusted p value, and genes are highlighted only when they pass both the effect-size and significance thresholds. The toy table makes the rule concrete: a very small p value is not enough if the logFC is small, and a large logFC is not enough if the p value is not significant. The final panel shows how these selected genes become the input to pathway enrichment and why different marker lists can lead to different pathway stories.",
    source:
      "Rich et al., Cell Systems 2026, Figure 4 legend and STAR Methods: Additional Biological Analysis.",
  },
  {
    after: 11,
    image: "formula_05_logfc_pvalues.png",
    title: "Formula explainer 5: logFC, pseudocounts, and adjusted p values",
    notes:
      "This explainer is for Figure 3 and the DE Methods. It reduces logFC to a tiny two-group example so the pseudocount issue is visible. When one group has zero expression, the pseudocount prevents division by zero, but the size of that pseudocount strongly changes the final logFC. The slide also introduces CCC as an agreement metric for comparing logFC vectors and shows why Bonferroni and Benjamini-Hochberg adjustment can move genes across the adjusted p < 0.05 threshold.",
    source:
      "Rich et al., Cell Systems 2026, Discussion: Marker selection and LogFC; STAR Methods: Differential Expression.",
  },
  {
    after: 10,
    image: "formula_04_downsampling_bootstrap.png",
    title: "Formula explainer 4: downsampling as a benchmark",
    notes:
      "This explainer belongs after Figure 2. The paper downsampled reads or cells, reran the workflow, and asked when the downsampled result became as different from the full result as Seurat was from Scanpy. The toy example shows the retained fraction f and compares a metric such as HVG Jaccard against the package-choice benchmark. The bootstrap panel summarizes the paper's null-distribution idea: deviations from the ideal metric in Seurat and Scanpy are combined to estimate the amount of noise expected from resampling alone.",
    source:
      "Rich et al., Cell Systems 2026, Figure 2 and STAR Methods: Data preparation; Quantification and statistical analysis.",
  },
  {
    after: 9,
    image: "formula_03_clustering_de_overlap.png",
    title: "Formula explainer 3: cluster agreement and marker overlap",
    notes:
      "This explainer follows the Figure 1 clustering and DE output panels. ARI is introduced as a cell-pair agreement score: if two cells are grouped together in one workflow, it checks whether the other workflow keeps them together, then corrects for agreement expected by chance. The marker example uses the same Jaccard logic as the earlier set slide, but now the objects are marker genes. The p-value flip table shows how a marker can be significant in only one workflow, which is why downstream interpretation can change.",
    source:
      "Rich et al., Cell Systems 2026, Figure 1D-E and STAR Methods: Clustering Analysis; Differential Expression.",
  },
  {
    after: 8,
    image: "formula_02_pca_snn.png",
    title: "Formula explainer 2: PCA angles and SNN neighborhoods",
    notes:
      "This explainer supports Figure 1B and 1C. The PCA part shows that the sine of the angle between PC vectors is near zero when the directions agree and near one when the directions are perpendicular. The graph part uses one cell's neighborhood to calculate Jaccard overlap and the SNN degree ratio. Together, these metrics explain why two workflows can begin from the same count matrix but pass different geometric and graph objects into clustering.",
    source:
      "Rich et al., Cell Systems 2026, Figure 1B-C and STAR Methods: Principal Component Analysis; K-Nearest Neighbors and Shared Nearest Neighbors.",
  },
  {
    after: 7,
    image: "formula_01_normalization_jaccard.png",
    title: "Formula explainer 1: normalization and set overlap",
    notes:
      "This explainer belongs after Figure 1A. It starts with a tiny count matrix and shows the log-normalization formula used to make cells with different library sizes comparable. The second panel then turns the HVG comparison into a Jaccard calculation: shared genes divided by the union of both selected sets. This is the same set-overlap logic used for cells, genes, HVGs, and later marker-gene sets in the paper.",
    source:
      "Rich et al., Cell Systems 2026, Figure 1A and STAR Methods: Filtering; Normalization; Highly Variable Gene Selection.",
  },
];

function contentType(filePath) {
  if (/\.(jpg|jpeg)$/i.test(filePath)) return "image/jpeg";
  return "image/png";
}

async function copyIfMissing(src, dst) {
  await fs.mkdir(path.dirname(dst), { recursive: true });
  try {
    await fs.stat(dst);
  } catch {
    await fs.copyFile(src, dst);
  }
}

function setNotes(slide, record) {
  slide.speakerNotes.setText(
    [
      "Ready-to-read note:",
      record.notes,
      "",
      "Source:",
      record.source,
    ].join("\n"),
  );
}

const presentation = await PresentationFile.importPptx(await FileBlob.load(SOURCE));
await copyIfMissing(SOURCE, BACKUP);

for (const record of INSERTIONS) {
  const afterSlide = presentation.slides.items[record.after - 1];
  if (!afterSlide) throw new Error(`Missing insertion anchor slide ${record.after}`);
  const { slide } = presentation.slides.insert({ after: afterSlide });
  slide.setViewportSize(SLIDE.width, SLIDE.height);
  const imagePath = path.join(ASSET_DIR, record.image);
  const data = await fs.readFile(imagePath);
  slide.images.add({
    data,
    contentType: contentType(imagePath),
    frame: { left: 0, top: 0, width: SLIDE.width, height: SLIDE.height },
    alt: record.title,
    fit: "fill",
  });
  setNotes(slide, record);
}

for (const slide of presentation.slides.items) {
  slide.setViewportSize(SLIDE.width, SLIDE.height);
}

const pptx = await PresentationFile.exportPptx(presentation);
await pptx.save(SOURCE);

console.log(
  JSON.stringify(
    {
      output: SOURCE,
      backup: BACKUP,
      inserted: INSERTIONS.length,
      slideCount: presentation.slides.count,
      insertionAnchors: INSERTIONS.map((item) => item.after).sort((a, b) => a - b),
    },
    null,
    2,
  ),
);
