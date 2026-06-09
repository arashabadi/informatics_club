import fs from "node:fs/promises";
import path from "node:path";

const ROOT = "/Users/ashi/github/informatics_club";
const ARTIFACT_TOOL =
  "/Users/ashi/.cache/codex-runtimes/codex-primary-runtime/dependencies/node/node_modules/@oai/artifact-tool/dist/artifact_tool.mjs";
const SOURCE = path.join(ROOT, "sctalk/06_seurat-vs-scanpy/06_sctalk0906.pptx");
const BACKUP = path.join(
  ROOT,
  "outputs/manual-20260609-formula-explainers/presentations/seurat-scanpy-formulas/backups/06_sctalk0906.before-extra-formula-explainers.pptx",
);
const ASSET_DIR = path.join(
  ROOT,
  "outputs/manual-20260609-formula-explainers/presentations/seurat-scanpy-formulas/assets",
);

const { PresentationFile, FileBlob } = await import(ARTIFACT_TOOL);

const SLIDE = { width: 1280, height: 720 };
const INSERTIONS = [
  {
    after: 10,
    image: "formula_08_graph_modularity.png",
    title: "Formula explainer 2b: graph limits and modularity",
    notes:
      "This method-math slide extends the Figure 1B-C graph explanation. The first panel shows the maximum possible Jaccard overlap when one SNN neighborhood is much larger than the other: even perfect nesting cannot exceed the smaller degree divided by the larger degree. The second panel gives the intuition behind exact KNN cost, where many distances must be checked as cell number grows, motivating approximate algorithms such as Annoy and PyNNDescent. The last panel simplifies modularity quality as observed within-community edges minus an expected random-edge penalty, which helps explain why implementation details in clustering optimizers can shift the final groups.",
    source:
      "Rich et al., Cell Systems 2026, STAR Methods: K-nearest neighbors and shared nearest neighbors; Clustering Analysis.",
  },
  {
    after: 8,
    image: "formula_07_hvg_scoring.png",
    title: "Formula explainer 1b: HVG scoring math",
    notes:
      "This method-math slide extends the Figure 1A HVG comparison. The tiny count matrix shows why raw variance alone is not enough: genes with larger means tend to have larger variance, so HVG methods score variability relative to the expected mean-variance relationship. The mean.var.plot style bins genes by mean expression and asks whether a gene's dispersion is unusually high for its bin. The vst/seurat_v3 style fits an expected variance curve, rescales each expression value, and ranks genes by normalized variance. This is the mathematical reason two workflows can start from the same normalized data but carry forward different HVG sets.",
    source:
      "Rich et al., Cell Systems 2026, STAR Methods: Highly Variable Gene Selection.",
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
