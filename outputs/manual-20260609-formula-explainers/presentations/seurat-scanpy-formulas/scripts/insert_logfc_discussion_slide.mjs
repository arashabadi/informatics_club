import fs from "node:fs/promises";
import path from "node:path";

const ROOT = "/Users/ashi/github/informatics_club";
const ARTIFACT_TOOL =
  "/Users/ashi/.cache/codex-runtimes/codex-primary-runtime/dependencies/node/node_modules/@oai/artifact-tool/dist/artifact_tool.mjs";
const SOURCE = path.join(ROOT, "sctalk/06_seurat-vs-scanpy/06_sctalk0906.pptx");
const BACKUP = path.join(
  ROOT,
  "outputs/manual-20260609-formula-explainers/presentations/seurat-scanpy-formulas/backups/06_sctalk0906.before-logfc-discussion-slide.pptx",
);
const IMAGE = path.join(
  ROOT,
  "outputs/manual-20260609-formula-explainers/presentations/seurat-scanpy-formulas/assets/logfc_discussion_summary.png",
);

const { PresentationFile, FileBlob } = await import(ARTIFACT_TOOL);

const presentation = await PresentationFile.importPptx(await FileBlob.load(SOURCE));
await fs.mkdir(path.dirname(BACKUP), { recursive: true });
try {
  await fs.stat(BACKUP);
} catch {
  await fs.copyFile(SOURCE, BACKUP);
}

const afterSlide = presentation.slides.items[17]; // after slide 18: existing logFC / p-value explainer
if (!afterSlide) throw new Error("Expected slide 18 to insert after.");

const { slide } = presentation.slides.insert({ after: afterSlide });
slide.setViewportSize(1280, 720);
slide.images.add({
  data: await fs.readFile(IMAGE),
  contentType: "image/png",
  frame: { left: 0, top: 0, width: 1280, height: 720 },
  alt: "LogFC discussion summary: pseudocounts, arithmetic versus geometric averaging, and consequences for Seurat and Scanpy",
  fit: "fill",
});

slide.speakerNotes.setText(
  [
    "Ready-to-read note:",
    "This discussion slide is the main warning behind the paper's logFC results. The authors are saying that Seurat and Scanpy do not simply compute the same fold-change statistic in two languages. Seurat reverses the log transform before averaging, so it behaves more like an arithmetic-mean calculation, while Scanpy averages logged values first and then reverses the transform, which behaves more like a geometric-mean calculation. The pseudocounts also differ: Seurat v5 uses a cluster-size-dependent pseudocount, Seurat v4 used a larger effective +1 pseudocount that compresses logFC toward zero, and Scanpy uses a very small epsilon. That is why a zero-expression edge case can become about -6.6 in a Seurat-like calculation but about -29.9 in Scanpy. The practical takeaway is that large logFC disagreements, especially extreme Scanpy values, can come from implementation choices rather than new biology.",
    "",
    "Legend:",
    "Discussion summary slide created from Rich et al.'s LogFC discussion. The slide condenses the intended logFC idea, the Seurat and Scanpy computational routes, pseudocount differences across Seurat v5, Seurat v4, and Scanpy, and the paper's minimal zero-expression example.",
    "",
    "Source:",
    "Rich et al., Cell Systems 2026, Discussion: Marker selection and LogFC; STAR Methods: Differential Expression.",
  ].join("\n"),
);

for (const item of presentation.slides.items) {
  item.setViewportSize(1280, 720);
}

const pptx = await PresentationFile.exportPptx(presentation);
await pptx.save(SOURCE);

console.log(
  JSON.stringify(
    {
      output: SOURCE,
      backup: BACKUP,
      insertedAfterSlide: 18,
      slideCount: presentation.slides.count,
    },
    null,
    2,
  ),
);
