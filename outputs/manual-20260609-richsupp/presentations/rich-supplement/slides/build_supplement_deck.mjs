import fs from "node:fs/promises";
import path from "node:path";

const ROOT = "/Users/ashi/github/informatics_club";
const ARTIFACT_TOOL =
  "/Users/ashi/.cache/codex-runtimes/codex-primary-runtime/dependencies/node/node_modules/@oai/artifact-tool/dist/artifact_tool.mjs";

const { PresentationFile, FileBlob } = await import(ARTIFACT_TOOL);

const SOURCE = path.join(ROOT, "sctalk/06_seurat-vs-scanpy/06_sctalk0906.pptx");
const OUT = path.join(ROOT, "sctalk/06_seurat-vs-scanpy/06_sctalk0906_with_supplemental.pptx");
const MANIFEST = path.join(
  ROOT,
  "outputs/manual-20260609-richsupp/presentations/rich-supplement/assets/supp_manifest.json",
);

const SLIDE = { width: 1280, height: 720 };
const IMAGE_BOX = { left: 300, top: 12, width: 860, height: 640 };

function shapeText(shape) {
  return shape?.text?.toString?.() ?? "";
}

function slideTexts(slide) {
  return slide.shapes.items.map(shapeText).filter(Boolean);
}

function isOldSupplementSlide(slide) {
  return slideTexts(slide).some((text) => /^FigS\d+\.\s*$/i.test(text.trim()));
}

function fitFrame(record) {
  const scale = Math.min(IMAGE_BOX.width / record.width, IMAGE_BOX.height / record.height);
  const width = record.width * scale;
  const height = record.height * scale;
  return {
    left: IMAGE_BOX.left + (IMAGE_BOX.width - width) / 2,
    top: IMAGE_BOX.top + (IMAGE_BOX.height - height) / 2,
    width,
    height,
  };
}

function setText(shape, text) {
  if (!shape) return;
  shape.text = text;
}

function styleSubtitle(shape) {
  try {
    shape.text.fontSize(18);
    shape.text.color("#365344");
    shape.text.typeface("Aptos");
  } catch {
    // Imported decks can expose limited text styling on some placeholder shapes.
  }
}

function addSubtitle(slide, text) {
  const subtitle = slide.shapes.add({
    geometry: "rect",
    position: { left: 58, top: 132, width: 215, height: 58 },
    fill: { color: "#FFFFFF", transparency: 100 },
    line: { color: "#FFFFFF", transparency: 100 },
  });
  subtitle.text = text;
  styleSubtitle(subtitle);
  return subtitle;
}

function findTitleShape(slide) {
  return slide.shapes.items.find((shape) => /^FigS\d+\.\s*$/i.test(shapeText(shape).trim()));
}

function findMarkerShape(slide) {
  return slide.shapes.items.find((shape) => /^\d+\s*$/.test(shapeText(shape).trim()));
}

function removeSlideById(presentation, slideId) {
  const index = presentation.slides.items.findIndex((slide) => slide.id === slideId);
  if (index >= 0) {
    presentation.slides.remove(index);
  }
}

const manifest = JSON.parse(await fs.readFile(MANIFEST, "utf8"));
const presentation = await PresentationFile.importPptx(await FileBlob.load(SOURCE));

const oldSupplementSlides = presentation.slides.items.filter(isOldSupplementSlide);
if (!oldSupplementSlides.length) {
  throw new Error("Could not find a FigS template slide in the source deck.");
}

const blankTemplate =
  oldSupplementSlides.find((slide) => slide.images.items.length === 0) ?? oldSupplementSlides[0];
const oldSupplementIds = oldSupplementSlides.map((slide) => slide.id);
const template = blankTemplate.duplicate();
template.moveTo(presentation.slides.items.length - 1);

for (const slideId of oldSupplementIds) {
  removeSlideById(presentation, slideId);
}

const appendedSlides = [];
for (const record of manifest) {
  const slide = template.duplicate();
  slide.moveTo(presentation.slides.items.length - 1);

  setText(findTitleShape(slide), `${record.label} `);
  addSubtitle(slide, record.title);

  const imageData = await fs.readFile(record.image);
  slide.images.add({
    data: imageData,
    contentType: "image/png",
    frame: fitFrame(record),
    alt: `${record.label} ${record.title}`,
    fit: "contain",
  });

  slide.speakerNotes.setText(
    [
      `${record.label} ${record.title}`,
      "",
      record.legend,
      "",
      `Source: Rich et al., Cell Systems 2026, supplemental PDF page ${record.sourcePage}.`,
    ].join("\n"),
  );
  appendedSlides.push(slide);
}

removeSlideById(presentation, template.id);

for (const slide of appendedSlides) {
  const marker = findMarkerShape(slide);
  const finalNumber = presentation.slides.items.findIndex((candidate) => candidate.id === slide.id) + 1;
  setText(marker, String(finalNumber));
}

for (const slide of presentation.slides.items) {
  slide.setViewportSize(SLIDE.width, SLIDE.height);
}

const pptx = await PresentationFile.exportPptx(presentation);
await pptx.save(OUT);

console.log(JSON.stringify({
  source: SOURCE,
  output: OUT,
  originalSlideCount: 22,
  removedPartialSupplementSlides: oldSupplementIds.length,
  appendedSupplementSlides: appendedSlides.length,
  finalSlideCount: presentation.slides.count,
}, null, 2));
