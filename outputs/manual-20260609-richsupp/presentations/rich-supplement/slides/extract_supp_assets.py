#!/usr/bin/env python3
import json
import re
import subprocess
from pathlib import Path

from PIL import Image, ImageChops


ROOT = Path("/Users/ashi/github/informatics_club")
PDF = ROOT / "sctalk/06_seurat-vs-scanpy/paper/supp.pdf"
WORK = ROOT / "outputs/manual-20260609-richsupp/presentations/rich-supplement"
RAW_DIR = WORK / "assets/supp_raw_450"
OUT_DIR = WORK / "assets/supp_crops"
TEXT_DIR = WORK / "assets/supp_page_text"
MANIFEST = WORK / "assets/supp_manifest.json"

DPI = 450

ITEMS = [
    ("fig", 1, 3, "Matrix workflow"),
    ("fig", 2, 5, "Aligned sequential"),
    ("fig", 3, 7, "Default controlled"),
    ("fig", 4, 9, "Aligned controlled"),
    ("fig", 5, 11, "Scanpy-aligned controlled"),
    ("fig", 6, 13, "Extended UMAP"),
    ("fig", 7, 15, "Read downsample: Seurat"),
    ("fig", 8, 17, "Read downsample: Scanpy"),
    ("fig", 9, 19, "Cell downsample: Seurat"),
    ("fig", 10, 21, "Cell downsample: Scanpy"),
    ("fig", 11, 23, "Read-downsampling metrics"),
    ("fig", 12, 25, "Cell-downsampling metrics"),
    ("fig", 13, 26, "Cell Ranger v7 vs v6"),
    ("fig", 14, 27, "UMI and feature counts"),
    ("fig", 15, 28, "Enrichr clusters 0-9"),
    ("fig", 16, 30, "Enrichr clusters 10-18"),
    ("fig", 17, 32, "Random seed variability"),
    ("fig", 18, 34, "KNN vs SNN graphs"),
    ("fig", 19, 35, "KNN/SNN workflow"),
    ("fig", 20, 37, "Mouse brain dataset"),
    ("fig", 21, 39, "NSCLC dataset"),
    ("fig", 22, 41, "Bootstrap seed metrics"),
    ("fig", 23, 43, "Bootstrap vs original"),
    ("table", 1, 45, "Seurat functions"),
    ("table", 2, 47, "Scanpy functions"),
]

CAPTION_CUT_FRACS = {
    25: 0.675,
    26: 0.715,
    27: 0.740,
    34: 0.395,
}


def run(cmd):
    subprocess.run(cmd, check=True)


def render_page(page):
    RAW_DIR.mkdir(parents=True, exist_ok=True)
    prefix = RAW_DIR / f"supp-{page:02d}"
    out_path = prefix.with_suffix(".png")
    if not out_path.exists():
        run([
            "pdftocairo",
            "-png",
            "-singlefile",
            "-r",
            str(DPI),
            "-f",
            str(page),
            "-l",
            str(page),
            str(PDF),
            str(prefix),
        ])
    return out_path


def crop_nonwhite(path, page, out_path):
    im = Image.open(path).convert("RGB")
    if page in CAPTION_CUT_FRACS:
        im = im.crop((0, 0, im.width, int(im.height * CAPTION_CUT_FRACS[page])))

    bg = Image.new("RGB", im.size, (255, 255, 255))
    diff = ImageChops.difference(im, bg).convert("L")
    mask = diff.point(lambda p: 255 if p > 12 else 0)
    bbox = mask.getbbox()
    if not bbox:
        bbox = (0, 0, im.width, im.height)

    pad = max(18, int(DPI * 0.04))
    left = max(0, bbox[0] - pad)
    top = max(0, bbox[1] - pad)
    right = min(im.width, bbox[2] + pad)
    bottom = min(im.height, bbox[3] + pad)
    cropped = im.crop((left, top, right, bottom))
    cropped.save(out_path, optimize=True)
    return {
        "rawWidth": im.width,
        "rawHeight": im.height,
        "cropBox": [left, top, right, bottom],
        "width": cropped.width,
        "height": cropped.height,
    }


def page_text(page):
    TEXT_DIR.mkdir(parents=True, exist_ok=True)
    out = TEXT_DIR / f"page-{page:02d}.txt"
    if not out.exists():
        run(["pdftotext", "-layout", "-f", str(page), "-l", str(page), str(PDF), str(out)])
    return out.read_text(errors="replace")


def clean_text(text):
    text = text.replace("\uFB00", "ff")
    text = text.replace("\uFB01", "fi")
    text = text.replace("\uFB02", "fl")
    text = text.replace("\u21B5", "ff")
    text = text.replace("↵", "ff")
    text = text.replace("−", "-")
    text = text.replace("≥", ">=")
    text = text.replace("≤", "<=")
    text = text.replace("θ", "theta")
    text = text.replace("×", "x")
    text = text.replace("α", "alpha")
    text = text.replace("β", "beta")
    text = text.replace("’", "'")
    text = text.replace("“", '"').replace("”", '"')
    text = re.sub(r"([A-Za-z])-\s+([A-Za-z])", r"\1\2", text)
    text = re.sub(r"\s+", " ", text)
    return text.strip()


def extract_legend(kind, number):
    label = "Figure" if kind == "fig" else "Table"
    pattern = re.compile(rf"Supplemental\s+{label}\s+{number}\s*:", re.IGNORECASE)
    pages = []
    for page in range(1, 49):
        text = page_text(page)
        match = pattern.search(text)
        if not match:
            continue
        end = len(text)
        next_match = re.search(r"Supplemental\s+(?:Figure|Table)\s+\d+\s*:", text[match.end():], re.IGNORECASE)
        if next_match:
            end = match.end() + next_match.start()
        pages.append(clean_text(text[match.start():end]))
    if not pages:
        return f"Supplemental {label} {number}: Legend not found in extracted PDF text."
    return max(pages, key=len)


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    records = []
    for kind, number, page, title in ITEMS:
        raw = render_page(page)
        out_name = f"{'figS' if kind == 'fig' else 'tableS'}{number:02d}.png"
        out_path = OUT_DIR / out_name
        crop_meta = crop_nonwhite(raw, page, out_path)
        label = f"FigS{number}." if kind == "fig" else f"Table S{number}."
        records.append({
            "kind": kind,
            "number": number,
            "label": label,
            "title": title,
            "sourcePage": page,
            "image": str(out_path),
            "legend": extract_legend(kind, number),
            **crop_meta,
        })

    MANIFEST.write_text(json.dumps(records, indent=2), encoding="utf-8")
    print(MANIFEST)
    print(f"wrote {len(records)} supplemental assets")


if __name__ == "__main__":
    main()
