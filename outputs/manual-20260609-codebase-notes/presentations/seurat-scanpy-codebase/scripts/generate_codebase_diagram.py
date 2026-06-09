#!/usr/bin/env python3
from pathlib import Path
from PIL import Image, ImageDraw, ImageFont

OUT = Path("/Users/ashi/github/informatics_club/outputs/manual-20260609-codebase-notes/presentations/seurat-scanpy-codebase/assets/codebase_structure.png")

W, H = 2200, 1040
GREEN = "#145A32"
BODY = "#29483A"
MUTED = "#66756C"
ORANGE = "#D66A1F"
BLUE = "#2B78B8"
YELLOW = "#F2B705"
RED = "#B81E14"
LIGHT = "#F5F8F6"
LINE = "#D8E2DC"
BLACK = "#111111"
WHITE = "#FFFFFF"


def font(size, bold=False):
    candidates = [
        "/System/Library/Fonts/Supplemental/Arial Bold.ttf" if bold else "/System/Library/Fonts/Supplemental/Arial.ttf",
        "/System/Library/Fonts/Helvetica.ttc",
        "/Library/Fonts/Arial.ttf",
    ]
    for candidate in candidates:
        try:
            return ImageFont.truetype(candidate, size)
        except OSError:
            continue
    return ImageFont.load_default()


F_TITLE = font(54, True)
F_SUB = font(31, False)
F_HEAD = font(34, True)
F_BODY = font(26, False)
F_SMALL = font(22, False)
F_MONO = font(22, False)


def rounded(draw, xy, radius=22, fill=WHITE, outline=LINE, width=3):
    draw.rounded_rectangle(xy, radius=radius, fill=fill, outline=outline, width=width)


def text(draw, xy, value, fill=BLACK, fnt=F_BODY, anchor=None):
    draw.text(xy, value, fill=fill, font=fnt, anchor=anchor)


def wrap(draw, value, max_width, fnt):
    words = value.split()
    lines, line = [], ""
    for word in words:
        test = word if not line else f"{line} {word}"
        if draw.textbbox((0, 0), test, font=fnt)[2] <= max_width:
            line = test
        else:
            if line:
                lines.append(line)
            line = word
    if line:
        lines.append(line)
    return lines


def bullet_block(draw, x, y, lines, color=BODY, fnt=F_BODY, leading=36):
    for item in lines:
        draw.ellipse((x, y + 9, x + 9, y + 18), fill=color)
        for j, line in enumerate(wrap(draw, item, 620, fnt)):
            text(draw, (x + 22, y + j * leading), line, fill=color, fnt=fnt)
        y += leading * max(1, len(wrap(draw, item, 620, fnt))) + 8
    return y


def arrow(draw, start, end, fill=GREEN, width=7):
    draw.line([start, end], fill=fill, width=width)
    sx, sy = start
    ex, ey = end
    if abs(ex - sx) >= abs(ey - sy):
        direction = 1 if ex > sx else -1
        pts = [(ex, ey), (ex - 26 * direction, ey - 15), (ex - 26 * direction, ey + 15)]
    else:
        direction = 1 if ey > sy else -1
        pts = [(ex, ey), (ex - 15, ey - 26 * direction), (ex + 15, ey - 26 * direction)]
    draw.polygon(pts, fill=fill)


img = Image.new("RGB", (W, H), WHITE)
draw = ImageDraw.Draw(img)

text(draw, (70, 56), "Reproducibility codebase: two linked workflows", fill=GREEN, fnt=F_TITLE)
text(draw, (72, 122), "Zenodo archive for Rich et al. organizes FASTQ-to-matrix generation separately from figure-level analysis.", fill=MUTED, fnt=F_SUB)

left = (80, 250, 1025, 790)
right = (1175, 250, 2120, 790)
middle = (820, 830, 1380, 958)

rounded(draw, left, fill="#F7FBF8", outline="#B8D7C2", width=4)
rounded(draw, right, fill="#F8FAFE", outline="#BCD1EB", width=4)
rounded(draw, middle, fill="#FFF9E8", outline="#E6C567", width=4)

text(draw, (125, 290), "matrix_generation/", fill=GREEN, fnt=F_HEAD)
text(draw, (125, 331), "FASTQ and reference files to count matrices", fill=MUTED, fnt=F_SMALL)

bullet_block(
    draw,
    130,
    380,
    [
        "main.py reads config.yaml, dispatches setup, downloads, downsampling, kb count, and Cell Ranger count",
        "fastq_processor.py stores project paths and routes methods to src/",
        "src/: project_setup.py, downsample_fastqs.py, kb.py, cellranger_count.py, read_counts.py",
        "config_files/ holds kb28, Cell Ranger v6/v7, and v7-like-old settings",
        "scripts/ bootstrap, align data, and organize outputs",
    ],
    color=BODY,
)

text(draw, (1220, 290), "analysis/", fill=BLUE, fnt=F_HEAD)
text(draw, (1220, 331), "Count matrices to figures, statistics, and version comparisons", fill=MUTED, fnt=F_SMALL)

bullet_block(
    draw,
    1225,
    380,
    [
        "rmd/ notebooks run Seurat-vs-Scanpy, version comparisons, downsampling summaries, and bootstrap plots",
        "yaml/ has one configuration per figure or supplemental figure family",
        "scripts/ contains plotting, stats, DE helper, scanpy14, and Ensembl conversion code",
        "env/ defines Docker plus Scanpy 1.9, Scanpy 1.4, and UMAP-version conda environments",
        "official_tutorials/ stores PBMC reference notebooks",
    ],
    color=BODY,
)

arrow(draw, (1025, 520), (1175, 520), fill=GREEN)
text(draw, (1100, 480), "matrices", fill=GREEN, fnt=F_SMALL, anchor="mm")

text(draw, (1100, 870), "outputs", fill="#745000", fnt=F_HEAD, anchor="mm")
text(draw, (1100, 918), "figures, metrics, DE tables, Enrichr/pathway summaries", fill="#745000", fnt=F_SMALL, anchor="mm")

chips = [
    ("43 YAML", YELLOW),
    ("22 Python", BLUE),
    ("9 Rmd", ORANGE),
    ("7 R scripts", RED),
    ("Docker + conda envs", GREEN),
]
x = 72
for label, color in chips:
    w = draw.textbbox((0, 0), label, font=F_MONO)[2] + 42
    draw.rounded_rectangle((x, 166, x + w, 226), radius=26, fill=color)
    text(draw, (x + 21, 184), label, fill=WHITE, fnt=F_MONO)
    x += w + 22

text(draw, (72, 985), "Reading the repo: choose a figure YAML -> run the matching Rmd -> helper scripts perform plotting/stats/DE; use matrix_generation only when regenerating matrices from FASTQs.", fill=MUTED, fnt=F_SMALL)

OUT.parent.mkdir(parents=True, exist_ok=True)
img.save(OUT)
print(OUT)
