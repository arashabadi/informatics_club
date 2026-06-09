from __future__ import annotations

import math
import textwrap
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont


OUT = Path("/Users/ashi/github/informatics_club/outputs/manual-20260609-formula-explainers/presentations/seurat-scanpy-formulas/assets")
W, H = 1280, 720
GREEN = "#145A32"
DARK = "#1F2A24"
BODY = "#344A3E"
MUTED = "#6E7F75"
LIGHT = "#F6F8F6"
LINE = "#C9D5CC"
ORANGE = "#D97A2B"
BLUE = "#2E73B8"
RED = "#B74242"
GRAY = "#EEF2EF"
WHITE = "#FFFFFF"


def font(size: int, bold: bool = False) -> ImageFont.FreeTypeFont:
    candidates = [
        "/System/Library/Fonts/Supplemental/Arial Bold.ttf" if bold else "/System/Library/Fonts/Supplemental/Arial.ttf",
        "/Library/Fonts/Arial Bold.ttf" if bold else "/Library/Fonts/Arial.ttf",
        "/System/Library/Fonts/Supplemental/DejaVu Sans Bold.ttf" if bold else "/System/Library/Fonts/Supplemental/DejaVu Sans.ttf",
    ]
    for candidate in candidates:
        try:
            return ImageFont.truetype(candidate, size)
        except OSError:
            pass
    return ImageFont.load_default()


F_TITLE = font(34, True)
F_SUB = font(18)
F_H = font(22, True)
F_BODY = font(18)
F_SMALL = font(14)
F_MONO = font(18)
F_MONO_B = font(19, True)


def draw_wrapped(draw: ImageDraw.ImageDraw, text: str, box, fill=BODY, fnt=F_BODY, width=48, spacing=6):
    x, y, w, _ = box
    lines = []
    for para in text.split("\n"):
        if not para:
            lines.append("")
        else:
            lines.extend(textwrap.wrap(para, width=width, break_long_words=False))
    line_h = fnt.getbbox("Ag")[3] - fnt.getbbox("Ag")[1] + spacing
    for line in lines:
        draw.text((x, y), line, fill=fill, font=fnt)
        y += line_h
    return y


def card(draw, xy, title: str | None = None):
    draw.rounded_rectangle(xy, radius=8, fill=WHITE, outline=LINE, width=2)
    if title:
        x1, y1, x2, _ = xy
        draw.text((x1 + 18, y1 + 14), title, fill=GREEN, font=F_H)


def footer(draw):
    draw.rectangle((0, 672, W, H), fill=GREEN)
    draw.text((42, 685), "Formula explainer for biologists", fill=WHITE, font=font(16, True))
    draw.text((510, 688), "Rich et al. single-cell workflow comparison", fill="#DCEBE0", font=font(13))
    draw.text((1038, 688), "Minimal count-matrix examples", fill="#DCEBE0", font=font(13))


def base(title: str, subtitle: str):
    img = Image.new("RGB", (W, H), "#FBFCFB")
    draw = ImageDraw.Draw(img)
    draw.text((42, 28), title, fill=DARK, font=F_TITLE)
    draw_wrapped(draw, subtitle, (42, 74, 1050, 50), fill=MUTED, fnt=F_SUB, width=110, spacing=4)
    draw.line((42, 116, 1238, 116), fill=LINE, width=2)
    footer(draw)
    return img, draw


def draw_table(draw, x, y, headers, rows, cell_w=82, cell_h=34, header_fill=GREEN, highlight=None):
    highlight = highlight or set()
    for j, h in enumerate([""] + headers):
        xx = x + j * cell_w
        draw.rectangle((xx, y, xx + cell_w, y + cell_h), fill=header_fill if j else "#E9F1EC", outline=LINE)
        draw.text((xx + 8, y + 8), h, fill=WHITE if j else BODY, font=F_SMALL)
    for i, (name, values) in enumerate(rows):
        yy = y + (i + 1) * cell_h
        draw.rectangle((x, yy, x + cell_w, yy + cell_h), fill="#E9F1EC", outline=LINE)
        draw.text((x + 8, yy + 8), name, fill=BODY, font=F_SMALL)
        for j, val in enumerate(values):
            xx = x + (j + 1) * cell_w
            fill = "#FFF0DF" if (i, j) in highlight else WHITE
            draw.rectangle((xx, yy, xx + cell_w, yy + cell_h), fill=fill, outline=LINE)
            draw.text((xx + 28, yy + 8), str(val), fill=DARK, font=F_SMALL)


def draw_set(draw, center, items, color, label):
    x, y = center
    draw.ellipse((x - 92, y - 62, x + 92, y + 62), fill=color, outline=LINE, width=2)
    draw.text((x - 72, y - 82), label, fill=DARK, font=F_SMALL)
    draw.text((x - 58, y - 14), "\n".join(items), fill=DARK, font=F_SMALL, spacing=2)


def slide_01():
    img, d = base(
        "Formula explainer 1: normalization and set overlap",
        "Figure 1A uses a count matrix, then asks how much the retained cell/gene/HVG sets overlap.",
    )
    card(d, (42, 140, 428, 628), "1) Start with a tiny count matrix")
    draw_table(
        d,
        66,
        190,
        ["GeneA", "GeneB", "GeneC"],
        [("Cell1", [10, 0, 5]), ("Cell2", [0, 8, 2]), ("Cell3", [4, 1, 0])],
        cell_w=78,
        cell_h=38,
        highlight={(0, 0)},
    )
    draw_wrapped(
        d,
        "For one entry: Cell1 has total counts 15. GeneA has 10 counts.",
        (66, 355, 325, 90),
        width=35,
    )
    d.text((66, 466), "Y_ig = log((X_ig / total_i) * 10000 + 1)", fill=DARK, font=F_MONO_B)
    d.text((66, 506), "GeneA in Cell1 = log((10/15)*10000 + 1) = 8.81", fill=BLUE, font=F_MONO)
    card(d, (462, 140, 830, 628), "2) Compare selected gene sets")
    draw_set(d, (590, 300), ["G1", "G2", "G3"], "#F6CFA8", "Seurat HVGs")
    draw_set(d, (700, 300), ["G2", "G3", "G4"], "#BFD5F0", "Scanpy HVGs")
    d.text((520, 420), "Shared = {G2, G3} = 2 genes", fill=DARK, font=F_BODY)
    d.text((520, 454), "Union = {G1, G2, G3, G4} = 4 genes", fill=DARK, font=F_BODY)
    d.text((520, 506), "Jaccard = shared / union = 2 / 4 = 0.50", fill=GREEN, font=F_MONO_B)
    card(d, (864, 140, 1238, 628), "Biologist takeaway")
    draw_wrapped(
        d,
        "Normalization rescales each cell so cells with different library sizes can be compared. Jaccard then asks a simple question: are the two workflows keeping the same objects? A Jaccard of 1 means identical sets; 0 means no overlap.",
        (890, 195, 312, 240),
        width=35,
    )
    d.rounded_rectangle((900, 470, 1205, 555), radius=8, fill="#EAF4EE", outline="#A4BCAF")
    d.text((922, 492), "Figure 1A: cells, genes, HVGs", fill=GREEN, font=F_H)
    d.text((922, 528), "same formula, different object set", fill=BODY, font=F_BODY)
    return img


def slide_02():
    img, d = base(
        "Formula explainer 2: PCA angles and SNN neighborhoods",
        "Figure 1B-C asks whether cells point in the same PC directions and keep the same graph neighbors.",
    )
    card(d, (42, 140, 424, 628), "1) PCA vector angle")
    d.line((108, 500, 108, 210), fill=LINE, width=3)
    d.line((108, 500, 350, 500), fill=LINE, width=3)
    d.line((108, 500, 308, 340), fill=ORANGE, width=7)
    d.line((108, 500, 275, 300), fill=BLUE, width=7)
    d.text((286, 350), "Seurat PC1", fill=ORANGE, font=F_SMALL)
    d.text((236, 280), "Scanpy PC1", fill=BLUE, font=F_SMALL)
    d.text((72, 178), "cos(theta) = dot(v1, v2)", fill=DARK, font=F_MONO_B)
    d.text((72, 218), "sin(theta) = sqrt(1 - cos(theta)^2)", fill=DARK, font=F_MONO_B)
    d.text((72, 555), "Example: cos = 0.96 -> sin = 0.28", fill=GREEN, font=F_MONO_B)
    card(d, (458, 140, 836, 628), "2) One cell's neighbors")
    d.text((490, 195), "Cell C1 neighbors", fill=DARK, font=F_H)
    d.text((490, 240), "Seurat: {C2, C3, C4}", fill=ORANGE, font=F_MONO_B)
    d.text((490, 278), "Scanpy: {C3, C4, C5, C6}", fill=BLUE, font=F_MONO_B)
    d.text((490, 350), "Shared = {C3, C4} = 2", fill=DARK, font=F_BODY)
    d.text((490, 384), "Union = {C2, C3, C4, C5, C6} = 5", fill=DARK, font=F_BODY)
    d.text((490, 435), "Jaccard = 2 / 5 = 0.40", fill=GREEN, font=F_MONO_B)
    d.text((490, 488), "Degree ratio = 3 / 4 = 0.75", fill=GREEN, font=F_MONO_B)
    d.text((490, 526), "|log(0.75)| = 0.29", fill=GREEN, font=F_MONO_B)
    card(d, (870, 140, 1238, 628), "How to read the figure")
    draw_wrapped(
        d,
        "For PCA, sin(theta) near 0 means the same direction; near 1 means almost perpendicular. For the SNN graph, Jaccard measures neighbor overlap, while degree ratio asks whether one graph gives a cell more edges than the other.",
        (896, 195, 310, 250),
        width=35,
    )
    d.rounded_rectangle((906, 486, 1198, 555), radius=8, fill="#EFF5FA", outline="#ABC3DD")
    d.text((928, 506), "Figure 1B-C", fill=BLUE, font=F_H)
    d.text((928, 538), "same cells, different geometry", fill=BODY, font=F_BODY)
    return img


def slide_03():
    img, d = base(
        "Formula explainer 3: cluster agreement and marker overlap",
        "Figure 1D-E turns geometry into cluster assignments and marker-gene sets.",
    )
    card(d, (42, 140, 420, 628), "1) ARI: pair agreement")
    draw_wrapped(
        d,
        "Adjusted Rand Index checks cell pairs. If two cells are together in Seurat, are they also together in Scanpy? Then it corrects for agreement expected by chance.",
        (70, 188, 300, 170),
        width=34,
    )
    d.text((70, 384), "ARI = 1.0  perfect same grouping", fill=GREEN, font=F_MONO_B)
    d.text((70, 426), "ARI = 0.0  chance-like agreement", fill=BODY, font=F_MONO_B)
    d.text((70, 468), "ARI < 0    worse than chance", fill=RED, font=F_MONO_B)
    d.rounded_rectangle((88, 530, 350, 585), radius=6, fill="#EAF4EE", outline="#A4BCAF")
    d.text((110, 548), "Figure 1D reports ARI = 0.71", fill=GREEN, font=F_BODY)
    card(d, (454, 140, 836, 628), "2) Marker set Jaccard")
    d.text((486, 196), "Seurat markers:", fill=ORANGE, font=F_H)
    d.text((486, 234), "{IL7R, CCR7, MS4A1}", fill=ORANGE, font=F_MONO_B)
    d.text((486, 292), "Scanpy markers:", fill=BLUE, font=F_H)
    d.text((486, 330), "{IL7R, MS4A1, LST1}", fill=BLUE, font=F_MONO_B)
    d.text((486, 410), "Shared = {IL7R, MS4A1} = 2", fill=DARK, font=F_BODY)
    d.text((486, 448), "Union = 4 genes", fill=DARK, font=F_BODY)
    d.text((486, 502), "Jaccard = 2 / 4 = 0.50", fill=GREEN, font=F_MONO_B)
    card(d, (870, 140, 1238, 628), "3) p-value flips")
    draw_table(
        d,
        900,
        208,
        ["Seurat", "Scanpy"],
        [("Gene A", [0.03, 0.20]), ("Gene B", [0.80, 0.02])],
        cell_w=95,
        cell_h=42,
    )
    d.text((900, 330), "threshold: adjusted p < 0.05", fill=DARK, font=F_BODY)
    d.text((900, 385), "Gene A is significant only in Seurat.", fill=ORANGE, font=F_BODY)
    d.text((900, 426), "Gene B is significant only in Scanpy.", fill=BLUE, font=F_BODY)
    draw_wrapped(
        d,
        "A flip means the biological story can change even for the same gene and same cluster.",
        (900, 490, 285, 90),
        width=32,
    )
    return img


def slide_04():
    img, d = base(
        "Formula explainer 4: downsampling as a benchmark",
        "Figure 2 asks: how much data loss inside one workflow matches the Seurat-vs-Scanpy difference?",
    )
    card(d, (42, 140, 432, 628), "1) Retain a fraction of data")
    draw_table(
        d,
        78,
        200,
        ["GeneA", "GeneB"],
        [("Cell1", [12, 4]), ("Cell2", [2, 10]), ("Cell3", [8, 0]), ("Cell4", [1, 3])],
        cell_w=80,
        cell_h=36,
    )
    d.text((78, 390), "Cell downsampling f = 0.50", fill=GREEN, font=F_MONO_B)
    d.text((78, 430), "Keep 2 of 4 cells -> rerun pipeline", fill=DARK, font=F_BODY)
    d.text((78, 482), "Read downsampling f = 0.04", fill=GREEN, font=F_MONO_B)
    d.text((78, 522), "Keep 4% of reads -> rebuild matrix", fill=DARK, font=F_BODY)
    card(d, (466, 140, 846, 628), "2) Compare to full data")
    d.text((492, 206), "Example metric: HVG Jaccard", fill=DARK, font=F_H)
    d.text((492, 260), "Full Seurat vs Scanpy = 0.60", fill=BLUE, font=F_MONO_B)
    d.text((492, 318), "Full vs downsampled at f=0.08 = 0.58", fill=ORANGE, font=F_MONO_B)
    d.text((492, 376), "Difference from benchmark = 0.02", fill=DARK, font=F_MONO_B)
    d.rounded_rectangle((492, 442, 808, 525), radius=8, fill="#EAF4EE", outline="#A4BCAF")
    d.text((516, 462), "Within 5% margin", fill=GREEN, font=F_H)
    d.text((516, 498), "so f=0.08 is enough for this metric", fill=BODY, font=F_SMALL)
    card(d, (880, 140, 1238, 628), "3) Bootstrap null idea")
    d.text((906, 200), "ideal metric = 1.00", fill=DARK, font=F_MONO_B)
    d.text((906, 248), "Seurat bootstrap: 0.92 -> eps_X = -0.08", fill=ORANGE, font=F_SMALL)
    d.text((906, 286), "Scanpy bootstrap: 0.95 -> eps_Y = -0.05", fill=BLUE, font=F_SMALL)
    d.text((906, 350), "combined noise:", fill=DARK, font=F_H)
    d.text((906, 390), "W = ideal + eps_X + eps_Y", fill=GREEN, font=F_MONO_B)
    d.text((906, 432), "W = 1 - 0.08 - 0.05 = 0.87", fill=GREEN, font=F_MONO_B)
    draw_wrapped(
        d,
        "This is the paper's way to judge whether an observed package difference is larger than expected resampling noise.",
        (906, 492, 286, 100),
        width=31,
    )
    return img


def slide_05():
    img, d = base(
        "Formula explainer 5: logFC, pseudocounts, and adjusted p values",
        "Figure 3 is mostly about DE math: small implementation choices can reshape marker calls.",
    )
    card(d, (42, 140, 432, 628), "1) Minimal logFC example")
    d.text((70, 194), "One gene, one cluster vs rest", fill=DARK, font=F_H)
    draw_table(
        d,
        70,
        240,
        ["Group 1", "Group 2"],
        [("mean count", [0, 1])],
        cell_w=105,
        cell_h=42,
    )
    d.text((70, 330), "logFC = log2((a1 + p1) / (a2 + p2))", fill=DARK, font=F_MONO_B)
    d.text((70, 388), "Seurat-like p: 0.01 / 0.0001", fill=ORANGE, font=F_SMALL)
    d.text((70, 422), "log2(0.01 / 1.0001) = -6.6", fill=ORANGE, font=F_MONO_B)
    d.text((70, 480), "Scanpy epsilon p: 1e-9", fill=BLUE, font=F_SMALL)
    d.text((70, 514), "log2(1e-9 / 1) = -29.9", fill=BLUE, font=F_MONO_B)
    card(d, (466, 140, 846, 628), "2) CCC: do logFC values agree?")
    d.text((492, 202), "rho_c = 2*rho*sx*sy /", fill=DARK, font=F_MONO_B)
    d.text((566, 232), "(sx^2 + sy^2 + (mux - muy)^2)", fill=DARK, font=F_MONO_B)
    d.text((492, 298), "High CCC needs:", fill=GREEN, font=F_H)
    d.text((516, 342), "1. correlation", fill=BODY, font=F_BODY)
    d.text((516, 380), "2. similar spread", fill=BODY, font=F_BODY)
    d.text((516, 418), "3. close to y = x", fill=BODY, font=F_BODY)
    d.rounded_rectangle((500, 492, 810, 552), radius=8, fill="#EFF5FA", outline="#ABC3DD")
    d.text((520, 510), "CCC near 1 = values agree", fill=BLUE, font=F_BODY)
    card(d, (880, 140, 1238, 628), "3) p-value correction")
    draw_table(
        d,
        908,
        202,
        ["raw p", "Bonf.", "BH"],
        [("Gene A", [0.02, 0.08, 0.04]), ("Gene B", [0.001, 0.004, 0.004])],
        cell_w=82,
        cell_h=40,
    )
    draw_wrapped(
        d,
        "Bonferroni multiplies by the number of tests. Benjamini-Hochberg uses the p-value rank. Different corrections and Wilcoxon tie handling can move genes across the 0.05 cutoff.",
        (908, 340, 278, 180),
        width=31,
    )
    d.text((908, 548), "threshold: adjusted p < 0.05", fill=GREEN, font=F_MONO_B)
    return img


def slide_06():
    img, d = base(
        "Formula explainer 6: volcano thresholds and pathway input",
        "Figure 4 uses DE outputs to decide which genes and pathways enter the biological story.",
    )
    card(d, (42, 140, 440, 628), "1) Volcano plot rules")
    d.text((72, 196), "x-axis: logFC", fill=DARK, font=F_H)
    d.text((72, 234), "y-axis: -log10(adjusted p)", fill=DARK, font=F_H)
    d.text((72, 292), "biological threshold:", fill=GREEN, font=F_BODY)
    d.text((72, 330), "|logFC| >= 1", fill=GREEN, font=F_MONO_B)
    d.text((72, 380), "statistical threshold:", fill=GREEN, font=F_BODY)
    d.text((72, 418), "p < 0.05  ->  -log10(p) > 1.30", fill=GREEN, font=F_MONO_B)
    d.rounded_rectangle((72, 492, 370, 558), radius=8, fill="#FFF0DF", outline="#E0B17E")
    d.text((92, 513), "Red or blue = passes both", fill=ORANGE, font=F_BODY)
    card(d, (474, 140, 844, 628), "2) Tiny gene examples")
    draw_table(
        d,
        500,
        204,
        ["logFC", "p", "call"],
        [("Gene A", [1.4, 0.01, "up"]), ("Gene B", [-1.2, 0.03, "down"]), ("Gene C", [0.5, 0.001, "not bio."]), ("Gene D", [2.1, 0.20, "not stat."])],
        cell_w=82,
        cell_h=40,
    )
    draw_wrapped(
        d,
        "A gene needs both enough effect size and enough statistical support. That is why a very small p-value alone is not automatically a biologically highlighted marker here.",
        (500, 408, 300, 120),
        width=34,
    )
    card(d, (878, 140, 1238, 628), "3) Pathway comparison")
    d.text((908, 200), "Input to Enrichr:", fill=DARK, font=F_H)
    d.text((908, 238), "genes passing both thresholds", fill=GREEN, font=F_BODY)
    d.text((908, 304), "Seurat genes -> pathways", fill=ORANGE, font=F_MONO_B)
    d.text((908, 346), "Scanpy genes -> pathways", fill=BLUE, font=F_MONO_B)
    d.text((908, 420), "Then compare pathway sets:", fill=DARK, font=F_BODY)
    d.text((908, 462), "shared / union", fill=GREEN, font=F_MONO_B)
    draw_wrapped(
        d,
        "Different marker lists can produce different pathway lists, even when many broad cell types look similar.",
        (908, 522, 280, 80),
        width=31,
    )
    return img


def slide_07():
    img, d = base(
        "Formula explainer 1b: HVG scoring math",
        "The Methods explain why different HVG algorithms can choose different genes after the same normalization.",
    )
    card(d, (42, 140, 432, 628), "1) Variance relative to mean")
    draw_table(
        d,
        70,
        198,
        ["Cell1", "Cell2", "Cell3"],
        [("GeneA", [10, 0, 5]), ("GeneB", [2, 3, 4])],
        cell_w=82,
        cell_h=42,
        highlight={(0, 0), (0, 1), (0, 2)},
    )
    d.text((70, 315), "GeneA mean = 5.0", fill=DARK, font=F_BODY)
    d.text((70, 350), "GeneA variance = 25.0", fill=DARK, font=F_BODY)
    d.text((70, 404), "dispersion d_g = variance / mean", fill=GREEN, font=F_MONO_B)
    d.text((70, 450), "d_A = 25 / 5 = 5.0", fill=ORANGE, font=F_MONO_B)
    draw_wrapped(
        d,
        "A gene can look variable just because it has a higher mean. HVG methods try to correct for that mean-variance relationship.",
        (70, 508, 310, 88),
        width=34,
    )
    card(d, (466, 140, 846, 628), "2) mean.var.plot / seurat")
    d.text((492, 202), "Put genes with similar mean", fill=DARK, font=F_BODY)
    d.text((492, 234), "into the same bin.", fill=DARK, font=F_BODY)
    d.text((492, 292), "z_d = (log(d_g) - mean_bin) / sd_bin", fill=GREEN, font=F_MONO_B)
    d.text((492, 352), "Toy bin:", fill=DARK, font=F_H)
    d.text((516, 392), "log(d_A) = 1.61", fill=ORANGE, font=F_MONO)
    d.text((516, 428), "mean_bin = 0.50", fill=BODY, font=F_MONO)
    d.text((516, 464), "sd_bin = 0.55", fill=BODY, font=F_MONO)
    d.text((516, 518), "z_d = 2.0 -> strong HVG", fill=GREEN, font=F_MONO_B)
    card(d, (880, 140, 1238, 628), "3) vst / seurat_v3")
    d.text((906, 202), "Fit expected variance", fill=DARK, font=F_BODY)
    d.text((906, 234), "from mean expression.", fill=DARK, font=F_BODY)
    d.text((906, 292), "z_ig = (X_ig - mean_g) / sqrt(var_pred)", fill=GREEN, font=F_MONO_B)
    d.text((906, 352), "Then rank genes by", fill=DARK, font=F_BODY)
    d.text((906, 384), "variance of the z values.", fill=DARK, font=F_BODY)
    d.text((906, 444), "Toy gene:", fill=DARK, font=F_H)
    d.text((930, 484), "observed var = 25", fill=ORANGE, font=F_MONO)
    d.text((930, 520), "predicted var = 5", fill=BODY, font=F_MONO)
    d.text((930, 562), "normalized var = 5.0", fill=GREEN, font=F_MONO_B)
    return img


def slide_08():
    img, d = base(
        "Formula explainer 2b: graph limits and modularity",
        "The graph formulas show why neighborhood size and clustering optimizers can change downstream groups.",
    )
    card(d, (42, 140, 432, 628), "1) Maximum possible Jaccard")
    d.text((70, 196), "Cell C1 degree:", fill=DARK, font=F_H)
    d.text((70, 244), "Seurat d1 = 3 neighbors", fill=ORANGE, font=F_MONO_B)
    d.text((70, 286), "Scanpy d2 = 8 neighbors", fill=BLUE, font=F_MONO_B)
    draw_wrapped(
        d,
        "Even if every Seurat neighbor is inside the larger Scanpy neighborhood, the union still has 8 cells.",
        (70, 342, 310, 100),
        width=34,
    )
    d.text((70, 470), "J_max = smaller degree / larger degree", fill=GREEN, font=F_MONO_B)
    d.text((70, 518), "J_max = 3 / 8 = 0.38", fill=GREEN, font=F_MONO_B)
    card(d, (466, 140, 846, 628), "2) Exact KNN cost")
    d.text((492, 202), "Exact KNN checks distances", fill=DARK, font=F_H)
    d.text((492, 238), "between many cell pairs.", fill=DARK, font=F_BODY)
    d.text((492, 304), "time ~ O(d * n^2)", fill=GREEN, font=F_MONO_B)
    d.text((492, 368), "Tiny example:", fill=DARK, font=F_H)
    d.text((516, 408), "n = 5 cells", fill=BODY, font=F_MONO)
    d.text((516, 444), "d = 2 PCs", fill=BODY, font=F_MONO)
    d.text((516, 480), "about 5 * 4 * 2 = 40", fill=ORANGE, font=F_MONO_B)
    d.text((516, 516), "coordinate comparisons", fill=ORANGE, font=F_MONO_B)
    draw_wrapped(
        d,
        "Large datasets use approximate searches such as Annoy or PyNNDescent.",
        (492, 560, 300, 50),
        width=34,
    )
    card(d, (880, 140, 1238, 628), "3) Modularity quality")
    d.text((906, 202), "Clustering rewards edges", fill=DARK, font=F_BODY)
    d.text((906, 234), "inside a community.", fill=DARK, font=F_BODY)
    d.text((906, 292), "quality = observed inside edges", fill=GREEN, font=F_MONO_B)
    d.text((986, 326), "- expected random edges", fill=GREEN, font=F_MONO_B)
    d.text((906, 390), "Partition A: 10 - 6 = 4", fill=ORANGE, font=F_MONO_B)
    d.text((906, 434), "Partition B: 10 - 8 = 2", fill=BLUE, font=F_MONO_B)
    draw_wrapped(
        d,
        "Seurat and Scanpy use different optimizer details, so the preferred partition can change even when the graph is similar.",
        (906, 500, 285, 95),
        width=31,
    )
    return img


SLIDES = [
    ("formula_01_normalization_jaccard.png", slide_01),
    ("formula_02_pca_snn.png", slide_02),
    ("formula_03_clustering_de_overlap.png", slide_03),
    ("formula_04_downsampling_bootstrap.png", slide_04),
    ("formula_05_logfc_pvalues.png", slide_05),
    ("formula_06_volcano_pathways.png", slide_06),
    ("formula_07_hvg_scoring.png", slide_07),
    ("formula_08_graph_modularity.png", slide_08),
]


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    for filename, maker in SLIDES:
        image = maker()
        image.save(OUT / filename, quality=95)
        print(OUT / filename)


if __name__ == "__main__":
    main()
