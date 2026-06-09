from __future__ import annotations

import textwrap
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont


OUT = Path(
    "/Users/ashi/github/informatics_club/outputs/manual-20260609-formula-explainers/"
    "presentations/seurat-scanpy-formulas/assets/logfc_discussion_summary.png"
)
W, H = 1280, 720
GREEN = "#145A32"
DARK = "#1F2A24"
BODY = "#33473B"
MUTED = "#6E7F75"
LINE = "#C9D5CC"
ORANGE = "#D97A2B"
BLUE = "#2E73B8"
RED = "#B74242"
WHITE = "#FFFFFF"
LIGHT = "#F6F8F6"


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
F_MONO = font(17)
F_MONO_B = font(18, True)


def wrapped(draw: ImageDraw.ImageDraw, text: str, box, fill=BODY, fnt=F_BODY, width=40, spacing=5):
    x, y, _, _ = box
    for para in text.split("\n"):
        lines = [""] if not para else textwrap.wrap(para, width=width, break_long_words=False)
        for line in lines:
            draw.text((x, y), line, fill=fill, font=fnt)
            y += fnt.getbbox("Ag")[3] - fnt.getbbox("Ag")[1] + spacing
    return y


def card(draw: ImageDraw.ImageDraw, xy, title: str):
    draw.rounded_rectangle(xy, radius=8, fill=WHITE, outline=LINE, width=2)
    draw.text((xy[0] + 18, xy[1] + 14), title, fill=GREEN, font=F_H)


def pill(draw: ImageDraw.ImageDraw, xy, text: str, fill: str, outline: str, color: str):
    draw.rounded_rectangle(xy, radius=8, fill=fill, outline=outline)
    tw = draw.textlength(text, font=F_BODY)
    draw.text((xy[0] + (xy[2] - xy[0] - tw) / 2, xy[1] + 15), text, fill=color, font=F_BODY)


def table(draw: ImageDraw.ImageDraw, x: int, y: int):
    widths = [92, 120, 126]
    headers = ["Method", "pseudocount", "effect"]
    rows = [
        ("Seurat v5", "1 / n_group", "cluster size affects logFC"),
        ("Seurat v4", "+1 after mean", "compresses logFC toward 0"),
        ("Scanpy", "epsilon = 1e-9", "can create +/-30 edge cases"),
    ]
    h = 42
    xx = x
    for w, header in zip(widths, headers):
        draw.rectangle((xx, y, xx + w, y + h), fill=GREEN, outline=LINE)
        draw.text((xx + 8, y + 10), header, fill=WHITE, font=F_SMALL)
        xx += w
    for i, row in enumerate(rows):
        yy = y + (i + 1) * h
        xx = x
        for j, (w, val) in enumerate(zip(widths, row)):
            draw.rectangle((xx, yy, xx + w, yy + h), fill=LIGHT if j == 0 else WHITE, outline=LINE)
            fill = ORANGE if "Seurat" in row[0] and j == 0 else BLUE if row[0] == "Scanpy" and j == 0 else BODY
            wrapped(draw, val, (xx + 8, yy + 8, w - 16, h - 10), fill=fill, fnt=F_SMALL, width=18, spacing=1)
            xx += w


def main() -> None:
    img = Image.new("RGB", (W, H), "#FBFCFB")
    d = ImageDraw.Draw(img)
    d.text((42, 28), "LogFC discussion: why effect sizes diverge", fill=DARK, font=F_TITLE)
    wrapped(
        d,
        "The discussion says Seurat and Scanpy are not just reporting the same log fold-change with different software wrappers.",
        (42, 76, 1100, 40),
        fill=MUTED,
        fnt=F_SUB,
        width=112,
    )
    d.line((42, 116, 1238, 116), fill=LINE, width=2)

    card(d, (42, 140, 430, 418), "1) Intended idea")
    d.text((70, 196), "logFC = log2(group 1 expression", fill=DARK, font=F_MONO_B)
    d.text((183, 226), "/ group 2 expression)", fill=DARK, font=F_MONO_B)
    wrapped(
        d,
        "The clean idea is simple: compare average expression between one cluster and the rest. But Y_ig is already normalized and logged, so the route back to expression matters.",
        (70, 282, 315, 90),
        width=37,
    )

    card(d, (462, 140, 846, 418), "2) Two computational routes")
    pill(d, (492, 196, 804, 252), "Seurat: undo log, then average", "#FFF0DF", "#E0B17E", ORANGE)
    d.text((520, 272), "mean(exp(Y) - 1)  ->  log ratio", fill=ORANGE, font=F_MONO_B)
    pill(d, (492, 318, 804, 374), "Scanpy: average Y, then undo log", "#EFF5FA", "#ABC3DD", BLUE)
    d.text((520, 394), "exp(mean(Y)) - 1  ->  log ratio", fill=BLUE, font=F_MONO_B)

    card(d, (878, 140, 1238, 418), "3) Pseudocounts differ")
    table(d, 890, 195)
    wrapped(
        d,
        "Cluster size and software version both enter the reported logFC.",
        (900, 374, 275, 35),
        fill=DARK,
        width=34,
    )

    card(d, (42, 438, 620, 648), "4) Minimal zero-expression example")
    d.text((70, 492), "One gene: group 1 mean = 0, group 2 mean = 1", fill=DARK, font=F_BODY)
    d.text((70, 532), "Seurat-like: log2(0.01 / 1.0001) = -6.6", fill=ORANGE, font=F_MONO_B)
    d.text((70, 574), "Scanpy:       log2(1e-9 / 1)      = -29.9", fill=BLUE, font=F_MONO_B)
    d.text((70, 614), "Same edge case, very different numeric output.", fill=RED, font=F_BODY)

    card(d, (654, 438, 1238, 648), "5) Discussion takeaways")
    wrapped(
        d,
        "Outliers near +/-30 in Scanpy versus near 0 in Seurat can come from the formula, not new biology. Seurat v4 compresses logFC toward 0; Seurat v5 is less compressed; Scanpy can be more extreme when one group is zero. The authors recommend implementing the desired logFC formula directly when effect-size comparability matters.",
        (682, 492, 510, 120),
        width=64,
    )

    d.rectangle((0, 672, W, H), fill=GREEN)
    d.text((42, 685), "Discussion summary for Figure 3 / DE", fill=WHITE, font=font(16, True))
    d.text((505, 688), "Rich et al. Cell Systems 2026", fill="#DCEBE0", font=font(13))
    d.text((963, 688), "logFC is a pipeline choice, not only a statistic", fill="#DCEBE0", font=font(13))

    OUT.parent.mkdir(parents=True, exist_ok=True)
    img.save(OUT, quality=95)
    print(OUT)


if __name__ == "__main__":
    main()
