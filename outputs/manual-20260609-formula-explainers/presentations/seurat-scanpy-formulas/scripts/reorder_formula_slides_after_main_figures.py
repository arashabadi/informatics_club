from __future__ import annotations

import argparse
import shutil
import tempfile
import zipfile
from pathlib import Path
from xml.etree import ElementTree as ET


P_NS = "http://schemas.openxmlformats.org/presentationml/2006/main"
R_NS = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
PKG_REL_NS = "http://schemas.openxmlformats.org/package/2006/relationships"

ET.register_namespace("p", P_NS)
ET.register_namespace("r", R_NS)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Move Figure 1 formula explainer slides after all Figure 1 panel slides."
    )
    parser.add_argument("pptx", type=Path)
    parser.add_argument("--backup", type=Path, required=True)
    return parser.parse_args()


def slide_targets(pptx: Path) -> tuple[list[str], ET.ElementTree, ET.Element, dict[str, str]]:
    with zipfile.ZipFile(pptx) as zf:
        presentation_xml = zf.read("ppt/presentation.xml")
        rels_xml = zf.read("ppt/_rels/presentation.xml.rels")
    tree = ET.ElementTree(ET.fromstring(presentation_xml))
    root = tree.getroot()
    sld_id_lst = root.find(f"{{{P_NS}}}sldIdLst")
    if sld_id_lst is None:
        raise RuntimeError("Could not find p:sldIdLst in presentation.xml")

    rel_root = ET.fromstring(rels_xml)
    rid_to_target = {
        rel.attrib["Id"]: rel.attrib["Target"].lstrip("/")
        for rel in rel_root.findall(f"{{{PKG_REL_NS}}}Relationship")
        if rel.attrib.get("Target", "").startswith("slides/")
    }
    ordered_targets = []
    for sld_id in list(sld_id_lst):
        rid = sld_id.attrib[f"{{{R_NS}}}id"]
        ordered_targets.append(rid_to_target[rid])
    return ordered_targets, tree, sld_id_lst, rid_to_target


def main() -> None:
    args = parse_args()
    pptx = args.pptx.resolve()
    backup = args.backup.resolve()
    backup.parent.mkdir(parents=True, exist_ok=True)
    if not backup.exists():
        shutil.copy2(pptx, backup)

    ordered_targets, tree, sld_id_lst, _ = slide_targets(pptx)
    if len(ordered_targets) < 21:
        raise RuntimeError(f"Expected at least 21 slides, found {len(ordered_targets)}")

    # Current order before this tidy-up:
    # 7 Fig1A, 8-9 Fig1A explainers, 10 Fig1B-C, 11-12 Fig1B-C explainers,
    # 13 Fig1D-E, 14 Fig1D-E explainer. Move all Figure 1 explainers after
    # the three Figure 1 panel slides so the presentation narrative is smoother.
    current_positions = list(range(1, len(ordered_targets) + 1))
    desired_positions = (
        list(range(1, 8))
        + [10, 13, 8, 9, 11, 12, 14]
        + list(range(15, len(ordered_targets) + 1))
    )
    if len(desired_positions) != len(current_positions):
        raise RuntimeError("Internal reorder map length mismatch")

    children = list(sld_id_lst)
    new_children = [children[i - 1] for i in desired_positions]
    for child in children:
        sld_id_lst.remove(child)
    for child in new_children:
        sld_id_lst.append(child)

    with tempfile.TemporaryDirectory() as tmpdir:
        tmp = Path(tmpdir) / pptx.name
        with zipfile.ZipFile(pptx, "r") as zin, zipfile.ZipFile(tmp, "w", zipfile.ZIP_DEFLATED) as zout:
            for item in zin.infolist():
                data = zin.read(item.filename)
                if item.filename == "ppt/presentation.xml":
                    data = ET.tostring(tree.getroot(), encoding="utf-8", xml_declaration=True)
                zout.writestr(item, data)
        shutil.move(tmp, pptx)

    print(
        {
            "pptx": str(pptx),
            "backup": str(backup),
            "slide_count": len(ordered_targets),
            "new_order_positions": desired_positions[:21],
        }
    )


if __name__ == "__main__":
    main()
