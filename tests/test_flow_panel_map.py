#!/usr/bin/env python3
"""Validate flow panel map and EV vs H filename parsing (stdlib only)."""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MAP_PATH = ROOT / "flow_panel_map.json"

FILENAME_RE_TECH = re.compile(
    r"^(?:ZZX[_-]?)?(EV|H)[-_ ]?([123])[-_ ]([12])[-_ ]+(?:PANEL[-_ ]?)?P?0?([123])[-_ ].*(unmixed|raw)\.fcs$",
    re.IGNORECASE,
)
FILENAME_RE_BIO = re.compile(
    r"^(?:ZZX[_-]?)?(EV|H)[-_ ]?([123])[-_ ]+(?:PANEL[-_ ]?)?P?0?([123])[-_ ].*(unmixed|raw)\.fcs$",
    re.IGNORECASE,
)


def parse_fcs_filename(name: str) -> dict | None:
    b = Path(name).name
    m = FILENAME_RE_TECH.match(b)
    if m:
        return {
            "group": m.group(1).upper(),
            "rep": m.group(2),
            "tech": m.group(3),
            "panel": "P" + m.group(4),
            "kind": m.group(5).lower(),
            "sample": f"{m.group(1).upper()}{m.group(2)}-{m.group(3)}",
            "bio_sample": f"{m.group(1).upper()}{m.group(2)}",
        }
    m = FILENAME_RE_BIO.match(b)
    if not m:
        return None
    return {
        "group": m.group(1).upper(),
        "rep": m.group(2),
        "tech": None,
        "panel": "P" + m.group(3),
        "kind": m.group(4).lower(),
        "sample": f"{m.group(1).upper()}{m.group(2)}",
        "bio_sample": f"{m.group(1).upper()}{m.group(2)}",
    }


def main() -> int:
    data = json.loads(MAP_PATH.read_text(encoding="utf-8"))
    errors: list[str] = []

    for panel, n_expect in (("P1", 25), ("P2", 11), ("P3", 26)):
        markers = data["panels"][panel]["markers"]
        if len(markers) != n_expect:
            errors.append(f"{panel} has {len(markers)} markers, expected {n_expect}")
        names = [m["marker"] for m in markers]
        if len(names) != len(set(names)):
            errors.append(f"{panel} has duplicate marker names")
        if names[0] != "L/D" or names[1] != "CD45":
            errors.append(f"{panel} must start with L/D then CD45")

    cases = {
        "EV-1_P1_unmixed.fcs": ("EV", "1", "P1", "unmixed"),
        "EV1_P2_unmixed.fcs": ("EV", "1", "P2", "unmixed"),
        "H3_P3_unmixed.fcs": ("H", "3", "P3", "unmixed"),
        "EV1-P3_unmixed.fcs": ("EV", "1", "P3", "unmixed"),
        "EV1_Panel3_unmixed.fcs": ("EV", "1", "P3", "unmixed"),
        "ZZX_EV1-1_P1_unmixed.fcs": ("EV", "1", "P1", "unmixed"),
        "ZZX_H2-2_P3_unmixed.fcs": ("H", "2", "P3", "unmixed"),
        "EV1-2_P1_unmixed.fcs": ("EV", "1", "P1", "unmixed"),
        "H1-1_P2_unmixed.fcs": ("H", "1", "P2", "unmixed"),
    }
    for fname, expected in cases.items():
        parsed = parse_fcs_filename(fname)
        if parsed is None:
            errors.append(f"failed to parse {fname}")
            continue
        got = (parsed["group"], parsed["rep"], parsed["panel"], parsed["kind"])
        if got != expected:
            errors.append(f"{fname} -> {got}, expected {expected}")

    if parse_fcs_filename("T-1_P1_unmixed.fcs") is not None:
        errors.append("legacy T-1 must not parse as EV/H")
    if parse_fcs_filename("T6-1_P1_unmixed.fcs") is not None:
        errors.append("legacy T6-1 must not parse as EV/H")
    if parse_fcs_filename("bad.fcs") is not None:
        errors.append("non-matching name should be None")
    if parse_fcs_filename("EV1_P1_unmixed.fcs")["group"] == "H":
        errors.append("EV1 must stay group EV")
    zzx = parse_fcs_filename("ZZX_EV1-1_P1_unmixed.fcs")
    if zzx is None or zzx["group"] != "EV" or zzx["sample"] != "EV1-1" or zzx["bio_sample"] != "EV1":
        errors.append("ZZX_EV1-1 must parse as EV tech EV1-1 / bio EV1")
    zzx_h = parse_fcs_filename("ZZX_H3-2_P2_unmixed.fcs")
    if zzx_h is None or zzx_h["group"] != "H" or zzx_h["sample"] != "H3-2":
        errors.append("ZZX_H3-2 must parse as H tech H3-2")
    if parse_fcs_filename("ZZX_EV1-1_P1_unmixed.fcs")["group"] == "H":
        errors.append("ZZX_EV must not be classified as H")

    if data.get("groups") != ["EV", "H"]:
        errors.append(f"groups should be [EV, H], got {data.get('groups')}")
    if data.get("comparison") != "H_vs_EV":
        errors.append(f"comparison should be H_vs_EV, got {data.get('comparison')}")
    if data.get("data_dir") is None or "fuction of cell" not in str(data.get("data_dir")):
        errors.append("data_dir should be E:/R/fuction of cell")

    af700 = data["fluorochrome_aliases"].get("AF700", [])
    if "AF700" not in data["fluorochrome_aliases"]:
        errors.append("fluorochrome_aliases must have AF700")
    if "R718" not in af700:
        errors.append("AF700 aliases must include R718 (Cytek detector name for the same ~700 nm channel)")
    p1_nk = next(m["fluorochrome"] for m in data["panels"]["P1"]["markers"] if m["marker"] == "NK1.1")
    if p1_nk != "AF700":
        errors.append(f"P1 NK1.1 should be AF700, got {p1_nk}")
    p3_nk = next(m["fluorochrome"] for m in data["panels"]["P3"]["markers"] if m["marker"] == "NK1.1")
    if p3_nk != "AF700":
        errors.append(f"P3 NK1.1 should be AF700, got {p3_nk}")
    p2_cd80 = next(m["fluorochrome"] for m in data["panels"]["P2"]["markers"] if m["marker"] == "CD80")
    if p2_cd80 != "APC":
        errors.append(f"P2 CD80 should be APC, got {p2_cd80}")
    p3_cd80 = next(m["fluorochrome"] for m in data["panels"]["P3"]["markers"] if m["marker"] == "CD80")
    if p3_cd80 != "APC":
        errors.append(f"P3 CD80 should be APC (synced with P2), got {p3_cd80}")
    p1_tnf = next(m["fluorochrome"] for m in data["panels"]["P1"]["markers"] if m["marker"] == "TNF-a")
    if p1_tnf != "APC":
        errors.append("P1 TNF-a stays APC; P2/P3 CD80 APC are different tubes")

    if parse_fcs_filename("flow_panel_map.json.txt") is not None:
        errors.append("json.txt is not an FCS name")

    if errors:
        print("FAIL")
        for e in errors:
            print(" -", e)
        return 1
    print("OK: panel map and filename parsing")
    return 0


if __name__ == "__main__":
    sys.exit(main())
