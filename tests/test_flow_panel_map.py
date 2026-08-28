#!/usr/bin/env python3
"""Validate flow panel map and T vs T6 filename parsing (stdlib only)."""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MAP_PATH = ROOT / "flow_panel_map.json"

FILENAME_RE = re.compile(r"^(T6|T)-([123])_(P[123])_(unmixed|raw)\.fcs$", re.IGNORECASE)


def parse_fcs_filename(name: str) -> dict | None:
    m = FILENAME_RE.match(Path(name).name)
    if not m:
        return None
    return {"group": m.group(1).upper().replace("T6", "T6"), "rep": m.group(2), "panel": m.group(3).upper(), "kind": m.group(4).lower()}


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

    # T6 must not parse as T
    cases = {
        "T-1_P1_unmixed.fcs": ("T", "1", "P1", "unmixed"),
        "T6-3_P2_unmixed.fcs": ("T6", "3", "P2", "unmixed"),
        "T-2_P3_raw.fcs": ("T", "2", "P3", "raw"),
        "T6-1_P3_unmixed.fcs": ("T6", "1", "P3", "unmixed"),
    }
    for fname, expected in cases.items():
        parsed = parse_fcs_filename(fname)
        if parsed is None:
            errors.append(f"failed to parse {fname}")
            continue
        got = (parsed["group"], parsed["rep"], parsed["panel"], parsed["kind"])
        if got != expected:
            errors.append(f"{fname} -> {got}, expected {expected}")

    assert parse_fcs_filename("T6-1_P1_unmixed.fcs")["group"] != "T"
    if parse_fcs_filename("bad.fcs") is not None:
        errors.append("non-matching name should be None")
    if parse_fcs_filename("T-1_P1_unmixed.fcs")["group"] == "T6":
        errors.append("T-1 must stay group T")

    if "AF700" not in data["fluorochrome_aliases"]["R718"]:
        errors.append("R718 aliases must include AF700 (Panel 1 NK1.1 naming)")

    # Windows 隐藏扩展名时，记事本另存为会变成 .json.txt
    if parse_fcs_filename("flow_panel_map.json.txt") is not None:
        errors.append("json.txt is not an FCS name")
    txt_alias = "flow_panel_map.json.txt"
    if not txt_alias.startswith("flow_panel_map"):
        errors.append("txt alias should still be the panel map")

    if errors:
        print("FAIL")
        for e in errors:
            print(" -", e)
        return 1
    print("OK: panel map and filename parsing")
    return 0


if __name__ == "__main__":
    sys.exit(main())
