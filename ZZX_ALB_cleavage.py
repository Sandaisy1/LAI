#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# ZZX DIA-NN：根据上清肽段丰度推断 ALB（P02768）的细胞剪切，不是胰酶酶切。
# 输入：report.pr_matrix（TSV；也接受 report.pr.matri）
# 元数据列：Protein.Group, Protein.Ids, Protein.Names, Genes, Proteotypic,
#           Stripped.Sequence, Modified.Sequence, Precursor.Charge, Precursor.Id
# 其余列：样品强度；默认最后两列为细胞上清。
# 用法：
#   python ZZX_ALB_cleavage.py
#   python ZZX_ALB_cleavage.py --data-dir C:/Users/Lenovo/Desktop/ZZX

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

# UniProt P02768 ALBU_HUMAN SV=2（人血清白蛋白前体，609 aa）
ALB_UNIPROT = "P02768"
ALB_GENE = "ALB"
ALB_NAME = "ALBU_HUMAN"
ALB_SEQUENCE = (
    "MKWVTFISLLFLFSSAYSRGVFRRDAHKSEVAHRFKDLGEENFKALVLIAFAQYLQQCPF"
    "EDHVKLVNEVTEFAKTCVADESAENCDKSLHTLFGDKLCTVATLRETYGEMADCCAKQEP"
    "ERNECFLQHKDDNPNLPRLVRPEVDVMCTAFHDNEETFLKKYLYEIARRHPYFYAPELLF"
    "FAKRYKAAFTECCQAADKAACLLPKLDELRDEGKASSAKQRLKCASLQKFGERAFKAWAV"
    "ARLSQRFPKAEFAEVSKLVTDLTKVHTECCHGDLLECADDRADLAKYICENQDSISSKLK"
    "ECCEKPLLEKSHCIAEVENDEMPADLPSLAADFVESKDVCKNYAEAKDVFLGMFLYEYAR"
    "RHPDYSVVLLLRLAKTYETTLEKCCAAADPHECYAKVFDEFKPLVEEPQNLIKQNCELFE"
    "QLGEYKFQNALLVRYTKKVPQVSTPTLVEVSRNLGKVGSKCCKHPEAKRMPCAEDYLSVV"
    "LNQLCVLHEKTPVSDRVTKCCTESLVNRRPCFSALEVDETYVPKEFNAETFTFHADICTL"
    "SEKERQIKKQTALVELVKHKPKATKEQLKAVMDDFAAFVEKCCKADDKETCFAEEGKKLV"
    "AASQAALGL"
)
assert len(ALB_SEQUENCE) == 609

# UniProt 加工：信号肽 1–18，propeptide 19–24，成熟链 25–609
UNIPROT_CUTS = [
    {
        "after_residue": 18,
        "type": "signal",
        "source": "uniprot_annotated",
        "protease_hint": "signal peptidase (ER)",
        "note_zh": "ER 信号肽酶切掉 1–18；成熟分泌蛋白不应再含这段",
    },
    {
        "after_residue": 24,
        "type": "propeptide",
        "source": "uniprot_annotated",
        "protease_hint": "furin-like proprotein convertase (Golgi, RXXR)",
        "note_zh": "高尔基体切除 RGVFRR（19–24），成熟 N 端为 D25",
    },
]

META_COLS = {
    "Protein.Group",
    "Protein.Ids",
    "Protein.Names",
    "Genes",
    "First.Protein.Description",
    "Proteotypic",
    "Stripped.Sequence",
    "Modified.Sequence",
    "Precursor.Charge",
    "Precursor.Id",
    "Precursor.Quantity",
    "Q.Value",
    "Protein.Q.Value",
    "PG.Q.Value",
    "Global.Q.Value",
    "Lib.Q.Value",
    "Proteotypic.Ids",
}

REPORT_NAMES = (
    "report.pr_matrix",
    "report.pr_matrix.tsv",
    "report.pr_matrix.txt",
    "report.pr.matrix",
    "report.pr.matri",
    "report.pr.matri.tsv",
)

DATA_DIR_CANDIDATES = (
    Path(r"C:/Users/Lenovo/Desktop/ZZX"),
    Path(r"C:/User/Lenovo/Desktop/ZZX"),
    Path("ZZX"),
    Path("."),
)


def log(lines: list[str], msg: str) -> None:
    print(msg)
    lines.append(msg)


def parse_num(value) -> float | None:
    if value is None:
        return None
    text = str(value).strip()
    if text in {"", "NA", "NaN", "nan", "None", "NULL", "#N/A"}:
        return None
    text = text.replace(",", "")
    try:
        return float(text)
    except ValueError:
        return None


def tokens(value: str) -> list[str]:
    text = str(value or "").strip()
    for sep in (";", "|", ",", "/"):
        text = text.replace(sep, " ")
    return [part for part in text.split() if part]


def is_alb_row(row: dict) -> bool:
    genes = {g.upper() for g in tokens(row.get("Genes", ""))}
    ids = {i.split("-")[0].upper() for i in tokens(row.get("Protein.Ids", ""))}
    names = str(row.get("Protein.Names", "")).upper()
    group = str(row.get("Protein.Group", "")).upper()
    if ALB_GENE in genes:
        return True
    if ALB_UNIPROT in ids:
        return True
    if ALB_UNIPROT in group:
        return True
    return ALB_NAME in names


def find_report(data_dir: Path | None) -> Path:
    dirs = []
    if data_dir:
        dirs.append(Path(data_dir))
    dirs.extend(DATA_DIR_CANDIDATES)
    seen = []
    for folder in dirs:
        folder = Path(folder)
        if folder in seen:
            continue
        seen.append(folder)
        if not folder.is_dir():
            continue
        for name in REPORT_NAMES:
            path = folder / name
            if path.is_file():
                return path
        for path in sorted(folder.glob("report.pr*")):
            if path.is_file() and "pg_matrix" not in path.name.lower():
                return path
    raise FileNotFoundError(
        "未找到 report.pr_matrix。请把 DIA-NN 肽段表放到 "
        "C:/Users/Lenovo/Desktop/ZZX 或用 --data-dir 指定。"
    )


def read_tsv(path: Path) -> tuple[list[str], list[dict]]:
    raw = path.read_bytes()
    for encoding in ("utf-8-sig", "utf-8", "gb18030", "latin-1"):
        try:
            text = raw.decode(encoding)
            break
        except UnicodeDecodeError:
            continue
    else:
        text = raw.decode("utf-8", errors="replace")
    first = text.splitlines()[0] if text else ""
    try:
        dialect = csv.Sniffer().sniff(first, delimiters="\t,;")
    except csv.Error:
        dialect = csv.excel_tab
    if first.count("\t") >= first.count(",") and first.count("\t") >= 3:
        dialect = csv.excel_tab
    reader = csv.DictReader(text.splitlines(), dialect=dialect)
    rows = [dict(r) for r in reader]
    if not reader.fieldnames:
        raise ValueError(f"{path} 没有表头")
    return list(reader.fieldnames), rows


_META_COL_SUFFIXES = (
    "Sequence",
    "Ids",
    "Names",
    "Genes",
    "Description",
    "Group",
)
_MS_FILE_EXTS = (".raw", ".mzml", ".dia", ".parquet", ".wiff2", ".wiff")


def sample_columns(fieldnames: list[str], n_last: int = 2) -> list[str]:
    extras = [c for c in fieldnames if c not in META_COLS]
    numeric_like = []
    for col in extras:
        if any(col.endswith(suffix) for suffix in _META_COL_SUFFIXES):
            continue
        numeric_like.append(col)
    if len(numeric_like) >= n_last:
        return numeric_like[-n_last:]
    if numeric_like:
        return numeric_like
    return extras[-n_last:] if extras else []


def short_sample_name(col: str, index: int) -> str:
    base = Path(str(col).replace(chr(92), "/")).name
    lower = base.lower()
    for ext in _MS_FILE_EXTS:
        if lower.endswith(ext):
            base = base[: -len(ext)]
            break
    return base.strip() or f"supernatant_{index}"


def map_peptide(sequence: str, peptide: str) -> list[tuple[int, int]]:
    """Return 1-based inclusive (start, end) matches."""
    hits = []
    start = 0
    pep = peptide.strip().upper()
    seq = sequence.upper()
    if not pep:
        return hits
    while True:
        i = seq.find(pep, start)
        if i < 0:
            break
        hits.append((i + 1, i + len(pep)))
        start = i + 1
    if hits:
        return hits
    # I/L 质谱无法区分时的回退
    seq_il = seq.replace("I", "L")
    pep_il = pep.replace("I", "L")
    start = 0
    while True:
        i = seq_il.find(pep_il, start)
        if i < 0:
            break
        hits.append((i + 1, i + len(pep)))
        start = i + 1
    return hits


def tryptic_n(seq: str, start: int, block_proline: bool = True) -> bool:
    if start <= 1:
        return True
    prev = seq[start - 2]
    this = seq[start - 1]
    if prev in "KR":
        return not (block_proline and this == "P")
    return False


def tryptic_c(seq: str, end: int, block_proline: bool = True) -> bool:
    if end >= len(seq):
        return True
    aa = seq[end - 1]
    nxt = seq[end] if end < len(seq) else ""
    if aa in "KR":
        return not (block_proline and nxt == "P")
    return False


def terminus_class(n_ok: bool, c_ok: bool) -> str:
    if n_ok and c_ok:
        return "fully_tryptic"
    if n_ok and not c_ok:
        return "semi_tryptic_neoC"
    if (not n_ok) and c_ok:
        return "semi_tryptic_neoN"
    return "non_tryptic"


def motif(seq: str, after: int) -> str:
    chars = []
    for k in range(after - 3, after + 5):
        chars.append(seq[k - 1] if 1 <= k <= len(seq) else "x")
    return "".join(chars[:4]) + "|" + "".join(chars[4:])


def protease_hint(seq: str, after: int) -> str:
    if after < 1 or after >= len(seq):
        return "unknown"
    p4, p3, p2, p1 = (
        seq[after - 4] if after >= 4 else "x",
        seq[after - 3] if after >= 3 else "x",
        seq[after - 2] if after >= 2 else "x",
        seq[after - 1],
    )
    p1p = seq[after]
    if after in (18, 24):
        return UNIPROT_CUTS[0]["protease_hint"] if after == 18 else UNIPROT_CUTS[1]["protease_hint"]
    if p1 == "R" and p4 == "R":
        return "furin-like (RXXR)"
    if p1 == "R" and p4 in "KR":
        return "proprotein convertase-like"
    if p1 == "D" and p1p in "GAS":
        return "caspase-like (D|x)"
    if p1 in "ASGC" and p1p not in "P":
        return "signal-peptidase-like small P1"
    if p1 in "KR":
        return "trypsin-like / K-or-R P1 (could be cellular or look tryptic after digest)"
    if p1 in "FLYM" and p1p not in "P":
        return "hydrophobic P1 (cathepsin/signal-like)"
    return "unassigned cellular protease"


def mean_or_none(values: list[float | None]) -> float | None:
    nums = [v for v in values if v is not None and not (isinstance(v, float) and math.isnan(v))]
    if not nums:
        return None
    return sum(nums) / len(nums)


def log10p(v: float | None) -> float | None:
    if v is None or v <= 0:
        return None
    return math.log10(v)


def region_mean(profile: list[float | None], start: int, end: int) -> float | None:
    vals = [profile[i - 1] for i in range(start, end + 1)]
    return mean_or_none(vals)


def write_csv(path: Path, rows: list[dict], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            out = {}
            for key in fieldnames:
                val = row.get(key, "")
                if val is None:
                    val = ""
                elif isinstance(val, float):
                    val = "" if math.isnan(val) else f"{val:.6g}"
                out[key] = val
            writer.writerow(out)


def svg_escape(text: str) -> str:
    return (
        str(text)
        .replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace(chr(34), "&quot;")
    )


def write_schematic_svg(path: Path, seq: str, candidates: list[dict], n_term_cov: dict) -> None:
    width, height = 980, 220
    n = len(seq)
    x0, y0, bar_w, bar_h = 40, 80, 900, 36

    def x_at(res: int) -> float:
        return x0 + bar_w * (res - 1) / n

    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}">',
        "<style>text{font-family:Arial,Helvetica,sans-serif;font-size:12px}</style>",
        f'<text x="40" y="28" font-size="16">ALB (P02768) cellular processing schematic</text>',
        f'<text x="40" y="48" fill="#555">signal 1–18 | propeptide 19–24 | mature 25–{n}; arrows = inferred cellular cuts</text>',
        f'<rect x="{x_at(1):.1f}" y="{y0}" width="{x_at(19)-x_at(1):.1f}" height="{bar_h}" fill="#f4c7c3" stroke="#333"/>',
        f'<rect x="{x_at(19):.1f}" y="{y0}" width="{x_at(25)-x_at(19):.1f}" height="{bar_h}" fill="#ffe199" stroke="#333"/>',
        f'<rect x="{x_at(25):.1f}" y="{y0}" width="{x_at(n)+bar_w/n-x_at(25):.1f}" height="{bar_h}" fill="#b7d7b0" stroke="#333"/>',
        f'<text x="{x_at(1)+4:.1f}" y="{y0 + 22}" font-size="11">SIGNAL</text>',
        f'<text x="{x_at(25)+8:.1f}" y="{y0 + 22}" font-size="11">MATURE ALBUMIN (secreted)</text>',
    ]
    shown = set()
    for cand in candidates:
        after = int(cand["after_residue"])
        if after in shown or after < 1 or after >= n:
            continue
        shown.add(after)
        x = x_at(after + 1)
        color = {"signal": "#c0392b", "propeptide": "#d35400"}.get(cand["type"], "#6c3483")
        parts.append(f'<line x1="{x:.1f}" y1="{y0-8}" x2="{x:.1f}" y2="{y0+bar_h+18}" stroke="{color}" stroke-width="2"/>')
        parts.append(
            f'<text x="{x:.1f}" y="{y0+bar_h+32}" font-size="10" text-anchor="middle" fill="{color}">'
            f'{after}|{after+1}</text>'
        )
    cov = (
        f"coverage 1-18={n_term_cov.get('signal')}  "
        f"19-24={n_term_cov.get('propeptide')}  "
        f"25-80={n_term_cov.get('mature_n')}  "
        f"81-{n}={n_term_cov.get('mature_rest')}"
    )
    parts.append(f'<text x="40" y="{height-16}" fill="#333">{svg_escape(cov)}</text>')
    parts.append("</svg>")
    path.write_text("\n".join(parts), encoding="utf-8")


def _fmt_cov(value) -> str:
    if value is None:
        return "supernatant: not detected"
    return f"supernatant: detected ({value:.3g})"


def write_full_length_cleavage_svg(
    path: Path,
    seq: str,
    candidates: list[dict],
    n_term_cov: dict,
) -> None:
    """Stepwise cartoon: full-length prepro-ALB -> signal cut -> propeptide cut -> mature (+ optional internal)."""
    n = len(seq)
    internals = [
        c
        for c in candidates
        if c.get("type") == "internal" and c.get("source") in {"neo_N", "neo_C"}
    ]
    internals.sort(key=lambda c: int(c["after_residue"]))
    # unique cut sites, strongest first then residue order
    uniq = []
    seen = set()
    for c in internals:
        after = int(c["after_residue"])
        if after in seen or after <= 24 or after >= n:
            continue
        seen.add(after)
        uniq.append(c)
        if len(uniq) >= 4:
            break

    n_steps = 3 + (1 if uniq else 0)
    width = 1120
    left, bar_w, bar_h = 130, 920, 42
    zoom_w = 300  # residues 1–24 expanded so SIGNAL/PRO are visible
    rest_w = bar_w - zoom_w
    step_gap = 152
    top = 96
    height = 86 + n_steps * step_gap + 56
    if uniq:
        height += 36

    def x_left(res: int) -> float:
        if res <= 24:
            return left + zoom_w * (res - 1) / 24.0
        return left + zoom_w + rest_w * (res - 25) / float(n - 24)

    def x_right(res: int) -> float:
        if res <= 24:
            return left + zoom_w * res / 24.0
        return left + zoom_w + rest_w * (res - 24) / float(n - 24)

    def x_cut(after: int) -> float:
        return x_right(after)

    domains = (
        (1, 18, "#e74c3c", "#fadbd8", "SIGNAL", "1–18"),
        (19, 24, "#e67e22", "#fdebd0", "PRO", "19–24"),
        (25, n, "#1e8449", "#d5f5e3", "MATURE", f"25–{n}"),
    )

    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        "<rect width='100%' height='100%' fill='#ffffff'/>",
        '<style>text{font-family:"Microsoft YaHei","SimHei","Noto Sans SC","PingFang SC",Arial,Helvetica,sans-serif}</style>',
        '<text x="40" y="28" font-size="18" font-weight="700">Full-length ALB cleavage  全长蛋白如何剪切</text>',
        '<text x="40" y="50" font-size="12" fill="#555">P02768 prepro-albumin (1–609). Cellular processing only (not trypsin). '
        "N-terminus (1–24) is expanded so SIGNAL / PROPEPTIDE cuts are visible.</text>",
    ]

    def draw_domains(y: float, keep_from: int, keep_to: int, ghost: bool = False) -> None:
        for a, b, stroke, fill, name, rng in domains:
            lo, hi = max(a, keep_from), min(b, keep_to)
            if lo > hi:
                continue
            x1, x2 = x_left(lo), x_right(hi)
            w = max(x2 - x1, 2)
            opacity = 0.22 if ghost else 1.0
            dash = ' stroke-dasharray="5 3"' if ghost else ""
            parts.append(
                f'<rect x="{x1:.1f}" y="{y}" width="{w:.1f}" height="{bar_h}" rx="6" '
                f'fill="{fill}" stroke="{stroke}" stroke-width="1.6"{dash} opacity="{opacity}"/>'
            )
            label_x = x1 + w / 2
            if w >= 36:
                parts.append(
                    f'<text x="{label_x:.1f}" y="{y + 17}" text-anchor="middle" font-size="11" font-weight="700" fill="{stroke}">{name}</text>'
                )
                parts.append(
                    f'<text x="{label_x:.1f}" y="{y + 32}" text-anchor="middle" font-size="10" fill="#333">{rng if (lo, hi) == (a, b) else f"{lo}–{hi}"}</text>'
                )

    def draw_scissors(after: int, y: float, color: str, label: str, motif_txt: str) -> None:
        x = x_cut(after)
        y0 = y - 16
        y1 = y + bar_h + 8
        parts.append(
            f'<line x1="{x:.1f}" y1="{y0}" x2="{x:.1f}" y2="{y1}" stroke="{color}" stroke-width="2.2"/>'
        )
        # scissors head
        parts.append(f'<polygon points="{x-7:.1f},{y0-2} {x+7:.1f},{y0-2} {x:.1f},{y0+12}" fill="{color}"/>')
        parts.append(
            f'<text x="{x:.1f}" y="{y + bar_h + 22}" text-anchor="middle" font-size="11" font-weight="700" fill="{color}">{svg_escape(label)}</text>'
        )
        parts.append(
            f'<text x="{x:.1f}" y="{y + bar_h + 36}" text-anchor="middle" font-size="10" fill="#666">{svg_escape(motif_txt)}</text>'
        )

    def step_badge(num: int, y: float, title: str, subtitle: str) -> None:
        cy = y + bar_h / 2
        parts.append(f'<circle cx="48" cy="{cy:.1f}" r="16" fill="#2c3e50"/>')
        parts.append(
            f'<text x="48" y="{cy + 5:.1f}" text-anchor="middle" font-size="14" font-weight="700" fill="#fff">{num}</text>'
        )
        parts.append(f'<text x="130" y="{y - 22}" font-size="13" font-weight="700">{svg_escape(title)}</text>')
        parts.append(f'<text x="130" y="{y - 6}" font-size="11" fill="#555">{svg_escape(subtitle)}</text>')

    def down_arrow(y_from: float, caption: str) -> None:
        x = left + bar_w / 2
        y1 = y_from + bar_h + 40
        y2 = y1 + 22
        parts.append(f'<line x1="{x:.1f}" y1="{y1}" x2="{x:.1f}" y2="{y2}" stroke="#2c3e50" stroke-width="2"/>')
        parts.append(f'<polygon points="{x-6:.1f},{y2} {x+6:.1f},{y2} {x:.1f},{y2+10}" fill="#2c3e50"/>')
        parts.append(
            f'<text x="{x + 12:.1f}" y="{y1 + 16}" font-size="11" fill="#2c3e50">{svg_escape(caption)}</text>'
        )

    # residue ticks on step 1
    def ticks(y: float) -> None:
        for res in (1, 18, 19, 24, 25, n):
            x = x_left(res) if res != n else x_right(n)
            parts.append(f'<line x1="{x:.1f}" y1="{y - 4}" x2="{x:.1f}" y2="{y}" stroke="#888" stroke-width="1"/>')
            parts.append(
                f'<text x="{x:.1f}" y="{y - 8}" text-anchor="middle" font-size="9" fill="#666">{res}</text>'
            )

    y1 = top
    step_badge(1, y1, "Full-length precursor  全长前体 prepro-ALB", f"residues 1–{n} (signal + propeptide + mature chain)")
    draw_domains(y1, 1, n)
    ticks(y1)
    draw_scissors(18, y1, "#c0392b", "cut 18|19", "SAYS|RGVF  signal peptidase (ER)")
    down_arrow(y1, "ER: remove SIGNAL 1–18")

    y2 = y1 + step_gap
    step_badge(
        2,
        y2,
        "After signal-peptide cleavage  切掉信号肽",
        f"released 1–18 ({_fmt_cov(n_term_cov.get('signal'))}); remaining proalbumin 19–{n}",
    )
    draw_domains(y2, 1, 18, ghost=True)
    draw_domains(y2, 19, n)
    draw_scissors(24, y2, "#d35400", "cut 24|25", "VFRR|DAHK  furin-like (Golgi, RXXR)")
    down_arrow(y2, "Golgi: remove PROPEPTIDE 19–24")

    y3 = y2 + step_gap
    step_badge(
        3,
        y3,
        "Mature albumin secreted  成熟链进入细胞上清",
        f"chain 25–{n}  ({_fmt_cov(n_term_cov.get('mature_n'))} at N-term 25–80; "
        f"{_fmt_cov(n_term_cov.get('mature_rest'))} for 81–{n})",
    )
    draw_domains(y3, 1, 24, ghost=True)
    draw_domains(y3, 25, n)
    # outline secreted product
    parts.append(
        f'<rect x="{x_left(25):.1f}" y="{y3 - 4}" width="{x_right(n) - x_left(25):.1f}" height="{bar_h + 8}" '
        f'rx="8" fill="none" stroke="#1e8449" stroke-width="2.4"/>'
    )

    if uniq:
        y4 = y3 + step_gap
        cuts = ", ".join(f"{int(c['after_residue'])}|{int(c['before_residue'])}" for c in uniq)
        step_badge(
            4,
            y4,
            "Additional cellular cuts on the mature chain  成熟链上的内部剪切",
            f"neo-N / neo-C peptide evidence (not trypsin K/R): {cuts}",
        )
        draw_domains(y4, 25, n)
        for c in uniq:
            after = int(c["after_residue"])
            motif_txt = str(c.get("motif_P4_P4prime") or "")
            hint = str(c.get("protease_hint") or "cellular protease")
            draw_scissors(after, y4, "#6c3483", f"cut {after}|{after + 1}", f"{motif_txt}  {hint}")
        # fragment labels under mature bar
        bounds = [25] + [int(c["after_residue"]) + 1 for c in uniq] + [n + 1]
        frag_y = y4 + bar_h + 48
        for i in range(len(bounds) - 1):
            a, b = bounds[i], bounds[i + 1] - 1
            if a > b:
                continue
            mx = (x_left(a) + x_right(b)) / 2
            parts.append(
                f'<text x="{mx:.1f}" y="{frag_y}" text-anchor="middle" font-size="11" fill="#6c3483">'
                f"fragment {a}–{b}</text>"
            )

    # legend
    ly = height - 28
    legend = [
        (left, "#e74c3c", "SIGNAL 1–18"),
        (left + 160, "#e67e22", "PROPEPTIDE 19–24"),
        (left + 360, "#1e8449", "MATURE 25–609 (secreted)"),
        (left + 620, "#888888", "dashed = released / not kept"),
    ]
    for x, color, lab in legend:
        parts.append(f'<rect x="{x}" y="{ly - 12}" width="16" height="12" fill="{color}" opacity="0.35" stroke="{color}"/>')
        parts.append(f'<text x="{x + 22}" y="{ly}" font-size="11" fill="#333">{lab}</text>')
    parts.append("</svg>")
    path.write_text("\n".join(parts), encoding="utf-8")


def write_coverage_svg(
    path: Path,
    seq: str,
    profiles: dict[str, list[float | None]],
    peptides: list[dict],
) -> None:
    n = len(seq)
    width = 1000
    track_h = 140
    pep_h = max(80, 8 + 7 * min(len(peptides), 40))
    height = 70 + track_h + pep_h
    x0, bar_w = 50, 920

    def x_at(res: int) -> float:
        return x0 + bar_w * (res - 0.5) / n

    ymax = 1.0
    for series in profiles.values():
        for v in series:
            lv = log10p(v)
            if lv is not None:
                ymax = max(ymax, lv)
    ymax *= 1.08
    colors = ["#1f77b4", "#d62728", "#2ca02c"]
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}">',
        "<style>text{font-family:Arial,Helvetica,sans-serif;font-size:11px}</style>",
        '<text x="50" y="24" font-size="16">ALB residue abundance (cell supernatant) and peptide map</text>',
        '<text x="50" y="42" fill="#555">Top: log10 intensity along prepro-ALB. Bottom: peptides '
        "(grey=trypsin digest; orange=cellular neo-N; blue=cellular neo-C; red=both).</text>",
        f'<rect x="{x0}" y="55" width="{bar_w}" height="{track_h}" fill="#fafafa" stroke="#ccc"/>',
    ]
    # domain shading
    for a, b, fill in ((1, 18, "#f4c7c3"), (19, 24, "#ffe199"), (25, n, "#eaf6e8")):
        xa, xb = x_at(a) - bar_w / n / 2, x_at(b) + bar_w / n / 2
        parts.append(
            f'<rect x="{xa:.1f}" y="55" width="{xb-xa:.1f}" height="{track_h}" fill="{fill}" opacity="0.45"/>'
        )
    for i, (name, series) in enumerate(profiles.items()):
        pts = []
        for res, val in enumerate(series, start=1):
            lv = log10p(val)
            if lv is None:
                continue
            y = 55 + track_h - (lv / ymax) * (track_h - 8)
            pts.append(f"{x_at(res):.1f},{y:.1f}")
        if pts:
            pt_str = " ".join(pts)
            color = colors[i % len(colors)]
            parts.append(
                f'<polyline fill="none" stroke="{color}" stroke-width="1.5" points="{pt_str}"/>'
            )
            parts.append(
                f'<text x="{x0 + 8}" y="{68 + i * 14}" fill="{color}">{svg_escape(name)}</text>'
            )
    y_pep0 = 55 + track_h + 16
    class_color = {
        "fully_tryptic": "#bbbbbb",
        "semi_tryptic_neoN": "#e67e22",
        "semi_tryptic_neoC": "#2980b9",
        "non_tryptic": "#c0392b",
    }
    shown = peptides[:80]
    for i, pep in enumerate(shown):
        y = y_pep0 + (i % 40) * 6
        x1 = x_at(int(pep["start"]))
        x2 = x_at(int(pep["end"]))
        color = class_color.get(pep["terminus_class"], "#888")
        parts.append(
            f'<rect x="{min(x1,x2):.1f}" y="{y:.1f}" width="{max(2, abs(x2-x1)):.1f}" height="4" '
            f'fill="{color}" opacity="0.85"/>'
        )
    parts.append(f'<text x="{x0}" y="{height-12}" fill="#555">residue 1 ... {n} (prepro-albumin)</text>')
    parts.append("</svg>")
    path.write_text("\n".join(parts), encoding="utf-8")


def infer_uniprot_evidence(profile: list[float | None], neo_n: dict, neo_c: dict, n: int) -> list[dict]:
    rows = []
    sig = region_mean(profile, 1, 18)
    pro = region_mean(profile, 19, 24)
    mat = region_mean(profile, 25, min(80, n))
    rest = region_mean(profile, 81, n) if n > 80 else mat

    def cov_label(v):
        return "detected" if v is not None and v > 0 else "not_detected"

    # signal 18|19
    n19 = neo_n.get(19, {})
    c18 = neo_c.get(18, {})
    signal_ok = cov_label(sig) == "not_detected" and cov_label(mat) == "detected"
    rows.append(
        {
            "after_residue": 18,
            "before_residue": 19,
            "type": "signal",
            "source": "uniprot_annotated",
            "motif_P4_P4prime": motif(ALB_SEQUENCE, 18),
            "protease_hint": "signal peptidase (ER)",
            "n_neo_N_peptides": n19.get("n", 0),
            "abundance_neo_N": n19.get("abundance"),
            "n_neo_C_peptides": c18.get("n", 0),
            "abundance_neo_C": c18.get("abundance"),
            "abundance_left_window": sig,
            "abundance_right_window": pro if pro is not None else mat,
            "supported": "yes" if signal_ok or n19.get("n") or c18.get("n") else "possible",
            "evidence_zh": (
                f"信号肽区 1–18 {cov_label(sig)}，成熟 N 端 25–80 {cov_label(mat)}。"
                f"{'上清几乎只有成熟链覆盖，符合 ER 切掉信号肽后分泌。' if signal_ok else '信号肽区仍有肽段，需检查是否未加工前体或错误匹配。'}"
                " 仅空洞不能单独定论；neo-N@19 才是切割的直接肽段证据。"
            ),
            "note_zh": UNIPROT_CUTS[0]["note_zh"],
        }
    )
    n25 = neo_n.get(25, {})
    c24 = neo_c.get(24, {})
    # 成熟 N 端 D25 前一位是 R24，胰酶消化后 DAHK 看起来像全胰酶肽，不能当成 neo-N
    pro_ok = cov_label(pro) == "not_detected" and cov_label(mat) == "detected"
    rows.append(
        {
            "after_residue": 24,
            "before_residue": 25,
            "type": "propeptide",
            "source": "uniprot_annotated",
            "motif_P4_P4prime": motif(ALB_SEQUENCE, 24),
            "protease_hint": "furin-like proprotein convertase (Golgi, RXXR)",
            "n_neo_N_peptides": n25.get("n", 0),
            "abundance_neo_N": n25.get("abundance"),
            "n_neo_C_peptides": c24.get("n", 0),
            "abundance_neo_C": c24.get("abundance"),
            "abundance_left_window": pro,
            "abundance_right_window": mat,
            "supported": "yes" if pro_ok or n25.get("n") or c24.get("n") else "possible",
            "evidence_zh": (
                f"propeptide 19–24 {cov_label(pro)}，成熟 N 端 25–80 {cov_label(mat)}。"
                "R24|D25 在胰酶消化后是合法胰酶位点，因此成熟 N 端肽段通常不是 neo-N；"
                "要用 19–24 覆盖缺失 + 25 之后高覆盖 来支持高尔基体切除 propeptide。"
            ),
            "note_zh": UNIPROT_CUTS[1]["note_zh"],
        }
    )
    return rows


def collect_neo(peptides: list[dict], which: str) -> dict[int, dict]:
    bucket: dict[int, dict] = {}
    for pep in peptides:
        if which == "N" and not pep["cellular_neo_N"]:
            continue
        if which == "C" and not pep["cellular_neo_C"]:
            continue
        pos = int(pep["start"] if which == "N" else pep["end"])
        rec = bucket.setdefault(pos, {"n": 0, "abundance": 0.0, "peptides": []})
        rec["n"] += 1
        rec["abundance"] += pep.get("abundance_mean") or 0.0
        rec["peptides"].append(pep["Stripped.Sequence"])
    return bucket


def merge_precursors(peptides: list[dict], sample_names: list[str]) -> list[dict]:
    """Sum charge states of the same stripped peptide and map coordinates."""
    merged: dict[tuple, dict] = {}
    for rec in peptides:
        key = (rec["Stripped.Sequence"], rec["start"], rec["end"])
        if key not in merged:
            item = dict(rec)
            item["_vals"] = {name: [] for name in sample_names}
            item["n_precursors"] = 0
            merged[key] = item
        item = merged[key]
        item["n_precursors"] += 1
        for name in sample_names:
            if rec.get(name) is not None:
                item["_vals"][name].append(rec[name])
    out = []
    for item in merged.values():
        for name in sample_names:
            vals = item["_vals"][name]
            item[name] = sum(vals) if vals else None
        item["abundance_mean"] = mean_or_none([item[name] for name in sample_names])
        item["Precursor.Charge"] = "merged"
        item["terminus_class"] = terminus_class(not item["cellular_neo_N"], not item["cellular_neo_C"])
        del item["_vals"]
        out.append(item)
    out.sort(key=lambda r: (r["start"], r["end"], r["Stripped.Sequence"]))
    return out


def abundance_steps(profile: list[float | None], window: int = 15) -> list[dict]:
    n = len(profile)
    out = []
    for after in range(window, n - window):
        left = mean_or_none(profile[after - window : after])
        right = mean_or_none(profile[after : after + window])
        if left is None or right is None:
            continue
        if left <= 0 and right <= 0:
            continue
        lstep = math.log2((right + 1.0) / (left + 1.0))
        if abs(lstep) < 1.5:
            continue
        out.append(
            {
                "after_residue": after,
                "before_residue": after + 1,
                "type": "internal",
                "source": "abundance_step",
                "motif_P4_P4prime": motif(ALB_SEQUENCE, after),
                "protease_hint": protease_hint(ALB_SEQUENCE, after),
                "n_neo_N_peptides": 0,
                "abundance_neo_N": None,
                "n_neo_C_peptides": 0,
                "abundance_neo_C": None,
                "abundance_left_window": left,
                "abundance_right_window": right,
                "log2_step": lstep,
                "supported": "weak",
                "evidence_zh": (
                    f"切割位点两侧 {window} aa 丰度 log2 变化 {lstep:.2f}。"
                    "空洞也可能来自肽段过短/疏水/漏切，不能单独当作细胞剪切。"
                ),
                "note_zh": "丰度阶跃只作提示，优先看 neo-N/neo-C",
            }
        )
    out.sort(key=lambda r: abs(r.get("log2_step") or 0), reverse=True)
    return out[:15]


def run(data_dir: Path | None, outdir: Path | None, n_last: int) -> Path:
    logs: list[str] = []
    report = find_report(data_dir)
    log(logs, f"[data] 读取 {report}")
    fieldnames, rows = read_tsv(report)
    qty_cols = sample_columns(fieldnames, n_last=n_last)
    if not qty_cols:
        raise ValueError("没有找到样品丰度列（表头 Precursor.Id 之后的上清强度列）")
    sample_names = [short_sample_name(c, i + 1) for i, c in enumerate(qty_cols)]
    log(logs, "[data] 上清丰度列: " + ", ".join(f"{n} <= {c}" for n, c in zip(sample_names, qty_cols)))
    log(logs, "[note] 只把非胰酶末端当细胞剪切候选；K/R 胰酶切口不是细胞剪切位点。")
    log(logs, "[note] 无生物学重复时不伪造 p 值，按覆盖与丰度描述。")

    alb_rows = [r for r in rows if is_alb_row(r)]
    if not alb_rows:
        genes = sorted({r.get("Genes", "") for r in rows if r.get("Genes")})
        raise SystemExit(
            "pr_matrix 里没有 ALB / P02768。表中基因包括: " + ", ".join(genes[:40])
        )
    log(logs, f"[filter] ALB 前体行 {len(alb_rows)} / 总行 {len(rows)}（含多电荷）")

    seq = ALB_SEQUENCE
    peptides = []
    unmapped = 0
    for row in alb_rows:
        pep = str(row.get("Stripped.Sequence") or "").strip().upper()
        if not pep:
            continue
        hits = map_peptide(seq, pep)
        if not hits:
            unmapped += 1
            continue
        abundances = [parse_num(row.get(c)) for c in qty_cols]
        if all(v is None for v in abundances):
            continue
        proteotypic = str(row.get("Proteotypic", "")).strip()
        for start, end in hits:
            n_ok = tryptic_n(seq, start)
            c_ok = tryptic_c(seq, end)
            rec = {
                "Protein.Group": row.get("Protein.Group", ""),
                "Protein.Ids": row.get("Protein.Ids", ""),
                "Genes": row.get("Genes", ""),
                "Proteotypic": proteotypic,
                "Stripped.Sequence": pep,
                "Modified.Sequence": row.get("Modified.Sequence", ""),
                "Precursor.Id": row.get("Precursor.Id", ""),
                "Precursor.Charge": row.get("Precursor.Charge", ""),
                "start": start,
                "end": end,
                "length": end - start + 1,
                "n_matches_on_ALB": len(hits),
                "prev_aa": seq[start - 2] if start > 1 else "",
                "next_aa": seq[end] if end < len(seq) else "",
                "tryptic_N": n_ok,
                "tryptic_C": c_ok,
                "cellular_neo_N": (not n_ok),
                "cellular_neo_C": (not c_ok),
                "terminus_class": terminus_class(n_ok, c_ok),
                "cleavage_after_if_neoN": start - 1 if not n_ok else "",
                "cleavage_after_if_neoC": end if not c_ok else "",
                "abundance_mean": mean_or_none(abundances),
            }
            for name, val in zip(sample_names, abundances):
                rec[name] = val
            peptides.append(rec)
    log(logs, f"[map] 定位到 P02768 的前体 {len(peptides)} 条；未映射 {unmapped} 条")
    peptide_level = merge_precursors(peptides, sample_names)

    n_neo = sum(1 for p in peptide_level if p["cellular_neo_N"] or p["cellular_neo_C"])
    if n_neo == 0:
        log(
            logs,
            "[warn] 没有任何 neo-N/neo-C 肽段。DIA 库很可能只搜了全胰酶肽，"
            "细胞剪切的直接末端证据会缺失；将主要依据 UniProt 加工位点 + 覆盖/丰度阶跃。",
        )
    else:
        log(logs, f"[cleavage] 细胞剪切相关肽段（neo-N 或 neo-C）{n_neo} 条")

    n = len(seq)
    profiles = {name: [None] * n for name in sample_names}
    profiles["mean"] = [None] * n
    counts = [0] * n
    for pep in peptide_level:
        for pos in range(pep["start"], pep["end"] + 1):
            counts[pos - 1] += 1
            for name in sample_names:
                val = pep.get(name)
                if val is None:
                    continue
                cur = profiles[name][pos - 1]
                profiles[name][pos - 1] = val if cur is None else cur + val
            val = pep.get("abundance_mean")
            if val is not None:
                cur = profiles["mean"][pos - 1]
                profiles["mean"][pos - 1] = val if cur is None else cur + val

    n_term_cov = {
        "signal": region_mean(profiles["mean"], 1, 18),
        "propeptide": region_mean(profiles["mean"], 19, 24),
        "mature_n": region_mean(profiles["mean"], 25, min(80, n)),
        "mature_rest": region_mean(profiles["mean"], 81, n) if n > 80 else None,
    }
    log(
        logs,
        "[coverage] 区段平均丰度 "
        f"信号肽1-18={n_term_cov['signal']}  pro19-24={n_term_cov['propeptide']}  "
        f"成熟N25-80={n_term_cov['mature_n']}  其余81-609={n_term_cov['mature_rest']}",
    )

    neo_n = collect_neo(peptide_level, "N")
    neo_c = collect_neo(peptide_level, "C")
    candidates = infer_uniprot_evidence(profiles["mean"], neo_n, neo_c, n)

    for pos, rec in sorted(neo_n.items(), key=lambda kv: -kv[1]["abundance"]):
        after = pos - 1
        if after in (18, 24):
            continue
        if after < 1:
            continue
        c_rec = neo_c.get(after, {})
        candidates.append(
            {
                "after_residue": after,
                "before_residue": pos,
                "type": "internal",
                "source": "neo_N",
                "motif_P4_P4prime": motif(seq, after),
                "protease_hint": protease_hint(seq, after),
                "n_neo_N_peptides": rec["n"],
                "abundance_neo_N": rec["abundance"],
                "n_neo_C_peptides": c_rec.get("n", 0),
                "abundance_neo_C": c_rec.get("abundance"),
                "abundance_left_window": region_mean(profiles["mean"], max(1, after - 14), after),
                "abundance_right_window": region_mean(profiles["mean"], pos, min(n, pos + 14)),
                "supported": "yes" if rec["n"] >= 1 else "possible",
                "evidence_zh": (
                    f"细胞 neo-N 肽段从残基 {pos} 开始（切割在 {after}|{pos}），"
                    f"{rec['n']} 条，Σ丰度 {rec['abundance']:.4g}。"
                    f"肽段: {', '.join(rec['peptides'][:8])}"
                ),
                "note_zh": "非胰酶 N 端，优先视为细胞内/分泌路径蛋白酶切割",
            }
        )
    for pos, rec in sorted(neo_c.items(), key=lambda kv: -kv[1]["abundance"]):
        after = pos
        if after in (18, 24) or after >= n:
            continue
        if any(int(c["after_residue"]) == after and c["source"] == "neo_N" for c in candidates):
            continue
        candidates.append(
            {
                "after_residue": after,
                "before_residue": after + 1,
                "type": "internal",
                "source": "neo_C",
                "motif_P4_P4prime": motif(seq, after),
                "protease_hint": protease_hint(seq, after),
                "n_neo_N_peptides": neo_n.get(after + 1, {}).get("n", 0),
                "abundance_neo_N": neo_n.get(after + 1, {}).get("abundance"),
                "n_neo_C_peptides": rec["n"],
                "abundance_neo_C": rec["abundance"],
                "abundance_left_window": region_mean(profiles["mean"], max(1, after - 14), after),
                "abundance_right_window": region_mean(profiles["mean"], after + 1, min(n, after + 15)),
                "supported": "yes",
                "evidence_zh": (
                    f"细胞 neo-C 肽段止于残基 {after}（切割在 {after}|{after+1}），"
                    f"{rec['n']} 条，Σ丰度 {rec['abundance']:.4g}。"
                    f"肽段: {', '.join(rec['peptides'][:8])}"
                ),
                "note_zh": "非胰酶 C 端，优先视为细胞内/分泌路径蛋白酶切割",
            }
        )

    candidates.extend(abundance_steps(profiles["mean"]))

    # 输出目录
    if outdir is None:
        outdir = report.parent / "cleavage_ALB"
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    pep_fields = [
        "Stripped.Sequence",
        "Modified.Sequence",
        "Precursor.Id",
        "Precursor.Charge",
        "start",
        "end",
        "length",
        "prev_aa",
        "next_aa",
        "tryptic_N",
        "tryptic_C",
        "cellular_neo_N",
        "cellular_neo_C",
        "terminus_class",
        "cleavage_after_if_neoN",
        "cleavage_after_if_neoC",
        "Proteotypic",
        "Protein.Ids",
        "Genes",
        "n_matches_on_ALB",
        "n_precursors",
        "abundance_mean",
        *sample_names,
    ]
    write_csv(outdir / "ALB_peptide_map.csv", peptide_level, pep_fields)

    cand_fields = [
        "after_residue",
        "before_residue",
        "type",
        "source",
        "supported",
        "motif_P4_P4prime",
        "protease_hint",
        "n_neo_N_peptides",
        "abundance_neo_N",
        "n_neo_C_peptides",
        "abundance_neo_C",
        "abundance_left_window",
        "abundance_right_window",
        "log2_step",
        "evidence_zh",
        "note_zh",
    ]
    write_csv(outdir / "ALB_candidate_cleavage_sites.csv", candidates, cand_fields)

    res_rows = []
    for i in range(n):
        row = {
            "residue": i + 1,
            "aa": seq[i],
            "region": (
                "signal"
                if i < 18
                else "propeptide"
                if i < 24
                else "mature"
            ),
            "n_peptides": counts[i],
            "abundance_mean": profiles["mean"][i],
        }
        for name in sample_names:
            row[name] = profiles[name][i]
        res_rows.append(row)
    write_csv(
        outdir / "ALB_residue_abundance.csv",
        res_rows,
        ["residue", "aa", "region", "n_peptides", "abundance_mean", *sample_names],
    )

    write_schematic_svg(outdir / "ALB_processing_schematic.svg", seq, candidates, n_term_cov)
    write_full_length_cleavage_svg(
        outdir / "ALB_full_length_cleavage.svg", seq, candidates, n_term_cov
    )
    write_coverage_svg(
        outdir / "ALB_peptide_coverage.svg",
        seq,
        {k: profiles[k] for k in list(sample_names) + ["mean"]},
        peptide_level,
    )

    # 中文结论
    signal_row = next(c for c in candidates if c["type"] == "signal")
    pro_row = next(c for c in candidates if c["type"] == "propeptide")
    internals = [c for c in candidates if c["type"] == "internal" and c["source"] in {"neo_N", "neo_C"}]
    log(logs, "")
    log(logs, "======== ALB 细胞剪切推断 ========")
    log(logs, f"信号肽 18|19 ({signal_row['motif_P4_P4prime']}): {signal_row['evidence_zh']}")
    log(logs, f"Propeptide 24|25 ({pro_row['motif_P4_P4prime']}): {pro_row['evidence_zh']}")
    if internals:
        log(logs, f"内部细胞剪切候选 {len(internals)} 个（neo 末端，已排除胰酶 K/R）:")
        for c in internals[:12]:
            log(logs, f"  {c['after_residue']}|{c['before_residue']}  {c['motif_P4_P4prime']}  {c['protease_hint']}")
    else:
        log(logs, "未发现内部 neo-N/neo-C。若搜库仅为胰酶特异性，内部剪切可能不可见。")
    log(logs, f"结果目录: {outdir}")
    log(logs, "全长剪切示意图: ALB_full_length_cleavage.svg")

    (outdir / "00_inference_log.txt").write_text("\n".join(logs) + "\n", encoding="utf-8")
    return outdir


def write_self_test_matrix(folder: Path) -> Path:
    """Tiny DIA-NN-like matrix: mature ALB peptides + one cellular neo-N + SERPINA5 decoy."""
    seq = ALB_SEQUENCE
    s1 = r"D:\astral\DATA\20260721_sup1.raw"
    s2 = r"D:\astral\DATA\20260721_sup2.raw"
    header = [
        "Protein.Group",
        "Protein.Ids",
        "Protein.Names",
        "Genes",
        "First.Protein.Description",
        "Proteotypic",
        "Stripped.Sequence",
        "Modified.Sequence",
        "Precursor.Charge",
        "Precursor.Id",
        s1,
        s2,
    ]
    # 成熟链胰酶肽（上清应有）；不含 1–24，模拟信号肽/propeptide 已切
    mature = [
        seq[24:34],   # 25–34 DAHKSEVAHR，R24|D25 看起来像胰酶位点
        seq[36:44],   # DLGEENFK
        seq[seq.find("LVNEVTEFAK"): seq.find("LVNEVTEFAK") + 10],
        seq[seq.find("SLHTLFGDK"): seq.find("SLHTLFGDK") + 9],
        seq[seq.find("AAFTECCQAADK"): seq.find("AAFTECCQAADK") + 12],
        seq[seq.find("LVAASQAALGL"): seq.find("LVAASQAALGL") + 11],  # C 端
    ]
    rows = []
    for i, pep in enumerate(mature):
        rows.append(
            ["P02768", "P02768", "ALBU_HUMAN", "ALB", "", "1", pep, pep, "2", pep + "2", str(80000 + i * 1000), str(70000 + i * 800)]
        )
    # 同一肽段电荷 3，测试合并
    pep0 = mature[0]
    rows.append(["P02768", "P02768", "ALBU_HUMAN", "ALB", "", "1", pep0, pep0, "3", pep0 + "3", "21000", "18000"])
    # 细胞剪切 neo-N：从胰酶肽内部开始，上一残基不是 K/R
    neo = seq[seq.find("LVNEVTEFAK") + 1 : seq.find("LVNEVTEFAK") + 10]  # VNEVTEFAK
    rows.append(["P02768", "P02768", "ALBU_HUMAN", "ALB", "", "1", neo, neo, "2", neo + "2", "15000", "12000"])
    # SERPINA5 不得被当成 ALB
    decoy = "AAAATGTIFTR"
    rows.append(["P05154", "P05154", "IPSP_HUMAN", "SERPINA5", "", "1", decoy, decoy, "2", decoy + "2", "390730", "96098.3"])
    path = folder / "report.pr_matrix"
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(header)
        writer.writerows(rows)
    return path


def self_test() -> None:
    import tempfile

    tmp = Path(tempfile.mkdtemp(prefix="zzx_alb_"))
    write_self_test_matrix(tmp)
    out = run(tmp, tmp / "cleavage_ALB", n_last=2)
    sites = list(csv.DictReader((out / "ALB_candidate_cleavage_sites.csv").open(encoding="utf-8")))
    peps = list(csv.DictReader((out / "ALB_peptide_map.csv").open(encoding="utf-8")))
    genes = {r["Genes"] for r in peps}
    assert genes == {"ALB"}, genes
    types = {r["type"] for r in sites}
    assert "signal" in types and "propeptide" in types, types
    assert any(r["source"] == "neo_N" for r in sites), "missing cellular neo-N candidate"
    assert any(r["Stripped.Sequence"] == "DAHKSEVAHR" for r in peps)
    schematic = (out / "ALB_full_length_cleavage.svg").read_text(encoding="utf-8")
    assert "Full-length ALB cleavage" in schematic
    assert "cut 18|19" in schematic and "cut 24|25" in schematic
    assert "cut 66|67" in schematic
    print(f"[self-test] OK  outputs in {out}")


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Infer cellular ALB cleavage from DIA-NN peptide matrix")
    p.add_argument("--data-dir", type=Path, default=None, help="含 report.pr_matrix 的目录")
    p.add_argument("--outdir", type=Path, default=None, help="输出目录，默认 <data>/cleavage_ALB")
    p.add_argument(
        "--n-last",
        type=int,
        default=2,
        help="用最后 N 个样品列作为细胞上清丰度（默认 2）",
    )
    p.add_argument("--self-test", action="store_true", help="用内置小矩阵检查脚本")
    return p


if __name__ == "__main__":
    args = build_parser().parse_args()
    if args.self_test:
        self_test()
    else:
        run(args.data_dir, args.outdir, args.n_last)
