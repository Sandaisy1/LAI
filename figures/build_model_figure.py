#!/usr/bin/env python3
"""Compose a publication model figure from the original Illustrator artwork.

Proteins, mitochondria, nucleus, lung, and the cell field are the author's
own vector drawings (lifted from Model 图 - S.ai). This script only rescales,
spaces, and labels them. New marks are limited to panels, arrows, and type.
"""

import base64
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent
ASSET = ROOT / "assets"
OUT = ROOT / "Model_red_protein_OMM_DRP1.svg"

FONT = "Arial, Helvetica, sans-serif"
INK = "#1C2834"
MUTED = "#5C6B76"
RED = "#C62828"


def load_sprite(name):
    text = (ASSET / name).read_text(encoding="utf-8")
    m = re.search(r'viewBox="([^"]+)"', text)
    x, y, w, h = [float(v) for v in m.group(1).split()]
    inner = re.search(r"<svg[^>]*>(.*)</svg>", text, re.S).group(1)
    inner = inner.replace("//>", "/>")
    return {"name": name, "x": x, "y": y, "w": w, "h": h, "inner": inner.strip()}


def place(sprite, tx, ty, target_w=None, target_h=None):
    if target_w is not None:
        scale = target_w / sprite["w"]
    else:
        scale = target_h / sprite["h"]
    ox = tx - sprite["x"] * scale
    oy = ty - sprite["y"] * scale
    pw, ph = sprite["w"] * scale, sprite["h"] * scale
    group = (
        f'<g transform="translate({ox:.3f},{oy:.3f}) scale({scale:.5f})">'
        f"{sprite['inner']}</g>"
    )
    return group, (tx, ty, pw, ph)


def T(x, y, text, size=24, fill=INK, weight="normal", anchor="start"):
    safe = (
        text.replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
    )
    return (
        f'<text x="{x:.1f}" y="{y:.1f}" font-family="{FONT}" font-size="{size}" '
        f'font-weight="{weight}" fill="{fill}" text-anchor="{anchor}">{safe}</text>'
    )


def lines(x, y, rows, size=22, fill=INK, weight="normal", anchor="start", leading=None):
    leading = leading or int(size * 1.35)
    chunks = [
        f'<text font-family="{FONT}" font-size="{size}" font-weight="{weight}" '
        f'fill="{fill}" text-anchor="{anchor}">'
    ]
    for i, row in enumerate(rows):
        safe = (
            row.replace("&", "&amp;")
            .replace("<", "&lt;")
            .replace(">", "&gt;")
        )
        dy = 0 if i == 0 else leading
        chunks.append(f'<tspan x="{x:.1f}" y="{y:.1f}" dy="{dy}">{safe}</tspan>')
    chunks.append("</text>")
    return "".join(chunks)


def build():
    mito = load_sprite("obj_00.svg")       # tan mitochondrion
    membrane = load_sprite("obj_01.svg")   # outer-membrane protein assembly
    blue = load_sprite("obj_12.svg")       # DRP1
    green = load_sprite("obj_07.svg")      # green partner
    red = load_sprite("obj_09.svg")        # red protein
    nucleus = load_sprite("obj_14.svg")
    cell_b64 = base64.b64encode((ASSET / "cell_migration.png").read_bytes()).decode()

    parts = []
    A = parts.append
    A(f'''<?xml version="1.0" encoding="UTF-8"?>
<svg xmlns="http://www.w3.org/2000/svg"
     xmlns:inkscape="http://www.inkscape.org/namespaces/inkscape"
     width="180mm" height="220mm" viewBox="0 0 1800 2200">
<title>Red protein suppresses breast-cancer lung metastasis</title>
<desc>Model figure assembled from the original shaded artwork. Red protein localizes to the outer mitochondrial membrane, binds DRP1, and blocks DRP1 binding to the green partner, oligomerization, fission, and fission-induced UPRmt. Without UPRmt, mitochondrial proteins do not enter the nucleus to induce EMT, and lung metastasis is suppressed.</desc>
<defs>
  <linearGradient id="amber" x1="0" y1="0" x2="0" y2="1">
    <stop offset="0" stop-color="#FFE0B2"/>
    <stop offset="1" stop-color="#EF6C00"/>
  </linearGradient>
  <radialGradient id="amberBall" cx="35%" cy="30%" r="70%">
    <stop offset="0" stop-color="#FFE8C2"/>
    <stop offset="0.55" stop-color="#FFA726"/>
    <stop offset="1" stop-color="#E65100"/>
  </radialGradient>
  <radialGradient id="lungFill" cx="32%" cy="28%" r="75%">
    <stop offset="0" stop-color="#FDE7EF"/>
    <stop offset="0.55" stop-color="#F48FB1"/>
    <stop offset="1" stop-color="#EC407A"/>
  </radialGradient>
  <radialGradient id="lungFillDeep" cx="30%" cy="30%" r="75%">
    <stop offset="0" stop-color="#F8BBD0"/>
    <stop offset="0.6" stop-color="#F06292"/>
    <stop offset="1" stop-color="#D81B60"/>
  </radialGradient>
</defs>
<g inkscape:groupmode="layer" inkscape:label="Background" id="Background">
  <rect width="1800" height="2200" fill="#FFFFFF"/>
</g>
''')

    # ----- title -----
    A(f'''<g inkscape:groupmode="layer" inkscape:label="Title" id="Title">
  <rect x="36" y="32" width="8" height="78" rx="4" fill="{RED}"/>
  {T(58, 54, "GRAPHICAL MODEL", size=15, fill=RED, weight="bold")}
  {T(58, 92, "Red protein suppresses breast-cancer lung metastasis", size=32, weight="bold")}
  {lines(58, 128, [
      "It docks on the outer mitochondrial membrane, holds DRP1 away from the green partner,",
      "and thereby shuts off fission-induced UPRmt, nuclear EMT signaling, and lung colonization.",
  ], size=18, fill=MUTED, leading=26)}
</g>''')

    # story path — type only, no stick-figure icons
    A(f'''<g inkscape:groupmode="layer" inkscape:label="Story path" id="StoryPath">
  <rect x="36" y="186" width="1728" height="92" rx="14" fill="#F7F4F1"/>
  <circle cx="70" cy="232" r="16" fill="{RED}"/>
  {T(70, 238, "1", size=16, fill="#fff", weight="bold", anchor="middle")}
  {T(96, 228, "Lung metastasis is suppressed", size=16, weight="bold")}
  {T(96, 250, "in breast cancer", size=14, fill=MUTED)}
  <circle cx="470" cy="232" r="16" fill="#1565C0"/>
  {T(470, 238, "2", size=16, fill="#fff", weight="bold", anchor="middle")}
  {T(496, 228, "Outer membrane: DRP1 is sequestered", size=16, weight="bold")}
  {T(496, 250, "no green-partner binding, no fission, no UPRmt", size=14, fill=MUTED)}
  <circle cx="1088" cy="232" r="16" fill="#5E35B1"/>
  {T(1088, 238, "3", size=16, fill="#fff", weight="bold", anchor="middle")}
  {T(1114, 228, "EMT genes stay off", size=16, weight="bold")}
  {T(1114, 250, "mitochondrial proteins do not enter the nucleus", size=14, fill=MUTED)}
</g>''')

    # panels
    A('''<g inkscape:groupmode="layer" inkscape:label="Panels" id="Panels">
  <rect x="36" y="300" width="860" height="980" rx="18" fill="#FFF9F5" stroke="#F0DFD4" stroke-width="1.4"/>
  <rect x="916" y="300" width="848" height="980" rx="18" fill="#F4F8F8" stroke="#D5E4E2" stroke-width="1.4"/>
  <rect x="36" y="1304" width="560" height="620" rx="18" fill="#FFFFFF" stroke="#E6E0D8" stroke-width="1.4"/>
  <rect x="620" y="1304" width="560" height="620" rx="18" fill="#FFFFFF" stroke="#E6E0D8" stroke-width="1.4"/>
  <rect x="1204" y="1304" width="560" height="620" rx="18" fill="#FFFFFF" stroke="#E6E0D8" stroke-width="1.4"/>
</g>''')

    A(f'''<g inkscape:groupmode="layer" inkscape:label="Panel titles" id="PanelTitles">
  <circle cx="68" cy="338" r="15" fill="#C45C26"/>
  {T(68, 344, "A", size=16, fill="#fff", weight="bold", anchor="middle")}
  {T(94, 344, "On the mitochondrion", size=22, weight="bold")}
  <circle cx="948" cy="338" r="15" fill="#1565C0"/>
  {T(948, 344, "B", size=16, fill="#fff", weight="bold", anchor="middle")}
  {T(974, 344, "At the outer membrane", size=22, weight="bold")}
  <circle cx="68" cy="1340" r="14" fill="#1565C0"/>
  {T(68, 1345, "C", size=14, fill="#fff", weight="bold", anchor="middle")}
  {T(92, 1346, "Fission and UPRmt", size=20, weight="bold")}
  <circle cx="652" cy="1340" r="14" fill="#5E35B1"/>
  {T(652, 1345, "D", size=14, fill="#fff", weight="bold", anchor="middle")}
  {T(676, 1346, "Nucleus and EMT", size=20, weight="bold")}
  <circle cx="1236" cy="1340" r="14" fill="#C2185B"/>
  {T(1236, 1345, "E", size=14, fill="#fff", weight="bold", anchor="middle")}
  {T(1260, 1346, "Lung metastasis", size=20, weight="bold")}
</g>''')

    # ----- author's mitochondrion -----
    mito_g, mito_box = place(mito, 150, 400, target_w=620)
    red_in, _ = place(red, 390, 620, target_w=54)
    red_edge, _ = place(red, 700, 520, target_w=64)
    A(f'''<g inkscape:groupmode="layer" inkscape:label="Mitochondrion" id="Mitochondrion">
  {mito_g}
  <g opacity="0.55">{red_in}</g>
  {red_edge}
  <path d="M430,640 C520,600 600,560 690,545" fill="none" stroke="{RED}" stroke-width="2.4" stroke-dasharray="7 6" stroke-linecap="round"/>
  <path d="M678,538 L706,548 L686,566 Z" fill="{RED}"/>
  {T(160, 382, "Mitochondrial pool", size=16, fill=MUTED)}
  {T(520, 382, "Outer-membrane form", size=16, fill=RED, weight="bold")}
  {lines(150, 1040, [
      "The red protein resides in mitochondria and can",
      "relocalize to the outer membrane. The organelle",
      "stays elongated when fission is blocked.",
  ], size=18, fill=INK, leading=26)}
</g>''')

    # ----- author's membrane assembly -----
    mem_g, mem_box = place(membrane, 980, 430, target_w=720)
    # key proteins under the assembly, author's own sprites
    red_k, _ = place(red, 1000, 1000, target_w=46)
    blu_k, _ = place(blue, 1240, 980, target_h=78)
    grn_k, _ = place(green, 1520, 990, target_w=56)
    A(f'''<g inkscape:groupmode="layer" inkscape:label="Membrane complex" id="MembraneComplex">
  {mem_g}
  {T(1340, 410, "OUTER MITOCHONDRIAL MEMBRANE", size=14, fill="#8D6E63", weight="bold", anchor="middle")}
  {red_k}
  {T(1060, 1030, "Red protein", size=16, fill=RED, weight="bold")}
  {blu_k}
  {T(1290, 1030, "DRP1", size=16, fill="#1565C0", weight="bold")}
  {grn_k}
  {T(1590, 1030, "Green protein", size=16, fill="#2E7D32", weight="bold")}
  <!-- inhibition between the key icons -->
  <line x1="1180" y1="1010" x2="1230" y2="1010" stroke="{RED}" stroke-width="3.2" stroke-linecap="round"/>
  <line x1="1222" y1="994" x2="1222" y2="1026" stroke="{RED}" stroke-width="3.2" stroke-linecap="round"/>
  {lines(980, 1100, [
      "Red protein binds DRP1 and blocks its contact",
      "with the green partner, so DRP1 cannot oligomerize.",
  ], size=18, fill=INK, leading=26)}
</g>''')

    # arrow from B down toward the cascade
    A(f'''<g inkscape:groupmode="layer" inkscape:label="Suppression arrow" id="Suppression">
  <line x1="900" y1="1248" x2="900" y2="1290" stroke="{RED}" stroke-width="3" stroke-linecap="round"/>
  <line x1="888" y1="1278" x2="912" y2="1290" stroke="{RED}" stroke-width="3" stroke-linecap="round"/>
  <line x1="912" y1="1278" x2="888" y2="1290" stroke="{RED}" stroke-width="3" stroke-linecap="round"/>
  {T(918, 1278, "suppresses the cascade below", size=16, fill=RED, weight="bold")}
</g>''')

    # ----- C fission -----
    mito2, _ = place(mito, 80, 1420, target_w=250)
    # a short row of the author's DRP1 along the mitochondrion
    drp_bits = []
    for i, x in enumerate((160, 230, 300)):
        g, _ = place(blue, x, 1468, target_h=46)
        drp_bits.append(g)
    grn_bit, _ = place(green, 230, 1548, target_w=34)
    A(f'''<g inkscape:groupmode="layer" inkscape:label="Fission" id="Fission">
  {mito2}
  {''.join(drp_bits)}
  {grn_bit}
  <circle cx="520" cy="1500" r="28" fill="none" stroke="{RED}" stroke-width="4"/>
  <line x1="500" y1="1480" x2="540" y2="1520" stroke="{RED}" stroke-width="4" stroke-linecap="round"/>
  <rect x="430" y="1560" width="132" height="28" rx="14" fill="#FFEBEE"/>
  {T(496, 1580, "BLOCKED", size=14, fill=RED, weight="bold", anchor="middle")}
  {lines(60, 1760, [
      "DRP1 would assemble on the green",
      "protein and split mitochondria.",
      "Fission-induced UPRmt stays off.",
  ], size=17, fill=INK, leading=24)}
</g>''')

    # ----- D nucleus -----
    nuc, nuc_box = place(nucleus, 700, 1400, target_w=300)
    ambers = []
    for x, y, r in ((650, 1520, 13), (675, 1560, 11), (640, 1575, 10)):
        ambers.append(
            f'<circle cx="{x}" cy="{y}" r="{r}" fill="url(#amberBall)" stroke="#E65100" stroke-width="1"/>'
            f'<ellipse cx="{x-4}" cy="{y-4}" rx="{r*0.35:.1f}" ry="{r*0.22:.1f}" fill="#fff" opacity="0.55"/>'
        )
    A(f'''<g inkscape:groupmode="layer" inkscape:label="Nucleus" id="Nucleus">
  {nuc}
  {''.join(ambers)}
  <line x1="690" y1="1540" x2="730" y2="1540" stroke="{RED}" stroke-width="3.2" stroke-linecap="round"/>
  <line x1="722" y1="1526" x2="722" y2="1554" stroke="{RED}" stroke-width="3.2" stroke-linecap="round"/>
  {T(640, 1490, "UPRmt", size=14, fill="#E65100", weight="bold")}
  {T(640, 1508, "proteins", size=14, fill="#E65100", weight="bold")}
  <rect x="1020" y="1688" width="132" height="28" rx="14" fill="#FFEBEE"/>
  {T(1086, 1708, "BLOCKED", size=14, fill=RED, weight="bold", anchor="middle")}
  {lines(650, 1820, [
      "UPRmt would carry mitochondrial",
      "proteins into the nucleus and",
      "turn on EMT genes. That relay stops.",
  ], size=17, fill=INK, leading=24)}
</g>''')

    # ----- E lung + the author's cell field -----
    A(f'''<g inkscape:groupmode="layer" inkscape:label="Metastasis" id="Metastasis">
  <g transform="translate(1360,1410)">
    <ellipse cx="70" cy="118" rx="78" ry="22" fill="#000" opacity="0.06"/>
    <path d="M78,8 C86,-8 98,-6 104,10 L104,36 C96,28 88,30 82,40 Z" fill="#F48FB1" stroke="#AD1457" stroke-width="2"/>
    <path d="M8,48 C-18,36 -24,78 -8,112 C8,146 48,158 78,140 C62,150 28,146 16,112 C6,84 10,58 8,48 Z" fill="url(#lungFill)" stroke="#C2185B" stroke-width="2.6" stroke-linejoin="round"/>
    <path d="M86,46 C128,18 188,28 196,70 C206,118 176,162 128,168 C108,170 86,150 84,124 C82,96 70,70 86,46 Z" fill="url(#lungFillDeep)" stroke="#C2185B" stroke-width="2.6" stroke-linejoin="round"/>
    <path d="M92,42 C70,70 62,100 70,128" fill="none" stroke="#AD1457" stroke-width="2.2" stroke-linecap="round"/>
    <path d="M96,58 C120,72 142,78 162,96" fill="none" stroke="#AD1457" stroke-width="2" stroke-linecap="round"/>
    <path d="M24,78 C40,86 52,108 48,126" fill="none" stroke="#AD1457" stroke-width="1.8" stroke-linecap="round"/>
  </g>
  {T(1590, 1490, "Lung", size=16, fill="#C2185B", weight="bold")}
  <line x1="1460" y1="1636" x2="1460" y2="1668" stroke="{RED}" stroke-width="3.2" stroke-linecap="round"/>
  <line x1="1444" y1="1660" x2="1476" y2="1660" stroke="{RED}" stroke-width="3.2" stroke-linecap="round"/>
  <image x="1240" y="1680" width="490" height="100" preserveAspectRatio="xMidYMid meet"
         href="data:image/png;base64,{cell_b64}"/>
  <rect x="1608" y="1604" width="132" height="28" rx="14" fill="#FFEBEE"/>
  {T(1674, 1624, "BLOCKED", size=14, fill=RED, weight="bold", anchor="middle")}
  {lines(1240, 1810, [
      "EMT makes these cells migratory.",
      "They would colonize the lung.",
      "That colonization is suppressed.",
  ], size=17, fill=INK, leading=24)}
</g>''')

    # legend
    lg_red, _ = place(red, 160, 1995, target_w=36)
    lg_blue, _ = place(blue, 430, 1978, target_h=52)
    lg_green, _ = place(green, 700, 1988, target_w=40)
    A(f'''<g inkscape:groupmode="layer" inkscape:label="Legend" id="Legend">
  <rect x="36" y="1950" width="1728" height="210" rx="16" fill="#FFFFFF" stroke="#E6E0D8" stroke-width="1.3"/>
  {T(56, 1982, "LEGEND  ·  drawings are from the original figure", size=14, fill=MUTED, weight="bold")}
  {lg_red}
  {T(210, 2020, "Red protein", size=18, weight="bold")}
  {lg_blue}
  {lines(500, 2004, ["DRP1", "dynamin-related protein 1"], size=16, leading=22)}
  {lg_green}
  {lines(760, 2004, ["Green protein", "DRP1-binding partner"], size=16, leading=22)}
  <circle cx="1080" cy="2016" r="12" fill="url(#amberBall)" stroke="#E65100" stroke-width="1"/>
  {lines(1104, 2004, ["UPRmt-associated protein", "moves to the nucleus if UPRmt is on"], size=16, leading=22)}
  <line x1="1480" y1="2016" x2="1536" y2="2016" stroke="{RED}" stroke-width="3.2" stroke-linecap="round"/>
  <line x1="1528" y1="2000" x2="1528" y2="2032" stroke="{RED}" stroke-width="3.2" stroke-linecap="round"/>
  {T(1550, 2022, "Inhibition", size=18, weight="bold")}
  {lines(56, 2088, [
      "Red protein is the brake. Blocking DRP1 stops fission, UPRmt, EMT transcription, and lung metastasis.",
  ], size=16, fill=MUTED)}
</g>''')

    A("</svg>\n")
    OUT.write_text("".join(parts), encoding="utf-8")
    print(OUT, "bytes", OUT.stat().st_size)
    print("mito box", tuple(round(v, 1) for v in mito_box))
    print("membrane box", tuple(round(v, 1) for v in mem_box))
    print("nucleus box", tuple(round(v, 1) for v in nuc_box))


if __name__ == "__main__":
    build()
