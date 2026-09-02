#!/usr/bin/env python3
"""两组病人血清蛋白质组：按蛋白丰度做排名气泡图。

输入优先为 DIA-NN report.pg_matrix；行过少时回退 report.pr_matrix
按 Protein.Group 取肽段强度中位数。分组必须来自 sample_annotation.csv。
"""
from __future__ import annotations

import csv
import math
import os
import re
import sys
from datetime import datetime
from pathlib import Path

# ---------------------------------------------------------------------------
# 依赖
# ---------------------------------------------------------------------------
try:
    import numpy as np
    import pandas as pd
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except ImportError:
    import subprocess

    subprocess.check_call(
        [sys.executable, "-m", "pip", "install", "pandas", "matplotlib", "numpy"],
        stdout=sys.stderr,
    )
    import numpy as np
    import pandas as pd
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

from matplotlib import font_manager

def _configure_cjk_font() -> None:
    for name in (
        "WenQuanYi Micro Hei",
        "Noto Sans CJK SC",
        "Source Han Sans SC",
        "Droid Sans Fallback",
        "Microsoft YaHei",
        "SimHei",
    ):
        if any(name.lower() in f.name.lower() for f in font_manager.fontManager.ttflist):
            plt.rcParams["font.sans-serif"] = [name, "DejaVu Sans"]
            break
    plt.rcParams["axes.unicode_minus"] = False


_configure_cjk_font()

META_COLS = {
    "Protein.Group",
    "Protein.Ids",
    "Protein.Names",
    "Genes",
    "First.Protein.Description",
    "N.Proteins",
    "Proteotypic",
    "Stripped.Sequence",
    "Modified.Sequence",
    "Precursor.Charge",
    "Precursor.Id",
    "Precursor.Lib.Index",
    "Q.Value",
}
META_PREFIX = re.compile(
    r"^(Q\.Value|PG\.Q|GG\.Q|Protein\.Q|Lib\.|Ms1\.|Normalisation)"
)
TOP_NS = (20, 30, 50)
MIN_DETECT_FRAC = 0.5


def log(msg: str, log_path: Path | None = None) -> None:
    line = f"{datetime.now():%H:%M:%S} | {msg}"
    print(line)
    if log_path is not None:
        with log_path.open("a", encoding="utf-8") as fh:
            fh.write(line + "\n")


DEFAULT_PROTEOMICS_DIRS = (
    Path(r"E:\天府\实验管理\课题\赵章寻\血清蛋白质组学"),
    Path("E:/天府/实验管理/课题/赵章寻/血清蛋白质组学"),
)


def resolve_proteomics_dir() -> Path:
    env = os.environ.get("SERUM_PROTEOMICS_DIR", "")
    candidates = [
        Path(env) if env else None,
        *DEFAULT_PROTEOMICS_DIRS,
        Path.cwd() / "serum_proteomics",
        Path.cwd(),
    ]
    markers = (
        "report.pg_matrix",
        "report.pr_matrix",
        "report.pg_matrix.tsv",
        "sample_annotation.csv",
    )
    for d in candidates:
        if d is None or not d.is_dir():
            continue
        if any((d / m).exists() for m in markers):
            return d.resolve()
    return Path.cwd().resolve()


def find_matrix(directory: Path, stem: str) -> Path | None:
    for name in (stem, f"{stem}.tsv", f"{stem}.txt"):
        p = directory / name
        if p.is_file():
            return p
    return None


def normalize_sample_token(name: str) -> str:
    name = str(name).replace("\\", "/")
    name = os.path.basename(name)
    return re.sub(r"\.(raw|wiff|mzML|d)$", "", name, flags=re.IGNORECASE)


def is_immunoglobulin_symbol(sym: object) -> bool:
    s = str(sym).strip().upper()
    s = re.sub(r"^(CON__|SP\||TR\|)", "", s)
    s = s.split("|")[-1]
    if s in {"", "-", ".", "NA", "NONE"}:
        return False
    if s in {"JCHAIN", "IGJ"}:
        return True
    # UniProt Ig 条目：P0DOX5 重链、P0DOY2 轻链等
    if re.match(r"^P0DO[XY]\d*$", s):
        return True
    return bool(re.match(r"^(IGH|IGK|IGL)", s))


def is_trypsin_symbol(sym: object) -> bool:
    s = str(sym).strip().upper().replace("-", "_")
    s = re.sub(r"^(CON__|SP\||TR\|)", "", s)
    s = s.split("|")[-1]
    if s in {
        "PRSS1", "PRSS2", "PRSS3", "TRY1", "TRY2", "TRY3", "TRYP", "TRYP_PIG",
        "TRY1_BOVIN", "TRY2_BOVIN", "TRY3_BOVIN",
        "P00761", "P00760", "P00763", "P07477", "P07478", "P35030",
    }:
        return True
    return bool(re.match(r"^(PRSS[123]|TRY[123]|TRYP)(_|$)", s))


def is_immunoglobulin_text(*parts: object) -> bool:
    blob = " ".join("" if p is None else str(p) for p in parts).lower()
    if not blob.strip():
        return False
    if "铁蛋白重链" in blob or "ferritin heavy" in blob:
        return False
    if "immunoglobulin superfamily" in blob and not re.search(
        r"\bigh[agmdvejk]|\bigk[cvlj]|\bigl[cvlj]", blob
    ):
        return False
    return bool(
        re.search(
            r"immunoglobulin\s+(heavy|kappa|lambda|alpha|gamma|mu|delta|epsilon)"
            r"|immunoglobulin heavy chain"
            r"|ig heavy chain"
            r"|免疫球蛋白重链"
            r"|免疫球蛋白"
            r"|\big\s*(heavy|kappa|lambda|gamma|alpha|mu)\b"
            r"|\big\s+(gamma|alpha|mu|kappa|lambda)\b"
            r"|heavy chain c region",
            blob,
        )
    )


def is_trypsin_text(*parts: object) -> bool:
    blob = " ".join("" if p is None else str(p) for p in parts).lower()
    if not blob.strip():
        return False
    if re.search(r"antitrypsin|anti-trypsin|trypsin inhibitor|抗胰蛋白酶", blob):
        return False
    if "tryptophan" in blob:
        return False
    return bool(re.search(r"\btrypsin(ogen)?\b|胰蛋白酶", blob))


def _meta_token_blob(meta: pd.DataFrame, i: int) -> tuple[list[str], tuple]:
    n_fields = []
    for col in ("Genes", "Protein.Names", "First.Protein.Description", "Protein.Ids", "Protein.Group"):
        n_fields.append(str(meta[col].iloc[i]) if col in meta.columns else "")
    parts = re.split(r"[,;|/ ]+", ",".join(n_fields))
    return [p for p in parts if p.strip()], tuple(n_fields)


def immunoglobulin_mask(meta: pd.DataFrame) -> np.ndarray:
    n = len(meta)
    out = np.zeros(n, dtype=bool)
    for i in range(n):
        parts, fields = _meta_token_blob(meta, i)
        out[i] = any(is_immunoglobulin_symbol(p) for p in parts) or is_immunoglobulin_text(*fields)
    return out


def trypsin_mask(meta: pd.DataFrame) -> np.ndarray:
    n = len(meta)
    out = np.zeros(n, dtype=bool)
    for i in range(n):
        parts, fields = _meta_token_blob(meta, i)
        out[i] = any(is_trypsin_symbol(p) for p in parts) or is_trypsin_text(*fields)
    return out


def exclusion_reasons(meta: pd.DataFrame) -> list[str]:
    """每行原因：immunoglobulin、trypsin，或空字符串。"""
    ig = immunoglobulin_mask(meta)
    tr = trypsin_mask(meta)
    out: list[str] = []
    for i in range(len(meta)):
        tags = []
        if ig[i]:
            tags.append("immunoglobulin")
        if tr[i]:
            tags.append("trypsin")
        out.append(";".join(tags))
    return out


def pick_official_symbol(value: object) -> str | None:
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return None
    text = str(value).strip()
    if text in {"", "-", ".", "NA"}:
        return None
    parts = [p.strip() for p in re.split(r"[,;|/]+", text)]
    parts = [p for p in parts if p and p not in {"-", ".", "NA"}]
    if not parts:
        return None

    def score(symbol: str) -> int:
        if re.match(r"^(XLOC|TCONS|CUFF)_", symbol, re.I):
            return 0
        if re.match(r"^LOC\d+$", symbol, re.I):
            return 1
        return 2

    return max(parts, key=score)


def is_intensity_col(name: str, series: pd.Series) -> bool:
    if name in META_COLS or META_PREFIX.match(name):
        return False
    if pd.api.types.is_numeric_dtype(series):
        return True
    if series.dtype == object:
        numeric = pd.to_numeric(series, errors="coerce")
        return numeric.notna().mean() >= 0.3
    return False


def read_diann_matrix(path: Path) -> tuple[pd.DataFrame, list[str]]:
    df = pd.read_csv(path, sep="\t")
    int_cols = [c for c in df.columns if is_intensity_col(c, df[c])]
    if not int_cols:
        raise SystemExit(f"没有找到样品强度列: {path}")
    for col in int_cols:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    return df, int_cols


def aggregate_pr_to_pg(pr: pd.DataFrame, int_cols: list[str]) -> tuple[pd.DataFrame, list[str]]:
    key = "Protein.Group" if "Protein.Group" in pr.columns else pr.columns[0]
    agg = pr.groupby(key, dropna=False)[int_cols].median()
    extra = [c for c in ("Protein.Ids", "Protein.Names", "Genes", "First.Protein.Description") if c in pr.columns]
    if extra:
        meta = pr.groupby(key, dropna=False)[extra].first()
        agg = meta.join(agg)
    agg = agg.reset_index()
    return agg, int_cols


def load_protein_matrix(directory: Path, log_path: Path) -> tuple[pd.DataFrame, list[str]]:
    pg_path = find_matrix(directory, "report.pg_matrix")
    pr_path = find_matrix(directory, "report.pr_matrix")
    if pg_path is not None:
        pg, int_cols = read_diann_matrix(pg_path)
        log(f"读取 pg_matrix: {pg_path} ({len(pg)} 行, {len(int_cols)} 个样品列)", log_path)
        if len(pg) >= 5:
            return pg, int_cols
        log("pg_matrix 蛋白行过少，回退 pr_matrix", log_path)
    if pr_path is None:
        raise SystemExit(f"找不到 report.pg_matrix 或 report.pr_matrix，目录: {directory}")
    pr, int_cols = read_diann_matrix(pr_path)
    log(f"读取 pr_matrix: {pr_path}", log_path)
    log("pr_matrix 按 Protein.Group 取肽段强度中位数，聚合到蛋白", log_path)
    return aggregate_pr_to_pg(pr, int_cols)


def write_annotation_template(path: Path, tokens: list[str]) -> None:
    uniq = list(dict.fromkeys(tokens))
    if len(uniq) == 2:
        df = pd.DataFrame({"sample": uniq, "group": uniq})
    else:
        df = pd.DataFrame({"sample": tokens, "group": [""] * len(tokens)})
    df.to_csv(path, index=False, encoding="utf-8-sig")


def load_annotation(
    directory: Path, intensity_cols: list[str], log_path: Path | None = None
) -> pd.DataFrame:
    anno_path = directory / "sample_annotation.csv"
    tokens = [normalize_sample_token(c) for c in intensity_cols]
    if not anno_path.is_file():
        write_annotation_template(anno_path, tokens)
        uniq = list(dict.fromkeys(tokens))
        if len(uniq) != 2:
            raise SystemExit(
                "缺少 sample_annotation.csv。已按样品列写出模板，请把 group 填成恰好两组后重跑。\n"
                f"文件: {anno_path}\n"
                f"当前强度列: {', '.join(tokens)}"
            )
        log(
            f"矩阵恰好两列样品，按 1-vs-1 分组（{uniq[0]} vs {uniq[1]}），不伪造 p 值",
            log_path,
        )
    anno = pd.read_csv(anno_path)
    anno.columns = [c.strip().lower() for c in anno.columns]
    if not {"sample", "group"}.issubset(anno.columns):
        raise SystemExit("sample_annotation.csv 必须包含 sample, group 两列")
    anno["sample_token"] = anno["sample"].map(normalize_sample_token)
    anno["group"] = anno["group"].astype(str).str.strip()
    groups = [g for g in dict.fromkeys(anno["group"]) if g not in {"", "nan", "None"}]
    if len(groups) != 2:
        raise SystemExit(f"必须恰好两组病人，当前: {', '.join(groups) if groups else '(空)'}")

    token_map = dict(zip(anno["sample_token"], anno["group"].astype(str)))
    raw_map = dict(zip(anno["sample"].astype(str), anno["group"].astype(str)))
    mapped = []
    missing = []
    for col, tok in zip(intensity_cols, tokens):
        group = token_map.get(tok, raw_map.get(tok) or raw_map.get(col))
        if group is None:
            missing.append(tok)
        mapped.append({"column": col, "sample": tok, "group": group})
    if missing:
        raise SystemExit(f"注释表对不上这些样品列: {', '.join(missing)}")
    info = pd.DataFrame(mapped)
    info["group"] = pd.Categorical(info["group"], categories=groups, ordered=True)
    return info


def preprocess(mat: pd.DataFrame, int_cols: list[str]) -> tuple[pd.DataFrame, np.ndarray]:
    raw = np.array(mat[int_cols].to_numpy(dtype=float), copy=True)
    raw[~np.isfinite(raw) | (raw < 0)] = np.nan
    detected = np.isfinite(raw) & (raw > 0)
    frac = np.nanmean(detected, axis=1)
    keep = frac >= MIN_DETECT_FRAC
    mat = mat.loc[keep].reset_index(drop=True)
    raw = raw[keep]
    log2x = np.log2(raw + 1.0)
    med = np.nanmedian(log2x, axis=0)
    log2x = log2x - (med - np.nanmedian(med))
    return mat, log2x


def welch_p(a: np.ndarray, b: np.ndarray) -> float:
    a = a[np.isfinite(a)]
    b = b[np.isfinite(b)]
    if a.size < 2 or b.size < 2:
        return np.nan
    va, vb = np.var(a, ddof=1), np.var(b, ddof=1)
    if va == 0 and vb == 0:
        return np.nan
    na, nb = a.size, b.size
    mean_diff = np.mean(b) - np.mean(a)
    se = math.sqrt(va / na + vb / nb)
    if se == 0:
        return np.nan
    t = mean_diff / se
    num = (va / na + vb / nb) ** 2
    den = (va / na) ** 2 / (na - 1) + (vb / nb) ** 2 / (nb - 1)
    df = num / den if den > 0 else np.nan
    if not np.isfinite(df) or df <= 0:
        return np.nan
    try:
        from scipy import stats

        return float(stats.t.sf(abs(t), df) * 2)
    except Exception:
        return float(math.erfc(abs(t) / math.sqrt(2)))


def rank_proteins(
    meta: pd.DataFrame, log2x: np.ndarray, sample_names: list[str] | None = None
) -> pd.DataFrame:
    mean_ab = np.nanmean(log2x, axis=1)
    genes_raw = meta["Genes"] if "Genes" in meta.columns else pd.Series([None] * len(meta))
    pgs = meta["Protein.Group"] if "Protein.Group" in meta.columns else pd.Series(range(len(meta)))
    genes = []
    seen: dict[str, int] = {}
    for raw, pg in zip(genes_raw, pgs):
        g = pick_official_symbol(raw) or str(pg)
        if g in seen:
            seen[g] += 1
            g = f"{g}_{seen[g]}"
        else:
            seen[g] = 0
        genes.append(g)
    out = pd.DataFrame(
        {
            "Protein.Group": pgs.astype(str).to_numpy(),
            "gene": genes,
            "description": meta["First.Protein.Description"].astype(str).to_numpy()
            if "First.Protein.Description" in meta.columns
            else [""] * len(meta),
            "mean_abundance": mean_ab,
        }
    )
    names = sample_names or [f"sample{i+1}" for i in range(log2x.shape[1])]
    for i, name in enumerate(names):
        out[f"sample_{normalize_sample_token(name)}"] = log2x[:, i]
    out["abundance_rank"] = out["mean_abundance"].rank(ascending=False, method="first").astype(int)
    return out.sort_values("abundance_rank").reset_index(drop=True)


def plot_rank_abundance_bubble(
    df: pd.DataFrame,
    n: int,
    outfile: Path,
    label_n: int = 12,
    title_prefix: str = "血清蛋白丰度排名气泡图",
    subtitle: str = "两样品平均丰度后排名；气泡大小一致",
) -> None:
    sub = df[df["abundance_rank"] <= n].sort_values("abundance_rank")
    fig, ax = plt.subplots(figsize=(8.2, 5.4))
    ax.plot(sub["abundance_rank"], sub["mean_abundance"], color="#9ECAE1", linewidth=0.9, zorder=2)
    ax.scatter(
        sub["abundance_rank"],
        sub["mean_abundance"],
        s=42,
        c="#2C7FB8",
        alpha=0.92,
        edgecolors="white",
        linewidths=0.4,
        zorder=3,
    )
    lab = sub[sub["abundance_rank"] <= min(label_n, n)]
    for _, row in lab.iterrows():
        ax.annotate(
            row["gene"],
            (row["abundance_rank"], row["mean_abundance"]),
            textcoords="offset points",
            xytext=(4, 6),
            fontsize=8,
            fontstyle="italic",
            color="0.2",
        )
    ax.set_xlabel("丰度排名")
    ax.set_ylabel("蛋白丰度值")
    ax.set_title(f"{title_prefix} (top{n})", loc="left")
    ax.text(
        0.0,
        1.02,
        subtitle,
        transform=ax.transAxes,
        fontsize=9,
        color="0.35",
    )
    ax.grid(linestyle=":", alpha=0.4)
    fig.tight_layout()
    fig.savefig(outfile, dpi=150, bbox_inches="tight")
    plt.close(fig)


def write_example_data(directory: Path) -> None:
    """写出小型 DIA-NN 风格示例，便于无真实矩阵时跑通流程。"""
    rng = np.random.default_rng(42)
    proteins = [
        ("P02768", "ALB", "Serum albumin"),
        ("P02647", "APOA1", "Apolipoprotein A-I"),
        ("P02787", "TF", "Serotransferrin"),
        ("P00738", "HP", "Haptoglobin"),
        ("P01023", "A2M", "Alpha-2-macroglobulin"),
        ("P01009", "SERPINA1", "Alpha-1-antitrypsin"),
        ("P02790", "HPX", "Hemopexin"),
        ("P02766", "TTR", "Transthyretin"),
        ("P01024", "C3", "Complement C3"),
        ("P00450", "CP", "Ceruloplasmin"),
        ("P02749", "APOH", "Beta-2-glycoprotein 1"),
        ("P02652", "APOA2", "Apolipoprotein A-II"),
        ("P04114", "APOB", "Apolipoprotein B-100"),
        ("P02751", "FN1", "Fibronectin"),
        ("P01008", "SERPINC1", "Antithrombin-III"),
        ("P00734", "F2", "Prothrombin"),
        ("P02774", "GC", "Vitamin D-binding protein"),
        ("P01857", "IGHG1", "Immunoglobulin heavy constant gamma 1"),
        ("P0DOX5", "IGHG1", "Immunoglobulin gamma-1 heavy chain"),
        ("P01834", "IGKC", "Immunoglobulin kappa constant"),
        ("P00761", "TRYP_PIG", "Trypsin"),
        ("P02741", "CRP", "C-reactive protein"),
        ("P05231", "IL6", "Interleukin-6"),
        ("P01375", "TNF", "Tumor necrosis factor"),
        ("P10145", "CXCL8", "Interleukin-8"),
        ("P05121", "SERPINE1", "Plasminogen activator inhibitor 1"),
        ("P01137", "TGFB1", "Transforming growth factor beta-1"),
        ("P01011", "SERPINA3", "Alpha-1-antichymotrypsin"),
        ("P02763", "ORM1", "Alpha-1-acid glycoprotein 1"),
        ("P19652", "ORM2", "Alpha-1-acid glycoprotein 2"),
        ("P00747", "PLG", "Plasminogen"),
        ("P02671", "FGA", "Fibrinogen alpha chain"),
        ("P02675", "FGB", "Fibrinogen beta chain"),
        ("P02679", "FGG", "Fibrinogen gamma chain"),
        ("P02743", "APCS", "Serum amyloid P-component"),
        ("P0DJI8", "SAA1", "Serum amyloid A-1 protein"),
        ("P02735", "SAA2", "Serum amyloid A-2 protein"),
        ("P06727", "APOA4", "Apolipoprotein A-IV"),
        ("P02654", "APOC1", "Apolipoprotein C-I"),
        ("P02656", "APOC3", "Apolipoprotein C-III"),
        ("P05090", "APOD", "Apolipoprotein D"),
        ("P02649", "APOE", "Apolipoprotein E"),
        ("P02753", "RBP4", "Retinol-binding protein 4"),
        ("P02765", "AHSG", "Alpha-2-HS-glycoprotein"),
        ("P04217", "A1BG", "Alpha-1B-glycoprotein"),
        ("P43652", "AFM", "Afamin"),
        ("P25311", "AZGP1", "Zinc-alpha-2-glycoprotein"),
    ]
    g1_samples = [f"GroupA_P{i:02d}" for i in range(1, 7)]
    g2_samples = [f"GroupB_P{i:02d}" for i in range(1, 7)]
    samples = g1_samples + g2_samples
    # 组 B 模拟炎症/急性期升高，白蛋白与部分载脂蛋白略降
    up = {"CRP", "SAA1", "SAA2", "HP", "ORM1", "ORM2", "SERPINA3", "C3", "FGA", "FGB", "FGG", "IL6", "CXCL8"}
    down = {"ALB", "APOA1", "APOA2", "TTR", "RBP4", "TF"}

    rows = []
    for acc, gene, desc in proteins:
        base = rng.uniform(14.0, 22.0)
        rec = {
            "Protein.Group": acc,
            "Protein.Ids": acc,
            "Protein.Names": desc,
            "Genes": gene,
            "First.Protein.Description": desc,
        }
        for s in samples:
            val = base + rng.normal(0, 0.18)
            if gene in up and s.startswith("GroupB"):
                val += rng.uniform(0.6, 1.4)
            if gene in down and s.startswith("GroupB"):
                val -= rng.uniform(0.35, 0.9)
            rec[s] = float(2 ** val)
        rows.append(rec)

    pg_path = directory / "report.pg_matrix"
    fieldnames = ["Protein.Group", "Protein.Ids", "Protein.Names", "Genes", "First.Protein.Description", *samples]
    with pg_path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    # 小型 pr_matrix：每蛋白 2 条肽段，供回退路径测试
    pr_path = directory / "report.pr_matrix"
    pr_fields = fieldnames + ["Stripped.Sequence", "Precursor.Id"]
    with pr_path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=pr_fields, delimiter="\t")
        writer.writeheader()
        for rec in rows[:8]:
            for k in range(1, 3):
                row = dict(rec)
                row["Stripped.Sequence"] = f"PEPTIDE{rec['Genes']}{k}"
                row["Precursor.Id"] = f"{rec['Protein.Group']}_{k}"
                for s in samples:
                    row[s] = rec[s] * rng.uniform(0.85, 1.15)
                writer.writerow(row)

    anno_path = directory / "sample_annotation.csv"
    with anno_path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=["sample", "group"])
        writer.writeheader()
        for s in g1_samples:
            writer.writerow({"sample": s, "group": "GroupA"})
        for s in g2_samples:
            writer.writerow({"sample": s, "group": "GroupB"})


def run_pipeline(
    result_subdir: str = "serum_proteomics_bubble",
    drop_immunoglobulin: bool = False,
    title_prefix: str = "血清蛋白丰度排名气泡图",
    subtitle: str = "两样品平均丰度后排名；气泡大小一致",
) -> Path:
    directory = resolve_proteomics_dir()
    if find_matrix(directory, "report.pg_matrix") is None and find_matrix(directory, "report.pr_matrix") is None:
        directory = Path.cwd() / "serum_proteomics"
        directory.mkdir(parents=True, exist_ok=True)
        write_example_data(directory)
    result_dir = directory / "results" / result_subdir
    result_dir.mkdir(parents=True, exist_ok=True)
    log_path = result_dir / f"log_{datetime.now():%Y%m%d_%H%M%S}.txt"
    log(f"工作目录: {directory}", log_path)

    mat, int_cols = load_protein_matrix(directory, log_path)
    log(f"样品列: {', '.join(int_cols)}；两列取平均后按丰度排名", log_path)
    meta, log2x = preprocess(mat, int_cols)
    log(f"过滤后蛋白数: {len(meta)}", log_path)
    if drop_immunoglobulin:
        reasons = exclusion_reasons(meta)
        drop = np.array([r != "" for r in reasons])
        drop_cols = [c for c in ("Protein.Group", "Protein.Ids", "Protein.Names", "Genes", "First.Protein.Description") if c in meta.columns]
        dropped = meta.loc[drop, drop_cols].copy()
        dropped.insert(len(dropped.columns), "reason", [reasons[i] for i, flag in enumerate(drop) if flag])
        drop_path = result_dir / "removed_Ig_trypsin.csv"
        dropped.to_csv(drop_path, index=False)
        log(
            f"去除免疫球蛋白/重链与胰蛋白酶 {int(drop.sum())} / {len(drop)} -> {drop_path}",
            log_path,
        )
        meta = meta.loc[~drop].reset_index(drop=True)
        log2x = log2x[~drop]
    ranked = rank_proteins(meta, log2x, int_cols)
    ranked.to_csv(result_dir / "protein_abundance_ranking.csv", index=False)
    log(f"写出排名表: {result_dir / 'protein_abundance_ranking.csv'}", log_path)

    ns = sorted({min(n, len(ranked)) for n in TOP_NS if n <= len(ranked)} | {len(ranked)})
    for n_use in ns:
        if n_use < 1:
            continue
        tag = "all" if n_use == len(ranked) else f"top{n_use}"
        # all 与 protein_abundance_ranking.csv 是同一张全表，不再重复写 all_ranked_proteins.csv
        if tag != "all":
            ranked[ranked["abundance_rank"] <= n_use].to_csv(
                result_dir / f"{tag}_ranked_proteins.csv", index=False
            )
        plot_rank_abundance_bubble(
            ranked,
            n_use,
            result_dir / f"{tag}_abundance_rank_bubble.png",
            title_prefix=title_prefix,
            subtitle=subtitle,
        )
        log(f"完成 {tag}", log_path)
    log(f"全部完成 -> {result_dir}", log_path)
    return result_dir


def main() -> None:
    run_pipeline()


if __name__ == "__main__":
    main()
