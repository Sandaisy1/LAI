#!/usr/bin/env python3
"""In silico knockdown of mitochondrial genes on a mammary EMT continuum.

Model
    GSE143607: MCF10A mammary epithelial cells, TGF-beta time points 0-7 days
    mixed in one 10x library (Watanabe et al.). This is the standard human
    breast epithelial system used to study migration and EMT. It is not a
    patient metastasis cohort.

Virtual knockdown
    The virtual cell is a shrinkage linear model fit on epithelial and
    intermediate cells only (the most mesenchymal 35% are held out).
    Technical covariates (log total counts and mitochondrial-transcript
    fraction) are residualized out. Each mitochondrial gene is then clamped
    to its 5th percentile and the rest of the transcriptome is predicted
    from the partial regression.

    A second fit residualizes log total counts only. A gene is called only
    when both fits, a within-training low-vs-high quintile contrast, and
    both random halves of the training cells agree on direction.

Readouts (gating markers CDH1, EPCAM, KRT18, VIM, FN1, CDH2 are excluded)
    EMT:        Hallmark EMT + GO epithelial to mesenchymal transition
    Migration:  GO positive / epithelial / mesenchymal cell migration,
                scored on genes that are not already in the EMT set
    Invasion:   GO extracellular-matrix disassembly plus a pre-specified
                invasion-effector panel, scored on genes outside EMT and
                the migration set
    Epithelial: tight-junction / keratin / epithelial-stability genes
                (these should fall if the state moves toward EMT)

Call rule, fixed before looking at rankings
    Nuclear-encoded mitochondrial gene, BH FDR < 0.05 on the combined
    score against expression-matched non-mitochondrial genes, favorable
    direction on all four readouts, the same directions in the
    count-only sensitivity fit and in the quintile contrast, and
    combined score > 0 in both training halves.
"""

from __future__ import annotations

import gzip
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.io import mmread
from scipy.sparse import csr_matrix
from scipy.stats import spearmanr
from sklearn.covariance import LedoitWolf
from statsmodels.stats.multitest import multipletests

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data"
GENESET_DIR = DATA / "genesets"
MTX = DATA / "GSE143607"
OUT = ROOT / "virtual_cell_mito_emt"
ACCESSION = "GSE143607"

GATE = ("CDH1", "EPCAM", "KRT18", "VIM", "FN1", "CDH2")
EPI_GATE = ("CDH1", "EPCAM", "KRT18")
MES_GATE = ("VIM", "FN1", "CDH2")

# Held-out epithelial identity. Gate genes are intentionally absent.
EPI_IDENTITY = (
    "OCLN", "CLDN1", "CLDN3", "CLDN4", "CLDN7", "KRT8", "KRT19", "KRT7",
    "MUC1", "CD24", "TJP1", "TJP2", "TJP3", "DSP", "JUP", "ESR1",
    "TACSTD2", "CDH3", "KRT5", "KRT14", "ESRP1", "ESRP2", "GRHL2",
    "OVOL1", "OVOL2", "ST14", "SPINT1", "SPINT2", "MAL2", "EPCAM",
)
# EPCAM is also a gate gene; it is removed from every readout below.

INVASION_EFFECTORS = (
    "MMP1", "MMP2", "MMP3", "MMP7", "MMP9", "MMP10", "MMP11", "MMP13",
    "MMP14", "MMP15", "MMP16", "PLAU", "PLAUR", "CTSB", "CTSD", "CTSL",
    "CTSS", "LOX", "LOXL2", "PLOD2", "P4HA1", "P4HA2", "POSTN", "SPP1",
    "TNC", "LAMC2", "LAMB3", "LAMA3", "ITGA2", "ITGA5", "ITGA6", "ITGAV",
    "ITGB1", "TGFBI", "SERPINE1", "ANGPTL4", "PTGS2", "CXCL8", "IL6",
    "MET", "CD44", "S100A4", "AXL", "RHOC", "RHOA", "RAC1", "CDC42",
    "ADAM12", "ADAMTS1", "FAP", "THBS1", "THBS2", "HAS2",
)

MARKERS = (
    "ZEB1", "SNAI1", "SNAI2", "TWIST1", "SPARC", "SERPINE1", "LOXL2",
    "LAMC2", "TGFBI", "COL1A1", "ITGAV", "MMP2", "MMP14", "AXL", "CD44",
    "OCLN", "CLDN4", "CLDN7", "KRT8", "KRT19", "CD24", "TJP1", "DSP",
    "GRHL2", "ESRP1",
)

REFERENCE_GENES = ("SNAI1", "ZEB1", "TGFB1", "VIM", "CDH1", "GAPDH", "ACTB")

CC_TERMS = (
    "Mitochondrial Envelope (GO:0005740)",
    "Mitochondrial Inner Membrane (GO:0005743)",
    "Mitochondrial Intermembrane Space (GO:0005758)",
    "Mitochondrial Matrix (GO:0005759)",
    "Mitochondrial Membrane (GO:0031966)",
    "Mitochondrial Outer Membrane (GO:0005741)",
    "Mitochondrial Proton-Transporting ATP Synthase Complex (GO:0005753)",
    "Mitochondrial Respiratory Chain Complex I (GO:0005747)",
    "Mitochondrial Respiratory Chain Complex III (GO:0005750)",
    "Mitochondrial Respiratory Chain Complex IV (GO:0005751)",
    "Respiratory Chain Complex I (GO:0045271)",
    "Respiratory Chain Complex III (GO:0045275)",
    "Respiratory Chain Complex IV (GO:0045277)",
    "Mitochondrial Ribosome (GO:0005761)",
    "TIM22 Mitochondrial Import Inner Membrane Insertion Complex (GO:0042721)",
    "TIM23 Mitochondrial Import Inner Membrane Translocase Complex (GO:0005744)",
    "Mitochondrial Outer Membrane Translocase Complex (GO:0005742)",
    "Mitochondrial Alpha-Ketoglutarate Dehydrogenase Complex (GO:0005947)",
)
BP_TERMS = (
    "Oxidative Phosphorylation (GO:0006119)",
    "Aerobic Electron Transport Chain (GO:0019646)",
    "Electron Transport Chain (GO:0022900)",
    "Respiratory Electron Transport Chain (GO:0022904)",
    "ATP Synthesis Coupled Electron Transport (GO:0042773)",
    "Mitochondrial ATP Synthesis Coupled Electron Transport (GO:0042775)",
    "Fatty Acid Beta-Oxidation (GO:0006635)",
    "Mitochondrial Translation (GO:0032543)",
    "Mitochondrial Transport (GO:0006839)",
    "Mitochondrion Organization (GO:0007005)",
    "Mitophagy (GO:0000423)",
    "Tricarboxylic Acid Metabolic Process (GO:0072350)",
)
KEGG_TERMS = (
    "Citrate cycle (TCA cycle)",
    "Oxidative phosphorylation",
    "Mitophagy",
    "Fatty acid degradation",
)
# 10x symbols for the 13 protein-coding mtDNA genes. Enrichr KEGG uses COX1/ND1/CYTB/ATP6.
MTDNA_GENES = (
    "MT-ND1", "MT-ND2", "MT-ND3", "MT-ND4", "MT-ND4L", "MT-ND5", "MT-ND6",
    "MT-CO1", "MT-CO2", "MT-CO3", "MT-CYB", "MT-ATP6", "MT-ATP8",
)
EMT_TERMS = (
    "Epithelial Mesenchymal Transition",
    "Epithelial To Mesenchymal Transition (GO:0001837)",
    "Positive Regulation Of Epithelial To Mesenchymal Transition (GO:0010718)",
)
MIGRATION_TERMS = (
    "Positive Regulation Of Cell Migration (GO:0030335)",
    "Epithelial Cell Migration (GO:0010631)",
    "Positive Regulation Of Epithelial Cell Migration (GO:0010634)",
    "Mesenchymal Cell Migration (GO:0090497)",
)
INVASION_TERMS = (
    "Extracellular Matrix Disassembly (GO:0022617)",
    "Positive Regulation Of Extracellular Matrix Disassembly (GO:0090091)",
)

MIN_SET = 15
FDR_CUT = 0.05
TRAIN_DROP_MESENCHYMAL_Q = 0.35
KD_QUANTILE = 0.05
HVG_N = 2000
RNG = np.random.default_rng(0)


def read_gmt(path: Path) -> dict[str, list[str]]:
    sets: dict[str, list[str]] = {}
    with path.open() as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            sets[parts[0]] = [gene for gene in parts[2:] if gene]
    return sets


def take_terms(library: dict[str, list[str]], names: tuple[str, ...], label: str) -> set[str]:
    missing = [name for name in names if name not in library]
    if missing:
        raise SystemExit(f"{label} is missing {missing}")
    genes: set[str] = set()
    for name in names:
        genes.update(library[name])
    return genes


def load_expression() -> tuple[csr_matrix, list[str]]:
    """Return cells x genes counts and unique gene symbols."""
    gene_path = MTX / "GSM4263710_MCF10Atimecourse_genes.tsv.gz"
    matrix_path = MTX / "GSM4263710_MCF10Atimecourse_matrix.mtx.gz"
    symbols: list[str] = []
    with gzip.open(gene_path, "rt") as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            symbols.append(parts[1] if len(parts) > 1 else parts[0])
    print("reading 10x matrix", flush=True)
    raw = csr_matrix(mmread(matrix_path)).astype(np.float32)
    totals = np.asarray(raw.sum(axis=1)).ravel()
    best: dict[str, int] = {}
    for i, gene in enumerate(symbols):
        if gene not in best or totals[i] > totals[best[gene]]:
            best[gene] = i
    order = [best[gene] for gene in best]
    genes = list(best.keys())
    # cells x genes
    matrix = raw[order].T.tocsr()
    return matrix, genes


def qc_mask(matrix: csr_matrix, genes: list[str]) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    n_counts = np.asarray(matrix.sum(axis=1)).ravel()
    n_genes = np.asarray((matrix > 0).sum(axis=1)).ravel()
    mt_cols = [i for i, gene in enumerate(genes) if gene.startswith("MT-")]
    mt_counts = np.asarray(matrix[:, mt_cols].sum(axis=1)).ravel()
    mt_frac = mt_counts / np.maximum(n_counts, 1.0)
    keep = (
        (n_genes >= 1000)
        & (n_genes <= 5000)
        & (n_counts >= 4000)
        & (n_counts <= 40000)
        & (mt_frac < 0.12)
    )
    return keep, n_counts, mt_frac


def lognorm_columns(matrix: csr_matrix, cell_idx: np.ndarray, gene_idx: np.ndarray, n_counts: np.ndarray) -> np.ndarray:
    block = matrix[cell_idx][:, gene_idx].toarray().astype(np.float64)
    scale = 10000.0 / np.maximum(n_counts[cell_idx], 1.0)
    return np.log1p(block * scale[:, None])


def zscore(values: np.ndarray) -> np.ndarray:
    sd = values.std()
    if sd < 1e-8:
        return np.zeros_like(values)
    return (values - values.mean()) / sd


def state_score(log_expr: np.ndarray, genes: list[str], gene_idx: np.ndarray) -> np.ndarray:
    lookup = {gene: i for i, gene in enumerate(genes)}
    cols = {gene: pos for pos, gene in enumerate(gene_idx_to_symbols(genes, gene_idx))}
    # log_expr columns follow gene_idx
    def mean_z(names: tuple[str, ...]) -> np.ndarray:
        present = [cols[name] for name in names if name in cols]
        if not present:
            raise SystemExit(f"gating genes missing: {names}")
        mat = np.column_stack([zscore(log_expr[:, j]) for j in present])
        return mat.mean(axis=1)

    return mean_z(EPI_GATE) - mean_z(MES_GATE)


def gene_idx_to_symbols(genes: list[str], gene_idx: np.ndarray) -> list[str]:
    return [genes[i] for i in gene_idx]


def residualize(expr: np.ndarray, covariates: np.ndarray) -> np.ndarray:
    design = np.column_stack([np.ones(expr.shape[0]), covariates])
    coef, _, _, _ = np.linalg.lstsq(design, expr, rcond=None)
    return expr - design @ coef


def virtual_shift(expr: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Clamp each gene to its 5th percentile and predict the other genes.

    Returns shift (genes x genes) and the clamp delta per gene.
    shift[j, g] is the predicted change of gene j when gene g is lowered.
    """
    work = expr.copy()
    stable = work.std(axis=0) >= 1e-6
    work[:, ~stable] = 0.0
    model = LedoitWolf().fit(work)
    cov = model.covariance_
    var = np.clip(np.diag(cov), 1e-8, None)
    beta = cov / var
    center = work.mean(axis=0)
    low = np.quantile(work, KD_QUANTILE, axis=0)
    delta = low - center
    delta[~stable] = 0.0
    shift = beta * delta
    shift[:, ~stable] = np.nan
    np.fill_diagonal(shift, np.nan)
    return shift, delta


def group_mean(shift: np.ndarray, rows: np.ndarray) -> np.ndarray:
    if rows.size == 0:
        return np.full(shift.shape[1], np.nan)
    return np.nanmean(shift[rows], axis=0)


def contrast(shift: np.ndarray, rows: np.ndarray, background: np.ndarray) -> np.ndarray:
    return group_mean(shift, rows) - group_mean(shift, background)


def empirical_greater(stat: np.ndarray, null: np.ndarray) -> np.ndarray:
    null = null[np.isfinite(null)]
    order = np.sort(null)
    n_ge = len(order) - np.searchsorted(order, stat, side="left")
    p = (1.0 + n_ge) / (1.0 + len(order))
    p = np.where(np.isfinite(stat), p, np.nan)
    return p


def z_against(stat: np.ndarray, null: np.ndarray) -> np.ndarray:
    sd = np.nanstd(null)
    if not np.isfinite(sd) or sd < 1e-8:
        return np.zeros_like(stat)
    return (stat - np.nanmean(null)) / sd


def index_of(symbols: list[str], wanted: set[str]) -> np.ndarray:
    return np.array([i for i, gene in enumerate(symbols) if gene in wanted], dtype=int)


def quintile_contrast(expr: np.ndarray, rows: np.ndarray, background: np.ndarray) -> np.ndarray:
    """Low-quintile minus high-quintile expression, program relative to background."""
    n_genes = expr.shape[1]
    out = np.full(n_genes, np.nan)
    for gene in range(n_genes):
        values = expr[:, gene]
        lo = np.quantile(values, 0.20)
        hi = np.quantile(values, 0.80)
        if not np.isfinite(lo) or hi <= lo:
            continue
        low = values <= lo
        high = values >= hi
        if low.sum() < 30 or high.sum() < 30:
            continue
        use_rows = rows[rows != gene] if rows.size else rows
        use_bg = background[background != gene] if background.size else background
        if use_rows.size < 5 or use_bg.size < 5:
            continue
        def rel(mask: np.ndarray) -> float:
            return float(expr[mask][:, use_rows].mean() - expr[mask][:, use_bg].mean())

        out[gene] = rel(low) - rel(high)
    return out


def program_bundle(shift: np.ndarray, sets: dict[str, np.ndarray], background: np.ndarray) -> dict[str, np.ndarray]:
    scores = {name: contrast(shift, rows, background) for name, rows in sets.items()}
    return scores


def combined_from_scores(scores: dict[str, np.ndarray], null_idx: np.ndarray) -> tuple[np.ndarray, dict[str, np.ndarray]]:
    pieces = {
        "emt": scores["emt"],
        "migration_independent": scores["migration_independent"],
        "invasion_independent": scores["invasion_independent"],
        "epithelial": scores["epithelial"],
    }
    z = {name: z_against(values, values[null_idx]) for name, values in pieces.items()}
    combined = z["emt"] + z["migration_independent"] + z["invasion_independent"] - z["epithelial"]
    return combined, z


def self_test() -> None:
    rng = np.random.default_rng(1)
    n = 500
    mito = rng.normal(size=n)
    emt = -1.8 * mito[:, None] + 0.4 * rng.normal(size=(n, 12))
    mig = -1.4 * mito[:, None] + 0.4 * rng.normal(size=(n, 10))
    inv = -1.2 * mito[:, None] + 0.4 * rng.normal(size=(n, 8))
    epi = 1.5 * mito[:, None] + 0.4 * rng.normal(size=(n, 8))
    bg = 0.05 * mito[:, None] + rng.normal(size=(n, 40))
    expr = np.hstack([mito[:, None], emt, mig, inv, epi, bg])
    shift, delta = virtual_shift(expr)
    assert delta[0] < 0
    sets = {
        "emt": np.arange(1, 13),
        "migration_independent": np.arange(13, 23),
        "invasion_independent": np.arange(23, 31),
        "epithelial": np.arange(31, 39),
    }
    background = np.arange(39, 79)
    scores = program_bundle(shift, sets, background)
    assert scores["emt"][0] > 0, scores["emt"][0]
    assert scores["migration_independent"][0] > 0
    assert scores["invasion_independent"][0] > 0
    assert scores["epithelial"][0] < 0, scores["epithelial"][0]
    q = quintile_contrast(expr, sets["emt"], background)
    assert q[0] > 0, q[0]
    print("self-test passed: lowering a gene anti-correlated with EMT raises the EMT program", flush=True)


def present(symbols_in_universe: list[str], genes: set[str] | tuple[str, ...]) -> set[str]:
    have = set(symbols_in_universe)
    return {gene for gene in genes if gene in have}


def main() -> None:
    self_test()
    OUT.mkdir(parents=True, exist_ok=True)
    hallmark = read_gmt(GENESET_DIR / "MSigDB_Hallmark_2020.gmt")
    go_cc = read_gmt(GENESET_DIR / "GO_CC_2023.gmt")
    go_bp = read_gmt(GENESET_DIR / "GO_BP_2023.gmt")
    kegg = read_gmt(GENESET_DIR / "KEGG_2021_Human.gmt")

    mito = (
        take_terms(go_cc, CC_TERMS, "GO CC")
        | take_terms(go_bp, BP_TERMS, "GO BP")
        | take_terms(kegg, KEGG_TERMS, "KEGG")
        | set(MTDNA_GENES)
    )
    emt = take_terms(hallmark, ("Epithelial Mesenchymal Transition",), "Hallmark") | take_terms(
        go_bp, EMT_TERMS[1:], "EMT GO"
    )
    migration = take_terms(go_bp, MIGRATION_TERMS, "migration")
    invasion = take_terms(go_bp, INVASION_TERMS, "invasion") | set(INVASION_EFFECTORS)
    epithelial = set(EPI_IDENTITY)

    drop = set(GATE)
    emt -= drop
    migration -= drop
    invasion -= drop
    epithelial -= drop
    migration_independent = set(migration) - emt
    invasion_independent = set(invasion) - emt - migration

    matrix, genes = load_expression()
    keep, n_counts, mt_frac = qc_mask(matrix, genes)
    print(f"cells {matrix.shape[0]}, after QC {int(keep.sum())}", flush=True)
    cell_idx = np.where(keep)[0]
    gate_idx = np.array([genes.index(gene) for gene in GATE])
    gate_log = lognorm_columns(matrix, cell_idx, gate_idx, n_counts)
    state = state_score(gate_log, genes, gate_idx)
    # state is aligned to cell_idx
    threshold = np.quantile(state, TRAIN_DROP_MESENCHYMAL_Q)
    train_local = state >= threshold
    train_cells = cell_idx[train_local]
    dropped_cells = cell_idx[~train_local]
    print(
        f"training cells {train_cells.size} (epithelial state >= q{TRAIN_DROP_MESENCHYMAL_Q:.2f}); "
        f"held-out mesenchymal cells {dropped_cells.size}",
        flush=True,
    )

    # Sanity: training cells are the epithelial side of the gate.
    gate_lookup = {gene: i for i, gene in enumerate(GATE)}
    cdh1_train = gate_log[train_local, gate_lookup["CDH1"]].mean()
    cdh1_drop = gate_log[~train_local, gate_lookup["CDH1"]].mean()
    vim_train = gate_log[train_local, gate_lookup["VIM"]].mean()
    vim_drop = gate_log[~train_local, gate_lookup["VIM"]].mean()
    if not (cdh1_train > cdh1_drop and vim_train < vim_drop):
        raise SystemExit(
            f"state gate failed: CDH1 train/drop {cdh1_train:.3f}/{cdh1_drop:.3f}, "
            f"VIM train/drop {vim_train:.3f}/{vim_drop:.3f}"
        )

    detected = np.asarray((matrix[train_cells] > 0).mean(axis=0)).ravel()
    hvg_pool = np.where(detected >= 0.10)[0]
    hvg_log = lognorm_columns(matrix, train_cells, hvg_pool, n_counts)
    variance = hvg_log.var(axis=0)
    order = np.argsort(variance)[::-1][:HVG_N]
    hvg_genes = {genes[hvg_pool[i]] for i in order}

    mito_tested = {gene for gene in mito if gene in genes and detected[genes.index(gene)] >= 0.10}
    readout_genes = emt | migration | invasion | epithelial
    readout_kept = {
        gene for gene in readout_genes if gene in genes and detected[genes.index(gene)] >= 0.05
    }
    forced = {
        gene for gene in (*REFERENCE_GENES, *MARKERS)
        if gene in genes and detected[genes.index(gene)] >= 0.05
    }
    universe = sorted(hvg_genes | mito_tested | readout_kept | forced)
    uni_idx = np.array([genes.index(gene) for gene in universe], dtype=int)
    print(f"universe {len(universe)}; mitochondrial genes tested {len(mito_tested)}", flush=True)

    emt_u = present(universe, emt)
    mig_u = present(universe, migration_independent)
    inv_u = present(universe, invasion_independent)
    epi_u = present(universe, epithelial)
    set_sizes = {
        "emt": len(emt_u),
        "migration_independent": len(mig_u),
        "invasion_independent": len(inv_u),
        "epithelial": len(epi_u),
    }
    print("readout sizes", set_sizes, flush=True)
    too_small = {name: n for name, n in set_sizes.items() if n < MIN_SET}
    if too_small:
        raise SystemExit(f"readout sets smaller than {MIN_SET}: {too_small}")

    expr = lognorm_columns(matrix, train_cells, uni_idx, n_counts)
    cov_primary = np.column_stack([
        np.log1p(n_counts[train_cells]),
        mt_frac[train_cells],
    ])
    cov_counts = np.log1p(n_counts[train_cells])[:, None]
    resid = residualize(expr, cov_primary)
    resid_counts = residualize(expr, cov_counts)
    # Drop genes with no residual variance from scoring later; keep the matrix aligned.
    print("fitting primary virtual cell", flush=True)
    shift, delta = virtual_shift(resid)
    print("fitting count-only sensitivity model", flush=True)
    shift_counts, _ = virtual_shift(resid_counts)

    symbol_index = {gene: i for i, gene in enumerate(universe)}
    sets = {
        "emt": index_of(universe, emt_u),
        "migration_independent": index_of(universe, mig_u),
        "invasion_independent": index_of(universe, inv_u),
        "epithelial": index_of(universe, epi_u),
    }
    program_union = set().union(emt_u, mig_u, inv_u, epi_u, GATE)
    null_genes = [
        gene for gene in universe
        if gene not in mito and gene not in program_union and not gene.startswith("MT-")
    ]
    null_idx = np.array([symbol_index[gene] for gene in null_genes], dtype=int)
    background = null_idx.copy()
    print(f"null genes {null_idx.size}", flush=True)

    scores = program_bundle(shift, sets, background)
    combined, z = combined_from_scores(scores, null_idx)
    scores_c = program_bundle(shift_counts, sets, background)
    combined_c, _ = combined_from_scores(scores_c, null_idx)

    print("quintile contrast", flush=True)
    q_scores = {name: quintile_contrast(resid, rows, background) for name, rows in sets.items()}
    q_combined = (
        z_against(q_scores["emt"], q_scores["emt"][null_idx])
        + z_against(q_scores["migration_independent"], q_scores["migration_independent"][null_idx])
        + z_against(q_scores["invasion_independent"], q_scores["invasion_independent"][null_idx])
        - z_against(q_scores["epithelial"], q_scores["epithelial"][null_idx])
    )

    print("split-half stability", flush=True)
    half_z = []
    perm = RNG.permutation(expr.shape[0])
    mid = expr.shape[0] // 2
    for part in (perm[:mid], perm[mid:]):
        part_resid = residualize(expr[part], cov_primary[part])
        part_shift, _ = virtual_shift(part_resid)
        part_scores = program_bundle(part_shift, sets, background)
        part_combined, _ = combined_from_scores(part_scores, null_idx)
        half_z.append(part_combined)
    half_a, half_b = half_z

    p_combined = empirical_greater(combined, combined[null_idx])
    # Hypothesis family is nuclear-encoded mitochondrial genes. mtDNA-encoded
    # transcripts are scored and reported, but they are not part of this FDR.
    mito_rows = [
        i for i, gene in enumerate(universe)
        if gene in mito_tested and not gene.startswith("MT-")
    ]
    mito_p = np.array([p_combined[i] for i in mito_rows])
    mito_fdr = np.full(len(universe), np.nan)
    if len(mito_rows):
        _, fdr, _, _ = multipletests(np.nan_to_num(mito_p, nan=1.0), method="fdr_bh")
        for row, value in zip(mito_rows, fdr):
            mito_fdr[row] = value

    # Descriptive association on all QC cells, not used for the call.
    state_by_cell = np.full(matrix.shape[0], np.nan)
    state_by_cell[cell_idx] = state
    spearman = np.full(len(universe), np.nan)
    # Compute Spearman against the epithelial state on QC cells for universe genes.
    # Use lognorm of universe on QC cells in chunks of genes to limit memory.
    qc_log = lognorm_columns(matrix, cell_idx, uni_idx, n_counts)
    for j in range(len(universe)):
        rho, _ = spearmanr(qc_log[:, j], state)
        spearman[j] = rho

    def direction_ok(score_map: dict[str, np.ndarray], i: int) -> bool:
        return bool(
            score_map["emt"][i] > 0
            and score_map["migration_independent"][i] > 0
            and score_map["invasion_independent"][i] > 0
            and score_map["epithelial"][i] < 0
        )

    records = []
    for i, gene in enumerate(universe):
        is_mito = gene in mito_tested
        is_ref = gene in REFERENCE_GENES
        if not is_mito and not is_ref:
            continue
        q_ok = bool(
            np.isfinite(q_scores["emt"][i])
            and q_scores["emt"][i] > 0
            and q_scores["migration_independent"][i] > 0
            and q_scores["invasion_independent"][i] > 0
            and q_scores["epithelial"][i] < 0
        )
        primary_ok = direction_ok(scores, i)
        sensitivity_ok = direction_ok(scores_c, i) and combined_c[i] > 0
        halves_ok = bool(half_a[i] > 0 and half_b[i] > 0)
        fdr = mito_fdr[i]
        called = bool(
            is_mito
            and (not gene.startswith("MT-"))
            and np.isfinite(fdr)
            and fdr < FDR_CUT
            and primary_ok
            and sensitivity_ok
            and q_ok
            and halves_ok
            and combined[i] > 0
        )
        records.append({
            "gene": gene,
            "mitochondrial": is_mito,
            "mtdna_encoded": gene.startswith("MT-"),
            "reference_only": is_ref and not is_mito,
            "mean_log_expr_training": float(expr[:, i].mean()),
            "clamp_delta": float(delta[i]),
            "score_emt": float(scores["emt"][i]),
            "score_migration_independent": float(scores["migration_independent"][i]),
            "score_invasion_independent": float(scores["invasion_independent"][i]),
            "score_epithelial": float(scores["epithelial"][i]),
            "z_emt": float(z["emt"][i]),
            "z_migration_independent": float(z["migration_independent"][i]),
            "z_invasion_independent": float(z["invasion_independent"][i]),
            "z_epithelial": float(z["epithelial"][i]),
            "z_combined": float(combined[i]),
            "z_combined_count_only": float(combined_c[i]),
            "z_combined_quintile": float(q_combined[i]),
            "z_combined_half_a": float(half_a[i]),
            "z_combined_half_b": float(half_b[i]),
            "p_empirical": float(p_combined[i]),
            "fdr_bh": float(fdr) if np.isfinite(fdr) else np.nan,
            "spearman_with_epithelial_state": float(spearman[i]),
            "direction_primary": primary_ok,
            "direction_count_only": bool(direction_ok(scores_c, i)),
            "direction_quintile": q_ok,
            "stable_halves": halves_ok,
            "called": called,
        })

    table = pd.DataFrame.from_records(records)
    mito_table = table[table["mitochondrial"]].sort_values(
        ["called", "fdr_bh", "z_combined"], ascending=[False, True, False]
    )
    called = mito_table[mito_table["called"]].copy()
    mito_table.to_csv(OUT / "mito_virtual_kd_all.csv", index=False)
    called.to_csv(OUT / "mito_virtual_kd_called.csv", index=False)
    table[table["reference_only"]].to_csv(OUT / "reference_gene_kd.csv", index=False)

    # Marker-level predicted shifts for the top nuclear genes.
    top = mito_table[~mito_table["mtdna_encoded"]].head(15)
    marker_rows = []
    for gene in MARKERS:
        if gene not in symbol_index:
            continue
        row = {"marker": gene}
        j = symbol_index[gene]
        for kd in top["gene"]:
            row[kd] = float(shift[j, symbol_index[kd]])
        marker_rows.append(row)
    markers = pd.DataFrame(marker_rows)
    markers.to_csv(OUT / "top15_predicted_marker_shift.csv", index=False)

    summary = {
        "accession": ACCESSION,
        "model": "MCF10A TGF-beta EMT continuum, 10x, mixed days 0-7",
        "n_cells_raw": int(matrix.shape[0]),
        "n_cells_qc": int(keep.sum()),
        "n_training_cells": int(train_cells.size),
        "n_heldout_mesenchymal_cells": int(dropped_cells.size),
        "gate_check_cdh1_train_vs_heldout": [float(cdh1_train), float(cdh1_drop)],
        "gate_check_vim_train_vs_heldout": [float(vim_train), float(vim_drop)],
        "n_universe": len(universe),
        "n_mito_tested": int(mito_table.shape[0]),
        "n_nuclear_mito_tested": int((~mito_table["mtdna_encoded"]).sum()),
        "n_null": int(null_idx.size),
        "readout_sizes": set_sizes,
        "n_called": int(called.shape[0]),
        "called_genes": called["gene"].tolist(),
        "fdr_cutoff": FDR_CUT,
        "kd_quantile": KD_QUANTILE,
        "training_rule": f"epithelial state score >= quantile {TRAIN_DROP_MESENCHYMAL_Q}",
    }
    (OUT / "run_summary.json").write_text(json.dumps(summary, indent=2))
    plot_results(state, threshold, mito_table, markers, summary)
    print(json.dumps(summary, indent=2), flush=True)


def plot_results(state: np.ndarray, threshold: float, mito_table: pd.DataFrame, markers: pd.DataFrame, summary: dict) -> None:
    nuclear = mito_table[~mito_table["mtdna_encoded"]].copy()
    fig, ax = plt.subplots(figsize=(6.2, 4.2))
    ax.hist(state, bins=40, color="#4C78A8", edgecolor="white")
    ax.axvline(threshold, color="#E45756", lw=1.5, label="training cutoff")
    ax.set_xlabel("Epithelial state score (gate markers only)")
    ax.set_ylabel("Cells")
    ax.set_title("MCF10A training cells sit on the epithelial side")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(OUT / "state_score_histogram.png", dpi=140)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6.4, 5.2))
    specific = (
        (nuclear["z_emt"] > 0)
        & (nuclear["z_migration_independent"] > 0)
        & (nuclear["z_invasion_independent"] > 0)
        & (nuclear["z_epithelial"] < 0)
    )
    colors = np.where(specific, "#E45756", "#4C78A8")
    ax.scatter(
        nuclear["z_emt"], nuclear["z_migration_independent"],
        c=colors, s=18, alpha=0.85, linewidths=0,
    )
    ax.axvline(0, color="0.6", lw=0.6)
    ax.axhline(0, color="0.6", lw=0.6)
    ax.set_xlabel("EMT program z (virtual knockdown)")
    ax.set_ylabel("Migration program z, genes outside EMT")
    ax.set_title("Red: all four readouts move toward EMT")
    fig.tight_layout()
    fig.savefig(OUT / "emt_vs_migration_z.png", dpi=140)
    plt.close(fig)

    show = nuclear.sort_values("z_combined", ascending=False).head(20).iloc[::-1]
    fig, ax = plt.subplots(figsize=(7.2, 6.2))
    colors = ["#E45756" if flag else "#9E9E9E" for flag in show["called"]]
    ax.barh(show["gene"], show["z_combined"], color=colors)
    ax.axvline(0, color="0.4", lw=0.6)
    ax.set_xlabel("Combined score (EMT + migration + invasion - epithelial)")
    ax.set_title("Largest scores; red would be FDR < 0.05 calls")
    fig.tight_layout()
    fig.savefig(OUT / "top20_combined_z.png", dpi=140)
    plt.close(fig)

    directional = nuclear[
        (nuclear["z_emt"] > 0)
        & (nuclear["z_migration_independent"] > 0)
        & (nuclear["z_invasion_independent"] > 0)
        & (nuclear["z_epithelial"] < 0)
    ].sort_values("z_combined", ascending=False)
    if not directional.empty:
        show_d = directional.head(15).iloc[::-1]
        fig, ax = plt.subplots(figsize=(7.2, 5.4))
        ax.barh(show_d["gene"], show_d["z_combined"], color="#4C78A8")
        ax.axvline(0, color="0.4", lw=0.6)
        ax.set_xlabel("Combined score")
        ax.set_title("Direction matches EMT, but none pass FDR < 0.05")
        fig.tight_layout()
        fig.savefig(OUT / "directional_but_not_significant.png", dpi=140)
        plt.close(fig)

    if not markers.empty and markers.shape[1] > 1:
        mat = markers.set_index("marker")
        fig_w = max(6.5, 0.46 * mat.shape[1] + 2)
        fig_h = max(4.5, 0.28 * mat.shape[0] + 1.5)
        fig, ax = plt.subplots(figsize=(fig_w, fig_h))
        limit = np.nanpercentile(np.abs(mat.to_numpy()), 98)
        limit = max(float(limit), 1e-3)
        image = ax.imshow(mat.to_numpy(), aspect="auto", cmap="RdBu_r", vmin=-limit, vmax=limit)
        ax.set_xticks(range(mat.shape[1]))
        ax.set_xticklabels(mat.columns, rotation=90)
        ax.set_yticks(range(mat.shape[0]))
        ax.set_yticklabels(mat.index)
        ax.set_title("Predicted log-expression change after virtual knockdown")
        fig.colorbar(image, ax=ax, fraction=0.03, pad=0.02, label="predicted delta")
        fig.tight_layout()
        fig.savefig(OUT / "top15_marker_shift_heatmap.png", dpi=140)
        plt.close(fig)

    _ = summary


if __name__ == "__main__":
    main()
