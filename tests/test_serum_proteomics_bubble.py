#!/usr/bin/env python3
"""血清蛋白质组排名气泡图：预处理、分组与回退路径。"""
from __future__ import annotations

import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import numpy as np
import pandas as pd

import serum_proteomics_ranked_bubble as sp


def test_pick_official_symbol() -> None:
    assert sp.pick_official_symbol("SAA2,SAA2-SAA4,SAA4") == "SAA2"
    assert sp.pick_official_symbol("ALB") == "ALB"
    assert sp.pick_official_symbol("-") is None


def test_annotation_requires_two_groups(tmp_path: Path) -> None:
    (tmp_path / "sample_annotation.csv").write_text(
        "sample,group\nA1,GroupA\nA2,GroupA\n", encoding="utf-8"
    )
    try:
        sp.load_annotation(tmp_path, ["A1", "A2"])
    except SystemExit as exc:
        assert "恰好两组" in str(exc)
    else:
        raise AssertionError("expected SystemExit for a single group")


def test_pg_matrix_fallback_to_pr(tmp_path: Path) -> None:
    pg = tmp_path / "report.pg_matrix"
    pg.write_text(
        "Protein.Group\tGenes\tS1\tS2\nP1\tALB\t100\t110\nP2\tAPOA1\t80\t90\n",
        encoding="utf-8",
    )
    pr = tmp_path / "report.pr_matrix"
    pr.write_text(
        "Protein.Group\tGenes\tStripped.Sequence\tS1\tS2\n"
        "P10\tHP\tAAA\t50\t60\n"
        "P10\tHP\tBBB\t70\t80\n"
        "P11\tTF\tCCC\t20\t30\n"
        "P12\tCRP\tDDD\t15\t25\n"
        "P13\tC3\tEEE\t40\t45\n"
        "P14\tTTR\tFFF\t10\t12\n",
        encoding="utf-8",
    )
    log_path = tmp_path / "log.txt"
    log_path.write_text("", encoding="utf-8")
    mat, cols = sp.load_protein_matrix(tmp_path, log_path)
    assert cols == ["S1", "S2"]
    assert len(mat) >= 5
    assert "P10" in set(mat["Protein.Group"].astype(str))
    hp = mat.loc[mat["Protein.Group"] == "P10", ["S1", "S2"]].iloc[0]
    assert float(hp["S1"]) == 60.0  # median of 50 and 70


def test_rank_by_abundance_not_raw_fc() -> None:
    meta = pd.DataFrame(
        {
            "Protein.Group": ["P1", "P2", "P3"],
            "Genes": ["HIGH", "MID", "LOW"],
            "First.Protein.Description": ["h", "m", "l"],
        }
    )
    # columns: A1 A2 B1 B2 ; HIGH has largest mean, LOW has largest |FC|
    log2x = np.array(
        [
            [10.0, 10.2, 10.1, 10.3],
            [8.0, 8.1, 8.4, 8.5],
            [4.0, 4.1, 7.0, 7.2],
        ]
    )
    sample_info = pd.DataFrame(
        {
            "column": ["A1", "A2", "B1", "B2"],
            "sample": ["A1", "A2", "B1", "B2"],
            "group": pd.Categorical(["GroupA", "GroupA", "GroupB", "GroupB"], ordered=True),
        }
    )
    ranked = sp.rank_proteins(meta, log2x, sample_info)
    assert list(ranked.sort_values("abundance_rank")["gene"]) == ["HIGH", "MID", "LOW"]
    assert ranked.loc[ranked["gene"] == "LOW", "fc_rank"].iloc[0] == 1
    assert ranked.loc[ranked["gene"] == "HIGH", "abundance_rank"].iloc[0] == 1


def test_preprocess_log2_and_filter() -> None:
    mat = pd.DataFrame(
        {
            "Protein.Group": ["P1", "P2"],
            "Genes": ["ALB", "ZZZ"],
            "A1": [100.0, np.nan],
            "A2": [110.0, np.nan],
            "B1": [90.0, np.nan],
            "B2": [95.0, np.nan],
        }
    )
    sample_info = pd.DataFrame(
        {
            "column": ["A1", "A2", "B1", "B2"],
            "sample": ["A1", "A2", "B1", "B2"],
            "group": pd.Categorical(["GroupA", "GroupA", "GroupB", "GroupB"], ordered=True),
        }
    )
    kept, log2x = sp.preprocess(mat, ["A1", "A2", "B1", "B2"], sample_info)
    assert list(kept["Genes"]) == ["ALB"]
    assert log2x.shape == (1, 4)
    # log2(x+1) then median-centered; all finite
    assert np.isfinite(log2x).all()


if __name__ == "__main__":
    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td)
        test_pick_official_symbol()
        test_annotation_requires_two_groups(tmp)
        test_pg_matrix_fallback_to_pr(tmp)
        test_rank_by_abundance_not_raw_fc()
        test_preprocess_log2_and_filter()
    print("all tests passed")
