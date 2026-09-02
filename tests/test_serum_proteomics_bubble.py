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


def test_immunoglobulin_filter() -> None:
    assert sp.is_immunoglobulin_symbol("IGHG1")
    assert sp.is_immunoglobulin_symbol("IGKC")
    assert sp.is_immunoglobulin_symbol("IGLV3-21")
    assert sp.is_immunoglobulin_symbol("JCHAIN")
    assert not sp.is_immunoglobulin_symbol("IGSF8")
    assert not sp.is_immunoglobulin_symbol("IGF1")
    assert not sp.is_immunoglobulin_symbol("ALB")
    meta = pd.DataFrame(
        {
            "Genes": ["ALB", "IGHG1", "IGSF8"],
            "Protein.Names": [
                "Serum albumin",
                "Immunoglobulin heavy constant gamma 1",
                "Immunoglobulin superfamily member 8",
            ],
            "First.Protein.Description": ["Serum albumin", "Ig gamma-1", "IgSF member"],
            "Protein.Group": ["P02768", "P01857", "Q969P0"],
        }
    )
    assert list(sp.immunoglobulin_mask(meta)) == [False, True, False]


def test_pick_official_symbol() -> None:
    assert sp.pick_official_symbol("SAA2,SAA2-SAA4,SAA4") == "SAA2"
    assert sp.pick_official_symbol("ALB") == "ALB"
    assert sp.pick_official_symbol("-") is None


def test_missing_annotation_two_samples_autowrite(tmp_path: Path) -> None:
    info = sp.load_annotation(tmp_path, ["GP_WJZ_11", "GP_WJZ_18"])
    csv_path = tmp_path / "sample_annotation.csv"
    assert csv_path.is_file()
    assert list(info["sample"]) == ["GP_WJZ_11", "GP_WJZ_18"]
    assert list(info["group"].cat.categories) == ["GP_WJZ_11", "GP_WJZ_18"]


def test_missing_annotation_three_samples_writes_template(tmp_path: Path) -> None:
    try:
        sp.load_annotation(tmp_path, ["A1", "A2", "B1"])
    except SystemExit as exc:
        assert "模板" in str(exc)
    else:
        raise AssertionError("expected SystemExit when more than two samples lack groups")
    text = (tmp_path / "sample_annotation.csv").read_text(encoding="utf-8-sig")
    assert "A1" in text and "B1" in text


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


def test_rank_by_mean_of_two_samples() -> None:
    meta = pd.DataFrame(
        {
            "Protein.Group": ["P1", "P2", "P3"],
            "Genes": ["HIGH", "MID", "LOW"],
            "First.Protein.Description": ["h", "m", "l"],
        }
    )
    log2x = np.array(
        [
            [10.0, 10.2],
            [8.0, 8.4],
            [4.0, 7.0],
        ]
    )
    ranked = sp.rank_proteins(meta, log2x, ["GP_WJZ_11", "GP_WJZ_18"])
    assert list(ranked.sort_values("abundance_rank")["gene"]) == ["HIGH", "MID", "LOW"]
    assert ranked.loc[ranked["gene"] == "HIGH", "abundance_rank"].iloc[0] == 1
    np.testing.assert_allclose(
        ranked.loc[ranked["gene"] == "MID", "mean_abundance"].iloc[0], 8.2
    )


def test_preprocess_log2_and_filter() -> None:
    mat = pd.DataFrame(
        {
            "Protein.Group": ["P1", "P2"],
            "Genes": ["ALB", "ZZZ"],
            "GP_WJZ_11": [100.0, np.nan],
            "GP_WJZ_18": [110.0, np.nan],
        }
    )
    kept, log2x = sp.preprocess(mat, ["GP_WJZ_11", "GP_WJZ_18"])
    assert list(kept["Genes"]) == ["ALB"]
    assert log2x.shape == (1, 2)
    assert np.isfinite(log2x).all()


def test_default_data_dir_is_zhao_serum() -> None:
    joined = " ".join(str(p) for p in sp.DEFAULT_PROTEOMICS_DIRS)
    assert r"天府" in joined or "天府" in joined
    assert "赵章寻" in joined
    assert "血清蛋白质组学" in joined
    assert "TG_BRCA" not in joined
    launcher = ROOT / "run_serum_proteomics_bubble.R"
    text = launcher.read_text(encoding="utf-8")
    assert "serum_proteomics_ranked_bubble.R" in text
    assert "天府" in text
    assert 'encoding = "UTF-8"' in text
    no_ig = (ROOT / "serum_proteomics_ranked_bubble_no_Ig.R").read_text(encoding="utf-8")
    assert "run_serum_abundance_bubble" not in no_ig
    assert "immunoglobulin" in no_ig.lower() or "免疫球蛋白" in no_ig
    assert "serum_proteomics_bubble_no_Ig" in no_ig


if __name__ == "__main__":
    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td)
        test_pick_official_symbol()
        test_immunoglobulin_filter()
        (tmp / "onegroup").mkdir()
        test_annotation_requires_two_groups(tmp / "onegroup")
        (tmp / "two").mkdir()
        test_missing_annotation_two_samples_autowrite(tmp / "two")
        (tmp / "three").mkdir()
        test_missing_annotation_three_samples_writes_template(tmp / "three")
        test_pg_matrix_fallback_to_pr(tmp)
        test_rank_by_mean_of_two_samples()
        test_preprocess_log2_and_filter()
        test_default_data_dir_is_zhao_serum()
    print("all tests passed")
