#!/usr/bin/env python3
"""去除各类免疫球蛋白后的血清蛋白丰度排名气泡图。不改原脚本出图。"""
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from serum_proteomics_ranked_bubble import run_pipeline


if __name__ == "__main__":
    run_pipeline(
        result_subdir="serum_proteomics_bubble_no_Ig",
        drop_immunoglobulin=True,
        title_prefix="去除免疫球蛋白后 血清蛋白丰度排名气泡图",
        subtitle="去掉 IGH/IGK/IGL/JCHAIN 等后，两样品平均丰度排名；气泡大小一致",
    )
