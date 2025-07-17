#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
common_intersect_with_ids_fast.py - 按 SNV/Indel 分别输出（高速向量化版本）

在比较 East_Asia.var.csv 与 Global.var.csv 中均为 Common 的变异后，
为 SNV 和 Indel 两种变异分别：
- 去重并排序赋唯一编号（var1, var2, …）
- 输出两列 CSV（East_Asia, Global）文件，用于后续韦恩图绘制。

直接使用 var.csv 中已存在的 Freq 和 Type 列，无需重新计算，完全向量化处理。
"""

import os
import numpy as np
import pandas as pd

def main():
    # 1. 读取 CSV，只取必要列
    usecols = ["CHROM", "POS", "REF", "ALT", "Freq", "Type"]
    ea = pd.read_csv(
        "/mnt/d/幽门螺旋杆菌/Script/分析结果/2-变异统计/output/比较东亚和全球/East_Asia.var.csv",
        usecols=usecols, encoding="utf-8"
    )
    gl = pd.read_csv(
        "/mnt/d/幽门螺旋杆菌/Script/分析结果/2-变异统计/output/比较东亚和全球/Global.var.csv",
        usecols=usecols, encoding="utf-8"
    )

    # 2. 筛选出 Common
    ea_c = ea[ea["Freq"] == "Common"].copy()
    gl_c = gl[gl["Freq"] == "Common"].copy()

    # 3. 构造 key 列
    for df in (ea_c, gl_c):
        df["key"] = (
            df["CHROM"].astype(str) + "_" +
            df["POS"].astype(str) + "_" +
            df["REF"] + "_" +
            df["ALT"]
        )

    # 4. 分类型处理：SNV 和 Indel
    for vtype in ["SNV", "Indel"]:
        # 提取各自 Type 对应的唯一 key
        ea_keys = ea_c.loc[ea_c["Type"] == vtype, "key"].unique()
        gl_keys = gl_c.loc[gl_c["Type"] == vtype, "key"].unique()

        # 计算并集并排序（np.union1d 会自动去重并返回排序后的结果）
        union_keys = np.union1d(ea_keys, gl_keys)

        # 生成 varID，如 var1, var2, …
        ids = np.arange(1, len(union_keys) + 1)
        ids_str = ["var{}".format(i) for i in ids]

        # 构建包含所有并集 key 的 DataFrame，并一次性填充 East_Asia / Global 列
        df_union = pd.DataFrame({
            "key": union_keys,
            "ID": ids_str
        })
        df_union["East_Asia"] = np.where(
            df_union["key"].isin(ea_keys),
            df_union["ID"],
            ""
        )
        df_union["Global"] = np.where(
            df_union["key"].isin(gl_keys),
            df_union["ID"],
            ""
        )

        # 只保留两列输出
        out_df = df_union[["East_Asia", "Global"]]

        # 写文件
        out_path = os.path.join(
            "/mnt/d/幽门螺旋杆菌/Script/分析结果/2-变异统计/",
            "output/比较东亚和全球",
            f"Common_Variants_{vtype}_IDs.csv"
        )
        out_df.to_csv(out_path, index=False, encoding="utf-8")
        print(f"Done. 已生成 {out_path}")

if __name__ == "__main__":
    main()
