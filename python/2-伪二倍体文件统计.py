#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
1-伪二倍体文件统计.py

统计伪装二倍体（真实单倍体）VCF.gz 中每个位点的变异：
- 按频率分类：Common/LowFreq/Rare/UltraRare
- 类型分类：SNV/Indel
- singleton/doubleton 统计

同时输出两份 CSV：
1) 汇总统计：Frequency/Type/Special 分类计数（含 Source 列）
2) 变异详情：CHROM, POS, REF, ALT, AC, Source
"""

import argparse
import sys
import csv
import os
from cyvcf2 import VCF

def main():
    parser = argparse.ArgumentParser(
        description="伪二倍体 VCF 统计脚本"
    )
    parser.add_argument("-i", "--vcf", required=True, help="输入 VCF.gz 文件")
    parser.add_argument("-o", "--out", required=True, help="输出统计结果 CSV")
    parser.add_argument("-v", "--var-out", required=True, help="输出变异详情 CSV")
    args = parser.parse_args()

    # 计算来源基础名
    base = os.path.basename(args.vcf)
    for ext in (".vcf.gz", ".vcf"):
        if base.endswith(ext):
            base = base[:-len(ext)]
            break

    try:
        vcf = VCF(args.vcf)
    except Exception as e:
        sys.exit(f"无法打开 VCF：{e}")

    # 初始化统计
    freq_counts    = {"Common":0, "LowFreq":0, "Rare":0, "UltraRare":0}
    type_counts    = {"SNV":0, "Indel":0}
    special_counts = {"Singleton":0, "Doubleton":0}

    # 打开变异详情输出文件
    try:
        var_f = open(args.var_out, "w", newline="", encoding="utf-8")
    except Exception as e:
        sys.exit(f"无法创建 {args.var_out}：{e}")
    var_writer = csv.writer(var_f)
    var_writer.writerow(["CHROM","POS","REF","ALT","AC","Source"])

    # 遍历每个位点
    for var in vcf:
        genos = var.genotypes
        n_alt = len(var.ALT)
        ac_list = [0] * n_alt
        an = 0

        for g in genos:
            a0, a1, _ = g
            # 缺失或杂合当缺失
            if a0 is None or a1 is None or a0 != a1:
                continue
            an += 1
            if a0 == 0:
                continue
            idx = a0 - 1
            if 0 <= idx < n_alt:
                ac_list[idx] += 1

        total_ac = sum(ac_list)
        if an == 0 or total_ac == 0:
            continue

        # 输出变异详情
        for alt, ac_val in zip(var.ALT, ac_list):
            if ac_val > 0:
                var_writer.writerow([var.CHROM, var.POS, var.REF, alt, ac_val, base])

        # singleton/doubleton 统计
        if total_ac == 1:
            special_counts["Singleton"] += 1
        elif total_ac == 2:
            special_counts["Doubleton"] += 1

        # 计算 MAF
        af = total_ac / an
        maf = af if af <= 0.5 else 1 - af

        # 频率分类
        if maf >= 0.05:
            freq_counts["Common"]    += 1
        elif maf >= 0.01:
            freq_counts["LowFreq"]   += 1
        elif maf >= 0.001:
            freq_counts["Rare"]      += 1
        else:
            freq_counts["UltraRare"] += 1

        # 类型分类
        if var.is_snp:
            type_counts["SNV"]   += 1
        else:
            type_counts["Indel"] += 1

    var_f.close()

    # 写入汇总统计 CSV（新增 Source 列）
    try:
        with open(args.out, "w", newline="", encoding="utf-8") as out_f:
            writer = csv.writer(out_f)
            writer.writerow(["Category","Class","Count","Source"])
            for cls, cnt in freq_counts.items():
                writer.writerow(["Frequency", cls, cnt, base])
            for cls, cnt in type_counts.items():
                writer.writerow(["Type", cls, cnt, base])
            for cls, cnt in special_counts.items():
                writer.writerow(["Special", cls, cnt, base])
    except Exception as e:
        sys.exit(f"无法写入 {args.out}：{e}")

    print(f"Done. 统计文件：{args.out}（含 Source 列）；变异详情：{args.var_out}")

if __name__ == "__main__":
    main()
