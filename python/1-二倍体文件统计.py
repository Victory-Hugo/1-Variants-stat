#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
1-二倍体文件统计.py

统计 VCF.gz 中每个位点的变异，按频率分类（Common/LowFreq/Rare/UltraRare）、类型分类（SNV/Indel）、
singleton/doubleton 统计，并将统计结果输出为一个 CSV 文件，同时将所有变异的位置信息、具体变异及其等位基因计数（AC）、
以及来源文件基础名写入另一个 CSV 文件，以便后续不同大陆 VCF 合并和韦恩图分析。
"""

import argparse
import sys
import csv
import os
from cyvcf2 import VCF

def main():
    parser = argparse.ArgumentParser(
        description="统计 VCF 文件中各类变异并输出统计 CSV 及变异详情 CSV（含 AC 和来源基础名）"
    )
    parser.add_argument(
        "--vcf", "-i", required=True,
        help="输入 VCF.gz 文件路径"
    )
    parser.add_argument(
        "--out", "-o", required=True,
        help="输出统计结果 CSV 文件路径"
    )
    parser.add_argument(
        "--var-out", "-v", required=True,
        help="输出变异详情 CSV 文件路径（含 AC 和来源基础名）"
    )
    args = parser.parse_args()

    # 计算来源文件基础名
    filename = os.path.basename(args.vcf)
    base = filename
    for ext in (".vcf.gz", ".vcf"):
        if base.endswith(ext):
            base = base[:-len(ext)]
            break

    # 打开 VCF
    try:
        vcf = VCF(args.vcf)
    except Exception as e:
        sys.exit(f"无法打开 VCF 文件：{e}")

    # 初始化统计字典
    freq_counts    = {"Common":0, "LowFreq":0, "Rare":0, "UltraRare":0}
    type_counts    = {"SNV":0, "Indel":0}
    special_counts = {"Singleton":0, "Doubleton":0}

    # 打开变异详情输出文件
    try:
        var_file = open(args.var_out, "w", newline="", encoding="utf-8")
    except Exception as e:
        sys.exit(f"无法创建变异详情文件：{e}")

    var_writer = csv.writer(var_file)
    # 增加一列 Source
    var_writer.writerow(["CHROM", "POS", "REF", "ALT", "AC", "Source"])

    # 遍历 VCF 条目，一次完成统计与详情输出
    for var in vcf:
        ac_info = var.INFO.get("AC")
        an = var.INFO.get("AN")
        # 跳过 AN=0 或 AC 全为0 的记录
        if an is None or an == 0 or ac_info is None:
            continue
        ac_list = ac_info if isinstance(ac_info, list) else [ac_info]
        if all(ac == 0 for ac in ac_list):
            continue

        # 写入每个 ALT 的变异详情和对应 AC 及来源基础名
        for alt, ac_val in zip(var.ALT, ac_list):
            var_writer.writerow([var.CHROM, var.POS, var.REF, alt, ac_val, base])

        # 统计 singleton/doubleton（基于总 AC）
        total_ac = sum(ac_list)
        if total_ac == 1:
            special_counts["Singleton"] += 1
        elif total_ac == 2:
            special_counts["Doubleton"] += 1

        # 计算次要等位基因频率
        af = total_ac / an
        maf = min(af, 1 - af)

        # 频率分类
        if maf >= 0.05:
            freq_counts["Common"]    += 1
        elif maf >= 0.01:
            freq_counts["LowFreq"]   += 1
        elif maf >= 0.001:
            freq_counts["Rare"]      += 1
        else:
            freq_counts["UltraRare"] += 1

        # 类型分类：SNV vs Indel
        if var.is_snp:
            type_counts["SNV"]   += 1
        else:
            type_counts["Indel"] += 1

    var_file.close()

    # 将统计结果写入 CSV
    try:
        with open(args.out, "w", newline="", encoding="utf-8") as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(["Category", "Class", "Count"])
            for cls, cnt in freq_counts.items():
                writer.writerow(["Frequency", cls, cnt])
            for cls, cnt in type_counts.items():
                writer.writerow(["Type", cls, cnt])
            for cls, cnt in special_counts.items():
                writer.writerow(["Special", cls, cnt])
    except Exception as e:
        sys.exit(f"无法写入统计结果文件：{e}")

    print(f"统计完成，统计结果已写入：{args.out}；变异详情（含 Source）已写入：{args.var_out}")

if __name__ == "__main__":
    main()
