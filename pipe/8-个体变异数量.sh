#!/usr/bin/env bash
# run_count_variants.sh

set -euo pipefail

# 输入参数
INPUT_DIR="/mnt/d/幽门螺旋杆菌/Script/分析结果/2-变异统计/conf/东亚低地和高地"   # VCF 文件所在目录
OUTPUT_DIR="/mnt/d/幽门螺旋杆菌/Script/分析结果/2-变异统计/output/东亚低地和高地"
PYTHON=/home/luolintao/miniconda3/envs/pyg/bin/python3
SCRIPT=/mnt/f/OneDrive/文档（科研）/脚本/Download/1-Variants-stat/python/8-个体变异数量.py

# 确保输出目录存在
mkdir -p "$OUTPUT_DIR"

# 找到所有 .vcf 或 .vcf.gz 文件，按文件名并行跑
find "$INPUT_DIR" -type f \( -name '*.vcf' -o -name '*.vcf.gz' \) | \
  parallel --eta '
    # 提取文件名，不带 .vcf 或 .vcf.gz 后缀
    f=$(basename {})
    base=${f%.vcf.gz}
    base=${base%.vcf}

    # 调用 Python 脚本
    '"$PYTHON"' '"$SCRIPT"' \
      --vcf {} \
      --out '"$OUTPUT_DIR"'/${base}_variants_per_genome.csv
  '
