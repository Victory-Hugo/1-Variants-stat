#!/usr/bin/env bash
set -euo pipefail

# 并行任务数，根据 CPU 核心数和 I/O 性能调节
JOBS=8

PYTHON="/home/luolintao/miniconda3/envs/pyg/bin/python3"
SCRIPT="/mnt/f/OneDrive/文档（科研）/脚本/Download/1-Variants-stat/python/2-伪二倍体文件统计.py"
DATA_DIR="/mnt/d/幽门螺旋杆菌/Script/分析结果/2-变异统计/conf/比较东亚和全球/"
OUT_DIR="/mnt/d/幽门螺旋杆菌/Script/分析结果/2-变异统计/output/比较东亚和全球/"

# 确保输出目录存在
mkdir -p "$OUT_DIR"

export PYTHON SCRIPT OUT_DIR

# 并行执行
find "$DATA_DIR" -maxdepth 1 -type f -name '*.vcf.gz' | \
parallel -j $JOBS '
  infile={}
  base=$(basename "$infile" .vcf.gz)
  "$PYTHON" "$SCRIPT" \
    --vcf "$infile" \
    --out  "$OUT_DIR/${base}.csv" \
    --var-out "$OUT_DIR/${base}.var.csv"
'

echo "All done."
