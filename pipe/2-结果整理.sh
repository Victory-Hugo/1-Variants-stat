#!/usr/bin/env bash
# 1. 定义输入/输出目录
INPUT_DIR="/mnt/d/幽门螺旋杆菌/Script/分析结果/2-变异统计/output/比较东亚和全球/var/"
OUTPUT_DIR="/mnt/d/幽门螺旋杆菌/Script/分析结果/2-变异统计/output/比较东亚和全球/csv/"

# 2. 确保输出目录存在
mkdir -p "${OUTPUT_DIR}"

# 3. 并行处理
#    - find 出所有 .var.csv
#    - parallel 每个文件执行一次脚本
find "${INPUT_DIR}" -type f -name '*.var.csv' | \
  parallel --bar \
  /home/luolintao/miniconda3/envs/pyg/bin/python3 \
  /mnt/f/OneDrive/文档（科研）/脚本/Download/1-Variants-stat/python/3-统计各项数量.py {} "${OUTPUT_DIR}"




# BASE_DIR='/mnt/d/幽门螺旋杆菌/Script/分析结果/2-变异统计/output/比较东亚和全球/'
# VAR_DIR="${BASE_DIR}/var/"
# CSV_DIR="${BASE_DIR}/csv/"

# mkdir -p "$VAR_DIR" "$CSV_DIR"

# for file in ${BASE_DIR}/*.var.csv;do
#  mv "$file" "${VAR_DIR}/";done

# for file in ${BASE_DIR}/*.csv;do
#  mv "${file}" "${CSV_DIR}/";done

# /home/luolintao/miniconda3/envs/pyg/bin/python3 \
#     /mnt/f/OneDrive/文档（科研）/脚本/Download/1-Variants-stat/python/4-结果整理.py \
#     --merge-dir ${CSV_DIR}/ \
#     --merge-out ${CSV_DIR}/merged.csv \
#     --var-dir ${VAR_DIR}/ \
#     --out ${VAR_DIR}/merged_all_sources.csv


# /home/luolintao/miniconda3/envs/pyg/bin/python3 \
#   /mnt/f/OneDrive/文档（科研）/脚本/Download/1-Variants-stat/python/5-韦恩数据.py \
#   --input  ${VAR_DIR}/merged_all_sources.csv \
#   --output ${VAR_DIR}/venn_counts.csv

