library(tidyverse)
library(tidyplots)
setwd('/mnt/c/Users/Administrator/Desktop')
variants_id_dir <- '/mnt/d/幽门螺旋杆菌/Script/分析结果/2-变异统计/output/东亚低地和高地/'
# 从变异文件目录下读取所有的csv文件
# 将其们合并为一个数据框
variants_id_files <- list.files(variants_id_dir, pattern = '\\.csv$', full.names = TRUE)
# 这里有很多csv文件，读取它们，只保留一次表头
variants_id_files <- variants_id_files[file.exists(variants_id_files)]
variants_id_list <- lapply(variants_id_files, function(file) {
  read.csv(file, header = TRUE, stringsAsFactors = FALSE)
})
# 合并所有数据框
variants_id_df <- do.call(rbind, variants_id_list)
variants_id_df |> as_tibble() -> df_variants_id
# 新增一列，将 variant_count统一除以1000
df_variants_id |>
   mutate(variant_count = variant_count / 1000) -> df_variants_id

df_variants_id |> 
    tidyplot(x = Source, y = variant_count, color = Source) |>
    add_boxplot(outlier.alpha = 0.25) |>
    add_median_value(hjust = 0.5, vjust = -3) |>  # 添加均值标签
    add_test_pvalue() |>
    adjust_y_axis_title('Variants per genome (K)') |>
    adjust_x_axis_title('Continent') |>
    adjust_x_axis(rotate_labels = 90) |>
    sort_x_axis_labels(.reverse = TRUE) |>
    save_plot('Boxplot.pdf')

df_variants_id |> select(Source) |> distinct() |> pull() 
# 仅保留North_Asia，South_Asia ，Southeast_Asia，
