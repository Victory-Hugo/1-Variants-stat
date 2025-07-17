library(tidyverse)
library(tidyplots)

setwd('/mnt/c/Users/Administrator/Desktop/')
# 读取数据文件
df_variants <- read.csv('/mnt/d/幽门螺旋杆菌/Script/分析结果/2-变异统计/output/比较东亚和全球/Global.csv', sep = ',')

# 转为tibble格式
df_variants |> as_tibble() -> df_variants
df_variants
#* 先绘制总体的柱状图
df_variants|> 
   filter(Category != 'Type') |>
   filter(Class != '') |> # 过滤掉 Class 列中 NA 的行
   mutate(Cat = as.factor(Category)) |>
   tidyplot(x = Class, y = Count, fill = Class) |>
   add_sum_bar() |> # 添加总和柱状图
   add_sum_value(fontsize = 4) |> # 添加数值标签
   adjust_x_axis(rotate_labels = 90) |> # 调整 x 轴标签
   sort_x_axis_labels(by = Cat) |> # 按照 Category 排序
   save_plot('总体变异统计.pdf') 

#* 再绘制分组的柱状图
df_variants |> 
   filter(Category != 'Type') |> 
   mutate(Cat = as.factor(Category)) |>
   tidyplot(x = Class, y = Count, fill = Class) |>
   add_sum_bar() |> # 添加总和柱状图
   add_sum_value(fontsize = 4) |> # 添加数值标签
   adjust_x_axis(rotate_labels = 90) |> # 调整 x 轴标签
   sort_x_axis_labels(by = Cat) |> # 按照 Cat 排序
   split_plot(by = Source, ncol = 2,widths =50,heights = 30) |>
   save_plot('变异统计.pdf') 

#* 绘制SNP和INDEL的柱状图
df_variants |> 
   filter(Category == 'Type') |>
   mutate(Cat = as.factor(Category)) |>
   tidyplot(x = Class, y = Count, fill = Class) |>
   add_sum_bar() |> # 添加总和柱状图
   add_sum_value(fontsize = 4) |> # 添加数值标签
   sort_x_axis_labels(by = Cat) |> # 按照 Category 排序
   save_plot('总体变异统计_SNP_INDEL.pdf') 

#* 再绘制分组的柱状图
df_variants |> 
   filter(Category == 'Type') |> 
   mutate(Cat = as.factor(Category)) |>
   tidyplot(x = Class, y = Count, fill = Class) |>
   add_sum_bar() |> # 添加总和柱状图
   add_sum_value(fontsize = 4) |> # 添加数值标签
   adjust_x_axis(rotate_labels = 90) |> # 调整 x 轴标签
   sort_x_axis_labels(by = Cat) |> # 按照 Cat 排序
   split_plot(by = Source, ncol = 3,widths =50,heights = 30) |>
   save_plot('变异统计_SNP_INDEL.pdf') 
