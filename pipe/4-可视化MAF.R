library(tidyverse)
library(tidyplots)
setwd('/mnt/c/Users/Administrator/Desktop/')

df_MAF <- read.csv('/mnt/f/OneDrive/文档（科研）/脚本/Download/1-Variants-stat/output/Bin_MAF_Comparison.csv', sep = ',')
df_MAF |> as_tibble() -> df_MAF
# 将宽数据转为长数据
df_MAF |> 
  pivot_longer(cols = -MAF, names_to = 'Source', values_to = 'global_count') |>
  filter(Source != 'global_count') -> df_long 
# 将MAF独特的种类转为一个向量
df_long |> select(MAF) |> distinct() |> pull(MAF) -> MAF
MAF
# MAF从小到大排序 - 自定义顺序而不是使用sort
# 定义MAF的正确顺序
maf_order <- c('<0.1', '0.1-0.3', '0.3-0.5', '0.5-0.7', '0.7-1.0',
               '1-2', '2-5', '5-10', '10-20', '20-30', '30-40', '40-50')

df_long |> 
  mutate(MAF = factor(MAF, levels = maf_order)) -> df_long
print(df_long,n = 50)

df_long  |> 
  tidyplot(x = MAF, y = global_count, fill = Source) |>
  add_barstack_absolute()  |>
  add_data_labels(label = global_count, position = "stack", size = 3) |>
  adjust_x_axis(rotate_labels = 90) |>
  save_plot("MAF分布.pdf")
