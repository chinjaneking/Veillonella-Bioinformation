library(ggplot2)
library(dplyr)
library(tidyr)

# 读取数据
df <- read.csv("pangenome_curve_data.csv")

# 转换为长格式用于ggplot
df_long <- df %>%
  select(sample_size, pan_mean, core_mean) %>%
  pivot_longer(cols = c(pan_mean, core_mean),
               names_to = "Type", values_to = "Genes")

# 添加标准差
df_long$std <- c(df$pan_std, df$core_std)

# 设置因子顺序
df_long$Type <- factor(df_long$Type,
                       levels = c("pan_mean", "core_mean"),
                       labels = c("Pan-genome", "Core genome"))

# 绘制泛基因组曲线
p <- ggplot(df_long, aes(x=sample_size, y=Genes, color=Type)) +
  geom_line(linewidth=1.2) +
  geom_ribbon(aes(ymin=Genes-std, ymax=Genes+std, fill=Type),
              alpha=0.2, color=NA) +
  scale_color_manual(values=c("Pan-genome"="#1976D2", "Core genome"="#388E3C")) +
  scale_fill_manual(values=c("Pan-genome"="#1976D2", "Core genome"="#388E3C")) +
  labs(x="Number of Genomes",
       y="Number of Gene Families",
       title="Pan-genome Analysis") +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    plot.title = element_text(hjust=0.5, face="bold", size=14),
    axis.title = element_text(size=12),
    panel.grid.minor = element_blank()
  )

ggsave("Figure_pangenome_curves_ggplot.pdf", p, width=8, height=6, dpi=300)

print("Pan-genome plot saved!")
