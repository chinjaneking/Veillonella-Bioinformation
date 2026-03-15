
library(ggplot2)
library(dplyr)

df <- read.csv("qc_report.csv")

# 基因组大小分布
p1 <- ggplot(df, aes(x=Genome_Size/1e6, fill=Pass_QC)) +
  geom_histogram(binwidth=0.1, alpha=0.7) +
  scale_fill_manual(values=c("TRUE"="#2E7D32", "FALSE"="#C62828"),
                    labels=c("Pass", "Fail")) +
  labs(x="Genome Size (Mb)", y="Count",
       title="Genome Size Distribution") +
  theme_minimal() +
  theme(legend.position="bottom")

ggsave("qc_genome_size.pdf", p1, width=6, height=4)

# GC含量分布
p2 <- ggplot(df, aes(x=GC_Content, fill=Pass_QC)) +
  geom_histogram(binwidth=1, alpha=0.7) +
  scale_fill_manual(values=c("TRUE"="#2E7D32", "FALSE"="#C62828")) +
  labs(x="GC Content (%)", y="Count",
       title="GC Content Distribution") +
  theme_minimal() +
  theme(legend.position="bottom")

ggsave("qc_gc_content.pdf", p2, width=6, height=4)

# N50 vs Contigs
p3 <- ggplot(df, aes(x=Contigs, y=N50/1e6, color=Pass_QC)) +
  geom_point(size=3, alpha=0.6) +
  scale_color_manual(values=c("TRUE"="#2E7D32", "FALSE"="#C62828"),
                     labels=c("Pass", "Fail")) +
  labs(x="Number of Contigs", y="N50 (Mb)",
       title="Assembly Quality") +
  theme_minimal() +
  theme(legend.position="bottom")

ggsave("qc_assembly.pdf", p3, width=6, height=4)

print("QC plots saved!")
