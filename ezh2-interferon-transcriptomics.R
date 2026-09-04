library(gemma.R)

# Dataset -----------------------------------------------------------------

gse <- "GSE183458.2"

# Download expression data and sample metadata
expression_data <- get_dataset_expression(gse)
metadata <- get_dataset_samples(gse)

# Define experimental condition
metadata$condition <- ifelse(
  grepl("shCtrl", metadata$sample.Name),
  "shCtrl",
  "shEZH2"
)


# Align expression data with metadata -------------------------------------

# First four columns contain gene annotation.
# Remaining columns contain expression values for the six samples.
sample_cols <- 5:ncol(expression_data)

# Extract GSM accessions from expression column names
expression_gsm <- sub(
  ".*___",
  "",
  colnames(expression_data)[sample_cols]
)

# Find the expression-column order corresponding to metadata order
match_order <- match(
  metadata$sample.Accession,
  expression_gsm
)

# Make sure every metadata sample was successfully matched
stopifnot(!anyNA(match_order))

# Reorder sample columns while retaining annotation columns
cols_to_keep <- c(1:4, sample_cols[match_order])

data_expr <- expression_data[
  ,
  cols_to_keep,
  with = FALSE
]

# Verify that expression samples and metadata are in identical order
stopifnot(
  identical(
    sub(
      ".*___",
      "",
      colnames(data_expr)[5:ncol(data_expr)]
    ),
    metadata$sample.Accession
  )
)

# Separate annotation and expression matrix -------------------------------

annotation <- data_expr[
  ,
  1:4,
  with = FALSE
]

expr <- as.matrix(
  data_expr[
    ,
    5:ncol(data_expr),
    with = FALSE
  ]
)

# Use readable sample names for downstream analysis
colnames(expr) <- metadata$sample.Name

# Final consistency check
stopifnot(
  ncol(expr) == nrow(metadata)
)

cor_expr <- cor(expr)

round(cor_expr, 3)
library(limma)

plotMDS(expr)

heatmap(
  cor_expr,
  scale = "none",
  symm = TRUE,
  margins = c(10, 10)
)
metadata$condition <- factor(
  metadata$condition,
  levels = c("shCtrl", "shEZH2")
)
design <- model.matrix(
  ~ 0 + condition,
  data = metadata
)
colnames(design) <- c("shCtrl", "shEZH2")
design
contrast_matrix <- makeContrasts(
  shEZH2_vs_shCtrl = shEZH2 - shCtrl,
  levels = design
)
contrast_matrix
fit <- lmFit(expr, design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

zero_var <- apply(expr, 1, var) == 0

sum(zero_var)
mean(zero_var)

expr_filt <- expr[!zero_var, ]
annotation_filt <- annotation[!zero_var, ]
dim(expr)
dim(expr_filt)
fit <- lmFit(expr_filt, design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

sum(fit$sigma == 0)
mean(fit$sigma == 0)

results <- topTable(
  fit2,
  coef = "shEZH2_vs_shCtrl",
  number = Inf,
  adjust.method = "BH"
)
head(results)
dim(results)
results$GeneSymbol <- annotation_filt$GeneSymbol[
  match(rownames(results), rownames(annotation_filt))
]

results$GeneName <- annotation_filt$GeneName[
  match(rownames(results), rownames(annotation_filt))
]

head(results[, c("GeneSymbol", "GeneName", "logFC", "P.Value", "adj.P.Val")])
results$significance <- "Not significant"

results$significance[
  results$adj.P.Val < 0.05 & results$logFC >= 1
] <- "Up"

results$significance[
  results$adj.P.Val < 0.05 & results$logFC <= -1
] <- "Down"
table(results$significance)

results$negLog10FDR <- -log10(results$adj.P.Val)
### volcano plot

library(ggplot2)

ggplot(
  results,
  aes(
    x = logFC,
    y = negLog10FDR,
    color = significance
  )
) +
  geom_point(alpha = 0.6, size = 1.5) +
  geom_vline(
    xintercept = c(-1, 1),
    linetype = "dashed"
  ) +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed"
  ) +
  labs(
    title = "Differential expression after EZH2 knockdown",
    x = "log2 fold change",
    y = "-log10 adjusted p-value"
  ) +
  theme_minimal()
results$significance <- factor(
  results$significance,
  levels = c("Down", "Not significant", "Up")
)

volcano_plot <- ggplot(
  results,
  aes(
    x = logFC,
    y = negLog10FDR,
    color = significance
  )
) +
  geom_point(alpha = 0.6, size = 1.4) +
  geom_vline(
    xintercept = c(-1, 1),
    linetype = "dashed"
  ) +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed"
  ) +
  labs(
    title = "Differential expression after EZH2 knockdown in MCF7 cells",
    x = "log2 fold change",
    y = "-log10 adjusted p-value",
    color = "Significance"
  ) +
  theme_minimal(base_size = 12)

volcano_plot

dir.create("figures", showWarnings = FALSE)

ggsave(
  "figures/volcano_plot.png",
  volcano_plot,
  width = 8,
  height = 6,
  dpi = 300
)
head(
  results[order(results$adj.P.Val), 
          c("GeneSymbol", "GeneName", "logFC", "adj.P.Val")],
  20
)

genes_to_check <- c(
  "CERS6-AS1",
  "CGA",
  "S100A8",
  "CD36",
  "CCN5",
  "LCN2",
  "S100A9"
)

expr[
  annotation$GeneSymbol %in% genes_to_check,
]
check_idx <- which(
  annotation$GeneSymbol %in% genes_to_check
)

cbind(
  GeneSymbol = annotation$GeneSymbol[check_idx],
  expr[check_idx, ]
)
up_genes <- results$GeneSymbol[
  results$significance == "Up" &
    !is.na(results$GeneSymbol) &
    results$GeneSymbol != ""
]

down_genes <- results$GeneSymbol[
  results$significance == "Down" &
    !is.na(results$GeneSymbol) &
    results$GeneSymbol != ""
]
up_genes <- unique(up_genes)
down_genes <- unique(down_genes)
length(up_genes)
length(down_genes)
head(up_genes)
head(down_genes)
background_genes <- results$GeneSymbol[
  !is.na(results$GeneSymbol) &
    results$GeneSymbol != ""
]

background_genes <- unique(background_genes)
length(background_genes)

#BiocManager::install("clusterProfiler", update = FALSE, ask = FALSE)
#BiocManager::install("org.Hs.eg.db", update = FALSE, ask = FALSE)

dir.create("results", showWarnings = FALSE)
writeLines(
  up_genes,
  "results/upregulated_genes.txt"
)

writeLines(
  down_genes,
  "results/downregulated_genes.txt"
)
writeLines(
  background_genes,
  "results/background_genes.txt"
)

gprof_up <- read.csv(
  "./Up_gProfiler_hsapiens_2026-09-04_18-25-12.csv"
)

names(gprof_up)


gprof_up_top <- gprof_up[
  order(gprof_up$adjusted_p_value),
][1:10, ]

gprof_up_top$term_name <- factor(
  gprof_up_top$term_name,
  levels = rev(gprof_up_top$term_name)
)


enrichment_plot <- ggplot(
  gprof_up_top,
  aes(
    x = negative_log10_of_adjusted_p_value,
    y = term_name,
    size = intersection_size
  )
) +
  geom_point() +
  labs(
    title = "Top enriched pathways among upregulated genes",
    x = "-log10 adjusted p-value",
    y = NULL,
    size = "Gene count"
  ) +
  theme_minimal(base_size = 12)

enrichment_plot

ggsave(
  "figures/enrichment_upregulated.png",
  enrichment_plot,
  width = 8,
  height = 6,
  dpi = 300
)
###########

gprof_down <- read.csv(
  "./Down_gProfiler_hsapiens_2026-09-04_18-28-13.csv"
)

names(gprof_down)

gprof_down_top <- gprof_down[
  order(gprof_down$adjusted_p_value),
][1:10, ]

gprof_down_top$term_name <- factor(
  gprof_down_top$term_name,
  levels = rev(gprof_down_top$term_name)
)

down_enrichment_plot <- ggplot(
  gprof_down_top,
  aes(
    x = negative_log10_of_adjusted_p_value,
    y = term_name,
    size = intersection_size
  )
) +
  geom_point() +
  labs(
    title = "Top enriched pathways among downregulated genes",
    x = "-log10 adjusted p-value",
    y = NULL,
    size = "Gene count"
  ) +
  theme_minimal(base_size = 12)

down_enrichment_plot


keep_terms <- c(
  "Interferon Signaling",
  "Interferon alpha/beta signaling",
  "response to virus",
  "innate immune response",
  "defense response",
  "response to biotic stimulus"
)
gprof_up_final <- gprof_up[
  gprof_up$term_name %in% keep_terms,
]
gprof_up_final$term_name <- factor(
  gprof_up_final$term_name,
  levels = gprof_up_final$term_name[
    order(gprof_up_final$adjusted_p_value, decreasing = TRUE)
  ]
)

final_enrichment_plot <- ggplot(
  gprof_up_final,
  aes(
    x = negative_log10_of_adjusted_p_value,
    y = term_name,
    size = intersection_size
  )
) +
  geom_point() +
  labs(
    title = "Interferon and antiviral pathways enriched after EZH2 knockdown",
    x = "-log10 adjusted p-value",
    y = NULL,
    size = "Gene count"
  ) +
  theme_minimal(base_size = 12)

final_enrichment_plot
ggsave(
  "figures/enrichment_upregulated_final.png",
  final_enrichment_plot,
  width = 8,
  height = 5,
  dpi = 300
)


genes_plot <- c("S100A8", "S100A9", "LCN2")

plot_idx <- which(
  annotation$GeneSymbol %in% genes_plot
)

gene_plot_data <- data.frame(
  GeneSymbol = rep(
    annotation$GeneSymbol[plot_idx],
    each = ncol(expr)
  ),
  Sample = rep(
    colnames(expr),
    times = length(plot_idx)
  ),
  Expression = as.vector(
    t(expr[plot_idx, ])
  )
)

gene_plot_data$condition <- metadata$condition[
  match(gene_plot_data$Sample, metadata$sample.Name)
]


head(gene_plot_data, 12)

###

gene_expression_plot <- ggplot(
  gene_plot_data,
  aes(
    x = condition,
    y = Expression,
    color = condition
  )
) +
  geom_jitter(
    width = 0.08,
    size = 2.5
  ) +
  facet_wrap(
    ~ GeneSymbol,
    scales = "free_y"
  ) +
  labs(
    title = "Expression of selected genes after EZH2 knockdown",
    x = NULL,
    y = "Expression"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none"
  )

gene_expression_plot

ggsave(
  "figures/gene_expression_selected.png",
  gene_expression_plot,
  width = 10,
  height = 6,
  dpi = 300
)
write.csv(
  results,
  "results/differential_expression_results.csv",
  row.names = FALSE
)

list.files("figures")
list.files("results")
png(
  "figures/mds_plot.png",
  width = 1200,
  height = 900,
  res = 150
)

plotMDS(expr)

dev.off()
png(
  "figures/correlation_heatmap.png",
  width = 1200,
  height = 1000,
  height = 1000,
  res = 150
)

heatmap(
  cor_expr,
  scale = "none",
  symm = TRUE,
  margins = c(10, 10)
)

dev.off()
