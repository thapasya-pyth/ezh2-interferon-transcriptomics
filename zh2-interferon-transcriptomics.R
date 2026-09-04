library(gemma.R)
grep("dataset", ls("package:gemma.R"), value = TRUE)
gse <- "GSE183458.2"

data <- get_dataset_expression(gse)
metadata <- get_dataset_samples(gse)
dim(metadata)
names(metadata)
head(metadata)
metadata$sample.FactorValues
lapply(metadata$sample.FactorValues, unlist)
metadata$sample.FactorValues
lapply(metadata$sample.FactorValues, unlist)
lapply(metadata$sample.Characteristics, unlist)
metadata$condition <- ifelse(
  grepl("shCtrl", metadata$sample.Name),
  "shCtrl",
  "shEZH2"
)
dim(data)
colnames(data)
metadata$sample.Name
all(colnames(data) == metadata$sample.Name)
colnames(data)
metadata$sample.Name
metadata[, c("sample.Name", "sample.Accession")]

sample_cols <- 5:ncol(data)

data_gsm <- sub(".*___", "", colnames(data)[sample_cols])

data_gsm
data_gsm %in% metadata$sample.Accession
sample_cols <- 5:ncol(data)

data_gsm <- sub(".*___", "", colnames(data)[sample_cols])

match_order <- match(metadata$sample.Accession, data_gsm)

data_expr <- data[, c(1:4, sample_cols[match_order])]

sub(".*___", "", colnames(data_expr)[5:ncol(data_expr)])
metadata$sample.Accession
cols_to_keep <- c(1:4, sample_cols[match_order])

data_expr <- data[, cols_to_keep, with = FALSE]
dim(data_expr)
annotation <- data_expr[, 1:4, with = FALSE]

expr <- as.matrix(
  data_expr[, 5:ncol(data_expr), with = FALSE]
)
colnames(expr) <- metadata$sample.Name
dim(expr)
head(expr[, 1:3])
