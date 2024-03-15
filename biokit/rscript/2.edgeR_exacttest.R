library(edgeR)
library(data.table)

args <- commandArgs(trailingOnly = TRUE)
exp_tpm_file = args[1]
group_file = args[2]

# 读取表达矩阵
exp_df <- read.table(exp_tpm_file,header=TRUE,row.names=1,sep=',')
samples = colnames(exp_df)

# 提取分组信息
group <- read.table(group_file,header=FALSE,sep=',')
group = as.factor(group$V1)

dge <- DGEList(counts = exp_df, group = group)

# 过滤低表达基因
# keep <- rowSums(cpm(dge) > 1 ) >= 1
# dge <- dge[keep, , keep.lib.sizes = FALSE]

# 进行差异分析
dge <- calcNormFactors(dge)
dge <- estimateCommonDisp(dge)
dge <- estimateTagwiseDisp(dge)
fit <- exactTest(dge)
fit$table$FDR <- p.adjust(fit$table$PValue, method = "fdr")
deg_table = fit$table

normalized_counts <- cpm(dge, normalized.lib.sizes = TRUE)
write.csv(normalized_counts, file = paste(dirname(exp_tpm_file), "exp_norm.csv",sep="/"))
write.csv(deg_table, file = paste(dirname(exp_tpm_file), "deg.csv",sep="/"))
