# coding='utf-8'
# Author:Tang Hongzhen
# Email: tanghongzhen34@gmail.com
# %%
# request CIBERSORT.R
# shell : export CIBERSORT_BINARY=/path/to/CIBERSORT.R
exp_df = read.csv('TCGA_PAAD_tpm_df.csv', row.names = 1)
all_result = data.frame()
indications = rep('paad', times = ncol(exp_df))

for (method in c('quantiseq', 'timer', 'cibersort', 'cibersort_abs', 'mcp_counter', 'xcell', 'epic', 'abis', 'consensus_tme', 'estimate')) {
  result <- immunedeconv::deconvolute(gene_expression = exp_df, method = method, tumor = TRUE, indications = indications)
  result$method = method
  all_result = rbind(all_result, result)
}
