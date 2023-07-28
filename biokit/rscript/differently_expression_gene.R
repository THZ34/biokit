library("edgeR")
library("stringi")

ncounter_parser <- function(ncounter_raw) {
  grp_line <- readLines(raw_exp_file, n = 2)[2]
  grp <- strsplit(grp_line, ',')
  grp <- grp[[1]]
  grp <- grp[5:length(grp)]
  grp <- as.factor(grp)

  exp_df <- read.csv(raw_exp_file)
  # 取[2:,4:]的表达量矩阵
  exp_df <- exp_df[3:dim(exp_df)[1], 5:(4 + length(grp))]
  exp_df <- as.exp_dfa.frame(lapply(exp_df, as.numeric))
  return(grp, exp_df)
}

deg_func1 <- function(exp_df, grp) {
  dge <- DGEList(counts = exp_df, group = grp)
  dge <- calcNormFactors(dge)
  design <- model.matrix(~grp)

  dge <- estimateGLMCommonDisp(dge, design)
  dge <- estimateGLMTagwiseDisp(dge, design)
  dge <- estimateGLMTrendedDisp(dge, design)

  fit <- glmFit(dge, design,robust = TRUE)
  lrt <- topTags(glmLRT(fit), n = nrow(dge$counts))

  return(lrt$table) }

deg_func2 <- function (exp_df,grp){
dge <- DGEList(counts=exp_df,genes=rownames(exp_df),grp = grp)
dge<- calcNormFactors(dge)
design <- model.matrix(~ grp)
dge <- estimateGLMCommonDisp(dge, design)
et_compare1 <- exactTest(dge)
de_compare1 <- decideTestsDGE(et_compare1)
tt_compare1 <- topTags(et_compare1, n=Inf)$table
}


