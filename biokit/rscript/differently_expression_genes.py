# coding='utf-8'
# Author:Tang Hongzhen
# Email: tanghongzhen34@gmail.com
# %%
from rpy2.robjects import pandas2ri
from rpy2.robjects import r

pandas2ri.activate()


# %%
def init_degfunc():
    r("""
    library('edgeR')
    deg_func <- function(exp_df, grp) {
        dge <- DGEList(counts = exp_df, group = grp)
        dge <- calcNormFactors(dge)
        
        design <- model.matrix(~grp)
        
        dge <- estimateGLMCommonDisp(dge, design)
        dge <- estimateGLMTagwiseDisp(dge, design)
        dge <- estimateGLMTrendedDisp(dge, design)
        
        fit <- glmFit(dge, design)
        lrt <- glmLRT(fit, contrast = c(-1, 1))
        return(lrt$table)
    }""")
