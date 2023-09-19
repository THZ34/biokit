#!/mnt/cfs/project/test_freshman/liaoyuwei/utitities/R/R-3.6.1/bin/Rscript
library(data.table)
library(xCell)

dt <- read.csv(rawfile, row.names = 1)
cell_score <- xCellAnalysis(dt)
cell_score <- as.data.table(cell_score, keep.rownames = T)
setnames(cell_score, "rn", "CellType")
fwrite(cell_score, file = paste0(gsub("\\/$", "", outdir), "/xCell.result.tsv"), sep = "\t", col.names = T)

