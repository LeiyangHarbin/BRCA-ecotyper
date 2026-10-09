# merge_recover.R
files <- c("GSE11121.txt","GSE12093.txt","GSE17705.txt","GSE2034.txt","GSE25066.txt","GSE2603.txt","GSE5327.txt","GSE45255.txt","GSE7390.txt","GSE1456-GPL96.txt")

# 读取函数
read_mat <- function(f) {
  m <- as.matrix(read.table(f, header=TRUE, row.names=1, sep="\t", check.names=FALSE))
  mode(m) <- "numeric"
  return(m)
}

mats <- lapply(files, read_mat)
common <- Reduce(intersect, lapply(mats, rownames))
if(length(common)==0) stop("没有共同基因，停止。")
mats <- lapply(mats, function(m) m[common, , drop=FALSE])
merged <- do.call(cbind, mats)
write.table(cbind(Gene=rownames(merged), merged), "merged_raw.txt", sep="\t", row.names=FALSE, quote=FALSE)

if(!requireNamespace("sva", quietly=TRUE)) stop("（BiocManager::install('sva')）")
library(sva)
batch <- rep(seq_along(mats), sapply(mats, ncol))
combat_mat <- ComBat(dat=merged, batch=batch, par.prior=TRUE, prior.plots=FALSE)
write.table(cbind(Gene=rownames(combat_mat), combat_mat), "merged_GPL96.txt", sep="\t", row.names=FALSE, quote=FALSE)
