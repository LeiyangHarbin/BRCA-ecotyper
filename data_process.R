

##################################################
rda_path <- "BRCA_fixed.rda"
out_prefix <- "BRCA_mnorm_combat_intersect"
GPL570_list <- c("GSE16446","GSE20685","GSE20711","GSE42568","GSE48390","GSE58812")
GPL96_list  <- c("GSE45255","GSE7390","GSE1456-GPL96")
scale_other_by_median <- FALSE  # MNORM 可选参数

required_pkgs <- c("Biobase", "sva")
for(p in required_pkgs){
  if(!requireNamespace(p, quietly = TRUE)){
    install.packages(p)
  }
  library(p, character.only = TRUE)
}

if(!exists("BRCA")){
  if(file.exists(rda_path) && nzchar(rda_path)){
    load(rda_path)
    if(!exists("BRCA")) stop("加载后未发现变量 BRCA，请确认 rda 文件里变量名")
  } else {
    stop("无法找到 BRCA_fixed.rda")
  }
} else {
  message("BRCA 已存在于 workspace，使用该变量。")
}

message("BRCA 中的 datasets：")
print(names(BRCA))

# ---- 检测 normal 的函数 ----
detect_normal_samples <- function(eset){
  ph <- pData(eset)
  candidate_cols <- intersect(colnames(ph), c("Tumor.or.normal","tumor_or_normal","characteristics_ch1",
                                              "source_name_ch1","tissue","sample_type","SampleType","sample_type"))
  normals <- rep(FALSE, nrow(ph))
  if(length(candidate_cols) > 0){
    for(col in candidate_cols){
      vals <- as.character(ph[[col]])
      vals[is.na(vals)] <- ""
      norm_flag <- grepl("normal|Normal|Normal tissue|healthy|adjacent|adj", vals, ignore.case = TRUE)
      normals <- normals | norm_flag
    }
  } else {
    for(col in colnames(ph)){
      vals <- as.character(ph[[col]])
      vals[is.na(vals)] <- ""
      normals <- normals | grepl("normal|Normal|Normal tissue|healthy|adjacent|adj", vals, ignore.case = TRUE)
    }
  }
  names(normals) <- rownames(ph)
  return(normals)
}

# ---- MNORM 函数 ----
mnorm_map_to_reference <- function(reference_mat, other_mat, ties.method = "average", scale_other_by_median = FALSE){
  if(!all(rownames(reference_mat) == rownames(other_mat))){
    stop("请先确保 reference_mat 与 other_mat 行名完全一致且顺序相同。")
  }
  G <- nrow(reference_mat)
  ref_medians <- apply(reference_mat, 1, median, na.rm = TRUE)
  ref_sorted <- sort(ref_medians, decreasing = FALSE)
  if(isTRUE(scale_other_by_median)){
    other_mat <- sweep(other_mat, 2, apply(other_mat, 2, median, na.rm = TRUE), "/")
  }
  rank_mat <- apply(other_mat, 2, function(col) rank(col, ties.method = ties.method))
  mapped <- matrix(NA, nrow = G, ncol = ncol(other_mat),
                   dimnames = list(rownames(other_mat), colnames(other_mat)))
  for(j in seq_len(ncol(rank_mat))){
    pos <- rank_mat[, j]
    pos[pos < 1] <- 1
    pos[pos > G] <- G
    low <- floor(pos)
    high <- ceiling(pos)
    frac <- pos - low
    mapped_vals <- ifelse(low == high,
                          ref_sorted[low],
                          ref_sorted[low] * (1 - frac) + ref_sorted[high] * frac)
    mapped[, j] <- mapped_vals
  }
  return(mapped)
}

# ---- 主流程 ----
report <- list()
removed_normals <- list()

get_expr_matrix_from_eset <- function(eset){
  as.matrix(Biobase::exprs(eset))
}

# 检查所选数据集是否存在
missing_in_brca <- setdiff(c(GPL570_list, GPL96_list), names(BRCA))
if(length(missing_in_brca) > 0){
  stop("下面这些数据集在 BRCA 中找不到，请核对名字：\n", paste(missing_in_brca, collapse=", "))
}

# 1) 处理 GPL570：剔除 normal，直接取表达矩阵（rownames 已为 Gene Symbol 且已处理重复）
gpl570_mats <- list()
for(ds in GPL570_list){
  eset <- BRCA[[ds]]
  message("处理（GPL570）: ", ds)
  normals_flag <- detect_normal_samples(eset)
  n_normals <- sum(normals_flag)
  removed_normals[[ds]] <- names(normals_flag)[which(normals_flag)]
  report[[ds]] <- list(n_samples = nrow(pData(eset)), n_normals = n_normals)
  eset_use <- if(n_normals > 0) eset[, !normals_flag, drop=FALSE] else eset
  expr <- get_expr_matrix_from_eset(eset_use)
  # 不做任何探针->symbol 或重复折叠
  gpl570_mats[[ds]] <- expr
  message("  基因数: ", nrow(expr), " 样本数: ", ncol(expr))
}

# 取 GPL570 的共同基因（交集）作为参考基因集合
common_genes_570 <- Reduce(intersect, lapply(gpl570_mats, rownames))
message("GPL570 6 个集合的交集基因数: ", length(common_genes_570))
if(length(common_genes_570) < 1000) warning("GPL570 交集基因数较少，请确认行名为 Gene Symbol 且一致。")
# 建立 reference 矩阵（按交集基因并按相同顺序）
ref_mat <- do.call(cbind, lapply(gpl570_mats, function(m) m[common_genes_570, , drop=FALSE]))
message("参考矩阵大小 (genes x samples): ", paste(dim(ref_mat), collapse = " x "))

# 2) 处理 GPL96：剔除 normal，直接取表达矩阵（rownames 为 Gene Symbol）
gpl96_raw <- list()
for(ds in GPL96_list){
  eset <- BRCA[[ds]]
  message("处理（GPL96）: ", ds)
  normals_flag <- detect_normal_samples(eset)
  n_normals <- sum(normals_flag)
  removed_normals[[ds]] <- names(normals_flag)[which(normals_flag)]
  report[[ds]] <- list(n_samples = nrow(pData(eset)), n_normals = n_normals)
  eset_use <- if(n_normals > 0) eset[, !normals_flag, drop=FALSE] else eset
  expr <- get_expr_matrix_from_eset(eset_use)
  gpl96_raw[[ds]] <- expr
  message("  基因数: ", nrow(expr), " 样本数: ", ncol(expr))
}

# 3) 计算最终用于合并的共同基因
#    final_genes = intersection(common_genes_570, each GPL96 rownames)
all_gene_sets <- c(list(common_genes_570), lapply(gpl96_raw, rownames))
final_genes <- Reduce(intersect, all_gene_sets)
message("所有选定数据集的共同基因数 (final_genes): ", length(final_genes))
if(length(final_genes) == 0) stop("没有共同基因，无法继续。请检查每个 ExpressionSet 的行名是否一致（Gene Symbol）。")

# 4) 用 final_genes 做 MNORM：ref_sub 与 each other_sub 行名和顺序都为 final_genes
ref_sub <- ref_mat[final_genes, , drop = FALSE]

gpl96_mats_mnormed <- list()
for(ds in names(gpl96_raw)){
  other_mat <- gpl96_raw[[ds]]
  # 确保其他矩阵包含 final_genes
  if(!all(final_genes %in% rownames(other_mat))){
    stop("在进行 MNORM 前，", ds, " 中缺少 final_genes 中的某些基因（应已计算为 intersection，但发现缺失）。")
  }
  other_sub <- other_mat[final_genes, , drop = FALSE]
  # 确保行顺序与 reference 一致
  other_sub <- other_sub[rownames(ref_sub), , drop = FALSE]
  mapped_other <- mnorm_map_to_reference(reference_mat = ref_sub, other_mat = other_sub,
                                         ties.method = "average", scale_other_by_median = scale_other_by_median)
  gpl96_mats_mnormed[[ds]] <- mapped_other
  message("MNORM: ", ds, " -> 完成 (", nrow(mapped_other), " genes x ", ncol(mapped_other), " samples).")
}

# 5) 合并：只使用 final_genes
combined_list <- list()
batch_vec <- c()
sample_meta <- data.frame(sample=character(), dataset=character(), stringsAsFactors = FALSE)

# GPL570 的样本（ref_sub 已包含所有 GPL570 的样本列）
# 但是 ref_sub 包含了所有 GPL570 数据的样本列（按我们之前拼接 ref_mat），
# 我们需要把这些样本的 dataset 标签分配回各自 dataset 名称：
# 先把每个 GPL570 的矩阵按 final_genes 子集化并加入 combined_list（保持 dataset 名称）
for(ds in GPL570_list){
  m <- gpl570_mats[[ds]]
  # 确保包含 final_genes
  if(!all(final_genes %in% rownames(m))) stop(paste0("GPL570 数据 ", ds, " 缺少 final_genes。"))
  m2 <- m[final_genes, , drop = FALSE]
  combined_list[[ds]] <- m2
  batch_vec <- c(batch_vec, rep(ds, ncol(m2)))
  sample_meta <- rbind(sample_meta, data.frame(sample = colnames(m2), dataset = ds, stringsAsFactors = FALSE))
}

# 加入 GPL96 的 MNORM 结果（按 final_genes）
for(ds in names(gpl96_mats_mnormed)){
  m <- gpl96_mats_mnormed[[ds]]
  combined_list[[ds]] <- m
  batch_vec <- c(batch_vec, rep(ds, ncol(m)))
  sample_meta <- rbind(sample_meta, data.frame(sample = colnames(m), dataset = ds, stringsAsFactors = FALSE))
}

combined_expr <- do.call(cbind, combined_list)
message("合并后的表达矩阵大小 (genes x samples): ", paste(dim(combined_expr), collapse = " x "))

# 6) 删除含 NA 的基因（不做任何填补）
na_per_gene <- rowMeans(is.na(combined_expr))
genes_with_na <- sum(na_per_gene > 0)
if(genes_with_na > 0){
  message("检测到含 NA 的基因数: ", genes_with_na, "（这些基因将被删除，脚本不进行填补）。")
}
gene_keep_idx <- which(na_per_gene == 0)
combined_expr_filt <- combined_expr[gene_keep_idx, , drop = FALSE]
message("删除含 NA 后的基因数: ", nrow(combined_expr_filt), "; 样本数: ", ncol(combined_expr_filt))

# 7) ComBat 批次校正（每个 dataset 为一批）
library(sva)
batch <- as.factor(batch_vec)
# 保证 sample_meta 顺序与矩阵列名一致
sample_meta <- sample_meta[match(colnames(combined_expr_filt), sample_meta$sample), , drop = FALSE]
if(any(is.na(sample_meta$sample))) stop("sample_meta 与表达矩阵样本名无法对齐，请检查。")
mod <- NULL
message("Running ComBat with ", length(unique(batch)), " batches ...")
combat_res <- ComBat(dat = as.matrix(combined_expr_filt), batch = batch, mod = mod, par.prior = TRUE, prior.plots = FALSE)
message("ComBat 完成。输出大小: ", paste(dim(combat_res), collapse = " x "))

# 8) 保存结果与报告
save(final_genes, combined_expr_filt, combat_res, sample_meta, removed_normals, report,
     file = paste0(out_prefix, "_results.rda"))
#write.csv(as.data.frame(combined_expr_filt), file = paste0(out_prefix, "_combined_expr_intersect.csv"))
#write.csv(as.data.frame(combat_res), file = paste0(out_prefix, "_combat_corrected_expr.csv"))
#write.csv(sample_meta, file = paste0(out_prefix, "_sample_meta.csv"))

message("===== 简要报告 =====")
for(ds in names(report)){
  message(ds, ": 总样本数=", report[[ds]]$n_samples, ", 被标记为 normal=", report[[ds]]$n_normals)
}
message("最终用于 ComBat 的样本数: ", ncol(combined_expr_filt))
message("最终用于 ComBat 的基因数: ", nrow(combined_expr_filt))  
message("结果已保存到: ", paste0(out_prefix, "_results.rda"))
message("脚本完成。")
