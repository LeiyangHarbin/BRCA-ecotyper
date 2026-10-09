# ==== 配置 ====
#GPL96+——"GSE45255","GSE7390","GSE1456-GPL96"
infile <- "BRCA_fixed.rda"
dataset_to_extract <- "GSE1456-GPL96"
out_txt <- "GSE1456-GPL96.txt"
out_rda <- "GSE1456-GPL96_matrix.rda"

# 反 log2 方法： "2pow" -> 2^x （默认）， "minus1" -> 2^x - 1
INVERSE_METHOD <- "2pow"

# NA 删除策略（可选）
# 可选值： "none" | "columns" | "rows" | "both"
# "none"    -> 不删除 NA（保留全部样本/基因，后续保留 NA）
# "columns" -> 删除含任意 NA 的样本（列）
# "rows"    -> 删除含任意 NA 的基因（行）
# "both"    -> 先删除含 NA 的列，再删除含 NA 的行（或反之，脚本实现为先列后行）
NA_REMOVE_METHOD <- "columns"

# 是否详细输出
VERBOSE <- TRUE
log_msg <- function(...) {
  if (!VERBOSE) return(invisible(NULL))
  txt <- paste(..., collapse = " ")
  message(sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), txt))
}

# ==== 开始 ====
log_msg("开始运行脚本。输入文件:", infile, "; 目标条目:", dataset_to_extract)
log_msg("设置：INVERSE_METHOD=", INVERSE_METHOD, "; NA_REMOVE_METHOD=", NA_REMOVE_METHOD)

# 依赖检测（仅提示）
if (!requireNamespace("Biobase", quietly = TRUE)) {
  log_msg("提示：未安装 Biobase（若处理 ExpressionSet，建议安装 Biobase）。")
} else {
  log_msg("检测到 Biobase；如可用将载入。")
  library(Biobase)
}

# 载入 rda
if (!file.exists(infile)) stop("找不到输入文件：", infile)
log_msg("正在载入 rda 文件：", infile)
before_objs <- ls()
load(infile)
after_objs <- ls()
new_objs <- setdiff(after_objs, before_objs)
log_msg("载入完成。新增对象数量：", length(new_objs), "（示例：", paste(head(new_objs, 5), collapse = ", "), "）")

# 找主对象（优先 BRCA）
main_obj <- NULL
if (exists("BRCA")) {
  main_obj <- "BRCA"
  log_msg("检测到对象 'BRCA'，将使用之。")
} else {
  objs <- ls()
  log_msg("rda 中存在对象：", paste(objs, collapse = ", "))
  for (o in objs) {
    val <- get(o)
    if (is.list(val)) {
      main_obj <- o
      log_msg("选择第一个 list 对象：", o, "作为主对象。")
      break
    }
  }
}
if (is.null(main_obj)) stop("rda 中未找到 BRCA，也未找到合适的 list 对象。请检查 infile。")
BRCA <- get(main_obj)
log_msg("使用 RDA 对象：", main_obj, "；对象类型：", paste(class(BRCA), collapse = "/"))

# 辅助函数：提取表达矩阵
get_expr <- function(x) {
  log_msg("get_expr: 输入类型：", paste(class(x), collapse = "/"))
  if (inherits(x, "ExpressionSet")) {
    log_msg("输入为 ExpressionSet，提取 exprs。")
    return(as.matrix(Biobase::exprs(x)))
  }
  if (is.list(x)) {
    for (nm in c("expr","exprs","expression","mat","matrix")) {
      if (!is.null(x[[nm]])) {
        log_msg("在 list 中找到字段 '", nm, "' 并提取为矩阵。")
        return(as.matrix(x[[nm]]))
      }
    }
    if (length(x) >= 1 && (is.matrix(x[[1]]) || is.data.frame(x[[1]]))) {
      log_msg("list 的第1个元素看起来像矩阵，尝试提取第1个元素。")
      return(as.matrix(x[[1]]))
    }
  }
  if (is.matrix(x) || is.data.frame(x)) {
    log_msg("输入直接为 matrix/data.frame，尝试做 as.matrix()。")
    return(as.matrix(x))
  }
  stop("无法识别该条目为表达矩阵的结构。")
}

# 辅助函数：提取 pheno（如有）
get_pheno <- function(x) {
  log_msg("get_pheno: 尝试从对象中提取 pheno（若存在）。")
  if (inherits(x, "ExpressionSet")) {
    pd <- as.data.frame(Biobase::pData(x))
    log_msg("从 ExpressionSet 提取 pheno，行数：", nrow(pd), "列数：", ncol(pd))
    return(pd)
  }
  if (is.list(x)) {
    for (nm in c("pheno","pdata","pData","meta","colData","phenoData")) {
      if (!is.null(x[[nm]])) {
        log_msg("在 list 中找到 pheno 字段 '", nm, "' 并返回。")
        return(as.data.frame(x[[nm]]))
      }
    }
  }
  log_msg("未找到 pheno（返回 NULL）。")
  return(NULL)
}

# 去重：保留 mean 最大的行
dedup_keep_maxmean <- function(mat) {
  if (is.null(rownames(mat))) stop("表达矩阵没有行名（基因）。")
  rn <- rownames(mat)
  n_before <- nrow(mat)
  if (!any(duplicated(rn))) {
    log_msg("dedup_keep_maxmean: 无重复基因名（行数=", n_before, "）。")
    return(mat)
  }
  dup_unique <- unique(rn[duplicated(rn)])
  log_msg("发现重复基因名（unique）数:", length(dup_unique), "（示例：", paste(head(dup_unique,5), collapse = ", "), "）")
  keep_idx <- tapply(seq_along(rn), rn, function(idxs) {
    if (length(idxs) == 1) return(idxs)
    rowmeans <- rowMeans(mat[idxs, , drop = FALSE], na.rm = TRUE)
    idx_chosen <- idxs[which.max(rowmeans)]
    idx_chosen
  })
  mat_dedup <- mat[unlist(keep_idx), , drop = FALSE]
  log_msg("重复基因处理：原行数=", n_before, " -> 处理后行数=", nrow(mat_dedup))
  mat_dedup
}

# log2 判断（启发式）
is_log2_matrix <- function(mat) {
  v <- as.numeric(mat)
  v <- v[is.finite(v)]
  if (length(v) == 0) return(FALSE)
  qs <- quantile(v, probs = c(0, .01, .25, .5, .75, .99, 1), na.rm = TRUE)
  log_msg("数值分位数（0,1%,25%,50%,75%,99%,100%）：", paste(round(qs, 4), collapse = ", "))
  if (any(v < 0)) {
    log_msg("检测到负值，倾向于认为是 log2 数据。")
    return(TRUE)
  }
  if (qs["99%"] <= 25) {
    log_msg("99% 分位数 <= 25，倾向于认为是 log2 数据。")
    return(TRUE)
  }
  log_msg("不认为是 log2 数据（99% 分位数 > 25 且无负值）。")
  return(FALSE)
}

# --- 找到目标数据条目 ---
hits <- names(BRCA)[names(BRCA) == dataset_to_extract]
if (length(hits) == 0) {
  hits <- names(BRCA)[grepl(dataset_to_extract, names(BRCA), ignore.case = TRUE, fixed = FALSE)]
}
if (length(hits) == 0) stop("在 BRCA 中未找到与 '", dataset_to_extract, "' 匹配的条目。请确认名称。")
if (length(hits) > 1) log_msg("模糊匹配到多个候选，默认使用第一个匹配：", paste(hits, collapse = ", "))
name <- hits[1]
log_msg("将提取条目：", name)

# 提取矩阵与 pheno
mat <- tryCatch(get_expr(BRCA[[name]]), error = function(e) stop("提取表达矩阵失败：", e$message))
pheno <- tryCatch(get_pheno(BRCA[[name]]), error = function(e) { log_msg("提取 pheno 失败：", e$message); NULL })

# 显示维度与样本名示例
if (!is.null(colnames(mat))) {
  log_msg("提取到的表达矩阵维度：", paste(dim(mat), collapse = " x "), "; 示例列名：", paste(head(colnames(mat),5), collapse = ", "))
} else {
  log_msg("注意：提取到的矩阵没有列名（样本名）。")
}
if (is.null(rownames(mat))) log_msg("警告：提取到的矩阵没有行名（基因）。")

# ---------- 简化的 normal 筛选（基于 pheno$Tumor.or.normal） ----------
original_colnames <- colnames(mat)
log_msg("保存原始样本名（用于 normal 检测）。示例：", paste(head(original_colnames,5), collapse = ", "))

normals_matched <- character(0)
normals_unmatched <- character(0)

if (!is.null(pheno) && "Tumor.or.normal" %in% colnames(pheno)) {
  log_msg("pheno 中存在列 'Tumor.or.normal'，开始检测包含 'normal' 的样本（不区分大小写）。")
  # 取该列值
  colvals <- as.character(pheno$Tumor.or.normal)
  # 找到包含 normal (case-insensitive) 的行索引
  idx <- grepl("normal", colvals, ignore.case = TRUE)
  if (any(idx, na.rm = TRUE)) {
    # 尝试优先用 pheno 的行名映射到表达矩阵样本名
    pheno_rows <- if (!is.null(rownames(pheno))) rownames(pheno) else NULL
    matched_samples <- character(0)
    for (i in which(idx)) {
      # 尝试按 pheno 行名匹配表达矩阵原始列名
      candidate <- NA
      if (!is.null(pheno_rows) && pheno_rows[i] %in% original_colnames) {
        candidate <- pheno_rows[i]
      } else {
        # 如果 pheno 中该字段本身恰好等于某个样本名，也映射
        if (!is.na(colvals[i]) && colvals[i] %in% original_colnames) candidate <- colvals[i]
      }
      if (!is.na(candidate)) {
        matched_samples <- c(matched_samples, candidate)
      } else {
        # 无法映射到样本名，则记录为 unmatched 候选（供人工检查）
        normals_unmatched <- c(normals_unmatched, paste0("pheno_row_", i, ":", colvals[i]))
      }
    }
    normals_matched <- unique(matched_samples)
    log_msg("在 pheno 的 'Tumor.or.normal' 中检测到含 'normal' 的条目数：", sum(idx, na.rm = TRUE),
            "；成功映射到表达矩阵样本数：", length(normals_matched),
            "；未映射到样本（需人工核查）数：", length(normals_unmatched))
    if (length(normals_unmatched) > 0) {
      log_msg("未映射示例（前50）：", paste(head(normals_unmatched,50), collapse = ", "))
    }
  } else {
    log_msg("pheno 的 'Tumor.or.normal' 列中未检测到包含 'normal' 的值，未移除样本。")
  }
} else {
  log_msg("pheno 为 NULL 或未包含 'Tumor.or.normal' 列，跳过基于该列的 normal 检测（保留全部样本）。")
}
# ---------- end normal 筛选 ----------

# 如果检测到 matched，就移除这些样本（基于 original_colnames）
if (length(normals_matched) > 0) {
  keep_cols <- setdiff(original_colnames, normals_matched)
  removed <- setdiff(original_colnames, keep_cols)
  log_msg("将移除的 normal/控制样本数（基于 pheno$Tumor.or.normal 匹配）：", length(removed),
          "; 示例移除列名前5：", paste(head(removed,5), collapse = ", "))
  if (length(keep_cols) == 0) stop("移除所有样本后无样本可用，请检查 pheno。")
  mat <- mat[, keep_cols, drop = FALSE]
  # 同步 pheno（若 pheno 的行名为样本名）
  if (!is.null(pheno) && !is.null(rownames(pheno)) && any(rownames(pheno) %in% colnames(mat))) {
    pheno <- pheno[rownames(pheno) %in% colnames(mat), , drop = FALSE]
    log_msg("已同步更新 pheno，pheno 行数现在：", nrow(pheno))
  } else {
    log_msg("pheno 未同步更新（pheno 行名与表达矩阵列名不匹配或 pheno 为 NULL）。")
  }
} else {
  log_msg("未检测到需移除的 normal 样本（或未能映射到表达矩阵），保留全部样本。")
}

# ========== 新增：删除 pheno$Tumor.or.normal 为 NA 的样本（若能映射） ==========
if (!is.null(pheno) && "Tumor.or.normal" %in% colnames(pheno)) {
  colvals2 <- pheno$Tumor.or.normal
  na_idx2 <- which(is.na(colvals2))
  if (length(na_idx2) > 0) {
    log_msg("发现 pheno$Tumor.or.normal 中存在 NA 的条目数：", length(na_idx2), "；尝试映射并删除对应样本。")
    pheno_rows2 <- if (!is.null(rownames(pheno))) rownames(pheno) else NULL
    # 能直接映射到表达矩阵列名的样本（优先使用 pheno 的行名）
    candidates_na <- character(0)
    if (!is.null(pheno_rows2)) {
      candidates_na <- intersect(pheno_rows2[na_idx2], colnames(mat))
    }
    # 记录未能映射的 NA（供人工检查）
    unmatched_na <- setdiff(pheno_rows2[na_idx2], candidates_na)
    if (length(candidates_na) > 0) {
      log_msg("将移除的样本（pheno 中 Tumor.or.normal 为 NA 并成功映射到表达矩阵列名）数：", length(candidates_na),
              "；示例前5：", paste(head(candidates_na,5), collapse = ", "))
      mat <- mat[, setdiff(colnames(mat), candidates_na), drop = FALSE]
      # 同步 pheno（若 pheno 的行名为样本名）
      if (!is.null(rownames(pheno))) {
        pheno <- pheno[rownames(pheno) %in% colnames(mat), , drop = FALSE]
        log_msg("已同步更新 pheno，pheno 行数现在：", nrow(pheno))
      } else {
        log_msg("pheno 未同步更新（pheno 无行名或行名与表达矩阵列名不匹配）。")
      }
    } else {
      log_msg("沒有样本能根据 pheno 行名映射到表达矩阵列名，因此未执行删除（请人工核查）。")
    }
    if (length(unmatched_na) > 0) {
      log_msg("pheno 中 Tumor.or.normal 为 NA 但未能映射到表达矩阵的 pheno 行名示例（前50）：",
              paste(head(unmatched_na,50), collapse = ", "))
    }
    # 检查是否删除后数据为空
    if (ncol(mat) == 0) stop("在删除 pheno 中 Tumor.or.normal 为 NA 的样本后无可用样本（列），已中止。")
  } else {
    log_msg("pheno$Tumor.or.normal 中无 NA 值，跳过此删除步骤。")
  }
} else {
  log_msg("pheno 为 NULL 或不包含 'Tumor.or.normal' 列，跳过删除 pheno NA 的步骤。")
}
# ======================================================

# 对列名做 make.names 以保证唯一（在移除 normal/NA 之后进行）
original_colnames_after <- colnames(mat)
colnames(mat) <- make.names(original_colnames_after, unique = TRUE)
if (!identical(original_colnames_after, colnames(mat))) {
  log_msg("样本名已标准化 / 唯一化（make.names）。示例：", paste(head(colnames(mat),5), collapse = ", "))
}

# 若无行名但第一列为 gene/symbol，尝试设为行名
if (is.null(rownames(mat))) {
  if (is.matrix(mat) || is.data.frame(mat)) {
    cn <- colnames(mat)
    if (!is.null(cn) && length(cn) >= 1 && grepl("gene|symbol", tolower(cn[1]))) {
      log_msg("发现第一列可能为基因名字段（列名:", cn[1], "），尝试将其设为行名。")
      rn <- as.character(mat[,1])
      mat <- as.matrix(mat[,-1,drop=FALSE])
      rownames(mat) <- rn
    } else stop("表达矩阵没有行名（基因），也无法自动识别基因列。请先在 BRCA 中用基因名作为行名。")
  }
}

# 强制数值（保留行/列名）
safe_coerce_numeric <- function(mat) {
  rn <- rownames(mat)
  cn <- colnames(mat)
  mat_num <- apply(mat, 2, function(col) suppressWarnings(as.numeric(as.character(col))))
  if (is.null(dim(mat_num))) { # single column
    mat_num <- matrix(mat_num, ncol = 1)
  }
  colnames(mat_num) <- cn
  rownames(mat_num) <- rn
  mat_num
}
v_check <- suppressWarnings(as.numeric(as.vector(mat)))
non_numeric_count <- sum(!is.finite(v_check), na.rm = TRUE)
if (non_numeric_count > 0) {
  log_msg("检测到非数值或非有限值元素数：", non_numeric_count, "；尝试强制转换为数值。")
  mat <- safe_coerce_numeric(mat)
}

log_msg("当前矩阵维度（行=基因, 列=样本）=", nrow(mat), "x", ncol(mat))
na_count_total <- sum(!is.finite(as.numeric(as.vector(mat))), na.rm = TRUE)
log_msg("表达矩阵中 NA / 非有限 值总数：", na_count_total)

# ========== NA 删除策略（根据 NA_REMOVE_METHOD） ==========
if (!(NA_REMOVE_METHOD %in% c("none", "columns", "rows", "both"))) {
  stop("未知的 NA_REMOVE_METHOD：", NA_REMOVE_METHOD, "；请设置为 'none'|'columns'|'rows'|'both'。")
}

if (NA_REMOVE_METHOD != "none") {
  # 计算含 NA 的行/列数
  na_per_sample_counts <- colSums(!is.finite(mat))
  na_per_gene_counts <- rowSums(!is.finite(mat))
  
  if (NA_REMOVE_METHOD %in% c("columns", "both")) {
    cols_with_any_na <- which(na_per_sample_counts > 0)
    if (length(cols_with_any_na) > 0) {
      removed_samples <- colnames(mat)[cols_with_any_na]
      log_msg("根据 NA_REMOVE_METHOD='", NA_REMOVE_METHOD, "'：检测到含 NA/非有限 值的样本数：", length(removed_samples),
              "；将移除这些样本（示例前5）：", paste(head(removed_samples,5), collapse = ", "))
      mat <- mat[, -cols_with_any_na, drop = FALSE]
      # 同步 pheno（若 pheno 的行名为样本名）
      if (!is.null(pheno) && !is.null(rownames(pheno))) {
        pheno <- pheno[rownames(pheno) %in% colnames(mat), , drop = FALSE]
        log_msg("已同步更新 pheno，pheno 行数现在：", nrow(pheno))
      } else {
        log_msg("pheno 未同步更新（pheno 行名与表达矩阵列名不匹配或 pheno 为 NULL）。")
      }
    } else {
      log_msg("根据 NA_REMOVE_METHOD='", NA_REMOVE_METHOD, "'：未检测到需移除的样本列。")
    }
  }
  
  if (NA_REMOVE_METHOD %in% c("rows", "both")) {
    # 刷新行 NA 计数（因为如果之前删除了列, 行的 NA 状况已经改变）
    na_per_gene_counts <- rowSums(!is.finite(mat))
    rows_with_any_na <- which(na_per_gene_counts > 0)
    if (length(rows_with_any_na) > 0) {
      removed_genes <- rownames(mat)[rows_with_any_na]
      log_msg("根据 NA_REMOVE_METHOD='", NA_REMOVE_METHOD, "'：检测到含 NA/非有限 值的基因数：", length(removed_genes),
              "；将移除这些基因（示例前5）：", paste(head(removed_genes,5), collapse = ", "))
      mat <- mat[-rows_with_any_na, , drop = FALSE]
    } else {
      log_msg("根据 NA_REMOVE_METHOD='", NA_REMOVE_METHOD, "'：未检测到需移除的基因行。")
    }
  }
  
  # 删除后检查
  if (ncol(mat) == 0) stop("按 NA_REMOVE_METHOD 删除后无可用样本（列），已中止。请改为 'none' 或改用阈值策略。")
  if (nrow(mat) == 0) stop("按 NA_REMOVE_METHOD 删除后无可用基因（行），已中止。请改为 'none' 或改用阈值策略。")
} else {
  log_msg("NA_REMOVE_METHOD='none'：不删除含 NA 的行/列，继续保留全部数据。")
}
# ======================================================

# 去重
n_genes_before <- nrow(mat)
mat <- dedup_keep_maxmean(mat)
n_genes_after <- nrow(mat)
log_msg("重复基因处理：基因数 ", n_genes_before, " -> ", n_genes_after, "（去重 ", n_genes_before - n_genes_after, "）")

# log2 检测
log2_flag <- is_log2_matrix(mat)
log_msg("检测到可能为 log2 空间：", log2_flag)

# 反变换
mat_nonlog <- NULL
if (log2_flag) {
  if (INVERSE_METHOD == "2pow") {
    log_msg("使用 2^x 反变换（INVERSE_METHOD = '2pow'）。")
    mat_nonlog <- tryCatch(2^mat, error = function(e) stop("反变换失败：", e$message))
  } else if (INVERSE_METHOD == "minus1") {
    log_msg("使用 2^x - 1 反变换（INVERSE_METHOD = 'minus1'），并对负值取 0。")
    mat_nonlog <- tryCatch(pmax(0, 2^mat - 1), error = function(e) stop("反变换失败：", e$message))
  } else stop("未知的 INVERSE_METHOD：", INVERSE_METHOD)
} else {
  log_msg("未检测到 log2 空间，直接按原值输出（请确保这就是你需要的非对数矩阵）。")
  mat_nonlog <- mat
}

# 反变换后汇总
safe_num_summary <- function(m) {
  v <- as.numeric(m)
  v <- v[is.finite(v)]
  if (length(v) == 0) return(list(n = 0))
  list(n = length(v),
       min = min(v),
       q1 = quantile(v, .25),
       median = median(v),
       mean = mean(v),
       q3 = quantile(v, .75),
       q99 = quantile(v, .99),
       max = max(v))
}
s_before <- safe_num_summary(mat)
s_after <- safe_num_summary(mat_nonlog)
log_msg("反变换前数值元素数:", s_before$n, "; 反变换后数值元素数:", s_after$n)
log_msg("反变换后 summary (min, q1, median, mean, q3, q99, max):",
        paste(round(c(s_after$min, s_after$q1, s_after$median, s_after$mean, s_after$q3, s_after$q99, s_after$max), 4), collapse = ", "))

# 最终 tidy 样本名
colnames(mat_nonlog) <- make.names(colnames(mat_nonlog), unique = TRUE)
log_msg("最终样本数（列）=", ncol(mat_nonlog), "；最终基因数（行）=", nrow(mat_nonlog))

# 写出 txt（首列 GeneSymbol）
df_out <- data.frame(GeneSymbol = rownames(mat_nonlog), mat_nonlog, check.names = FALSE, stringsAsFactors = FALSE)
safe_write <- function(df, file) {
  ok <- tryCatch({
    write.table(df, file = file, sep = "\t", quote = FALSE, row.names = FALSE)
    TRUE
  }, error = function(e) {
    log_msg("写出文件失败：", e$message)
    FALSE
  })
  ok
}
log_msg("写出非对数矩阵 txt（首列 GeneSymbol），路径：", out_txt)
ok <- safe_write(df_out, out_txt)
if (ok) {
  log_msg("已成功写出非对数矩阵 txt：", out_txt, "; 文件大小：",
          if (file.exists(out_txt)) format(file.size(out_txt), big.mark = ",") else "未知")
} else log_msg("写出 txt 文件失败，请检查磁盘权限或路径。")

# 保存 rda 备份
save_obj <- list(
  dataset_name = name,
  original_extracted = mat,
  pheno = pheno,
  detected_log2 = log2_flag,
  inverse_method = INVERSE_METHOD,
  na_remove_method = NA_REMOVE_METHOD,
  mat_nonlog = mat_nonlog,
  extracted_gene_count = nrow(mat),
  extracted_sample_count = ncol(mat),
  timestamp = Sys.time(),
  sessionInfo = sessionInfo()
)
log_msg("保存备份 rda 到：", out_rda)
save_ok <- tryCatch({ save(save_obj, file = out_rda); TRUE }, error = function(e) { log_msg("保存 rda 失败：", e$message); FALSE })
if (save_ok) log_msg("已保存备份 rda：", out_rda) else stop("保存备份 rda 失败，已中止。")

# 结束
log_msg("最终输出（mat_nonlog）维度： 行（基因）=", nrow(mat_nonlog), "；列（样本）=", ncol(mat_nonlog))
log_msg("完成。脚本结束时间：", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
