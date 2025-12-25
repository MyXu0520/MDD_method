###############################   计算腺癌PMD  ###################################
suppressPackageStartupMessages({
  library(MethylSeekR)
  library(GenomicRanges)
  library(data.table)
})

# 防止某些函数内部调用 png() 还需要 X11
options(bitmapType = "cairo")

# ========== 路径设置 ==========
input_dir  <- "/data3/xumy_PMD/data"
output_dir <- file.path(input_dir, "EAC_GEJ")

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# 只使用这 12 个样本：
sample_names <- c(
  paste0("EAC_", c(1:4,6)),
  paste0("GEJ_", 1:7)
)
# 对应的输入文件名：EAC_1methyCov.bed ... GEJ_7methyCov.bed
files <- file.path(input_dir, paste0(sample_names, "methyCov.bed"))

# 检查文件是否存在
missing <- !file.exists(files)
if (any(missing)) {
  stop("以下文件未找到，请检查文件名或路径：\n",
       paste(files[missing], collapse = "\n"))
}

# ========== 参数，可以按需修改 ==========
minCov    <- 5      # 至少覆盖 reads 数
num.cores <- 4      # 用几个核跑 HMM（没有并行就设为 1）

# ========== 辅助函数：从文件读取并构建 GRanges ==========
# 当前格式：V1=chr, V2=pos, V3=beta, V4=meth, V5=unmeth
load_methyl_gr <- function(f, minCov = 5) {
  message("读取文件: ", f)
  dt <- fread(f)
  
  if (ncol(dt) < 5) {
    stop("文件列数 < 5，检查 methyCov.bed 格式: ", f)
  }
  setnames(dt, 1:5, c("chr", "pos", "beta", "meth", "unmeth"))
  
  dt[, T := meth + unmeth]
  dt[, M := meth]
  
  dt <- dt[T >= minCov]
  if (nrow(dt) == 0) {
    stop("过滤后没有 CpG，检查 minCov 或原始数据: ", f)
  }
  
  gr <- GRanges(
    seqnames = dt$chr,
    ranges   = IRanges(start = dt$pos, end = dt$pos),
    strand   = "*",
    T        = as.integer(dt$T),
    M        = as.integer(dt$M)
  )
  
  sort(gr)
}

# ========== 主循环：对每个样本调用 segmentPMDs ==========
for (i in seq_along(sample_names)) {
  sample_name <- sample_names[i]
  f           <- files[i]
  
  message("===== 处理样本: ", sample_name, " =====")
  
  meth.gr <- load_methyl_gr(f, minCov = minCov)
  
  # 用当前数据估一个 seqLengths
  seqLengths_vec <- tapply(end(meth.gr), as.character(seqnames(meth.gr)), max)
  seqLengths     <- as.numeric(seqLengths_vec)
  names(seqLengths) <- names(seqLengths_vec)
  
  # 选一个染色体来训练 HMM
  chr.sel <- names(seqLengths)[1]
  message("用于训练 HMM 的染色体: ", chr.sel)
  
  # ⭐ 关键改动：强制画到 pdf，不往屏幕上画
  pdf_file <- file.path(output_dir, paste0(sample_name, "_PMD_HMM.pdf"))
  
  PMDsegments.gr <- segmentPMDs(
    m           = meth.gr,
    chr.sel     = chr.sel,
    pdfFilename = pdf_file,   # 不再是 NULL
    seqLengths  = seqLengths,
    num.cores   = num.cores
  )
  
  # 只保留 PMD 区域（type == "PMD"）
  if (!"type" %in% colnames(mcols(PMDsegments.gr))) {
    stop("PMDsegments.gr 中没有 'type' 列，请用 str(PMDsegments.gr) 看一下实际列名。")
  }
  PMD_only <- PMDsegments.gr[mcols(PMDsegments.gr)$type == "PMD"]
  
  if (length(PMD_only) == 0) {
    warning("样本 ", sample_name, " 没有检测到 PMD 区域。")
  }
  
  # 输出 BED：chr, start(0-based), end, nCG
  pmd_dt <- data.table(
    chr   = as.character(seqnames(PMD_only)),
    start = start(PMD_only) - 1L,   # BED 0-based
    end   = end(PMD_only),         # BED 右闭
    nCG   = mcols(PMD_only)$nCG
  )
  
  out_bed <- file.path(output_dir, paste0(sample_name, "_PMD.bed"))
  fwrite(pmd_dt, out_bed, sep = "\t", col.names = TRUE)
  message("PMD 结果已保存到: ", out_bed)
}

message("全部 12 个样本 PMD 计算完成。")


###############################   计算核心腺癌PMD  ###################################
suppressPackageStartupMessages({
  library(GenomicRanges)
  library(data.table)
})

# ========== 路径和样本名 ==========
pmd_dir   <- "/data3/xumy_PMD/data/EAC_GEJ"
out_file  <- file.path(pmd_dir, "core_PMD_all12.bed")

sample_names <- c(
  paste0("EAC_", c(1:4,6)),
  paste0("GEJ_", 1:7)
)

pmd_files <- file.path(pmd_dir, paste0(sample_names, "_PMD.bed"))

# 检查文件是否存在
missing <- !file.exists(pmd_files)
if (any(missing)) {
  stop("以下 PMD 文件未找到，请检查路径或名字：\n",
       paste(pmd_files[missing], collapse = "\n"))
}

# ========== 读入 BED -> GRanges ==========
# 你的 PMD.bed: chr, start(0-based), end, nCG
read_pmd_gr <- function(f) {
  message("读取 PMD 文件: ", f)
  dt <- fread(f)
  if (ncol(dt) < 3) {
    stop("PMD bed 至少应有 3 列: chr, start, end: ", f)
  }
  setnames(dt, 1:3, c("chr", "start", "end"))
  
  GRanges(
    seqnames = dt$chr,
    ranges   = IRanges(
      start = dt$start + 1L,   # BED 0-based -> GRanges 1-based
      end   = dt$end
    ),
    strand   = "*"
  )
}

gr_list <- lapply(pmd_files, read_pmd_gr)

# ========== 计算 12 个样本的交集（核心 PMD） ==========
message("计算 12 个样本的 PMD 交集（核心 PMD）...")

core_pmd <- Reduce(
  f = function(x, y) GenomicRanges::intersect(x, y, ignore.strand = TRUE),
  x = gr_list
)

# 可选：再 reduce 一下避免碎片重复（一般 intersect 已经是非重叠的，可留可不留）
core_pmd <- reduce(core_pmd, ignore.strand = TRUE)

if (length(core_pmd) == 0L) {
  warning("注意：12 个样本之间没有任何共同 PMD 区域，core_pmd 为空。")
}

# ========== 输出成 BED ==========
core_dt <- data.table(
  chr   = as.character(seqnames(core_pmd)),
  start = start(core_pmd) - 1L,   # 再转回 BED 0-based
  end   = end(core_pmd)
)

fwrite(core_dt, out_file, sep = "\t", col.names = TRUE)
message("核心 PMD 已保存到: ", out_file)
