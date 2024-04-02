globalVariables(c("marker_id", "real_emd", "result", "."))

myEMD <-  function(A, B, binSize = NULL) {
  stopifnot(is.numeric(A) & is.numeric(B))
  if (is.null(binSize)) 
    binSize <- 2 * stats::IQR(c(A[A != 0], B[B != 0]))/length(c(A[A != 
                                                                    0], B[B != 0]))^(1/3)
  emd <- tryCatch({
    bins <- seq(floor(min(c(A, B))), ceiling(max(c(A, B))), by = binSize)
    if (max(bins) < max(A, B)) 
      bins <- c(bins, bins[length(bins)] + binSize)
    histA <- graphics::hist(A, breaks = bins, plot = FALSE)
    histB <- graphics::hist(B, breaks = bins, plot = FALSE)
    densA <- histA$density
    densA <- densA/sum(densA)
    densB <- histB$density
    densB <- densB/sum(densB)
    return(CyEMD:::emdC(densA, densB))
  }, 
  error = function(e) {
    message(paste("Error in cluster", cluster_id, ":", e))
    return(0)
  })
  return(emd)
}

rowwiseEMD <- function(mat, condition, binSize = NULL) {
  stopifnot(is.matrix(mat), is.numeric(mat), nlevels(as.factor(condition)) == 
              2, ncol(mat) == length(condition))
  condition <- as.factor(condition)
  result <- apply(mat, 1, function(marker) {
    grouped <- split(marker, condition)
    myEMD(grouped[[1]], grouped[[2]], cluster_id=cluster_id)
  })
  out_dt <- data.table::as.data.table(result, keep.rownames = "marker_id")
  out_dt[, `:=`(marker_id, as.factor(marker_id))]
  out_dt
}

#' CyEMD
#'
#' Differential analysis method using the Earth Mover´s Distance to compare normalized distributions.
#' @param sce SingleCellExperiment containing expression data and metadata.
#' @param condition Name of vector specifying the comparison of interest.
#' @param binSize Bin width of histograms. If NULL (default), a flexible bin width is estimated by the Freedman-Diaconis rule evaluated on all non-zero values.
#' @param nperm Number of permutations. The default is 100.
#' @param assay Name of assay containing relevant data. The default is "exprs".
#' @param seed Seed initialization for reproducibility. The default is 1.
#' @param parallel Logical value indicating whether calculations should be performed in parallel. The default is FALSE.
#' @param replace Logical value indicating whether permutations should be performed with or without replacement. The default is FALSE.
#' @returns data.table containing marker IDs, EMDs and empirical p-values
#' @examples
#' path <- system.file("extdata", "pbmc/sce.rds", package = "CyEMD")
#' sce <- readRDS(path)
#' cyEMD(sce, condition = "condition")
#' @importFrom data.table :=
#' @export
cyEMD <- function(sce, condition, binSize=NULL, nperm=100, assay="exprs", seed=1, parallel=FALSE, replace=FALSE) {
  sceEI <- CATALYST::ei(sce)
  if(length(unique(sceEI[[condition]])) == 1){
    message("Only one condition detected. Returning NA.")
    empty_dt <- data.table::data.table(
      marker_id = rownames(rowData(sce)),
      emd = rep(NA, nrow(rowData(sce))),
      p_val = rep(NA, nrow(rowData(sce)))
    )
    return(empty_dt)
  }
  bppar <- BiocParallel::bpparam()
  if (!parallel) 
    bppar <- BiocParallel::SerialParam(progressbar = TRUE)
  set.seed(1)
  assay <- match.arg(assay, names(SummarizedExperiment::assays(sce)))
  data <- SummarizedExperiment::assay(sce, assay)
  emd_real <- rowwiseEMD(mat = data, condition = sce[[condition]], 
                         binSize = binSize, cluster_id=cluster_id)
  data.table::setnames(emd_real, "result", "real_emd")
  data.table::setkey(emd_real, marker_id)
  if(replace){
    allowed_perms <- RcppAlgos::permuteCount(sceEI[[condition]], n = nperm, 
                                              seed = seed, repetition = FALSE)
    if(allowed_perms < nperm){
      message(paste0("Allowed permutations (", 
                  allowed_perms, "=", length(sceEI[[condition]]), 
                  "!) <", nperm, ". Returning NA."))
      empty_dt <- data.table::data.table(
        marker_id = rownames(rowData(sce)),
        emd = rep(NA, nrow(rowData(sce))),
        p_val = rep(NA, nrow(rowData(sce)))
      )
      return(empty_dt)
    }
    perms <- RcppAlgos::permuteSample(sceEI[[condition]], n = nperm, 
                                      seed = seed, repetition = FALSE)
  }else{
    allowed_perms <- RcppAlgos::permuteCount(v = unique(sceEI[[condition]]), 
                                             m = length(sceEI[[condition]]),
                                             freqs = table(sceEI[[condition]]),
                                             n = nperm, 
                                             seed = seed)
    if(allowed_perms < nperm){
      message(paste0("Allowed permutations (",
                  allowed_perms, "=", length(sceEI[[condition]]), "!/(", 
                  sum(sceEI[[condition]] == unique(sceEI[[condition]])[1]), "! * ", sum(sceEI[[condition]] == unique(sceEI[[condition]])[2]), "!)",
                  ") <", nperm, ". Returning NA"))
      empty_dt <- data.table::data.table(
        marker_id = rownames(rowData(sce)),
        emd = rep(NA, nrow(rowData(sce))),
        p_val = rep(NA, nrow(rowData(sce)))
      )
      return(empty_dt)
    }
    perms <- RcppAlgos::permuteSample(v = unique(sceEI[[condition]]), 
                                      m = length(sceEI[[condition]]),
                                      freqs = table(sceEI[[condition]]),
                                      n = nperm, 
                                      seed = seed)
  }
  perm_res <- BiocParallel::bplapply(as.data.frame(t(unclass(perms))), 
                                     function(perm, sceEI, data, binSize) {
                                       condition_permutation_cells <- rep(perm, times = sceEI$n_cells)
                                       rowwiseEMD(mat = data, condition = condition_permutation_cells, 
                                                  binSize = binSize)
                                     }, sceEI, data, binSize, BPPARAM = bppar)
  all_perms <- data.table::rbindlist(perm_res, idcol = "permutation")
  data.table::setkey(all_perms, marker_id)
  res_agg <- all_perms[emd_real][, .(p_val = (sum(result >= 
                                                    real_emd) + 1)/(nperm + 1)), by = c("marker_id", "real_emd")]
  data.table::setnames(res_agg, "real_emd", "emd")
  return(res_agg)
}

