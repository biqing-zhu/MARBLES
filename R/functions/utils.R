.sample_gene_inds <- function(gs, ns) {
  cluster_ids <- colnames(ns)
  vapply(cluster_ids, function(k)
    split(gs, factor(ns[, k], levels = c('ee', 'ep', 'de', 'dp', 'dm', 'db'))),
    vector("list", length(cats)))
}

.check_sce <- function(x, req_group = TRUE) {
    stopifnot(is(x, "SingleCellExperiment"))
    stopifnot(c("cluster_id", "sample_id") %in% colnames(colData(x)))
    if (req_group)
        stopifnot("group_id" %in% colnames(colData(x)))
}

.check_args_simData <- function(u) {
    if (!u$force && u$ng != nrow(u$x))
        stop("Number of simulated genes should match with reference,\n", 
            "  but 'ng != nrow(x)'; please specify 'force = TRUE' if\n", 
            "  simulation should be forced regardlessly (see '?simData').")
    if (!is.null(u$phylo_tree) && u$p_type != 0)
        stop("Only one of arguments 'p_type' or 'phylo_tree'\n",
            "  can be specified; see '?simData' for 'Details'.")
    # assure number of simulated clusters matches with specified phylogeny
    if (!is.null(u$phylo_tree)) {
        kids_phylo <- .get_clusters_from_phylo(u$phylo_tree)
        nk_phylo <- length(kids_phylo)
        ns_phylo <- as.numeric(gsub("[a-z]", "", kids_phylo))
        if (!all(sort(ns_phylo) == seq_len(nk_phylo)))
            stop("Some clusters appear to be missing from 'phylo_tree';\n",
                "  please make sure all clusters up to ", 
                dQuote(kids_phylo[which.max(ns_phylo)]), " are present.")
        # possibly update number of clusters 'nk'
        if (nk_phylo != u$nk) u$nk <- nk_phylo
    }
    stopifnot(
        is.numeric(u$nc), length(u$nc) == 1, u$nc > 0, as.integer(u$nc) == u$nc,
        is.numeric(u$nk), length(u$nk) == 1, u$nk > 0, as.integer(u$nk) == u$nk,
        is.numeric(u$ns), length(u$ns) %in% c(1, 2), u$ns > 0, as.integer(u$ns) == u$ns,
        is.numeric(u$de_perc), length(u$de_perc) == 1, u$de_perc <= 1, u$de_perc >= 0,
        is.logical(u$paired), length(u$paired) == 1,
        is.numeric(u$p_ep), length(u$p_ep) == 1, u$p_ep > 0, u$p_ep < 1,
        is.numeric(u$p_dp), length(u$p_dp) == 1, u$p_dp > 0, u$p_dp < 1,
        is.numeric(u$p_dm), length(u$p_dm) == 1, u$p_dm > 0, u$p_dm < 1,
        is.numeric(u$p_type), length(u$p_type) == 1, u$p_type >= 0, u$p_type <= 1,
        is.numeric(u$lfc), is.numeric(u$lfc), length(u$lfc) == 1, u$lfc >= 0,
        is.numeric(u$ng), length(u$ng) == 1, u$ng > 0, as.integer(u$ng) == u$ng,
        is.logical(u$force), length(u$force) == 1,
        is.numeric(u$phylo_pars), length(u$phylo_pars) == 2, u$phylo_pars >= 0)
    if (!is.null(u$rel_lfc))
        stopifnot(is.numeric(u$rel_lfc), 
            length(u$rel_lfc) == u$nk, u$rel_lfc >= 0)
    return(list(nk = u$nk))
}

.check_arg_assay <- function(x, y) {
    stopifnot(is.character(y), length(y) == 1, y %in% assayNames(x))
    if (sum(assayNames(x) == y) > 1)
        stop("Argument 'assay' was matched to multiple times.\n ", 
            " Please assure that the input SCE has unique 'assayNames'.")
}

.check_args_aggData <- function(u) {
    stopifnot(is.character(u$by), length(u$by) <= 2, 
        u$by %in% colnames(colData(u$x)))
    stopifnot(is.logical(u$scale), length(u$scale) == 1)
    if (u$scale & (!u$assay %in% c("cpm", "CPM") | u$fun != "sum"))
        stop("Option 'scale = TRUE' only valid for", 
            " 'assay = \"cpm/CPM\"' and 'fun = \"sum\"'.")
}

.check_args_pbDS <- function(u) {
    if (!is.null(u$design))
        stopifnot(is.matrix(u$design),
            !is.null(rownames(u$design)),
            !is.null(colnames(u$design)))
    stopifnot(
        is.null(u$contrast) | is.matrix(u$contrast),
        is.null(u$coef) | is.numeric(unlist(u$coef)),
        is.numeric(u$min_cells), length(u$min_cells) == 1,
        is.logical(u$verbose), length(u$verbose) == 1,
        is.logical(u$treat), length(u$treat) == 1)
}

.sample_cell_md <- function(n, ids, probs = NULL) {
    ns <- vapply(ids, length, numeric(1))
    if (is.null(probs)) 
        probs <- vector("list", 3)
    probs <- lapply(seq_along(probs), function(i) {
        if (!is.null(probs[[i]])) {
            return(probs[[i]])
        } else {
            rep(1 / ns[i], ns[i])
        }
    })
    cd <- vapply(seq_along(probs), function(i) 
        sample(ids[[i]], n, TRUE, probs[[i]]), 
        character(n))
    cd <- data.frame(cd, row.names = NULL)
    colnames(cd) <- c("cluster_id", "sample_id", "group_id")
    cd$group_id <- factor(cd$group_id, levels = ids[[3]])
    return(cd)
}

.split_cells <- function(x, 
    by = c("cluster_id", "sample_id")) {
    if (is(x, "SingleCellExperiment"))
        x <- colData(x)
    cd <- data.frame(x[by], check.names = FALSE)
    cd <- data.table(cd, cell = rownames(x)) %>% 
        split(by = by, sorted = TRUE, flatten = FALSE)
    map_depth(cd, length(by), "cell")
}

.sim <- function(
    cat = c("ee", "ep", "de", "dp", "dm", "db"),
    cs_g1, cs_g2, m_g1, m_g2, d, lfc, ep, dp, dm) {
    
    cat <- match.arg(cat)
    ng1 <- length(cs_g1)
    ng2 <- length(cs_g2)
    
    re <- switch(cat,
        "ee" = {
            list(
                .nb(cs_g1, d, m_g1),
                .nb(cs_g2, d, m_g2))
        },
        "ep" = {
            g1_hi <- sample(ng1, round(ng1 * ep))
            g2_hi <- sample(ng2, round(ng2 * ep))
            list(
                .nb(cs_g1[-g1_hi], d, m_g1),
                .nb(cs_g1[ g1_hi], d, m_g1, lfc), # 50% g1 hi
                .nb(cs_g2[-g2_hi], d, m_g2),
                .nb(cs_g2[ g2_hi], d, m_g2, lfc)) # 50% g2 hi
        },
        "de" = {
            list(
                .nb(cs_g1, d, m_g1, -lfc), # lfc < 0 => all g1 hi
                .nb(cs_g2, d, m_g2,  lfc)) # lfc > 0 => all g2 hi
        },
        "dp" = {
            props <- sample(c(dp, 1 - dp), 2)
            g1_hi <- sample(ng1, round(ng1 * props[1]))
            g2_hi <- sample(ng2, round(ng2 * props[2]))
            list(                           
                .nb(cs_g1[-g1_hi], d, m_g1), 
                .nb(cs_g1[ g1_hi], d, m_g1,  lfc), # lfc > 0 => dp/(1-dp)% up
                .nb(cs_g2[-g2_hi], d, m_g2), 
                .nb(cs_g2[ g2_hi], d, m_g2, -lfc)) # lfc < 0 => (1-dp)/dp% up
        },
        "dm" = {
            g1_hi <- sample(ng1, round(ng1 * dm))
            g2_hi <- sample(ng2, round(ng2 * dm))
            list(
                .nb(cs_g1[-g1_hi], d, m_g1),
                .nb(cs_g1[ g1_hi], d, m_g1, -lfc), # lfc < 0 => 50% g1 hi
                .nb(cs_g2[-g2_hi], d, m_g2),
                .nb(cs_g2[ g2_hi], d, m_g2,  lfc)) # lfc > 0 => 50% g2 hi
        }, 
        "db" = {
            if (sample(c(TRUE, FALSE), 1)) {
                # all g1 mi, 50% g2 hi
                g2_hi <- sample(ng2, round(ng2 * 0.5))
                list(
                    .nb(cs_g1, d, m_g1, abs(lfc), 0.5),
                    .nb(cs_g2[-g2_hi], d, m_g2, -lfc), 
                    .nb(cs_g2[ g2_hi], d, m_g2,  lfc)) 
            } else {
                # all g2 mi, 50% g1 hi
                g1_hi <- sample(ng1, round(ng1 * 0.5))
                list(
                    .nb(cs_g2, d, m_g2, abs(lfc), 0.5), 
                    .nb(cs_g1[-g1_hi], d, m_g1, -lfc),  
                    .nb(cs_g1[ g1_hi], d, m_g1,  lfc))  
            }
        }
    )
    cs <- map(re, "counts")
    cs <- do.call("cbind", cs)
    ms <- map(re, "means")
    rmv <- vapply(ms, is.null, logical(1))
    ms <- ms[!rmv] %>% 
        map_depth(2, mean) %>% 
        map_depth(1, unlist) %>% 
        data.frame %>% 
        as.matrix
    ms <- switch(cat, 
        ee = ms,
        de = ms,
        db = if (ng2 == 0) {
            as.matrix(
                ms[, 1])
        } else {
            cbind(
                ms[, 1],
                rowMeans(ms[, c(2, 3)]))
        }, if (ng2 == 0) {
            as.matrix(
                rowMeans(ms[, c(1, 2)]))
        } else {
            cbind(
                rowMeans(ms[, c(1, 2)]),
                rowMeans(ms[, c(3, 4)]))
        })
    ms <- split(ms, col(ms))
    names(ms) <- c("A", "B")[c(ng1, ng2) != 0]
    list(cs = cs, ms = ms)
}

.nb <- function(cs, d, m, lfc = NULL, f = 1) {
    n_gs <- length(d)
    n_cs <- length(cs)
    if (is.null(lfc)) {
        lfc <- rep(0, n_gs)
    } else {
        lfc[lfc < 0] <- 0
    }
    fc <- f * (2 ^ lfc)
    fc <- rep(fc, each = n_cs)
    ds <- rep(1/d, each = n_cs)
    ms <- c(t(m[, cs])) * fc 
    y <- rnbinom(n_gs * n_cs, size = ds, mu = ms)
    y <- matrix(y, byrow = TRUE, 
        nrow = n_gs, ncol = n_cs, 
        dimnames = list(names(d), cs))
    ms <- split(ms, rep(seq_len(nrow(m)), each = n_cs))
    list(counts = y, means = ms)
}

.run_method <- function(m) {
#   start_time <- Sys.time()
  res <- pbDS(pb, method = m, verbose = FALSE)
#   end_time <- Sys.time()
  tbl <- resDS(sim, res)
  return(right_join(tbl, gi, by = c("gene", "cluster_id")))
}

aggregateData <- function(x, 
    assay = NULL, by = c("cluster_id", "sample_id"), 
    fun = c("sum", "mean", "median", "prop.detected", "num.detected"), 
    scale = FALSE, verbose = TRUE, BPPARAM = SerialParam(progressbar = verbose)) {
    
    # check validity of input arguments
    fun <- match.arg(fun)
    if (is.null(assay)) 
        assay <- assayNames(x)[1] 
    .check_arg_assay(x, assay)
    .check_args_aggData(as.list(environment()))
    stopifnot(is(BPPARAM, "BiocParallelParam"))
    
    # assure 'by' colData columns are factors
    # so that missing combinations aren't dropped
    for (i in by) 
        if (!is.factor(x[[i]])) 
            x[[i]] <- factor(x[[i]])
    
    # compute pseudo-bulks
    pb <- .pb(x, by, assay, fun, BPPARAM)
    if (scale & length(by) == 2) {
        # compute library sizes
        cs <- if (assay == "counts" && fun == "sum")
            pb else .pb(x, by, "counts", "sum", BPPARAM)
        ls <- lapply(cs, colSums)
        # scale pseudobulks by CPM
        pb <- lapply(seq_along(pb), function(i) pb[[i]] / 1e6 * ls[[i]])
        names(pb) <- names(ls)
    }
    
    # construct SCE
    md <- metadata(x)
    md$agg_pars <- list(assay = assay, by = by, fun = fun, scale = scale)
    pb <- SingleCellExperiment(pb, rowData = rowData(x), metadata = md)
    
    # tabulate number of cells
    cd <- data.frame(colData(x)[, by])
    for (i in names(cd))
        if (is.factor(cd[[i]]))
            cd[[i]] <- droplevels(cd[[i]])
    ns <- table(cd)
    if (length(by) == 2) {
        ns <- asplit(ns, 2)
        ns <- map(ns, ~c(unclass(.)))
    } else ns <- c(unclass(ns))
    int_colData(pb)$n_cells <- ns
    
    # propagate 'colData' columns that are unique across 2nd 'by'
    if (length(by) == 2) {
        cd <- colData(x)
        ids <- colnames(pb)
        counts <- vapply(ids, function(u) {
            m <- as.logical(match(cd[, by[2]], u, nomatch = 0))
            vapply(cd[m, ], function(u) length(unique(u)), numeric(1))
        }, numeric(ncol(colData(x))))
        cd_keep <- apply(counts, 1, function(u) all(u == 1))
        cd_keep <- setdiff(names(which(cd_keep)), by)
        if (length(cd_keep) != 0) {
            m <- match(ids, cd[, by[2]], nomatch = 0)
            cd <- cd[m, cd_keep, drop = FALSE]
            rownames(cd) <- ids
            colData(pb) <- cd
        }
    }
    return(pb)
}
                         
Execute_BPSC <- function(object, coef = 2, BPSC.parallel = TRUE,...){

  object_BPSC <- as.matrix(SummarizedExperiment::assay(object, "normcounts"))

  controllds <- which(object$label == levels(factor(object$label)[1]))
  design <- model.matrix(~ object$label)
  resbp <- BPSC::BPglm(data = object_BPSC, controlIds = controllds,
                       design = design, coef = coef, estIntPar = FALSE, useParallel = BPSC.parallel,...)
  FDR <- p.adjust(resbp$PVAL, method = "BH")
  result_BPSC <- list(gene_names = names(resbp$PVAL), pvalue = resbp$PVAL, FDR = FDR)
  return(result_BPSC)
}      
                         
Execute_DEsingle <- function(object, DEsingle.parallel = TRUE){

  object_DEsingle <- as.matrix(SummarizedExperiment::assay(object, "counts"))
  results_DEsingle <- DEsingle::DEsingle(counts = object_DEsingle,
                                         group = factor(object$label),
                                         parallel = DEsingle.parallel)
  result_DEsingle_DE <- list(gene_names = row.names(results_DEsingle), pvalue = results_DEsingle$pvalue,
                             FDR = results_DEsingle$pvalue.adj.FDR)
  return(result_DEsingle_DE)
}

Execute_monocle <- function(object, monocle.cores = 1,...){
  options(warn = -1)
  objects_monocle <- as.matrix(SummarizedExperiment::assay(object, "counts"))


  counts <-  floor(objects_monocle)
  object_monocle <- monocle::newCellDataSet(as.matrix(counts),
                                            phenoData = new("AnnotatedDataFrame",
                                                            data = data.frame(condition = object$label,
                                                                              row.names = colnames(counts))),
                                            lowerDetectionLimit = 0.5,  expressionFamily = VGAM::negbinomial.size())


  object_monocle <- BiocGenerics::estimateSizeFactors(object_monocle)
  object_monocle <- BiocGenerics::estimateDispersions(object_monocle)

  result_monocle <- monocle::differentialGeneTest(object_monocle, fullModelFormulaStr = "~ condition", cores = monocle.cores,...)


  result_monocle <- list(gene_names = rownames(result_monocle),
                         pvalue = result_monocle$pval,
                         FDR = result_monocle$qval)
  return(result_monocle)
}
                         
Execute_scDD <- function(object, alpha1 = 0.01, mu0 = 0, s0 = 0.01, a0 = 0.01, b0 = 0.01, scDD.permutation = 0, ...){

  object_scDD <- as.matrix(SummarizedExperiment::assay(object, "normcounts"))

  ### creat scDD input
  condition <- object$label
  names(condition) <- colnames(object)

  SDSumExp <- SingleCellExperiment::SingleCellExperiment(assays = list(normcounts = object_scDD),
                                   colData = data.frame(condition))

  prior_param= list(alpha = alpha1, mu0 = mu0, s0 = s0, a0 = a0, b0 = b0)
  scDatExSim <- scDD::scDD(SDSumExp, prior_param=prior_param, condition = "condition", min.nonzero = NULL, min.size = 3, permutations = scDD.permutation, ...)
  res <- scDD::results(scDatExSim)
  result_scDD <- list(gene_names = as.vector(res$gene),
                      pvalue = res$combined.pvalue,
                      FDR = res$combined.pvalue.adj)
  return(result_scDD)
}    
                         
Execute_Ttest <- function(object){

  cpm <- as.matrix(SummarizedExperiment::assay(object, "normcounts"))
  logcpm <- log2(cpm + 1)
  groups <- factor(object$label)
  idx <- seq_len(nrow(logcpm))
  names(idx) <- rownames(logcpm)
  ttest_p <- sapply(idx, function(i){
    t.test(logcpm[i,]~ groups)$p.value
  })
  FDR <- p.adjust(ttest_p, method = "BH")
  result_Ttest <- list(gene_names = names(ttest_p),
                       pvalue = ttest_p,
                       FDR = FDR)
  return(result_Ttest)
}  
                         
Execute_Wilcoxon <- function(object){

  cpm <- as.matrix(SummarizedExperiment::assay(object, "normcounts"))
  groups <- factor(object$label)
  idx <- 1:nrow(cpm)
  names(idx) <- rownames(cpm)
  wilcox_p <- sapply(idx, function(i){
    wilcox.test(cpm[i,]~ groups)$p.value
  })
  FDR <- p.adjust(wilcox_p, method = "BH")
  result_Wilcoxon <- list(gene_names = names(wilcox_p),
                          pvalue = wilcox_p,
                          FDR = FDR)
  return(result_Wilcoxon)
}
                         
Execute_Seurat <- function(object,  method.Seurat = "bimod"){

  object_Seurat <- as.matrix(SummarizedExperiment::assay(object, "counts"))

  tmp <- object$label
  names(tmp) <- colnames(object_Seurat)
  meta.data <- data.frame(groups = tmp)

  Seurat.input <- Seurat::CreateSeuratObject(counts = object_Seurat, meta.data = meta.data)
  Seurat.input <- Seurat::NormalizeData(object = Seurat.input)
  res <- Seurat::FindMarkers(object = Seurat.input, ident.1 = levels(as.factor(tmp))[1],
                             ident.2 = levels(as.factor(tmp))[2], group.by = 'groups',
                             logfc.threshold = -Inf, test.use = method.Seurat,
                             only.pos = FALSE, verbose = FALSE)
  results_Seurat <- list(gene_names = row.names(res),
                         pvalue = res$p_val,
                         FDR = res$p_val_adj)

  return(results_Seurat)
}
                         
Execute_zingeR.edgeR <- function(object, maxit.EM = 200){

  object_zingerR.edgeR <- as.matrix(SummarizedExperiment::assay(object, "counts"))

  dge <- edgeR::DGEList(round(object_zingerR.edgeR), group = factor(object$label))
  dge <- edgeR::calcNormFactors(dge)
  groups <- object$label
  design <- stats::model.matrix(~groups)
  weights <- zingeR::zeroWeightsLS(counts=dge$counts, design=design, maxit = maxit.EM, normalization="TMM")
  dge$weights <- weights

  dge = edgeR::estimateDisp(dge, design)
  fit = edgeR::glmFit(dge, design)
  lrt = zingeR::glmWeightedF(fit, coef=2, independentFiltering = TRUE)

  result_zingerR.edgeR <- list(gene_names = rownames(lrt$table),
                               pvalue = lrt$table$PValue,
                               FDR = lrt$table$padjFilter)
}                         