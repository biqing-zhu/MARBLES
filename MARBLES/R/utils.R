.sample_gene_inds <- function(gs, ns) {
  cluster_ids <- colnames(ns)
  vapply(cluster_ids, function(k)
    split(gs, factor(ns[, k], levels = c('ee', 'de'))),
    vector("list", length(cats)))
}

.check_sce <- function(x, req_group = TRUE) {
  stopifnot(is(x, "SingleCellExperiment"))
  stopifnot(c("cluster_id", "sample_id") %in% colnames(colData(x)))
  if (req_group)
    stopifnot("group_id" %in% colnames(colData(x)))
}

.check_args_simData <- function(u) {
  stopifnot(
    is.numeric(u$nc), length(u$nc) == 1, u$nc > 0, as.integer(u$nc) == u$nc,
    is.numeric(u$nk), length(u$nk) == 1, u$nk > 0, as.integer(u$nk) == u$nk,
    is.numeric(u$ns), length(u$ns) %in% c(1, 2), u$ns > 0, as.integer(u$ns) == u$ns,
    is.numeric(u$de_perc), length(u$de_perc) == 1, u$de_perc <= 1, u$de_perc >= 0,
    is.numeric(u$lfc), is.numeric(u$lfc), length(u$lfc) == 1, u$lfc >= 0,
    is.numeric(u$ng), length(u$ng) == 1, u$ng > 0, as.integer(u$ng) == u$ng
  )
  return(list(nk = u$nk))
}

.sample_cell_md <- function(n, ids) {
  ns <- vapply(ids, length, numeric(1))
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
  cat, cs_g1, cs_g2, m_g1, m_g2, d, lfc) {
  
  cat <- match.arg(cat)
  ng1 <- length(cs_g1)
  ng2 <- length(cs_g2)
  
  re <- switch(cat,
               "ee" = {
                 list(
                   .nb(cs_g1, d, m_g1),
                   .nb(cs_g2, d, m_g2))
               },
               "de" = {
                 list(
                   .nb(cs_g1, d, m_g1, -lfc), # lfc < 0 => all g1 hi
                   .nb(cs_g2, d, m_g2,  lfc)) # lfc > 0 => all g2 hi
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
  start_time <- Sys.time()
  res <- pbDS(pb, method = m, verbose = FALSE)
  end_time <- Sys.time()
  tbl <- resDS(sim, res)
  return(list(res = right_join(tbl, gi, by = c("gene", "cluster_id")), time = as.numeric(difftime(end_time, start_time), units="secs")))
}
