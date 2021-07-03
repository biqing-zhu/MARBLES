#' Run data simulation
#'
#' Simulation of complex scRNA-seq data based on the sampled underlying truth X
#'
#' @param x a \code{\link[SingleCellExperiment]{SingleCellExperiment}}.
#' @param nc number of cells to simulate.
#' @param ns number of samples to simulate.
#' @param nk number of cell types to simulate.  
#' @param ng number of genes to simulate. 
#' @param de_perc percentage of cell types to be DE when initializing X0.
#' @param lfc numeric value to use as mean logFC 
#'   (logarithm base 2) for DE genes.
#' @param c_c cell type relationship network matrix. 1 means connected, and 0 means not connected.
#' @param \eqn{\gamma},\eqn{\beta} \eqn{\Phi} used for sampling the state matrix X.
#' @param iter number of iterations when sampling X.
#' @return The estimated model parameters and the DE status
#' @return a \code{\link[SingleCellExperiment]{SingleCellExperiment}}
#'   containing multiple clusters & samples across two groups 
#'   as well as the following metadata: \describe{
#'   \item{cell metadata (\code{colData(.)})}{a \code{DataFrame} containing,
#'   for each cell, it's cluster, sample, and group ID.}
#'   \item{experiment metadata (\code{metadata(.)})}{
#'   \describe{
#'   \item{\code{experiment_info}}{a \code{data.frame} 
#'   summarizing the experimental design.}
#'   \item{\code{n_cells}}{the number of cells for each sample.}
#'   \item{\code{gene_info}}{a \code{data.frame} containing, for each gene
#'   in each cluster, it's differential distribution \code{category},
#'   mean \code{logFC} (NA for genes for categories "ee"),
#'   gene used as reference (\code{sim_gene}), dispersion \code{sim_disp},
#'   and simulation means for each group \code{sim_mean.A/B}.}
#'   \item{\code{ref_sids/kids}}{the sample/cluster IDs used as reference.}
#'   \item{\code{args}}{a list of the function call's input arguments.}}}}
#' @export
#' 
run_simulation <- function(x, nc = 2e3, ns = 3, nk, ng = nrow(x), 
                                  de_perc, lfc = 2, c_c, gamma, beta, iter = 5){
  
  ## Store all input arguments to be returned in final output
  args <- c(as.list(environment()))
  
  ## Check validity of input arguments
  .check_sce(x, req_group = FALSE)
  args_tmp <- .check_args_simData(as.list(environment()))
  nk <- args$nk <- args_tmp$nk
  
  ## Reference IDs
  nk0 <- length(kids0 <- purrr::set_names(levels(x$cluster_id)))
  ns0 <- length(sids0 <- purrr::set_names(levels(x$sample_id)))
  
  ## Simulation IDs
  nk <- length(kids <- purrr::set_names(paste0("cluster", seq_len(nk))))
  sids <- purrr::set_names(paste0("sample", seq_len(ns)))
  gids <- purrr::set_names(c("A", "B"))
  
  ## Sample reference clusters & samples
  ref_kids <- stats::setNames(sample(kids0, nk, nk > nk0), kids)
  ref_sids <- replicate(length(gids), 
                          sample(sids0, ns, ns > ns0))
  dimnames(ref_sids) <- list(sids, gids)
  
  rel_lfc <- rep(1, nk)
  names(rel_lfc) <- kids

  ## Initialize count matrix
  gs <- paste0("gene", seq_len(ng))
  cs <- paste0("cell", seq_len(nc))
  y <- matrix(0, ng, nc, dimnames = list(gs, cs))
  
  ## Sample cell metadata
  cd <- .sample_cell_md(
    n = nc, ids = list(kids, sids, gids))
  rownames(cd) <- cs
  cs_idx <- .split_cells(cd, by = colnames(cd))
  n_cs <- purrr::modify_depth(cs_idx, -1, length)
  
  ## Split input cells by cluster-sample
  cs_by_ks <- .split_cells(x)
  
  ## Sample nb. of genes to simulate per category & gene indices
  state_mat <- matrix(, 0, nk)  
  for (i in 1:ng) {
    
    state_mat <- rbind(state_mat, initialization(c_c, paraMRF = c(gamma, beta), de_perc, iter = iter)[[1]])
    
  }  
  state_mat[state_mat == 1] <- 'de'
  state_mat[state_mat == '0'] <- 'ee'
  colnames(state_mat) <- paste0('cluster', 1:nk)
  gs_idx <- .sample_gene_inds(gs, state_mat)
  n_dd <- table(factor(as.vector(state_mat), levels = c('ee', 'de')), rep(paste0('cluster', 1:nk), each = ng))
  n_dd <- n_dd[, paste0('cluster', 1:nk)]
  
  ## For ea. cluster, sample set of genes to simulate from
  gs_by_k <- stats::setNames(sample(rownames(x), ng, ng > nrow(x)), gs)
  gs_by_k <- replicate(nk, gs_by_k)
  colnames(gs_by_k) <- kids
  
  ## Split by cluster & categroy
  gs_by_k <- S4Vectors::split(gs_by_k, col(gs_by_k))
  gs_by_k <- stats::setNames(map(gs_by_k, purrr::set_names, gs), kids)
  gs_by_kc <- lapply(kids, function(k) 
    lapply(S4Vectors::unfactor(cats), function(c) 
      gs_by_k[[k]][gs_idx[[c, k]]])) 
  
  ## Sample logFCs
  lfc <- vapply(kids, function(k) 
    lapply(S4Vectors::unfactor(cats), function(c) { 
      
      n <- n_dd[c, k]
      if (c == "ee") return(rep(NA, n))
      signs <- sample(c(-1, 1), n, TRUE)
      lfcs <- stats::rgamma(n, 4, 4/lfc) * signs
      names(lfcs) <- gs_by_kc[[k]][[c]]
      lfcs * rel_lfc[k]
      
    }), vector("list", length(cats)))
  
  ## Compute NB parameters
  m <- lapply(sids0, function(s) {
    
    b <- paste0("beta.", s)
    b <- exp(rowData(x)[[b]])
    m <- outer(b, exp(x$offset), "*")
    dimnames(m) <- dimnames(x); m
    
  })
  d <- rowData(x)$dispersion 
  names(d) <- rownames(x)
  
  ## Initialize list of depth two to store 
  ## Simulation means in each cluster & group
  sim_mean <- lapply(kids, function(k) 
    lapply(gids, function(g) 
      stats::setNames(numeric(ng), gs)))
  
  ## Run simulation -----------------------------------------------------------
  for (k in kids) {
    
    for (s in sids) {
      
      ## Get reference samples, clusters & cells
      s0 <- ref_sids[s, ]
      k0 <- ref_kids[k]
      cs0 <- cs_by_ks[[k0]][s0]
      
      ## Get output cell indices
      ci <- cs_idx[[k]][[s]]
      
      for (c in cats[n_dd[, k] != 0]) {
        
        ## Sample cells to simulate from
        cs_g1 <- sample(cs0[[1]], n_cs[[k]][[s]][[1]], TRUE)
        cs_g2 <- sample(cs0[[2]], n_cs[[k]][[s]][[2]], TRUE)
        
        ## Get reference genes & output gene indices
        gs0 <- gs_by_kc[[k]][[c]] 
        gi <- gs_idx[[c, k]]
        
        ## Get NB parameters
        m_g1 <- m[[s0[[1]]]][gs0, cs_g1, drop = FALSE]
        m_g2 <- m[[s0[[2]]]][gs0, cs_g2, drop = FALSE]
        d_kc <- d[gs0]
        lfc_kc <- lfc[[c, k]]
        
        re <- .sim(c, cs_g1, cs_g2, m_g1, m_g2, d_kc, lfc_kc)
        y[gi, unlist(ci)] <- re$cs
        
        for (g in gids) sim_mean[[k]][[g]][gi] <- ifelse(
          is.null(re$ms[[g]]), NA, list(re$ms[[g]]))[[1]]
        
      }
    }
  }
  ## Construct gene metadata table storing ------------------------------------
  ## gene | cluster_id | category | logFC, gene, disp, mean used for sim.
  gi <- data.frame(
    gene = unlist(gs_idx),
    cluster_id = rep.int(rep(kids, each = length(cats)), c(n_dd)),
    category = rep.int(rep(cats, nk), c(n_dd)),
    logFC = unlist(lfc),
    sim_gene = unlist(gs_by_kc),
    sim_disp = d[unlist(gs_by_kc)]) %>% 
    dplyr::mutate_at("gene", as.character)
  ## Add true simulation means
  sim_mean <- sim_mean %>%
    map(bind_cols) %>% 
    bind_rows(.id = "cluster_id") %>% 
    mutate(gene = rep(gs, nk))
  
  gi <- full_join(gi, sim_mean, by = c("gene", "cluster_id")) %>% 
    dplyr::rename("sim_mean.A" = "A", "sim_mean.B" = "B") %>% 
    dplyr::mutate_at("cluster_id", factor)
  
  ## Reorder
  o <- order(as.numeric(gsub("[a-z]", "", gi$gene)))
  gi <- gi[o, ]; rownames(gi) <- NULL
  
  ## Construct SCE ------------------------------------------------------------
  ## Cell metadata including group, sample, cluster IDs
  cd$group_id <- droplevels(cd$group_id)
  cd$sample_id <- factor(paste(cd$sample_id, cd$group_id, sep = "."))
  m <- match(levels(cd$sample_id), cd$sample_id)
  gids <- cd$group_id[m]
  o <- order(gids)
  sids <- levels(cd$sample_id)[o]
  cd <- cd %>% 
    dplyr::mutate_at("cluster_id", factor, levels = kids) %>% 
    dplyr::mutate_at("sample_id", factor, levels = sids) 
  ## Simulation metadata including used reference samples/cluster, 
  ## List of input arguments, and simulated genes' metadata
  ei <- data.frame(sample_id = sids, group_id = gids[o])
  md <- list(
    experiment_info = ei,
    n_cells = table(cd$sample_id),
    gene_info = gi,
    ref_sids = ref_sids,
    ref_kids = ref_kids, 
    args = args)
  ## Return SCE 
  SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = as.matrix(y)),
    colData = cd, metadata = md)
}

## Get hidden states
initialization <- function(c_c, paraMRF, de_perc, iter = 5) {
  n_cell <- nrow(c_c)
  x <- rbinom(n_cell, 1, de_perc)
  x_trace <- matrix(, n_cell, iter + 1)
  x_trace[,1] <- x
  gamma <- paraMRF[1]; beta <- paraMRF[2]
  
  for(it in 1:iter) {
    
    for (c in 1:n_cell) {
      
      u0 <- sum(c_c[c, -c]) - c_c[c, -c] %*% x[-c]
      u1 <- c_c[c, -c] %*% x[-c]
      a <- gamma - beta * u0
      b <- -beta * u1
      prob <- exp(a) / (exp(a) + exp(b))
      x_trace[c, it + 1] <- x[c] <- (prob >= runif(1)) + 0
      
    }
  }
  return(list(x, x_trace))
}