source('./functions/utils.R')

library(extraDistr)
library(pscl)
library(MixedPoisson)
library(data.table)
library(purrr)
library(dplyr)
set.seed(2020)

# # Find proper paraMRF
# for(beta in seq(0.1, 1, 0.1)) {
#   # for(gamma in seq(-1,0.6, 0.2)) {
#     print('beta:')
#     print(beta)
#     # print('gamma:')
#     # print(gamma)
#     print(initialization(c_c, paraMRF = c(.3, 0.8))[[2]])
#   # }
# }
# x <- rbinom(n_cell, 1, .3)
# optim(c(0.2, 0.8), fn = initialization_optim, c_c = c_c_10, x = x, lower = c(-Inf, 1e-3), method = "L-BFGS-B")$par

# initialization_optim <- function(paraMRF, c_c, x) {
#   n_cell <- nrow(c_c)
#   # x <- rep(1, n_cell)
#   ## X0
#   # x <- rep(c(1, 0), each = 5)
#   # phi <- optim(phi_init, fn = lx, x = x, c_c = c_c, lower = c(-Inf, -Inf, 0), method = "L-BFGS-B")$par
#   iter <- 5
#   x_trace <- matrix(, n_cell, iter + 1)
#   x_trace[,1] <- x
#   gamma <- paraMRF[1]; beta <- paraMRF[2]
#   # gamma0 <- phi[1]; gamma1 <- phi[2]; beta <- phi[3]

#   for(it in 1:iter) {
#     for (c in 1:n_cell) {
#         u0 <- sum(c_c[c, -c]) - c_c[c, -c] %*% x[-c]
#         u1 <- c_c[c, -c] %*% x[-c]
#         a <- gamma - beta * u0
#         b <- -beta * u1
#         prob <- exp(a) / (exp(a) + exp(b))
#         x_trace[c, it + 1] <- x[c] <- (prob >= runif(1)) + 0
#     }
#   }
#   u0 <- u1 <- 0
#   for (c in 1:n_cell) {
#     if (x[c] == 0) {
#        u0 <- u0 + sum(c_c[c, -c]) - c_c[c, -c] %*% x[-c]
#     } else {
#        u1 <- u1 + c_c[c, -c] %*% x[-c]
#     }
#   }
#   return((u1 - u0)^2)
# }

names(cats) <- cats <- c("ee", "ep", "de", "dp", "dm", "db")
cats <- factor(cats, levels = cats)

## Get hidden states
initialization <- function(c_c, paraMRF, de_perc, iter = 5) {
  n_cell <- nrow(c_c)
  # x <- rep(1, n_cell)
  ## X0
  x <- rbinom(n_cell, 1, de_perc)
  # x <- rep(c(1, 0), each = 5)
  # phi <- optim(phi_init, fn = lx, x = x, c_c = c_c, lower = c(-Inf, -Inf, 0), method = "L-BFGS-B")$par
  x_trace <- matrix(, n_cell, iter + 1)
  x_trace[,1] <- x
  gamma <- paraMRF[1]; beta <- paraMRF[2]
  # gamma0 <- phi[1]; gamma1 <- phi[2]; beta <- phi[3]

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

simulation <- function(alpha_inf, beta_inf, alpha_disp, beta_disp, alpha_mean, beta_mean, taus, c_c, de_perc){
  state <- initialization(c_c, paraMRF = c(0, 2), de_perc)[[1]]
  types <- length(state)

  ## Sample ZINB parameters for simulation
  mu_hc_vec <- mu_pd_vec <- rgamma(types, alpha_mean, beta_mean)
  disp_vec <- rgamma(types, alpha_disp, beta_disp)
  inf_vec <- rbeta(types, alpha_inf, beta_inf)

  ## Get taus (up/down)
  expr <- list()
  tau <- rep(NA, types)
  results <- list()
  for (i in 1:types) {
    if (state[i] == 1) {
      tau[i] <- sample(taus, 1)
      mu_pd_vec[i] <- mu_hc_vec[i] * tau[i]
    }
  }

  ### Simulation
  ## HC
  for(i in 1:types) {
    repeat{
      attempt <- rzinb(300, size = disp_vec[i] , prob = disp_vec[i] / (disp_vec[i]  + mu_hc_vec[i]), pi = inf_vec[i])
      if(sum(attempt) != 0 & min(attempt) == 0) {
        expr[[i]] <- attempt
        break
      }
    }
  }
  ## PD
  for(i in 1:types) {
    repeat{
      attempt <- rzinb(300, size = disp_vec[i] , prob = disp_vec[i] / (disp_vec[i]  + mu_pd_vec[i]), pi = inf_vec[i])
      if(sum(attempt) != 0 & min(attempt) == 0) {
        expr[[i + types]] <- attempt
        break
      }
    } 
  }
  names(expr) <- c(paste0('HC_',1:types), paste0('PD_',1:types))
  results$expr <- expr
  results$tau <- tau
  results$state <- state
  return(results)
}

simulation_nb <- function(alpha_size, beta_size, alpha_mu, beta_mu, taus, c_c){
  state <- initialization(c_c, paraMRF = c(1, 1, 2))[[1]]
  types <- length(state)
  mu_hc_vec <- mu_pd_vec <- rgamma(types, alpha_mu, beta_mu)
  size_vec <- rgamma(types, alpha_size, beta_size)
  expr <- list()
  tau <- rep(NA, types)
  results <- list()
  for (i in 1:types) {
    if (state[i] == 1) {
      tau[i] <- sample(taus, 1)
      mu_pd_vec[i] <- mu_hc_vec[i] * tau[i]
    }
  }
  ## HC
  for(i in 1:types) {
    repeat{
      attempt <- rnbinom(1000, size = size_vec[i], mu = mu_hc_vec[i])
      if(sum(attempt) != 0 & min(attempt) == 0) {
        expr[[i]] <- attempt
        break
      }
    }
  }
  ## PD
  for(i in 1:types) {
    repeat{
      attempt <- rnbinom(1000, size = size_vec[i], mu = mu_pd_vec[i])
      if(sum(attempt) != 0 & min(attempt) == 0) {
        expr[[i + types]] <- attempt
        break
      }
    } 
  }
  names(expr) <- c(paste0('HC_',1:types), paste0('PD_',1:types))
  results$expr <- expr
  results$tau <- tau
  results$state <- state
  return(results)
}

simulation_pseudobulk <- function(x, nc = 2e3, ns = 3, nk,
    probs = NULL, de_perc, paired = FALSE,
    p_ep = 0.5, p_dp = 0.3, p_dm = 0.5,
    p_type = 0, lfc = 2, rel_lfc = NULL, 
    phylo_tree = NULL, phylo_pars = c(ifelse(is.null(phylo_tree), 0, 0.1), 3),
    ng = nrow(x), force = FALSE, c_c, gamma, beta, iter = 5){

  # throughout this code...
    # k: cluster ID
    # s: sample ID
    # g: group ID
    # c: DD category
    # 0: reference
    
    # store all input arguments to be returned in final output
    args <- c(as.list(environment()))
    
    # check validity of input arguments
    .check_sce(x, req_group = FALSE)
    args_tmp <- .check_args_simData(as.list(environment()))
    nk <- args$nk <- args_tmp$nk
    
    # reference IDs
    nk0 <- length(kids0 <- set_names(levels(x$cluster_id)))
    ns0 <- length(sids0 <- set_names(levels(x$sample_id)))
    
    # simulation IDs
    nk <- length(kids <- set_names(paste0("cluster", seq_len(nk))))
    sids <- set_names(paste0("sample", seq_len(ns)))
    gids <- set_names(c("A", "B"))

    # sample reference clusters & samples
    ref_kids <- setNames(sample(kids0, nk, nk > nk0), kids)
    if (paired) { 
        # use same set of reference samples for both groups
        ref_sids <- sample(sids0, ns, ns > ns0)
        ref_sids <- replicate(length(gids), ref_sids)
    } else {
        # draw reference samples at random for each group
        ref_sids <- replicate(length(gids), 
            sample(sids0, ns, ns > ns0))
    }
    dimnames(ref_sids) <- list(sids, gids)
    
    if (is.null(rel_lfc)) {
        rel_lfc <- rep(1, nk)
    }
    if (is.null(names(rel_lfc))) {
        names(rel_lfc) <- kids
    } else {
        stopifnot(names(rel_lfc) %in% kids0)
    }
    
    # initialize count matrix
    gs <- paste0("gene", seq_len(ng))
    cs <- paste0("cell", seq_len(nc))
    y <- matrix(0, ng, nc, dimnames = list(gs, cs))
    
    # sample cell metadata
    cd <- .sample_cell_md(
        n = nc, probs = probs,
        ids = list(kids, sids, gids))
    rownames(cd) <- cs
    cs_idx <- .split_cells(cd, by = colnames(cd))
    n_cs <- modify_depth(cs_idx, -1, length)
    
    # split input cells by cluster-sample
    cs_by_ks <- .split_cells(x)

    # sample nb. of genes to simulate per category & gene indices
    state_mat <- matrix(, 0, nk)  
    for (i in 1:ng) {
      state_mat <- rbind(state_mat, initialization(c_c, paraMRF = c(gamma, beta), de_perc, iter = 5)[[1]])
    }  
    state_mat[state_mat == 1] <- 'de'
    state_mat[state_mat == '0'] <- 'ee'
    colnames(state_mat) <- paste0('cluster', 1:nk)
    gs_idx <- .sample_gene_inds(gs, state_mat)
    n_dd <- table(factor(as.vector(state_mat), levels = c('ee', 'ep', 'de', 'dp', 'dm', 'db')), rep(paste0('cluster', 1:nk), each = ng))
    n_dd <- n_dd[, paste0('cluster', 1:nk)]
    
    # for ea. cluster, sample set of genes to simulate from
    gs_by_k <- setNames(sample(rownames(x), ng, ng > nrow(x)), gs)
    gs_by_k <- replicate(nk, gs_by_k)
    colnames(gs_by_k) <- kids

    # when 'phylo_tree' is specified, induce hierarchical cluster structure
    if (!is.null(phylo_tree)) {                                  
        res <- .impute_shared_type_genes(x, gs_by_k, gs_idx, phylo_tree, phylo_pars)
        gs_by_k <- res$gs_by_k
        class  <- res$class
        specs <- res$specs
    # otherwise, simply impute type-genes w/o phylogeny
    } else if (p_type != 0) {
        res <- .impute_type_genes(x, gs_by_k, gs_idx, p_type)
        stopifnot(!any(res$class == "shared"))
        gs_by_k <- res$gs_by_k
        class <- res$class
        specs <- res$specs
    } else {
        class <- rep("state", ng)
        specs <- rep(NA, ng)
        names(class) <- names(specs) <- gs
    }

    # split by cluster & categroy
    gs_by_k <- split(gs_by_k, col(gs_by_k))
    gs_by_k <- setNames(map(gs_by_k, set_names, gs), kids)
    gs_by_kc <- lapply(kids, function(k) 
        lapply(unfactor(cats), function(c) 
            gs_by_k[[k]][gs_idx[[c, k]]])) 
    
    # sample logFCs
    lfc <- vapply(kids, function(k) 
        lapply(unfactor(cats), function(c) { 
            n <- n_dd[c, k]
            if (c == "ee") return(rep(NA, n))
            signs <- sample(c(-1, 1), n, TRUE)
            lfcs <- rgamma(n, 4, 4/lfc) * signs
            names(lfcs) <- gs_by_kc[[k]][[c]]
            lfcs * rel_lfc[k]
        }), vector("list", length(cats)))

    # compute NB parameters
    m <- lapply(sids0, function(s) {
        b <- paste0("beta.", s)
        b <- exp(rowData(x)[[b]])
        m <- outer(b, exp(x$offset), "*")
        dimnames(m) <- dimnames(x); m
    })
    d <- rowData(x)$dispersion 
    names(d) <- rownames(x)
    
    # initialize list of depth two to store 
    # simulation means in each cluster & group
    sim_mean <- lapply(kids, function(k) 
        lapply(gids, function(g) 
            setNames(numeric(ng), gs)))
   
    # run simulation -----------------------------------------------------------
    for (k in kids) {
        for (s in sids) {
            # get reference samples, clusters & cells
            s0 <- ref_sids[s, ]
            k0 <- ref_kids[k]
            cs0 <- cs_by_ks[[k0]][s0]
            
            # get output cell indices
            ci <- cs_idx[[k]][[s]]
            
            for (c in cats[n_dd[, k] != 0]) {
                # sample cells to simulate from
                cs_g1 <- sample(cs0[[1]], n_cs[[k]][[s]][[1]], TRUE)
                cs_g2 <- sample(cs0[[2]], n_cs[[k]][[s]][[2]], TRUE)
                
                # get reference genes & output gene indices
                gs0 <- gs_by_kc[[k]][[c]] 
                gi <- gs_idx[[c, k]]
                
                # get NB parameters
                m_g1 <- m[[s0[[1]]]][gs0, cs_g1, drop = FALSE]
                m_g2 <- m[[s0[[2]]]][gs0, cs_g2, drop = FALSE]
                d_kc <- d[gs0]
                lfc_kc <- lfc[[c, k]]
                
                re <- .sim(c, cs_g1, cs_g2, m_g1, m_g2, d_kc, lfc_kc, p_ep, p_dp, p_dm)
                y[gi, unlist(ci)] <- re$cs
                
                for (g in gids) sim_mean[[k]][[g]][gi] <- ifelse(
                    is.null(re$ms[[g]]), NA, list(re$ms[[g]]))[[1]]
            }
        }
    }
    # construct gene metadata table storing ------------------------------------
    # gene | cluster_id | category | logFC, gene, disp, mean used for sim.
    gi <- data.frame(
        gene = unlist(gs_idx),
        cluster_id = rep.int(rep(kids, each = length(cats)), c(n_dd)),
        category = rep.int(rep(cats, nk), c(n_dd)),
        logFC = unlist(lfc),
        sim_gene = unlist(gs_by_kc),
        sim_disp = d[unlist(gs_by_kc)]) %>% 
        mutate_at("gene", as.character)
    # add true simulation means
    sim_mean <- sim_mean %>%
        map(bind_cols) %>% 
        bind_rows(.id = "cluster_id") %>% 
        mutate(gene = rep(gs, nk))
  
    gi <- full_join(gi, sim_mean, by = c("gene", "cluster_id")) %>% 
        dplyr::rename("sim_mean.A" = "A", "sim_mean.B" = "B") %>% 
        mutate_at("cluster_id", factor)

    # reorder
    o <- order(as.numeric(gsub("[a-z]", "", gi$gene)))
    gi <- gi[o, ]; rownames(gi) <- NULL
    
    # construct SCE ------------------------------------------------------------
    # cell metadata including group, sample, cluster IDs
    cd$group_id <- droplevels(cd$group_id)
    cd$sample_id <- factor(paste(cd$sample_id, cd$group_id, sep = "."))
    m <- match(levels(cd$sample_id), cd$sample_id)
    gids <- cd$group_id[m]
    o <- order(gids)
    sids <- levels(cd$sample_id)[o]
    cd <- cd %>% 
        mutate_at("cluster_id", factor, levels = kids) %>% 
        mutate_at("sample_id", factor, levels = sids) 
    # gene metadata storing gene classes & specificities
    rd <- DataFrame(class = factor(class, 
        levels = c("state", "shared", "type")))
    rd$specs <- as.list(specs)
    # simulation metadata including used reference samples/cluster, 
    # list of input arguments, and simulated genes' metadata
    ei <- data.frame(sample_id = sids, group_id = gids[o])
    md <- list(
        experiment_info = ei,
        n_cells = table(cd$sample_id),
        gene_info = gi,
        ref_sids = ref_sids,
        ref_kids = ref_kids, 
        args = args)
    # return SCE 
    SingleCellExperiment(
        assays = list(counts = as.matrix(y)),
        colData = cd, rowData = rd, metadata = md)
}