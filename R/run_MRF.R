#' Run MRF
#'
#' Run the MRF model to estimate the DE state for each gene across each cell type.
#'
#' @param paraMRF Starting value of the model parameter \\Phi.
#' @param expr Pseudo-bulk expression list for a gene. Each item is a vector of the pseudo-bulk expression across all individuals for a cell type in a condition. The list names should start with 'Cond1_' or 'Cond2_', and followed by the name of the cell type. 
#' @param c_c Binary cell type relationship network matrix. 1 means connected, and 0 means not connected.
#' @param x_init Initial DE status. Should be a binary vector of length nk. 1 means DE, and 0 means EE.
#' @return The estimated model parameters and the DE status: \describe{
#' \item{x_mat}{Trace of the DE status for each cell type in the ICM algorithm. A 2-dimensional array: (num of cell types)\\cdot(num of iterations). The first column represents the initialization, and the last column represents the final DE status.}
#' \item{theta_mat}{Trace of the estimated model parameters \\Theta in the ICM algorithm. A 2-dimensional array: 2\\cdot(num of iterations). The first column represents the initialization, and the last column represents the final parameters.}
#' \item{phi_mat}{Trace of the estimated model parameters \\Phi in the ICM algorithm. A 2-dimensional array: 2\\cdot(num of iterations). The first column represents the initialization, and the last column represents the final parameters.}}
#' @export
#' 
run_MRF <- function(paraMRF, expr, c_c, x_init) {
  
  types_num <- length(expr) / 2
  types <- unique(substr(names(expr), 7, nchar(names(expr))))

  phi_mat <- matrix(paraMRF, nrow = 2, ncol = 1)
  rownames(phi_mat) <- c('gamma', 'beta')
  theta_mat <- matrix(c(1, 1), nrow = 2, ncol = 1)
  rownames(theta_mat) <- c('alpha', 'beta')
  x_mat <- matrix(x_init, nrow = types_num, ncol = 1)
  rownames(x_mat) <- types
  
  ## Number of patients
  mc <- lengths(expr)[grep('^Cond1', names(expr))]
  nc <- lengths(expr)[grep('^Cond2', names(expr))]
  
  ## Conditional density of y
  log_ly_sum <- function(par, x) {   
    
    alpha <- par[1]
    beta <- par[2]
    log_ly <- rep(NA, types_num)
      
    for(j in 1:types_num) {
        
      cond1 <- expr[[paste0('Cond1_', types[j])]]
      cond2 <- expr[[paste0('Cond2_', types[j])]]
      log_ly[j] <- x[j] * (2 * alpha * log(beta) + lgamma(sum(cond1) + alpha) + lgamma(sum(cond2) + alpha) - 
                             (2 * lgamma(alpha) + sum(lfactorial(cond1)) + sum(lfactorial(cond2)) + 
                                (sum(cond1) + alpha) * log(mc[j] + beta) + (sum(cond2) + alpha) * log(nc[j] + beta))) + 
                   (1 - x[j]) * (alpha * log(beta) + lgamma(sum(cond1) + sum(cond2) + alpha) - 
                             (lgamma(alpha) + sum(lfactorial(cond1)) + sum(lfactorial(cond2)) + 
                                (sum(cond2) + sum(cond1) + alpha) * log(mc[j] + nc[j] + beta)))
    }
    
    -sum(log_ly)
    
  }
  
  ## Iteratively update parameters  
  iter <- 0
  repeat{
      
    ## Estimate theta
    theta_mat <- cbind(theta_mat, optim(theta_mat[,ncol(theta_mat)], fn = log_ly_sum, x = x_mat[, ncol(x_mat)], lower = 1e-8, method = "L-BFGS-B")$par)
      
    ## Estimate phi
    phi_new <- tryCatch(optim(phi_mat[,ncol(phi_mat)], fn = lx, x = x_mat[, ncol(x_mat)], c_c = c_c, lower = c(-Inf, 1e-8), method = "L-BFGS-B")$par, error=function(e) NULL)
    if (is.null(phi_new)) {
      
      x_mat <- cbind(x_mat, 0)
      break
      
    } else {
      
      phi_mat <- cbind(phi_mat, phi_new)
      x_new <- x_mat[, ncol(x_mat)]
        
      ## Update x
      for(c in 1:types_num){

        cond1 <- expr[[paste0('Cond1_', types[c])]]
        cond2 <- expr[[paste0('Cond2_', types[c])]]
          
        log_fy_0 <- theta_mat[1, ncol(theta_mat)] * log(theta_mat[2, ncol(theta_mat)]) + lgamma(sum(cond1) + sum(cond2) + theta_mat[1, ncol(theta_mat)]) - 
                        (lgamma(theta_mat[1, ncol(theta_mat)]) + sum(lfactorial(cond1)) + sum(lfactorial(cond2)) + 
                           (sum(cond2) + sum(cond1) + theta_mat[1, ncol(theta_mat)]) * log(mc[c] + nc[c] + theta_mat[2, ncol(theta_mat)]))
        log_fy_1 <- 2 * theta_mat[1, ncol(theta_mat)] * log(theta_mat[2, ncol(theta_mat)]) + lgamma(sum(cond1) + theta_mat[1, ncol(theta_mat)]) + 
                        lgamma(sum(cond2) + theta_mat[1, ncol(theta_mat)]) - 
                          (2 * lgamma(theta_mat[1, ncol(theta_mat)]) + sum(lfactorial(cond1)) + sum(lfactorial(cond2)) + (sum(cond1) + theta_mat[1, ncol(theta_mat)]) * 
                             log(mc[c] + theta_mat[2, ncol(theta_mat)]) + (sum(cond2) + theta_mat[1, ncol(theta_mat)]) * log(nc[c] + theta_mat[2, ncol(theta_mat)]))
          
        u0 <- sum(c_c[c, -c]) - c_c[c, -c] %*% x_new[-c]
        u1 <- c_c[c, -c] %*% x_new[-c]
        a <- phi_mat[1, ncol(phi_mat)] - phi_mat[2, ncol(phi_mat)] * u0
        b <- -phi_mat[2, ncol(phi_mat)] * u1
        lprob <- a + log_fy_1 - matrixStats::logSumExp(c(a + log_fy_1, b + log_fy_0))
        x_new[c] <- (exp(lprob) >= runif(1)) + 0
        
      }
      
      x_mat <- cbind(x_mat, x_new)
       iter <- iter + 1
      if(((sum(x_mat[,ncol(x_mat)] != x_mat[,(ncol(x_mat) - 1)]) == 0) & (sum(theta_mat[,ncol(theta_mat)] - theta_mat[,(ncol(theta_mat) - 1)]) < 1e-4) & 
          (sum(phi_mat[,ncol(phi_mat)] - phi_mat[,(ncol(phi_mat) - 1)]) < 1e-4)) | iter > 100) {
        
        break
        
      }
    }  
  }
  
  results_list <- list()
  results_list$x_mat <- x_mat
  results_list$theta_mat <- theta_mat
  results_list$phi_mat <- phi_mat

  return(results_list)        
}


## x conditional likelihood
lx <- function(par, x, c_c) {  
  
  gamma <- par[1]
  beta <- par[2]
  log_lx <- c()
  for(j in 1:length(x)) {
    
    u1 <- c_c[j, -j] %*% x[-j]
    u0 <- sum(c_c[j, -j]) - c_c[j, -j] %*% x[-j]
    log_lx[j] <- (1 - x[j]) * (-beta * u1) + x[j] * (gamma - beta * u0) - matrixStats::logSumExp(c(-beta * u1, gamma - beta * u0))
    
  }
  
  -sum(log_lx)
}