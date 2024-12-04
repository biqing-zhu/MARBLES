library(matrixStats)

## x conditional likelihood
lx <- function(par, x, c_c) {   
  gamma <- par[1]
  beta <- par[2]
  log_lx <- c()
  for(j in 1:length(x)) {
    # total_neigh <- sum(c_c[j, -j])
    u1 <- c_c[j, -j] %*% x[-j]
    u0 <- sum(c_c[j, -j]) - c_c[j, -j] %*% x[-j]
    log_lx[j] <- (1 - x[j]) * (-beta * u1) + x[j] * (gamma - beta * u0) - logSumExp(c(-beta * u1, gamma - beta * u0))
  }
  -sum(log_lx)
}

get_pst_prob <- function(lfy1, lfy0, paraMRF, states, iterGibbs = 20000, br = 10000, skip = 5, c_c){
  gamma <- paraMRF[1]; beta = paraMRF[2]
  count <- 0
  statessum <- states * 0
  # pfdr1_mat <- matrix(, nrow = 0, ncol = length(states))
  cat("\n")
  for (iter in 1:iterGibbs){
    if (iter%%200 == 0){
      cat("\r", "Posterior sampling,", floor(iter/iterGibbs*100), "%", "completed")
    }
    for (c in 1:dim(c_c)[1]) {
        u0 <- sum(c_c[c, -c]) - c_c[c, -c] %*% states[-c]
        u1 <- c_c[c, -c] %*% states[-c]
        a <- gamma - beta * u0
        b <- -beta * u1
        lprob <- a + lfy1[c] - logSumExp(c(a + lfy1[c], b + lfy0[c]))
        states[c] <- (exp(lprob) >= runif(1)) + 0
    }
    if (iter >= br && (iter%%skip)==0) {
      count <- count + 1
      statessum <- statessum + states
    }
    # pfdr1_mat <- rbind(pfdr1_mat, statessum / count)
  }
  pfdr1 <- statessum / count
  return(pfdr1)
}

icm <- function(paraMRF, expr, distr, c_c, x_init) {
    types_num <- length(expr) / 2
    types <- unique(substr(names(expr), 4, nchar(names(expr))))
    # iter <- 10
    phi_mat <- matrix(paraMRF, nrow = 2, ncol = 1)
    rownames(phi_mat) <- c('gamma', 'beta')
    x_mat <- matrix(x_init, nrow = types_num, ncol = 1)
    rownames(x_mat) <- types

    if(distr == 'ZINB') {
      ### ZINB
      theta_hc <- list()
      theta_pd <- list()
      theta <- list()

      for(i in 1:types_num) {
        theta_hc[[i]] <- rep(NA, 3)
        theta_pd[[i]] <- rep(NA, 3)
        theta[[i]] <- rep(NA, 3)

        # hc <- expr[paste0('HC', i),]
        # pd <- expr[paste0('PD', i),]
        hc <- expr[[paste0('HC_', i)]]
        pd <- expr[[paste0('PD_', i)]]

        # if(t.test(hc, pd)$p.value < threshold) {
        #   x_mat[i, 1] <- 1
        # }

        ## Estimate theta
        m_hc <- zeroinfl(hc ~ 1 | 1, dist = "negbin")
        m_pd <- zeroinfl(pd ~ 1 | 1, dist = "negbin")
        m <- zeroinfl(c(hc, pd) ~ 1 | 1, dist = "negbin")
        theta_hc[[i]][1] <- 1 / (exp(-m_hc$coefficients$zero) + 1)
        theta_hc[[i]][2] <- m_hc$theta
        theta_hc[[i]][3] <- m_hc$fitted.values[1]
        theta_pd[[i]][1] <- 1 / (exp(-m_pd$coefficients$zero) + 1)
        theta_pd[[i]][2] <- m_pd$theta
        theta_pd[[i]][3] <- m_pd$fitted.values[1]
        theta[[i]][1] <- 1 / (exp(-m$coefficients$zero) + 1)
        theta[[i]][2] <- m$theta
        theta[[i]][3] <- m$fitted.values[1]
      }

      iter <- 0
      repeat{
    
        ## Estimate phi
        phi_mat <- cbind(phi_mat, optim(phi_mat[,ncol(phi_mat)], fn = lx, x = x_mat[, ncol(x_mat)], c_c = c_c, lower = c(-Inf, -Inf, 1e-8), method = "L-BFGS-B")$par)

        x_new <- x_mat[, ncol(x_mat)]
        ## Update x
        for(j in 1:types_num){
          # hc <- expr[paste0('HC', j),]
          # pd <- expr[paste0('PD', j),]
          hc <- expr[[paste0('HC_', types[j])]]
          pd <- expr[[paste0('PD_', types[j])]]
          hc <- hc[hc < quantile(hc, 0.9)]
          pd <- pd[pd < quantile(pd, 0.9)]
          log_fy_x1 <- sum(dzinb(pd, size = theta_pd[[j]][2] , prob = theta_pd[[j]][2] / (theta_pd[[j]][2]  + theta_pd[[j]][3]), pi = theta_pd[[j]][1], log = TRUE) + dzinb(hc, size = theta_hc[[j]][2] , prob = theta_hc[[j]][2] / (theta_hc[[j]][2]  + theta_hc[[j]][3]), pi = theta_hc[[j]][1], log = TRUE))
          log_px_x1 <- phi_mat[2, ncol(phi_mat)] - phi_mat[3, ncol(phi_mat)] * (sum(c_c[j, -j]) - c_c[j, -j] %*% x_new[-j])
          sum_x1 <- log_fy_x1 + log_px_x1
          log_fy_x0 <- sum(dzinb(c(hc, pd), size = theta[[j]][2] , prob = theta[[j]][2] / (theta[[j]][2]  + theta[[j]][3]), pi = theta[[j]][1], log = TRUE))
          log_px_x0 <- phi_mat[1, ncol(phi_mat)] - phi_mat[3, ncol(phi_mat)] * (c_c[j, -j] %*% x_new[-j])
          sum_x0 <- log_fy_x0 + log_px_x0
          if(sum_x1 > sum_x0) {
            x_new[j] <- 1
          } else {
             x_new[j] <- 0
          }
          x_mat <- cbind(x_mat, x_new)
          iter <- iter + 1
         }
        if(((sum(x_mat[,ncol(x_mat)] != x_mat[,(ncol(x_mat) - 1)]) == 0) & (sum(phi_mat[,ncol(phi_mat)] - phi_mat[,(ncol(phi_mat) - 1)]) < 1e-4)) | iter > 100) {
              break  
        }
        print(phi_mat)
        print(x_mat) 
      }
    }

    ### Poisson Gamma
    else if (distr == 'PG') {

      # expr <- lapply(expr, function(vec) {vec[vec %in% as.numeric(names(table(vec))[table(vec) >= 4])]})
      mc <- lengths(expr)[grep('^HC', names(expr))]
      nc <- lengths(expr)[grep('^PD', names(expr))]
      theta_mat <- matrix(c(1, 1), nrow = 2, ncol = 1)
      rownames(theta_mat) <- c('alpha', 'beta')
      
      ## Conditional density of y
      log_ly_sum <- function(par, x) {   
        alpha <- par[1]
        beta <- par[2]
        log_ly <- rep(NA, types_num)
        
        for(j in 1:types_num) {

          hc <- expr[[paste0('HC_', types[j])]]
          pd <- expr[[paste0('PD_', types[j])]]
          log_ly[j] <- x[j] * (2 * alpha * log(beta) + lgamma(sum(hc) + alpha) + lgamma(sum(pd) + alpha) - (2 * lgamma(alpha) + sum(lfactorial(hc)) + sum(lfactorial(pd)) + (sum(hc) + alpha) * log(mc[j] + beta) + (sum(pd) + alpha) * log(nc[j] + beta))) + 
                      (1 - x[j]) * (alpha * log(beta) + lgamma(sum(hc) + sum(pd) + alpha) - (lgamma(alpha) + sum(lfactorial(hc)) + sum(lfactorial(pd)) + (sum(pd) + sum(hc) + alpha) * log(mc[j] + nc[j] + beta)))
        }
        # print('log_ly:')
        # print(log_ly)
        -sum(log_ly)
      }

      iter <- 0
      repeat{

        ## Estimate theta
        theta_mat <- cbind(theta_mat, optim(theta_mat[,ncol(theta_mat)], fn = log_ly_sum, x = x_mat[, ncol(x_mat)], lower = 1e-8, method = "L-BFGS-B")$par)

        ## Estimate phi
        # phi_mat <- cbind(phi_mat, optim(phi_mat[,ncol(phi_mat)], fn = lx, x = x_mat[, ncol(x_mat)], c_c = c_c, lower = c(-Inf, 1e-8), method = "L-BFGS-B")$par)
        phi_new <- tryCatch(optim(phi_mat[,ncol(phi_mat)], fn = lx, x = x_mat[, ncol(x_mat)], c_c = c_c, lower = c(-Inf, 1e-8), method = "L-BFGS-B")$par, error=function(e) NULL)
        if (is.null(phi_new)) {
           x_mat <- cbind(x_mat, 0)
           break
        } else {
           phi_mat <- cbind(phi_mat, phi_new)
           x_new <- x_mat[, ncol(x_mat)]
        
          ## Update x
          for(c in 1:types_num){
            # hc <- expr[paste0('HC', j),]
            # pd <- expr[paste0('PD', j),]
            hc <- expr[[paste0('HC_', types[c])]]
            pd <- expr[[paste0('PD_', types[c])]]

            # if (mean(hc == 0) > zero_perc & mean(pd == 0) > zero_perc) {
            #    x_new[j] <- 0
            # } else{
            log_fy_0 <- theta_mat[1, ncol(theta_mat)] * log(theta_mat[2, ncol(theta_mat)]) + lgamma(sum(hc) + sum(pd) + theta_mat[1, ncol(theta_mat)]) - (lgamma(theta_mat[1, ncol(theta_mat)]) + sum(lfactorial(hc)) + sum(lfactorial(pd)) + (sum(pd) + sum(hc) + theta_mat[1, ncol(theta_mat)]) * log(mc[c] + nc[c] + theta_mat[2, ncol(theta_mat)]))
            # log_px_0 <- -phi_mat[2, ncol(phi_mat)] * (c_c[j, -j] %*% x_new[-j])
            # sum_0 <- log_fy_0 + log_px_0
            log_fy_1 <- 2 * theta_mat[1, ncol(theta_mat)] * log(theta_mat[2, ncol(theta_mat)]) + lgamma(sum(hc) + theta_mat[1, ncol(theta_mat)]) + lgamma(sum(pd) + theta_mat[1, ncol(theta_mat)]) - (2 * lgamma(theta_mat[1, ncol(theta_mat)]) + sum(lfactorial(hc)) + sum(lfactorial(pd)) + (sum(hc) + theta_mat[1, ncol(theta_mat)]) * log(mc[c] + theta_mat[2, ncol(theta_mat)]) + (sum(pd) + theta_mat[1, ncol(theta_mat)]) * log(nc[c] + theta_mat[2, ncol(theta_mat)]))
          
            u0 <- sum(c_c[c, -c]) - c_c[c, -c] %*% x_new[-c]
            u1 <- c_c[c, -c] %*% x_new[-c]
            a <- phi_mat[1, ncol(phi_mat)] - phi_mat[2, ncol(phi_mat)] * u0
            b <- -phi_mat[2, ncol(phi_mat)] * u1
            lprob <- a + log_fy_1 - logSumExp(c(a + log_fy_1, b + log_fy_0))
            x_new[c] <- (exp(lprob) >= runif(1)) + 0

            # log_px_1 <- phi_mat[1, ncol(phi_mat)] - phi_mat[2, ncol(phi_mat)] * (sum(c_c[j, -j]) - c_c[j, -j] %*% x_new[-j])
            # sum_1 <- log_fy_1 + log_px_1
            # if(sum_1 > sum_0) {
            #   x_new[j] <- 1
            # } else {
            #   x_new[j] <- 0
            # }
            # }
          }
          x_mat <- cbind(x_mat, x_new)
          iter <- iter + 1
          if(((sum(x_mat[,ncol(x_mat)] != x_mat[,(ncol(x_mat) - 1)]) == 0) & (sum(theta_mat[,ncol(theta_mat)] - theta_mat[,(ncol(theta_mat) - 1)]) < 1e-4) & (sum(phi_mat[,ncol(phi_mat)] - phi_mat[,(ncol(phi_mat) - 1)]) < 1e-4)) | iter > 100) {
              break
          }
        }       
      }
      print(x_mat[,ncol(x_mat)])
      print(phi_mat[,ncol(phi_mat)]) 
      
      # ## Get final probabilities
      # lfy0 <- lfy1 <- c()
      # for(j in 1:types_num){
      #     hc <- expr[[paste0('HC_', types[j])]]
      #     pd <- expr[[paste0('PD_', types[j])]]
      #     lfy0 <- c(lfy0, theta_mat[1, ncol(theta_mat)] * log(theta_mat[2, ncol(theta_mat)]) + lgamma(sum(hc) + sum(pd) + theta_mat[1, ncol(theta_mat)]) - (lgamma(theta_mat[1, ncol(theta_mat)]) + sum(lfactorial(hc)) + sum(lfactorial(pd)) + (sum(pd) + sum(hc) + theta_mat[1, ncol(theta_mat)]) * log(mc[j] + nc[j] + theta_mat[2, ncol(theta_mat)])))
      #     lfy1 <- c(lfy1, 2 * theta_mat[1, ncol(theta_mat)] * log(theta_mat[2, ncol(theta_mat)]) + lgamma(sum(hc) + theta_mat[1, ncol(theta_mat)]) + lgamma(sum(pd) + theta_mat[1, ncol(theta_mat)]) - (2 * lgamma(theta_mat[1, ncol(theta_mat)]) + sum(lfactorial(hc)) + sum(lfactorial(pd)) + (sum(hc) + theta_mat[1, ncol(theta_mat)]) * log(mc[j] + theta_mat[2, ncol(theta_mat)]) + (sum(pd) + theta_mat[1, ncol(theta_mat)]) * log(nc[j] + theta_mat[2, ncol(theta_mat)])))
      # }
      # pst_prob_1 <- get_pst_prob(lfy1, lfy0, paraMRF = phi_mat[, ncol(phi_mat)], x_mat[,ncol(x_mat)], c_c = c_c)
    }
  results_list <- list()
  results_list$x_mat <- x_mat
  results_list$theta_mat <- theta_mat
  results_list$phi_mat <- phi_mat
  # results_list$pst_prob_1 <- pst_prob_1

  return(results_list)           
}           
