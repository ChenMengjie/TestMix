initialization_discrete <- function(YZ, K){

  logYZ <- log(YZ+1)
  ini.est <- kmeans(logYZ, centers = K)
  alpha <- min(ini.est$centers)
  aa <- which.min(ini.est$centers)
  mu <- ini.est$centers[-aa]
  return(mu)

}

combinatorial_hypothesis <- function(id1, id2, K){

  n1 <- length(id1)
  n2 <- length(id2)
  n12 <- n1 + n2
  coding <- NULL
  dis <- 2^K/2
  for(i in 1:K){
    coding <- cbind(coding, rep(c(0, 1), each = dis))
    dis <- dis/2
  }

  take_nonzero <- function(x){
    if(is.na(x) | is.infinite(x)) return(0)
    else return(x)
  }

  likelihoods <- NULL
  dis <- 2^K/2
  for(i in 1:K){
    n1_prime <- length(id1[id1==i])
    n2_prime <- length(id2[id2==i])
    likelihood1 <- take_nonzero(log(n1_prime/n1))*n1_prime + take_nonzero(log(n2_prime/n2))*n2_prime
    likelihood2 <- take_nonzero(log((n1_prime + n2_prime)/n12))*(n1_prime + n2_prime)
    likelihoods <- cbind(likelihoods, rep(c(likelihood1, likelihood2), each = dis))
    dis <- dis/2
  }
  likelihoods_sum  <- apply(likelihoods, 1, sum)
  hypo <- list(coding = coding, likelihoods = likelihoods_sum)
  return(hypo)
}

test_one_para <- function(hypo_coding, likelihoods, ID){
  flag <- hypo_coding[, ID] == 1
  test_statistics <- max(likelihoods[!flag]) - max(likelihoods[flag])
  p_value <- 1 - pchisq(2*test_statistics, 1)
  return(p_value)
}

test_at_least_one_nonzero_para <- function(hypo_coding, likelihoods){
  null.hypo <- apply(hypo_coding[, -1], 1, function(y){all( y == 1)})
  test_statistics <- max(likelihoods[!null.hypo]) - max(likelihoods[null.hypo])
  p_value <- 1 - pchisq(2*test_statistics, 1)
  return(p_value)
}

test_at_least_one_para <- function(hypo_coding, likelihoods){
  null.hypo <- apply(hypo_coding, 1, function(y){all(y == 1)})
  test_statistics <-  max(likelihoods[!null.hypo]) - max(likelihoods[null.hypo])
  p_value <- 1 - pchisq(2*test_statistics, 1)
  return(p_value)
}


estimate_parameters_discrete <- function(Y, Z, K, psi = 10, iter = 50, steps = 50, gamma = 0.6, down = 0.05){

  YZ <- c(Y, Z)
  psi_vec <- rep(psi, K-1)
  n <- length(YZ)
  EM_res <- EM_discrete_mix(YZ, psi_vec, K, n, steps, iter, gamma, down)

  n1 <- length(Y)
  posterior_mat <- EM_res$posterior_mat
  posterior_y_mat <- posterior_mat[1:n1, ]
  posterior_z_mat <- posterior_mat[(n1+1):n, ]

  idy <- apply(posterior_y_mat, 1, which.max)
  idz <- apply(posterior_z_mat, 1, which.max)

  delta_y <- table(idy)
  delta_z <- table(idz)
  delta1 <- delta_y/sum(delta_y)
  delta2 <- delta_z/sum(delta_z)

  est <- EM_res$est
  mu_vec <- est[seq(1, 2*K-2, by=2)]
  psi_vec <- est[seq(2, 2*K-2, by=2)]
  res <- list(mu = mu_vec, psi = psi_vec, idy = idy, idz = idz, delta1 = delta1, delta2 = delta2)

  return(res)
}


joint_inference <- function(Y, Z, K, psi = 10, iter = 10, steps = 25, gamma = 0.6, down = 0.05){

  est <- estimate_parameters_discrete(Y, Z, K, psi, iter, steps, gamma, down)
  idy <- est$idy
  idz <- est$idz

  hypo <- combinatorial_hypothesis(idy, idz, K)
  hypo_coding <- hypo$coding
  likelihoods <- hypo$likelihoods

  pvalue_vector <- NULL
  pvalue_vector <- c(pvalue_vector, test_at_least_one_para(hypo_coding, likelihoods))
  pvalue_vector <- c(pvalue_vector, test_at_least_one_nonzero_para(hypo_coding, likelihoods))

  for(i in 1:K){
    pvalue_vector <- c(pvalue_vector, test_one_para(hypo_coding, likelihoods, i))
  }

  names(pvalue_vector) <- c("AtleastOne", "AtleastOneNonzero", "Zero", paste0("Prop", 1:(K-1)))
  res <- list(mu = est$mu, psi = est$psi, idy = est$idy, idz = est$idz, delta1 = est$delta1, delta2 = est$delta2,
              pvalue_vector = pvalue_vector, K = K)

  return(res)
}



initialization <- function(Y, Z, K){

  YZ <- c(Y, Z)
  ini.est <- kmeans(YZ, centers = K)
  alpha <- min(ini.est$centers)
  aa <- which.min(ini.est$centers)
  mu <- ini.est$centers[-aa]
  groups <- c(1:K)[-aa]

  omega <- NULL
  for(i in groups){
    flag <- which(ini.est$cluster == i)
    omega <- c(omega, mean((YZ[flag]-ini.est$centers[i])^2))
  }

  n1 <- length(Y)
  clusterY <- ini.est$cluster[1:length(Y)]
  posterior_y <- clusterY == aa
  for(i in groups){
    posterior_y <- rbind(posterior_y, clusterY == i)
  }
  posterior_y <- apply(posterior_y, 2, as.numeric)

  clusterZ <- ini.est$cluster[(n1+1):length(YZ)]
  posterior_z <- clusterZ == aa
  for(i in groups){
    posterior_z <- rbind(posterior_z, clusterZ == i)
  }
  posterior_z <- apply(posterior_z, 2, as.numeric)

  res <- list(alpha = alpha, mu = mu, omega = omega, posterior_y = posterior_y, posterior_z = posterior_z)

  return(res)
}

estimate_parameters_continuous <- function(K, Y, Z, iter = 25){

  res <- initialization(Y, Z, K)
  alpha <- res$alpha
  mu <- res$mu
  omega <- res$omega
  posterior_y <- res$posterior_y
  posterior_z <- res$posterior_z
  beta <- 1

  logY <- log(Y+0.1)
  logZ <- log(Z+0.1)
  smallY <- ifelse(Y == 0, 0.001, Y)
  smallZ <- ifelse(Z == 0, 0.001, Z)

  for(j in 1:iter){

    for(i in 2:K){
      mu[i-1] <- (sum(posterior_y[i, ]*Y) + sum(posterior_z[i, ]*Z))/(sum(posterior_y[i, ]) + sum(posterior_z[i, ]))
      omega[i-1] <- (sum(posterior_y[i, ]*(Y - mu[i-1])^2) + sum(posterior_z[i, ]*(Z - mu[i-1])^2))/(sum(posterior_y[i, ]) + sum(posterior_z[i, ]))
    }

    digmma_alpha <- (sum(posterior_y[1, ]*(logY + log(beta))) + sum(posterior_z[1, ]*(logZ + log(beta))))/(sum(posterior_y[1, ]) + sum(posterior_z[1, ]))

    fr <- function(x){
      abs(digamma(x) - digmma_alpha)
    }

    if(alpha <= 2){
      alpha_prime <- optim(0.1, fr, method = "Brent", lower = 0, upper = 2)$par
      beta_prime <- alpha_prime*(sum(posterior_y[1, ]) + sum(posterior_z[1, ]))/(sum(posterior_y[1, ]*Y) + sum(posterior_z[1, ]*Z))
      if(beta_prime <= 500) {
        beta <- beta_prime
        alpha <- alpha_prime
      }
    }

    #alpha <- optim(0.1, fr, method = "Brent", lower = 0, upper = 2)$par

    posterior_y[1, ] <- dgamma(smallY, alpha, beta)
    posterior_z[1, ] <- dgamma(smallZ, alpha, beta)

    for(i in 2:K){
      posterior_y[i, ] <- dnorm(Y, mu[i-1], sqrt(omega[i-1]))
      posterior_z[i, ] <- dnorm(Z, mu[i-1], sqrt(omega[i-1]))
    }

    posterior_y <- apply(posterior_y, 2, function(x){x/sum(x)})
    posterior_z <- apply(posterior_z, 2, function(x){x/sum(x)})

  }

  idy <- apply(posterior_y, 2, which.max)
  idz <- apply(posterior_z, 2, which.max)
  delta_y <- table(idy)
  delta_z <- table(idz)
  delta1 <- delta_y/sum(delta_y)
  delta2 <- delta_z/sum(delta_z)

  res <- list(alpha = alpha, beta = beta, mu = mu, omega = omega,
              idy = idy, idz = idz, delta1 = delta1, delta2 = delta2)

  return(res)
}


loglikelihood_given_mixture_continuous <- function(Y, Z, idy, idz, alpha, beta, mu, omega){

  n1 <- length(Y)
  n2 <- length(Z)
  smallY <- ifelse(Y == 0, 0.001, Y)
  smallZ <- ifelse(Z == 0, 0.001, Z)

  all_likelihood <- 0
  for(i in 1:n1){
    id <- idy[i]
    if(id == 1){
      all_likelihood <- all_likelihood + dgamma(smallY[i], alpha, beta, log = TRUE)
    } else {
      all_likelihood <- all_likelihood + dnorm(Y[i], mu[id-1], sqrt(omega[id-1]), log = TRUE)
    }
  }

  for(j in 1:n2){
    id <- idz[j]
    if(id == 1){
      all_likelihood <- all_likelihood + dgamma(smallZ[j], alpha, beta, log = TRUE)
    } else {
      all_likelihood <- all_likelihood + dnorm(Z[j], mu[id-1], sqrt(omega[id-1]), log = TRUE)
    }
  }

  totalK <- length(mu)*2 + 2
  AIC <- 2*totalK - 2*all_likelihood
  BIC <- log(n1 + n2)*totalK - 2*all_likelihood
  est <- list(AIC = AIC, BIC = BIC, likelihood = all_likelihood)
  return(est)
}

model_selection <- function(logY, logZ, krange = c(2, 6), iter = 50, mindis = 0.2, min.prop = 0.1){

  fromK <- min(krange)
  toK <- max(krange)

  BIC.list <- NULL
  AIC.list <- NULL
  li.list <- NULL
  obs.prop.list <- NULL
  obs.min.mu.list <- NULL
  gamma.list <- NULL
  for(k in fromK:toK){
    est <- estimate_parameters_continuous(k, logY, logZ, iter)
    res <- loglikelihood_given_mixture_continuous(logY, logZ, est$idy, est$idz, est$alpha, est$beta, est$mu, est$omega)
    if(k == 2){
      obs.min.mu.list <- c(obs.min.mu.list, 10)
    } else {
      obs.min.mu.list <- c(obs.min.mu.list, min(diff(sort(est$mu))))
    }
    aa <- table(c(est$idy, est$idz))
    obs.prop.list <- c(obs.prop.list, min(aa/sum(aa)))
    gamma.list <- c(gamma.list, est$alpha/est$beta)
    BIC.list <- c(BIC.list, res$BIC)
    AIC.list <- c(AIC.list, res$AIC)
    li.list <- c(li.list, res$likelihood)
  }

  flag <- which(obs.prop.list >= min.prop & obs.min.mu.list >= mindis)
  K <- max(flag)

  return(K)
}


TestMix_counts <- function(Y, Z, krange = c(2, 6), iter = 50, iter.poisson = 20, mindis = 0.2, min.prop = 0.1, psi = 10, steps = 25, gamma = 0.6, down = 0.05){
  Y <- round(Y)
  Z <- round(Z)
  logY <- log(Y+1)
  logZ <- log(Z+1)
  testK <- model_selection(logY, logZ, krange = krange, iter = iter, mindis = mindis, min.prop = min.prop)
  res <- joint_inference(Y, Z, testK, psi = psi, iter = iter.poisson, steps = steps, gamma = gamma, down = down)
  return(res)
}


joint_inference_continuous <- function(Y, Z, K, iter = 50){

  est <- estimate_parameters_continuous(K, Y, Z, iter = iter)
  idy <- est$idy
  idz <- est$idz
  hypo <- combinatorial_hypothesis(idy, idz, K)
  hypo_coding <- hypo$coding
  likelihoods <- hypo$likelihoods

  pvalue_vector <- NULL
  pvalue_vector <- c(pvalue_vector, test_at_least_one_para(hypo_coding, likelihoods))
  pvalue_vector <- c(pvalue_vector, test_at_least_one_nonzero_para(hypo_coding, likelihoods))

  for(i in 1:K){
    pvalue_vector <- c(pvalue_vector, test_one_para(hypo_coding, likelihoods, i))
  }

  names(pvalue_vector) <- c("AtleastOne", "AtleastOneNonzero", "Zero", paste0("Prop", 1:(K-1)))

  res <- list(alpha = est$alpha, beta = est$beta, mu = est$mu, omega = est$omega,
              idy = est$idy, idz = est$idz, delta1 = est$delta1, delta2 = est$delta2,
              pvalue_vector = pvalue_vector, K = K)

  return(res)
}


TestMix_continuous <- function(Y, Z, krange = c(2, 6), iter = 50, mindis = 0.2, min.prop = 0.1){
  testK <- model_selection(Y, Z, krange = krange, iter = iter, mindis = mindis, min.prop = min.prop)
  res <- joint_inference_continuous(Y, Z, testK, iter = iter)
  return(res)
}


