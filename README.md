
### Installation

**TestMix** relies on the following R packages: **Rcpp**, **RcppArmadillo**. All packagess are hosted on CRAN. 
  ```R
  install.packages ("Rcpp")
  install.packages ("RcppArmadillo")
  ```

**scPoissonGamma** can be installed from github directly as follows:

  ```R
  install.packages ("devtools")
  library(devtools)
  install_github("ChenMengjie/TestMix")
  ```
  
  
### Simulation example 1: counts

```R
mu <- log(c(5, 50, 200))
psi <- c(10, 5, 10)


prop1 <- c(0.2, 0.4, 0, 0.4)
prop2 <- c(0.3, 0, 0.3, 0.4)
n1 <- 200
n2 <- 300

set.seed(500)
id1 <- sample(c(1:4), n1, prob = prop1, replace = T)
id2 <- sample(c(1:4), n2, prob = prop2, replace = T)

Y <- NULL
for(i in 1:length(id1)){
	if(id1[i] == 1){
		Y <- c(Y, 0)
	} else {
	  E <- rgamma(1, psi[id1[i]-1], psi[id1[i]-1])
	 mu_prime <- exp(mu[id1[i]-1])
	  Y <- c(Y, rpois(1, E*mu_prime))
	}
}

Z <- NULL
for(j in 1:length(id2)){
	if(id2[j] == 1){
		Z <- c(Z, 0)
	} else {
	  E <- rgamma(1, psi[id2[j]-1], psi[id2[j]-1])
	  mu_prime <- exp(mu[id2[j]-1])
	  Z <- c(Z, rpois(1, E*mu_prime))
	}
}


K <- 4
Y <- round(Y)
Z <- round(Z)

library(TestMix)

logY <- log(Y+1)
logZ <- log(Z+1)


test1 <- TestMix_counts(Y, Z)  #### input count data 
test2 <- TestMix_continuous(logY, logZ)  #### input log transformed data 

 ```

### Simulation example 2: continuous

```R
alpha <- 0.3
beta <- 1
mu <- c(1, 3, 5)
omega <- c(0.5, 0.5, 0.5)


prop1 <- c(0.2, 0.4, 0, 0.4)
prop2 <- c(0.3, 0, 0.3, 0.4)
n1 <- 200
n2 <- 300

set.seed(500)
id1 <- sample(c(1:4), n1, prob = prop1, replace = T)
id2 <- sample(c(1:4), n2, prob = prop2, replace = T)

Y <- NULL
for(i in 1:length(id1)){
	if(id1[i] == 1){
		Y <- c(Y, rgamma(1, 0.3))
	} else {
		Y <- c(Y, rnorm(1, mu[id1[i] - 1], omega[id1[i] - 1]))
	}
}

Z <- NULL
for(i in 1:length(id2)){
	if(id2[i] == 1){
		Z <- c(Z, rgamma(1, 0.3))
	} else {
		Z <- c(Z, rnorm(1, mu[id2[i] - 1], omega[id2[i] - 1]))
	}
}


library(TestMix)
Y <- Y[Y>0]
Z <- Z[Z>0]
test3 <- TestMix_continuous(Y, Z)  #### input log transformed data 

 ```


  
### Author

**Mengjie Chen** (UChicago)
