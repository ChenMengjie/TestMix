#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
using namespace Rcpp;
using namespace std;


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

arma::vec initialization_discrete(arma::vec YZ, int K){
  Environment myEnv("package:TestMix");
  Function initialization_discrete = myEnv["initialization_discrete"];
  Rcpp::NumericVector initialization_discrete_res = wrap(initialization_discrete(YZ, K));
  return initialization_discrete_res;
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

arma::vec log_poisson_likelihood_mix(arma::vec Y, double psi, double mu, int n, arma::vec log_factorial_Y){

  arma::vec likelihood = arma::zeros<arma::vec>(n);
  double common_term = -lgamma(psi) + psi*log(psi);

  for(int i = 0; i < n; ++i){
    double e_term = exp(mu);
    double psi_Yi = psi + Y(i);
    likelihood(i) = lgamma(psi_Yi) - log_factorial_Y(i) - psi_Yi*log(psi + e_term) + Y(i)*mu;
  }

  likelihood = likelihood + common_term;
  return(likelihood);
}


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

arma::vec gradient_all_mix(arma::vec Y, arma::vec Z,  double psi, double mu, arma::vec posterior_y, arma::vec posterior_z, int n1, int n2){

  arma::vec gradient = arma::zeros<arma::vec>(2);
  double exp_mu = exp(mu);
  double exp_mu_psi_inv = 1/(exp_mu + psi);
  double frac = exp_mu*exp_mu_psi_inv;

  arma::vec Y_psi = Y + psi;
  arma::vec Z_psi = Z + psi;
  arma::vec a1 = Y - Y_psi*frac;
  arma::vec a2 = Z - Z_psi*frac;
  gradient(0) = sum(a1%posterior_y) + sum(a2%posterior_z);

  double com_term = log(psi) + 1 + log(exp_mu_psi_inv) - R::digamma(psi);

  arma::vec b1 = arma::zeros<arma::vec>(n1);
  arma::vec b2 = arma::zeros<arma::vec>(n2);
  for(int i = 0; i < n1; ++i){
    b1(i) = R::digamma(Y_psi(i)) - Y_psi(i)*exp_mu_psi_inv + com_term;
  }
  for(int i = 0; i < n2; ++i){
    b2(i) = R::digamma(Z_psi(i)) - Z_psi(i)*exp_mu_psi_inv + com_term;
  }

  gradient(1) = sum(b1%posterior_y) + sum(b2%posterior_z);

  return gradient;
}


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

arma::vec gradient_all_combined_mix(arma::vec Y, double psi, double mu, arma::vec posterior, int n){

  arma::vec gradient = arma::zeros<arma::vec>(2);
  double exp_mu = exp(mu);
  double exp_mu_psi_inv = 1/(exp_mu + psi);
  double frac = exp_mu*exp_mu_psi_inv;

  arma::vec Y_psi = Y + psi;
  arma::vec a1 = Y - Y_psi*frac;
  gradient(0) = sum(a1%posterior);

  double com_term = log(psi) + 1 + log(exp_mu_psi_inv) - R::digamma(psi);

  arma::vec b1 = arma::zeros<arma::vec>(n);
  for(int i = 0; i < n; ++i){
    b1(i) = R::digamma(Y_psi(i)) - Y_psi(i)*exp_mu_psi_inv + com_term;
  }

  gradient(1) = sum(b1%posterior);

  return gradient;
}


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

double test_stepsize_for_psi_mix(arma::vec Y, double gra_psi, double ll, double psi, double mu, arma::vec posterior, int n, double gamma, double down, arma::vec log_factorial_Y){

  double gra_psi2 = gra_psi*gra_psi*gamma;
  double start = sqrt(abs(psi/gra_psi))/2;

  double aa = start;
  double selected = psi;
  while(aa > 0){
    double aa2 = aa*aa;
    double psi_prime = psi + aa2*gra_psi;
    if(psi_prime > 0){
      double lpsi_prime = sum(log_poisson_likelihood_mix(Y, psi, mu, n, log_factorial_Y));
      if(lpsi_prime - ll - aa2*gra_psi2 > 0){
        selected = psi_prime;
        break;
      }
    }
    aa = aa - start*down;
  }

  return selected;
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

double test_stepsize_for_mu_mix(arma::vec Y, double gra_mu, double ll, double psi, double mu, arma::vec posterior, int n, double gamma, double down, arma::vec log_factorial_Y){

  double gra_mu_2 = gra_mu*gra_mu*gamma;
  double start = sqrt(abs(mu/gra_mu))/2;

  double aa = start;
  double selected = mu;
  while(aa > 0){
    double aa2 = aa*aa;
    double mu1_prime = mu + aa2*gra_mu;
    double lmu1_prime = sum(log_poisson_likelihood_mix(Y, psi, mu, n, log_factorial_Y));
    if(lmu1_prime - ll - aa2*gra_mu_2 > 0){
      selected = mu1_prime;
      break;
    }

    aa = aa - start*down;
  }

  return selected;
}


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

arma::vec gradient_descent_mix(arma::vec Y, arma::mat posterior_mat, arma::vec psi_vec, arma::vec mu_vec,
                               arma::vec log_factorial_Y,int n, int K, int steps, double gamma, double down){

  arma::vec est = arma::zeros<arma::vec>(2*K - 2);

  for(int k = 0; k < K-1; ++k){

    double psi = psi_vec(k);
    double mu = mu_vec(k);
    arma::vec posterior = posterior_mat.col(k + 1);
    arma::vec gradient = gradient_all_combined_mix(Y, psi, mu, posterior, n);
    double ll = sum(log_poisson_likelihood_mix(Y, psi, mu, n, log_factorial_Y));

    double mu_prime = 0; double psi_prime = 0;

    for(int i = 0; i < steps; ++i){

      if(abs(gradient(0)) >= 0.00001){
        mu_prime = test_stepsize_for_mu_mix(Y, gradient(0), ll, psi, mu, posterior, n, gamma, down, log_factorial_Y);
      } else {
        mu_prime = mu;
      }
      if(abs(gradient(1)) >= 0.00001){
        psi_prime = test_stepsize_for_psi_mix(Y, gradient(1), ll, psi, mu, posterior, n, gamma, down, log_factorial_Y);
      } else {
        psi_prime = psi;
      }
      mu = mu_prime; psi = psi_prime;
      gradient = gradient_all_combined_mix(Y, psi, mu, posterior, n);
      ll = sum(log_poisson_likelihood_mix(Y, psi, mu, n, log_factorial_Y));

    }

    est(2*k) = mu;
    est(2*k + 1) = psi;
  }

   return est;
}


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

arma::rowvec take_exp_weight(arma::rowvec x){

  arma::rowvec lx = x;
  if(x(0) == 1){
    x(0) = 0;
    arma::rowvec exp_x = exp(x);
    lx = exp_x/sum(exp_x);
  } else {
    arma::rowvec exp_x = exp(x);
    exp_x(0) = 0;
    lx = exp_x/sum(exp_x);
  }
  return lx;
}


// [[Rcpp::export]]
double log_factorial(int Y){
  double res = 0;
  for(int kk = 1; kk <= Y; ++kk){
    res += log(kk);
  }
  return res;
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

arma::vec log_factorial_calculated(int N){

  arma::vec values = arma::zeros<arma::vec>(N+1);

  for(int kk = 1; kk <= N; ++kk){
    values(kk) = values(kk-1) + log(kk);
  }

  return values;
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

Rcpp::List EM_discrete_mix(arma::vec YZ, arma::vec psi_vec, int K, int n, int steps, int iter, double gamma, double down){

  arma::vec mu_vec = initialization_discrete(YZ, K);

  arma::vec calculated_values = log_factorial_calculated(YZ.max());
  arma::vec log_factorial_YZ = arma::zeros<arma::vec>(n);
  for(int i = 0; i < n; ++i){
    log_factorial_YZ(i) = calculated_values(YZ(i));
  }


  arma::mat posterior_mat = arma::zeros<arma::mat>(n, K);
  arma::uvec YZ_ind_zero = arma::find(YZ == 0);

  int n_0 = YZ_ind_zero.n_elem;
  for(int j = 0; j < n_0; ++j){
    int id = YZ_ind_zero(j);
    posterior_mat(id, 0) = 1;
  }

  for(int k = 0; k < K-1; ++k){
    double psi = psi_vec(k);
    double mu = mu_vec(k);
    posterior_mat.col(k + 1) = log_poisson_likelihood_mix(YZ, psi, mu, n, log_factorial_YZ);
  }

  for(int i = 0; i < n; ++i){
    arma::rowvec prob_vec = take_exp_weight(posterior_mat.row(i));
    posterior_mat.row(i) = prob_vec;
  }

  arma::vec est = arma::zeros<arma::vec>(2*K - 2);

  for(int ll = 0; ll < iter; ++ll){
    est = gradient_descent_mix(YZ, posterior_mat, psi_vec, mu_vec, log_factorial_YZ, n, K, steps, gamma, down);
    for(int j = 0; j < n_0; ++j){
      int id = YZ_ind_zero(j);
      posterior_mat(id, 0) = 1;
    }

    for(int k = 0; k < K-1; ++k){
      double mu = est(2*k);
      double psi = est(2*k + 1);
      posterior_mat.col(k + 1) = log_poisson_likelihood_mix(YZ, psi, mu, n, log_factorial_YZ);
    }

    for(int i = 0; i < n; ++i){
      arma::rowvec prob_vec = take_exp_weight(posterior_mat.row(i));
      posterior_mat.row(i) = prob_vec;
    }

  }

  return Rcpp::List::create(Rcpp::Named("posterior_mat") = posterior_mat,
                            Rcpp::Named("est") = est);
}






