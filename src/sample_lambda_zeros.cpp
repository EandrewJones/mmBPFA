// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

// Log likelihood ratio of adding features with loadings `lambda` to a variable with
// residuals `e`, integrating out the new factor scores g_i ~ N(0, I):
// -N/2 log(1 + |lambda|^2) + |lambda|^2 / (2 (1 + |lambda|^2)) * sum(e^2)
double feature_loglik_ratio(const arma::vec& e, const arma::vec& lambda) {
    if (lambda.n_elem == 0) return 0;
    const double s = arma::dot(lambda, lambda);
    return -0.5 * e.n_elem * std::log1p(s) + 0.5 * s / (1 + s) * arma::dot(e, e);
}

// [[Rcpp::export(name = ".feature_loglik_ratio")]]
double feature_loglik_ratio_r(const arma::vec& e, const arma::vec& lambda) {
    return feature_loglik_ratio(e, lambda);
}

//' Sampling step for lambda and binary, sparsity-inducing zero matrix
//'
//' @keywords internal
//'
//' @param x A multivariate Gaussian copula matrix of dimenion N x P.
//' @param lambda A standard multivariate normal matrix of factor loadings of dimension P x K.
//' @param zeros A binary, sparsity inducing matrix of dimenion P x K.
//' @param omega A standard multivariate normal matrix of factor scores of dimenion N x K.
//' @param alpha A numeric vector of item-intercepts of length P.
//' @param gamma_k A numeric vector of factor precisions of length K.
//' @param dc A integer vector of counts for the number of times each dimensions has been sampled.
//' @param ibp_a A scalar double, Indian Buffet Process alpha parameter.
//' @param ibp_b A scalar double, Indian Buffet Process beta parameter.
//' @param tau A tunable hyperparameter for the poisson draws of new dishes to be sampled.
//' @param sparse A boolean indicating whether to include sparsity-inducing prior.
//' @param infinite A boolen indicating whether to sample new potential dishes at each step.
//'
//' @return zeros A binary, sparsity inducing matrix of dimenion P x K.
//' @return lambda A standard multivariate normal matrix of factor loadings of dimension P x K.
//' @return dc A integer vector of counts for the number of times each dimensions has been sampled.
//'
// [[Rcpp::export(name = "sample_lambda_and_zeros")]]
List sample_lambda_and_zeros(
    arma::mat x,
    arma::mat lambda,
    arma::mat zeros,
    arma::mat omega,
    arma::vec alpha,
    arma::vec gamma_k,
    std::vector<int> dc,
    const double ibp_a,
    const double ibp_b,
    const double tau,
    const bool sparse,
    const bool infinite
    ) {
        // get dimensions
        int n = omega.n_rows, p = lambda.n_rows, K = lambda.n_cols;

        // initialize variables
        arma::mat new_lambda(p, K), new_zeros(p, K), alpha_mat(n, p), resids(n, p);
        std::copy(lambda.begin(), lambda.end(), new_lambda.begin()); // duplicate lambda
        std::copy(zeros.begin(), zeros.end(), new_zeros.begin()); // duplicate zeros
        alpha_mat.each_row() = alpha.t();
        arma::vec gam_ks(K);
        arma::rowvec mu_k(p);
        double log_rp = 0, log_p1 = 0, p1 = 0, lambda_jk = 0, zero_jk = 0;

        // calculate gamma and colsums of sparsity matrix
        for (int k = 0; k < K; k++) {
            gam_ks(k) = arma::as_scalar(trans(omega.col(k)) * omega.col(k) + gamma_k(k));
        }

        // Loop col-row wise, faster
        for (int k = 0; k < K; k++)
        {
            // calculate residuals with current features set to 0
            arma::uvec idx(K, arma::fill::ones);
            idx(k) = 0;
            if (K == 1)
            {
                resids = x + alpha_mat;
            }
            else
            {
                resids = x - omega.cols(arma::find(idx == 1)) * new_lambda.cols(arma::find(idx == 1)).t() + alpha_mat;
            }

            // calculate mu for dimension k
            mu_k = (1 / gam_ks(k)) * (omega.col(k).t() * resids);

            // number of variables currently loading on feature k
            double m_k = sum(new_zeros.col(k));

            // sample current features
            for (int j = 0; j < p; j++) {
                // without the beta-process prior every loading is active
                if (sparse == false) {
                    new_lambda(j, k) = R::rnorm(mu_k(j), sqrt(1 / gam_ks(k)));
                    new_zeros(j, k) = 1;
                    continue;
                }

                // number of other variables loading on feature k
                int m_jk = m_k - new_zeros(j, k);

                // singletons are resampled by the new-feature step instead
                if (m_jk == 0) continue;

                // two-parameter IBP prior odds
                log_rp = log(m_jk) - log(p + ibp_b - 1 - m_jk);

                // log posterior odds of z_jk = 1
                log_p1 = 0.5 * (log(gamma_k(k)) - log(gam_ks(k))) +
                         (0.5 * gam_ks(k) * pow(mu_k(j), 2.0)) +
                         log_rp;
                // numerically stable inverse logit
                if (log_p1 >= 0) {
                    p1 = 1 / (1 + exp(-log_p1));
                } else {
                    p1 = exp(log_p1) / (1 + exp(log_p1));
                }

                // Bernoulli prior
                zero_jk = R::rbinom(1, p1);

                // sample lambda_jk according to zero_jk
                lambda_jk = zero_jk * R::rnorm(mu_k(j), sqrt(1 / gam_ks(k)));
                new_lambda(j, k) = lambda_jk;
                m_k += zero_jk - new_zeros(j, k);
                new_zeros(j, k) = zero_jk;
            }
        }

        // propose new features (Knowles and Ghahramani 2011): for each variable j,
        // replace its current singleton features with kappa ~ Pois(tau * nu) new ones
        if (infinite == true) {
            const double nu = (ibp_a * ibp_b) / (ibp_b + p - 1);
            const arma::rowvec m_k = arma::sum(new_zeros, 0);
            std::vector<int> new_rows;
            std::vector<double> new_vals;

            for (int j = 0; j < p; j++) {
                // current singletons of variable j
                arma::uvec single = arma::find((new_zeros.row(j) == 1) % (m_k == 1));
                arma::rowvec shared_lambda = new_lambda.row(j);
                arma::vec old_lambda = shared_lambda.elem(single);

                // residuals excluding the singletons
                shared_lambda.elem(single).zeros();
                arma::vec resids_j = x.col(j) - omega * shared_lambda.t() + alpha_mat.col(j);

                // proposal
                const int kappa_old = single.n_elem;
                const int kappa_new = R::rpois(tau * nu);
                arma::vec prop_lambda = Rcpp::as<arma::vec>(Rcpp::rnorm(kappa_new, 0, 1));

                double log_acc = feature_loglik_ratio(resids_j, prop_lambda) -
                                 feature_loglik_ratio(resids_j, old_lambda) +
                                 R::dpois(kappa_new, nu, true) - R::dpois(kappa_old, nu, true) +
                                 R::dpois(kappa_old, tau * nu, true) - R::dpois(kappa_new, tau * nu, true);

                if (log(R::runif(0, 1)) < log_acc) {
                    for (arma::uword c = 0; c < single.n_elem; c++) {
                        new_zeros(j, single(c)) = 0;
                        new_lambda(j, single(c)) = 0;
                    }
                    for (int c = 0; c < kappa_new; c++) {
                        new_rows.push_back(j);
                        new_vals.push_back(prop_lambda(c));
                    }
                }
            }

            // each accepted new feature gets its own column
            const int n_new = new_rows.size();
            arma::mat add_lambda(p, n_new, arma::fill::zeros);
            arma::mat add_zeros(p, n_new, arma::fill::zeros);
            for (int c = 0; c < n_new; c++) {
                add_zeros(new_rows[c], c) = 1;
                add_lambda(new_rows[c], c) = new_vals[c];
            }
            new_lambda = arma::join_rows(new_lambda, add_lambda);
            new_zeros = arma::join_rows(new_zeros, add_zeros);

            // extend dc length by new features
            dc.insert(dc.end(), n_new, 0);
        }

        // update dc
         for (int& i : dc)
             i += 1;

        // output values
        List output_list(3);
        output_list[0] = new_zeros;
        output_list[1] = new_lambda;
        output_list[2] = dc;

        return output_list;
}