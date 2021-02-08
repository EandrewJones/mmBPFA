// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

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

            // sample current features
            for (int j = 0; j < p; j++) {
                // number of features
                int m_jk = sum(new_zeros.col(k)) - new_zeros(j, k);

                // switch for sparsity inducing prior
                if (sparse == true) {
                    log_rp = log(m_jk) - log(p - m_jk - 1);
                }

                // Beta prior (constrain large values to avoid NaNs)
                log_p1 = 0.5 * (log(gamma_k(k)) - log(gam_ks(k))) +
                         (0.5 * gam_ks(k) * pow(mu_k(j), 2.0)) +
                         log_rp;
                if (log_p1 < -10) {
                    log_p1 = -10;
                } else if (log_p1 > 10) {
                    log_p1 = 10;
                }
                p1 = 1 / (1 + exp(-log_p1));

                // Bernoulli prior
                zero_jk = R::rbinom(1, p1);

                // sample lambda_jk according to zero_jk
                lambda_jk = zero_jk * R::rnorm(mu_k(j), sqrt(1 / gam_ks(k)));
                new_lambda(j, k) = lambda_jk;
                new_zeros(j, k) = zero_jk;
            }
        }

        // sample new features
        if (infinite == true) {
            // sample new possible new dimensions
            double prior_exp = (ibp_a * ibp_b) / (ibp_b + p - 1);
            arma::vec pois_draws = rpois(p, tau * prior_exp), K_new(p);
            for (int j = 0; j < p; j++)
            {
                if (R::runif(0, 1) < 0.9)
                {
                    K_new(j) = pois_draws(j);
                }
                else
                {
                    K_new(j) = 1;
                }
            }

            // initialize new feature matrices
            arma::mat add_lambda(p, max(K_new), arma::fill::zeros);
            arma::mat add_zeros(p, max(K_new), arma::fill::zeros);

            // get number of new samples for row j
            for (int j = 0; j < p; j++) {
                if (K_new(j) > 0)
                {
                    // calculate residuals for feature j
                    arma::vec resids_j = x.col(j) - trans(new_lambda.row(j) * omega.t()) +
                                         alpha_mat.col(j);

                    // calculate MH log likelihood
                    arma::mat I(K_new(j), K_new(j), arma::fill::eye);
                    arma::vec prop_lambda = Rcpp::rnorm(K_new(j), 0, 1);
                    arma::mat M = prop_lambda * prop_lambda.t() + I;
                    arma::mat mn = M.i() * prop_lambda * resids_j.t();
                    double log_al = (0.5 * n * K_new(j) * log(2 * M_PI)) +
                                    (0.5 * n * log(det(M))) -
                                    (0.5 * sum(diagvec((mn.t() * M) * mn)));


                    // calculate MH prior
                    double log_ap = R::dpois(K_new(j), prior_exp, TRUE) -
                                    R::dpois(K_new(j), tau * prior_exp, TRUE);

                    // MH acceptance rate
                    double log_acc = log_al + log_ap;

                    // MH step
                    if (log_acc > 0) {
                        // always accept
                        add_zeros(j, arma::span(0, K_new(j) - 1)) = arma::ones<arma::rowvec>(K_new(j));
                        add_lambda(j, arma::span(0, K_new(j) - 1)) = prop_lambda.t();
                    } else if (log_acc > -10) {
                        // MH draw if acceptance rate is not unreasonably low
                        if (R::runif(0, 1) < exp(log_acc)) {
                            add_zeros(j, arma::span(0, K_new(j) - 1)) = arma::ones<arma::rowvec>(K_new(j));
                            add_lambda(j, arma::span(0, K_new(j) - 1)) = prop_lambda.t();
                        }
                    }
                }
            }
            // concatenate new features to matrices
            new_lambda = arma::join_rows(new_lambda, add_lambda);
            new_zeros = arma::join_rows(new_zeros, add_zeros);

            // extend dc length by new features
            std::vector<int> add_dc(max(K_new), 0);
            dc.reserve(dc.size() + add_dc.size());
            dc.insert(dc.end(), add_dc.begin(), add_dc.end());
        }


        // Normalize lambda
        K = new_lambda.n_cols;
        arma::vec lambda_sum(p);
        for (int j = 0; j < p; j++) {
            double total = 0;
            for (int k = 0; k < K; k++) {
                total += pow(new_lambda(j, k), 2.0);
            }
            lambda_sum(j) = total;
        }
        arma::mat lambda_denom(p, K);
        lambda_denom.each_col() = sqrt(1 + lambda_sum);
        new_lambda = new_lambda / lambda_denom;

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