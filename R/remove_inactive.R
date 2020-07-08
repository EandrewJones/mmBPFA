# Remove inactive features
# 
# @param zeros P x K sparsity-inducing binary matrix.
# @param lambda P x K matrix of factor loadings.
# @param omega N x K matrix of factor scores.
# @param gamma_k A numeric vector factor precisions of length K.
# @param dc A numeric vector counting number of times each dimension has been sampled.
# 
# @return List of zeros, lambda, omega, gamma_k, dc.
# 
remove_inactive <- function(
    zeros,
    lambda,
    omega,
    gamma_k,
    dc
    ) {
    # get colsumns of sparsity-inducing matrix
    m_k <- colSums(zeros)

    # drop inactive features
    drop <- which(m_k == 0)

    zeros <- zeros[, -drop]
    lambda <- lambda[, -drop]
    omega <- omega[, -drop]
    gamma_k <- gamma_k[-drop]
    dc <- dc[-drop]

    # store output
    output_list <- list()
    output_list[["zeros"]] <- as.matrix(zeros)
    output_list[["lambda"]] <- as.matrix(lambda)
    output_list[["omega"]] <- as.matrix(omega)
    output_list[["gamma_k"]] <- gamma_k
    output_list[["dc"]] <- dc
    return(output_list)
}


# When K is very weakly loaded on, remove
remove_weakly_active <- function(
    zeros,
    lambda,
    omega,
    gamma_k,
    dc
    ) {
    m_k <- colSums(zeros)
    keep <- which(dc >= 10 & m_k > (0.01 * dim(zeros)[1]))

    if (length(keep) == 0) {
        t <- 0
        c <- 10
        while (t == 0) {
            c <- c - 1
            keep <- which(dc >= c & m_k > (0.01 * dim(zeros)[1]))
            if (length(keep) != 0) {
                t <- t + 1
            }
        }
    }

    # store outputs
    new_zeros <- zeros[, keep]
    new_lambda <- lambda[, keep]
    new_omega <- omega[, keep]
    new_gamma_k <- gamma_k[keep]
    new_dc <- dc[keep]

    output_list <- list()
    output_list[["zeros"]] <- as.matrix(new_zeros)
    output_list[["lambda"]] <- as.matrix(new_lambda)
    output_list[["omega"]] <- as.matrix(new_omega)
    output_list[["gamma_k"]] <- new_gamma_k
    output_list[["dc"]] <- new_dc
    return(output_list)
}
