# Fits the model and finds the ordinary least squares estimator.
# Parameters: X (nxp design matrix), y (n-dim vector)
# Returns a list with the following information:
# - est_beta: the estimate for beta
# - cov_unscaled: (X^TX)^{-1}, multiply by an estimate of sigma^2 to get estimated
#     variance of est_beta
# - hat_matrix: H = X(X^TX)^{-1}X^T, the hat matrix
# - leverages: diagonal elements of H, the leverages
# - residuals: vector of fitted residuals \hat{e}_i, i=1,2,...,n
# - y_fitted: fitted values \hat{y}_i, i=1,2,...,n
fit_ols_model = function(X,y) {
  cov_unscaled = solve(t(X) %*% X) # can reuse this calculation of (X^TX)^{-1}
  est_beta = c(cov_unscaled %*% t(X) %*% y)
  y_fitted = X %*% est_beta
  residuals = y - y_fitted
  hat_matrix = X %*% cov_unscaled %*% t(X)
  leverages = diag(hat_matrix)
  return(list(est_beta=est_beta, 
              y_fitted=y_fitted,
              cov_unscaled=cov_unscaled,
              hat_matrix=hat_matrix,
              leverages=leverages,
              residuals=residuals))
}


