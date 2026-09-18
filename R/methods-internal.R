## Internal helpers shared by the coef() and predict() methods.

## Coefficient matrix (intercept row followed by the slopes) of a fitted
## path, evaluated at the requested lambda values `s` by linear
## interpolation between the fitted values (`lambda` is the sequence the
## path was fitted on; for nonconvex fits it is the penalty sequence).
coef_path <- function(object, s = NULL, lambda = object$lambda) {
  b0 <- t(as.matrix(object$b0))
  rownames(b0) <- "(Intercept)"
  nbeta <- rbind2(b0, object$beta)
  if (!is.null(s)) {
    if (!is.numeric(s) || any(s < 0))
      stop("s must be a vector of non-negative lambda values")
    vnames <- dimnames(nbeta)[[1]]
    dimnames(nbeta) <- list(NULL, NULL)
    lamlist <- lambda.interp(lambda, s)
    nbeta <- nbeta[, lamlist$left, drop = FALSE] %*%
      Diagonal(x = lamlist$frac) +
      nbeta[, lamlist$right, drop = FALSE] %*%
      Diagonal(x = 1 - lamlist$frac)
    dimnames(nbeta) <- list(vnames, paste(seq(along = s)))
  }
  nbeta
}

## coef() with the "coefficients"/"nonzero" switch
coef_path_type <- function(object, s, type, lambda = object$lambda) {
  nbeta <- coef_path(object, s, lambda)
  if (type == "coefficients") return(nbeta)
  nonzero(nbeta[-1, , drop = FALSE], bystep = TRUE)
}

## linear predictor for new data
predict_path <- function(object, newx, s, lambda = object$lambda) {
  if (missing(newx))
    stop("newx (a matrix of new predictor values) is required")
  newx <- as.matrix(newx)
  if (ncol(newx) != nrow(object$beta))
    stop("newx must have ", nrow(object$beta), " columns")
  nbeta <- coef_path(object, s, lambda)
  as.matrix(as.matrix(cbind2(1, newx)) %*% nbeta)
}

## resolve the `s` argument of the cv methods to lambda values
cv_lambda <- function(object, s) {
  if (is.numeric(s)) {
    s
  } else if (is.character(s)) {
    object[[match.arg(s, c("lambda.1se", "lambda.min"))]]
  } else {
    stop("Invalid form for s")
  }
}
