#' Fit an Overlapping Group Lasso Model
#'
#' @description
#' This function fits an overlapping group lasso model with hierarchical regularization, allowing both sparsity within groups and across groups.
#'
#' @param x A numeric matrix of predictor variables (no missing values allowed).
#' @param y A numeric vector of response variable values.
#' @param group An integer vector defining the group membership for each predictor. Default is NULL (each predictor forms its own group).
#' @param family A character string specifying the model family. Options are "gaussian" (default) and "binomial".
#' @param nlambda Number of lambda values. Default is 100.
#' @param lambda.factor Factor determining minimum lambda as a fraction of maximum lambda. Default depends on dimensionality.
#' @param lambda Numeric vector of lambda values. If provided, overrides `nlambda` and `lambda.factor`.
#' @param pf_group Penalty factor for groups. Default is square root of group size.
#' @param pf_sparse Penalty factor for individual predictors. Default is 1 for each predictor.
#' @param intercept Logical; whether to include an intercept. Default is TRUE.
#' @param asparse1 Sparsity penalty factor controlling group-level sparsity. Default is 1.
#' @param asparse2 Sparsity penalty factor controlling within-group sparsity. Default is 0.05.
#' @param standardize Logical; if TRUE, standardizes predictors before fitting. Default is TRUE.
#' @param lower_bnd Numeric vector specifying lower bounds for coefficients. Default is -Inf.
#' @param upper_bnd Numeric vector specifying upper bounds for coefficients. Default is Inf.
#' @param weights Optional numeric vector of observation weights. Currently limited functionality.
#' @param offset Optional numeric vector specifying a known component to be included in the linear predictor.
#' @param warm Optional initial values for optimization.
#' @param trace_it Integer indicating the verbosity level. Default is 0 (no output).
#' @param dfmax Maximum number of groups allowed in the model. Default derived from groups.
#' @param pmax Maximum number of predictors allowed in the model. Default derived from groups.
#' @param eps Numeric convergence threshold for optimization. Default is 1e-08.
#' @param maxit Maximum number of iterations for optimization. Default is 3e+06.
#' @param cn Additional internal numeric parameter for optimization.
#' @param drgix,drgiy Numeric vectors specifying indices for specific group and predictor structures.
#' @param cn_s,cn_e Numeric vectors specifying starting and ending indices for substructures.
#' @param random_asparse Logical; if TRUE, randomly selects sparsity parameters. Default is FALSE.
#'
#' @return An object of class `sparsegl` containing:
#' \item{call}{The matched function call.}
#' \item{lambda}{The lambda values used for fitting.}
#' \item{asparse1, asparse2}{Sparsity parameters used.}
#' \item{nobs}{Number of observations.}
#' \item{pf_group, pf_sparse}{Penalty factors used.}
#' \item{coefficients}{Estimated coefficients.}
#' Additional components relevant to model diagnostics and fitting.
#'
#'
overlapping_gl <- function(
  x, y, group = NULL, family = c("gaussian", "binomial"),
  nlambda = 100, lambda.factor = ifelse(nobs < nvars, 0.01, 1e-04),
  lambda = NULL, 
  pf_group = sqrt(bs), 
  pf_sparse = rep(1, nvars),
  intercept = TRUE, 
  asparse1=1, # 
  asparse2=0.05, 
  standardize = TRUE,
  lower_bnd = -Inf, upper_bnd = Inf,
  weights = NULL, offset = NULL, warm = NULL,
  trace_it = 0,
  dfmax = as.integer(max(group)) + 1L,
  pmax = min(dfmax * 1.2, as.integer(max(group))),
  eps = 1e-08, maxit = 3e+06,cn,drgix,drgiy,cn_s,cn_e,random_asparse=FALSE) {

  this.call <- match.call()
  if (!is.matrix(x) && !inherits(x, "sparseMatrix")) {
    cli::cli_abort("`x` must be a matrix.")
  }

  if (any(is.na(x))) cli::cli_abort("Missing values in `x` are not supported.")

  y <- drop(y)
  if (!is.null(dim(y))) cli::cli_abort("`y` must be a vector or 1-column matrix.")
  np <- dim(x)
  nobs <- as.integer(np[1])
  nvars <- as.integer(np[2])
  vnames <- colnames(x)

  if (is.null(vnames)) vnames <- paste("V", seq(nvars), sep = "")

  if (length(y) != nobs) {
    cli::cli_abort("`x` has {nobs} rows while `y` has {length(y)}.")
  }

  #    group setup
  if (is.null(group)) {
    group <- 1:nvars
  } else {
    if (length(group) != nvars) {
      cli::cli_abort(c(
        "The length of `group` is {length(group)}.",
        "It must match the number of columns in `x`: {nvars}"
      ))
    }
  }

  bn <- as.integer(max(group))  # number of groups
  bs <- as.integer(as.numeric(table(group)))  # number of elements in each group

  if (!identical(as.integer(sort(unique(group))), as.integer(1:bn))) {
    cli::cli_abort("Groups must be consecutively numbered 1, 2, 3, ...")
  }

  # if (asparse1 > 1) {
  #   cli::cli_abort(c(
  #     "`asparse` must be less than or equal to 1.",
  #     i = "You may want {.fn glmnet::glmnet} instead."
  #   ))
  # }

  # if ((asparse1 < 0)|(asparse2 < 0)) {
  #   asparse1 <- 0
  #   asparse2 <- 0
  #   cli::cli_warn("`asparse` must be in {.val [0, 1]}, running ordinary group lasso.")
  # }
  if (any(pf_sparse < 0)) cli::cli_abort("`pf_sparse` must be non-negative.")
  if (any(is.infinite(pf_sparse))) {
    cli::cli_abort(
      "`pf_sparse` may not be infinite. Simply remove the column from `x`."
    )
  }
  if (any(pf_group < 0)) cli::cli_abort("`pf_group` must be non-negative.")
  if (any(is.infinite(pf_group))) {
    cli::cli_abort(c(
      "`pf_group` must be finite.",
      i = "Simply remove the group from `x`."
    ))
  }
  if (all(pf_sparse == 0)) {
    # if (asparse1 > 0) {
    #   cli::cli_abort(
    #     "`pf_sparse` is identically 0 but `asparse` suggests some L1 penalty is desired."
    #   )
    # } else {
    #   cli::cli_warn("`pf_sparse` was set to 1 because `asparse` = {.val {0}}.")
    #   pf_sparse = rep(1, nvars)
    # }
  }

  ## Note: should add checks to see if any columns are completely unpenalized
  ## This is not currently expected.

  iy <- cumsum(bs) # last column of x in each group
  ix <- c(0, iy[-bn]) + 1 # first column of x in each group
  ix <- as.integer(ix)
  iy <- as.integer(iy)
  group <- as.integer(group)

  #parameter setup
  if (length(pf_group) != bn) {
    cli::cli_abort(
      "The length of `pf_group` must be the same as the number of groups: {.val {bn}}."
    )
  }
  if (length(pf_sparse) != nvars) {
    cli::cli_abort(
      "The length of `pf_sparse` must be equal to the number of predictors: {.val {nvars}}."
    )
  }

  pf_sparse <- pf_sparse / sum(pf_sparse) * nvars
  maxit <- as.integer(maxit)
  pf_group <- as.double(pf_group)
  pf_sparse <- as.double(pf_sparse)
  eps <- as.double(eps)
  dfmax <- as.integer(dfmax)
  pmax <- as.integer(pmax)

  #lambda setup
  nlam <- as.integer(nlambda)
  if (is.null(lambda)) {
    if (lambda.factor >= 1) {
      cli::cli_abort("`lambda.factor` must be less than {.val {1}}.")
    }
    flmin <- as.double(lambda.factor)
    ulam <- double(1)
  } else {
    #flmin = 1 if user define lambda
    flmin <- as.double(1)
    if (any(lambda < 0)) cli::cli_abort("`lambda` must be non-negative.")
    ulam <- as.double(rev(sort(lambda)))
    nlam <- as.integer(length(lambda))
  }
  intr <- as.integer(intercept)

  ### check on upper/lower bounds
  if (any(lower_bnd > 0)) cli::cli_abort("`lower_bnd` must be non-positive.")
  if (any(upper_bnd < 0)) cli::cli_abort("`upper_bnd` must be non-negative.")
  lower_bnd[lower_bnd == -Inf] <- -9.9e30
  upper_bnd[upper_bnd == Inf] <- 9.9e30
  if (length(lower_bnd) < bn) {
    if (length(lower_bnd) == 1) lower_bnd <- rep(lower_bnd, bn)
    else cli::cli_abort("`lower_bnd` must be length {.val {1}} or length {.val {bn}}.")
  } else {
    lower_bnd <- lower_bnd[seq_len(bn)]
  }
  if (length(upper_bnd) < bn) {
    if (length(upper_bnd) == 1) upper_bnd <- rep(upper_bnd, bn)
    else cli::cli_abort("`upper_bnd` must be length {.val {1}} or length {.val {bn}}.")
  } else {
    upper_bnd <- upper_bnd[seq_len(bn)]
  }
  storage.mode(upper_bnd) <- "double"
  storage.mode(lower_bnd) <- "double"

  # call R sub-function
  fam <- validate_family(family)
  if(random_asparse){
    cli::cli_warn(c(
      "Randomly select penalty factor 1"
    ))
    sample_log_alpha1=runif(nlam,min = -1,max = 3)
    asparse1=exp(sample_log_alpha1)
    
  }
  
  if(random_asparse){
    cli::cli_warn(c(
      "Randomly select penalty factor 2"
    ))
    sample_alpha2=runif(nlam,min = 0.01,max = 0.20)
    asparse2=sample_alpha2
    
  }
  if (fam$check == "char") {
    family <- match.arg(family)
    if (!is.null(weights)) {
      cli::cli_warn(c(
        "Currently, `weights` are only supported when `family` has class {.cls family}.",
        i = "Estimating unweighted sparse group lasso. See {.fn sparsegl::sparsegl}."
      ))
    }
    if (!is.null(offset)) {
      cli::cli_warn(c(
        "Currently, `offset` is only supported when `family` has class {.cls family}.",
        i = "Estimating sparse group lasso without any offset. See {.fn sparsegl::sparsegl}."
      ))
    }
    
    fit <- switch(
      family,
      gaussian = sgl_ls_my(
        bn, bs, ix, iy, nobs, nvars, x, y, pf_group, pf_sparse,
        dfmax, pmax, nlam, flmin, ulam, eps, maxit, vnames, group, intr,
        asparse1=as.double(asparse1),asparse2=as.double(asparse2), standardize, lower_bnd, upper_bnd,cn=cn,drgix=drgix,drgiy=drgiy,cn_s=cn_s,cn_e=cn_e),
      binomial = sgl_logit(
        bn, bs, ix, iy, nobs, nvars, x, y, pf_group, pf_sparse,
        dfmax, pmax, nlam, flmin, ulam, eps, maxit, vnames, group, intr,
        asparse1=as.double(asparse1),asparse2=as.double(asparse2), standardize, lower_bnd, upper_bnd,cn=cn,drgix=drgix,drgiy=drgiy,cn_s=cn_s,cn_e=cn_e)
    )
  }
  # if (fam$check == "fam") {
  #   fit <- sgl_irwls(
  #     bn, bs, ix, iy, nobs, nvars, x, y, pf_group, pf_sparse,
  #     dfmax, pmax, nlam, flmin, ulam, eps, maxit, vnames, group, intr,
  #     as.double(asparse), standardize, lower_bnd, upper_bnd, weights,
  #     offset, fam$family, trace_it, warm
  #   )
  # }

  # output
  if (is.null(lambda)) fit$lambda <- lamfix(fit$lambda)
  fit$call <- this.call
  fit$asparse1 <- asparse1
  fit$asparse2 <- asparse2
  fit$nobs <- nobs
  fit$pf_group <- pf_group
  fit$pf_sparse <- pf_sparse
  class(fit) <- c(class(fit), "sparsegl")
  fit
}

