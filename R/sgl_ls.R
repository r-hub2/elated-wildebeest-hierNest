sgl_ls_my <- function(
  bn, bs, ix, iy, nobs, nvars, x, y, pf, pfl1, dfmax, pmax, nlam,
  flmin, ulam, eps, maxit, vnames, group, intr, asparse1,asparse2, standardize,
  lower_bnd, upper_bnd,
  cn,drgix,drgiy,cn_s,cn_e) {
  # call Fortran core
  is.sparse <- FALSE
  #print("sgl_ls_my")
  if (!is.numeric(y)) rlang::abort("For family = 'gaussian', y must be numeric.")
  if (inherits(x,"sparseMatrix")) {
    is.sparse <- TRUE
    x <- as_dgCMatrix(x)
  }
  ym <- mean(y)
  if (intr) {
    y <- y - ym
    nulldev <- mean(y^2)
  } else {
    nulldev <- mean((y - ym)^2)
  }
  if (standardize) {
    sx <- sqrt(Matrix::colSums(x^2))
    sx[sx < sqrt(.Machine$double.eps)] <- 1 # Don't divide by zero!]
    xs <- 1 / sx
    x <- x %*% Matrix::Diagonal(x = xs)
  }
  if (is.sparse) {
    xidx <- as.integer(x@i + 1)
    xcptr <- as.integer(x@p + 1)
    xval <- as.double(x@x)
    nnz <- as.integer(utils::tail(x@p, 1))
  }

  gamma <- calc_gamma(x, ix, iy, bn)

  if (!is.sparse) {
    fit <- dotCall64::.C64(
      "sparse_four",
      SIGNATURE = c("integer", "integer", "integer", "integer", "double","integer", "integer", "double", "double", "double", "double",
                    "integer", "integer", "integer", "double", "double","double", "integer", "integer", 
                    "integer", "double", "double", "integer","integer", "double", "integer", "integer", "double",
                    "double", "double","double", "double","integer","integer","integer","integer","integer"),
      # Read only
      bn = bn, bs = bs, ix = ix, iy = iy, gam = gamma, nobs = nobs,
      nvars = nvars, x = as.double(x), y = as.double(y), pf = pf,
      pfl1 = pfl1,
      # Read / write
      dfmax = dfmax, pmax = pmax, nlam = nlam, flmin = flmin, ulam = ulam,
      eps = eps, maxit = maxit, intr = as.integer(intr),
      # Write only
      nalam = integer_dc(1), b0 = numeric_dc(nlam),
      beta = numeric_dc(nvars * nlam),
      activeGroup = integer_dc(pmax), nbeta = integer_dc(nlam),
      alam = numeric_dc(nlam), npass = integer_dc(1),
      jerr = integer_dc(1), mse = numeric_dc(nlam),
      # read only
      alsparse1 = asparse1,alsparse2 = asparse2, lb = lower_bnd, ub = upper_bnd,
      cn=cn,drgix=drgix,drgiy=drgiy,cn_s=cn_s,cn_e=cn_e,
      INTENT = c(rep("r", 11), rep("rw", 8), rep("w", 9), rep("r", 9)),
      NAOK = TRUE,
      PACKAGE = "hierNest")
  } else { # sparse design matrix
    fit <- dotCall64::.C64(
      "spmat_four",
      SIGNATURE = c("integer", "integer", "integer", "integer", "double",
                    "integer", "integer", "double", "integer", "integer",
                    "integer", "double", "double", "double", "integer", "integer",
                    "integer", "double", "double", "double", "integer",
                    "integer", "integer", "double", "double", "integer",
                    "integer", "double", "integer", "integer", "double",
                    "double", "double", "double","double","integer","integer","integer","integer","integer"),
      # Read only
      bn = bn, bs = bs, ix = ix, iy = iy, 
      
      
      gam = gamma, nobs = nobs,
      nvars = nvars, x = as.double(xval), xidx = xidx, xcptr = xcptr,
      nnz = nnz, y = as.double(y), pf = pf, pfl1 = pfl1,
      # Read write
      dfmax = dfmax, pmax = pmax, nlam = nlam, flmin = flmin,
      ulam = ulam, eps = eps, maxit = maxit, intr = as.integer(intr),
      # Write only
      nalam = integer_dc(1), b0 = numeric_dc(nlam),
      beta = numeric_dc(nvars * nlam), activeGroup = integer_dc(pmax),
      nbeta = integer_dc(nlam), alam = numeric_dc(nlam),
      npass = integer_dc(1), jerr = integer_dc(1), mse = numeric_dc(nlam),
      # Read only
      alsparse1 = as.double(asparse1), alsparse2 = as.double(asparse2), lb = lower_bnd, ub = upper_bnd,
      cn=cn,drgix=drgix,drgiy=drgiy,cn_s=cn_s,cn_e=cn_e,
      INTENT = c(rep("r", 14), rep("rw", 8), rep("w", 9), rep("r", 9)),
      NAOK = TRUE,
      PACKAGE = "hierNest")
  }
  #print(as.double(asparse1))
  # output
  outlist <- getoutput(x, group, fit, maxit, pmax, nvars, vnames, eps)
  if (standardize) outlist$beta <- outlist$beta * xs

  if (intr) {
    outlist$b0 <- outlist$b0 + ym
  } else {
    outlist$beta[1,1]=mean(y)
    #outlist$b0 <- outlist$b0
    outlist$b0 <- rep(0, dim(outlist$beta)[2])
  }
  outlist$npasses <- fit$npass
  outlist$jerr <- fit$jerr
  outlist$group <- group
  outlist$mse <- fit$mse[seq(fit$nalam)]
  outlist$dev.ratio <- 1 - outlist$mse / nulldev
  outlist$nulldev <- nulldev
  outlist$gamma=gamma
  class(outlist) <- c("lsspgl")
  outlist
}
