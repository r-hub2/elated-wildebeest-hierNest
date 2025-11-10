
cv.sparsegl <- function(
    x, y, group = NULL, family = c("gaussian", "binomial"),
    lambda = NULL,
    pred.loss = c("default", "mse", "deviance", "mae", "misclass","ROC"),
    nfolds = 10, foldid = NULL,
    weights = NULL, offset = NULL,
    cn,drgix,drgiy,cn_s,cn_e, partition=c("general","subgroup"), subgroupnumber,
    #asparse1=NULL,asparse2=NULL,
    ...) {

  fam <- validate_family(family)


  # not allowed for some families
  pred.loss <- match.arg(pred.loss)
  if (pred.loss == "misclass") {
    bugger <- FALSE
    if (fam$check == "char") if (fam$family != "binomial") bugger <- TRUE
    if (fam$check == "fam") if (fam$family$family != "binomial") bugger <- TRUE
    if (bugger) {
      cli::cli_abort(c(
        "When `pred.loss` is {.val {pred.loss}}, `family` must be either:",
        `!` = "{.val {'binomial'}}, or {.fn stats::binomial}."
      ))
    }
  }

  N <- nrow(x)
  ###Fit the model once to get dimensions etc of output
  y <- drop(y)
  # random_asparse=FALSE
  # if(is.null(asparse1)){
  #   random_asparse=TRUE
  # }
  sparsegl.object <- overlapping_gl(x, y, group, lambda = lambda, family = family,
                                 # asparse1=ifelse(is.null(asparse1),NULL,asparse1),
                                 # asparse2=ifelse(is.null(asparse2),NULL,asparse2),
                                 weights = weights, offset = offset,cn=cn,drgix=drgix,drgiy=drgiy,cn_s=cn_s,cn_e=cn_e, ...)
  lambda <- sparsegl.object$lambda
  # predict -> coef
  if (is.null(foldid)){
    if(partition=="general"){
      foldid <- sample(rep(seq(nfolds), length = N))
    }
    if(partition=="subgroup"){
      if(family=="binomial"){
        #print("binomial leave out")
        foldid=c()
        startinx=1
        endinx=subgroupnumber[1]
        for(i in 1:length(subgroupnumber)){
          foldid_temp=1:subgroupnumber[i]
          foldid_temp[y[startinx:endinx]==0]=sample(rep(seq(nfolds), length = sum(y[startinx:endinx]==0)))
          foldid_temp[y[startinx:endinx]==1]=sample(rep(seq(nfolds), length = sum(y[startinx:endinx]==1)))
          foldid=c(foldid,foldid_temp)
        }
      }else{
        foldid=c()
        for(i in 1:length(subgroupnumber)){
          foldid=c(foldid,sample(rep(seq(nfolds), length = subgroupnumber[i])))
        }
      }
    }
    
  }else{
    nfolds <- max(foldid)
  } 
  if (nfolds < 2) {
    cli::cli_abort(
      "`nfolds` must be at least {.val {2}}. `nfolds` = {.val {10}} is recommended."
    )
  }
  outlist <- as.list(seq(nfolds))
  ###Now fit the nfold models and store them
  for (i in seq(nfolds)) {
    test_fold <- foldid == i
    outlist[[i]] <- overlapping_gl(
      x = x[!test_fold, , drop = FALSE],
      y = y[!test_fold], group = group, lambda = lambda, family = family, 
      # sparsegl.object$asparse1,
      # sparsegl.object$asparse2,
      weights = weights[!test_fold], offset = offset[!test_fold],cn=cn,drgix=drgix,drgiy=drgiy,cn_s=cn_s,cn_e=cn_e, ...)
  }
  ###What to do depends on the pred.loss and the model fit
  cvstuff <- cverror(sparsegl.object, outlist, lambda, x, y, foldid,
                     pred.loss, weights)
  cvm <- cvstuff$cvm
  cvsd <- cvstuff$cvsd
  cvname <- cvstuff$name
  nz <- predict(sparsegl.object, type = "nonzero")
  nnzero <- sapply(nz, length)
  active_grps <- sapply(nz, function(x) length(unique(group[x])))
  out <- list(lambda = lambda, cvm = cvm, cvsd = cvsd, cvupper = cvm + cvsd,
              cvlo = cvm - cvsd, name = cvname,
              nnzero = nnzero, active_grps = active_grps,
              hierNest.fit = sparsegl.object,
              call = match.call())
  lamin <- getmin(lambda, cvm, cvsd)
  cv.inx=order(cvm)[1]
  obj <- c(out, as.list(lamin))
  obj[["cv.inx"]]=cv.inx
  obj[["foldid"]]=foldid
  obj[["minvalue"]]=min(cvm)
  class(obj) <- "cv.sparsegl"
  obj
}



cverror <- function(fullfit, outlist, lambda, x, y, foldid, pred.loss, ...) {
  UseMethod("cverror")
}

#' @export
cverror.lsspgl <- function(
    fullfit, outlist, lambda, x, y, foldid,
    pred.loss = c("default", "mse", "deviance", "mae"),
    ...) {

  typenames <- c(default = "Mean squared error", mse = "Mean squared error",
                 deviance = "Mean squared error", mae = "Mean absolute error")
  pred.loss <- match.arg(pred.loss)
  predmat <- matrix(NA, length(y), length(lambda))
  nfolds <- max(foldid)
  nlams <- double(nfolds)
  for (i in seq_len(nfolds)) {
    test_fold <- foldid == i
    fitobj <- outlist[[i]]
    preds <- predict(fitobj, x[test_fold, , drop = FALSE], type = "link")
    nlami <- length(outlist[[i]]$lambda)
    predmat[test_fold, seq_len(nlami)] <- preds
    nlams[i] <- nlami
  }
  cvraw <- switch(pred.loss, mae = abs(y - predmat), (y - predmat)^2)
  N <- length(y) - apply(is.na(predmat), 2, sum)
  cvm <- apply(cvraw, 2, mean, na.rm = TRUE)
  scaled <- scale(cvraw, cvm, FALSE)^2
  cvsd <- sqrt(apply(scaled, 2, mean, na.rm = TRUE) / (N - 1))
  list(cvm = cvm, cvsd = cvsd, name = typenames[pred.loss])
}


#' @export
cverror.logitspgl <- function(
    fullfit, outlist, lambda, x, y, foldid,
    pred.loss = c("default", "mse", "deviance", "mae", "misclass","ROC"),
    ...) {
  #print("cverror_binomial")
  typenames <- c(default = "Binomial deviance", mse = "Mean squared error",
                 deviance = "Binomial deviance", mae = "Mean absolute error",
                 misclass = "Missclassification error")
  pred.loss <- match.arg(pred.loss)
  prob_min <- 1e-05
  fmax <- log(1 / prob_min - 1)
  fmin <- -fmax
  y <- as.factor(y)
  y <- as.numeric(y) - 1 # 0 / 1 valued
  nfolds <- max(foldid)
  predmat <- matrix(NA, length(y), length(lambda))
  nlams <- double(nfolds)
  for (i in seq_len(nfolds)) {
    test_fold <- foldid == i
    fitobj <- outlist[[i]]
    preds <- predict(fitobj, x[test_fold, , drop = FALSE], type = "response")
    nlami <- length(outlist[[i]]$lambda)
    predmat[test_fold, seq_len(nlami)] <- preds
    nlams[i] <- nlami
  }
  predmat <- pmin(pmax(predmat, fmin), fmax)
  binom_deviance <- function(m) stats::binomial()$dev.resids(y, m, 1)
  if(pred.loss=="ROC"){
    cvm=1:dim(predmat)[2]
    for(i in 1:dim(predmat)[2]){
      roc_obj=tryCatch({
        (pROC::roc(y ~ predmat[,i]))
      },error=function(x){})
      
      cvm[i]=ifelse(is.null(roc_obj$auc),NA,1-roc_obj$auc)
    }
    
    N <- length(y) - apply(is.na(predmat), 2, sum)
    cvsd=NULL
    
  }else{
    cvraw <- switch(
      pred.loss,
      mse = (y - predmat)^2,
      mae = abs(y - predmat),
      misclass = y != ifelse(predmat > 0.5, 1, 0),
      apply(predmat, 2, binom_deviance)
    )
    N <- length(y) - apply(is.na(predmat), 2, sum)
    cvm <- apply(cvraw, 2, mean, na.rm = TRUE)
    cvsd <- sqrt(apply(scale(cvraw, cvm, FALSE)^2, 2, mean, NA.RM = TRUE) /
                   (N - 1))
  }
  
  list(cvm = cvm, cvsd = cvsd, name = typenames[pred.loss])
}

#' @export
cverror.irlsspgl <- function(
    fullfit, outlist, lambda, x, y, foldid,
    pred.loss = c("default", "mse", "deviance", "mae", "misclass"),
    weights = NULL, ...) {
  typenames <- c(default = "Deviance", mse = "Mean squared error",
                 deviance = "Deviance", mae = "Mean absolute error")
  pred.loss <- match.arg(pred.loss)

  nfolds <- max(foldid)
  predmat <- matrix(NA, length(y), length(lambda))
  nlams <- double(nfolds)
  for (i in seq_len(nfolds)) {
    test_fold <- foldid == i
    fitobj <- outlist[[i]]
    preds <- predict(fitobj, x[test_fold, , drop = FALSE], type = "response")
    nlami <- length(outlist[[i]]$lambda)
    predmat[test_fold, seq_len(nlami)] <- preds
    nlams[i] <- nlami
  }

  dev_fun <- function(m) fullfit$family$dev.resids(y, m, 1)
  cvraw <- switch(
    pred.loss,
    mse = (y - predmat)^2,
    mae = abs(y - predmat),
    misclass = y != ifelse(predmat > 0.5, 1, 0),
    apply(predmat, 2, dev_fun)
  )

  N <- length(y) - apply(is.na(predmat), 2, sum)
  if (is.null(weights)) weights <- rep(1, nrow(cvraw))
  cvm <- apply(cvraw, 2, stats::weighted.mean, na.rm = TRUE, w = weights)
  cvsd <- sqrt(
    apply(
      scale(cvraw, cvm, FALSE)^2,
      2, stats::weighted.mean, w = weights, NA.RM = TRUE
    ) / (N - 1)
  )
  list(cvm = cvm, cvsd = cvsd, name = typenames[pred.loss])
}
