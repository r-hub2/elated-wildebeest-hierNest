#' Cross-validated hierarchical nested regularization for subgroup models
#'
#' @description
#' Fits regularization paths for hierarchical subgroup-specific penalized learning problems,
#' leveraging nested group structure (such as Major Diagnostic Categories [MDC] and Diagnosis-Related Groups [DRG])
#' with options for lasso or overlapping group lasso penalties. Performs cross-validation
#' to select tuning parameters for penalization and subgroup structure.
#'
#' This function enables information sharing across related subgroups by reparameterizing 
#' covariate effects into overall, group-specific, and subgroup-specific components and
#' supports structured shrinkage through hierarchical regularization. Users may select
#' between the lasso (hierNest-Lasso) and overlapping group lasso (hierNest-OGLasso)
#' frameworks, as described in Jiang et al. (2024, submitted, see details below).
#'
#' @param x Matrix of predictors, of dimension \eqn{n \times p}{n * p}; each row is an observation.
#'   Can be a dense or sparse matrix.
#' @param y Response variable. For `family="gaussian"`, should be numeric. For `family="binomial"`, should be a
#'   factor with two levels or a numeric vector with two unique values.
#' @param group Optional vector or factor indicating group assignments for variables. Used for custom grouping.
#' @param family Character string specifying the model family. Options are `"gaussian"` (default) for least-squares regression,
#'   or `"binomial"` for logistic regression.
#' @param nlambda Number of lambda values to use for regularization path. Default is 100.
#' @param lambda.factor Factor determining the minimal value of lambda in the sequence,
#'   where `min(lambda) = lambda.factor * max(lambda)`. See Details.
#' @param pred.loss Character string indicating loss to minimize during cross-validation. Options include `"default"`, `"mse"`,
#'   `"deviance"`, `"mae"`, `"misclass"`, and `"ROC"`.
#' @param lambda Optional user-supplied sequence of lambda values (overrides `nlambda`/`lambda.factor`).
#' @param pf_group Optional penalty factors on the groups, as a numeric vector. Default adjusts for group size.
#' @param pf_sparse Optional penalty factors on the l1-norm (for sparsity), as a numeric vector.
#' @param intercept Logical; whether to include an intercept in the model. Default is TRUE.
#' @param asparse1 Relative weight(s) for the first (e.g., group) layer of the overlapping group lasso penalty. Default is c(0.5, 20).
#' @param asparse2 Relative weight(s) for the second (e.g., subgroup) layer. Default is c(0.01, 0.2).
#' @param asparse1_num Number of values in asparse1 grid (for grid search). Default is 4.
#' @param asparse2_num Number of values in asparse2 grid (for grid search). Default is 4.
#' @param standardize Logical; whether to standardize predictors prior to model fitting. Default is TRUE.
#' @param lower_bnd Lower bound(s) for coefficient values. Default is \code{-Inf}.
#' @param upper_bnd Upper bound(s) for coefficient values. Default is \code{Inf}.
#' @param eps Convergence tolerance for optimization. Default is 1e-8.
#' @param maxit Maximum number of optimization iterations. Default is 3e6.
#' @param hier_info Required for `method = "overlapping"`; a matrix describing the hierarchical structure of the subgroups
#'   (see Details).
#' @param method Character; either `"overlapping"` for overlapping group lasso, `"sparsegl"` for sparse group lasso,
#'   or `"general"` for other hierarchical regularization. Default is `"overlapping"`.
#' @param partition Character string; determines subgroup partitioning. Default is `"subgroup"`.
#' @param cvmethod Cross-validation method. Options include `"general"` (default), `"grid_search"`, or `"user_supply"` (for user-supplied grid).
#'
#' @details
#' The hierarchical nested framework decomposes covariate effects into overall, group, and subgroup-specific components,
#' with regularization encouraging fusion or sparsity across these hierarchical levels. The function can fit both
#' the lasso penalty (allowing arbitrary zero/non-zero patterns) and the overlapping group lasso penalty (enforcing
#' hierarchical selection structure), as described in Jiang et al. (2024, submitted).
#'
#' The argument `hier_info` must be supplied for `"overlapping"` method, and encodes the hierarchical
#' relationship between groups and subgroups (e.g., MDCs and DRGs).
#'
#' Cross-validation is used to select tuning parameters, optionally over a grid for hierarchical penalty weights
#' (asparse1, asparse2), and the regularization parameter lambda.
#'
#' @return An object containing the fitted hierarchical model and cross-validation results, including:
#'   \item{fit}{Fitted model object.}
#'   \item{lambda}{Sequence of lambda values considered.}
#'   \item{cv_error}{Cross-validation error/loss for each combination of tuning parameters.}
#'   \item{best_params}{Best tuning parameters selected.}
#'   \item{...}{Additional diagnostic and output fields.}
#'
#' @references
#' Jiang, Z., Huo, L., Hou, J., & Huling, J. D.
#' "Heterogeneous readmission prediction with hierarchical effect decomposition and regularization".
#' 
#'
#' @seealso
#' [glmnet::glmnet()], [hierNest::cv.sparsegl()]
#'
#' @examples
#' # Simulated data:
#' set.seed(123)
#' n <- 60; p <- 4
#' x <- matrix(rnorm(n * p), n, p)
#' y <- rbinom(n, 1, plogis(x[,1] - 0.5*x[,2]))
#' # Create toy hier_info: 2 groups, 3 subgroups in each
#' hier_info <- cbind(rep(1:2, each=n/2), rep(1:4, each = n/4))
#' library(rTensor)
#' cv.fit=hierNest::cv.hierNest(x,y,method="overlapping", hier_info=hier_info,  cvmethod = "general",intercept = FALSE,family="binomial")
#' @export
#' 
#' 
cv.hierNest = function(x, y, 
                    group = NULL,
                    family = c("gaussian", "binomial"),
                    nlambda=100,
                    lambda.factor= NULL,
                    pred.loss = c("default", "mse", "deviance", "mae", "misclass","ROC"),
                    lambda=NULL,
                    pf_group=NULL,
                    pf_sparse=NULL,
                    intercept=FALSE,
                    asparse1=c(0.5,20),
                    asparse2=c(0.01,0.2),
                    asparse1_num=4,
                    asparse2_num=4,
                    standardize=TRUE,
                    lower_bnd = -Inf,
                    upper_bnd = Inf,
                    eps = 1e-08, 
                    maxit=3e+06,
                    hier_info=NULL, 
                    method="overlapping",
                    partition="subgroup",
                    cvmethod="general"){
 
  if(method=="overlapping"){
    
    if (is.null(hier_info)) {
      cli::cli_abort("must provide the hierarchy information through `hier_info`.")
    }else{
      iden_matrix=matrix(nrow = NROW(x),ncol = (1+max(hier_info[,1])+max(hier_info[,2])))
      
      iden_matrix[,1]=1
      
      curr_ix=2
      
      
      drgix_single=1:max(hier_info[,1])
      drgiy_single=1:max(hier_info[,1])

      for(i in 1:max(hier_info[,1])){
        iden_matrix[,curr_ix]=ifelse(hier_info[,1]==i,1,0)
        drgix_single[i]=ifelse(i==1,2,curr_ix)
        curr_ix=curr_ix+1
        
        hier_curr=hier_info[hier_info[,1]==i,2]
        
        for(j in (min(hier_curr)):(max(hier_curr))){
          iden_matrix[,curr_ix]=ifelse(hier_info[,2]==j,1,0)
          curr_ix=curr_ix+1
        }
        drgiy_single[i]=curr_ix-1
      }
      
      p=NCOL(x)
      
      design=(t(khatri_rao(t(iden_matrix),t(cbind(matrix(rep(1,NROW(x)),ncol = 1),x)))))

      trans_design.inx=1:NCOL(design)
      
      for(i in 1:(p+1)){
        for(j in 1:NCOL(iden_matrix)){
          trans_design.inx[(i-1)*NCOL(iden_matrix)+j]=(j-1)*(p+1)+i
        }
      }
      x.design=design[,trans_design.inx]
      cn=max(hier_info[,1])*(p+1)
      drgix=1:cn
      drgiy=1:cn
      ncol_single=NCOL(iden_matrix)
      
      
      drgix[1:length(drgix_single)]=drgix_single
      drgiy[1:length(drgix_single)]=drgiy_single
      
      for(i in 1:(p)){
        drgix[(i*length(drgix_single)+1):((i+1)*length(drgix_single))]=drgix_single+i*ncol_single
        drgiy[(i*length(drgiy_single)+1):((i+1)*length(drgiy_single))]=drgiy_single+i*ncol_single
      }
      
      # drgix=drgix-1
      # drgiy=drgiy-1
      # 
      
      cn_s=(1:(p+1))*max(hier_info[,1])-max(hier_info[,1])+1
      cn_e=(1:(p+1))*max(hier_info[,1])
      
      
      group_use=rep(1:(p+1), each=ncol_single)
      #group_use=group_use[-1]
      
      
      bs <- as.integer(as.numeric(table(group_use)))
      
      
      if(is.null(pf_group)){
        pf_group=sqrt(bs)
      }
      
      
      np <- dim(x.design)
      nobs <- as.integer(np[1])
      nvars <- as.integer(np[2])
      
      if(is.null(pf_sparse)){
        pf_sparse=rep(1, nvars)
      }
      if(is.null(lambda.factor)){
        lambda.factor= ifelse(nobs < nvars, 0.01, 1e-04)
      }
      
      
      x.design.spars=as(x.design,"sparseMatrix")
      subgroupnumber=1:max(hier_info[,2])
      
      for(i in 1:max(hier_info[,2])){
        subgroupnumber[i]=sum(hier_info[,2]==i)
      }
      
      if(cvmethod=="general"){
        
        res= cv.sparsegl(x.design.spars,y,
                                  group =group_use,family=family,
                                  cn=cn,
                                  drgix=drgix,
                                  drgiy=drgiy,
                                  cn_s=cn_s,
                                  cn_e=cn_e,
                                  intercept = intercept,
                                  pred.loss =pred.loss,
                                  subgroupnumber=subgroupnumber,partition = partition,
                                  asparse1=asparse1,asparse2=asparse2,
                                  nlambda=nlambda,lambda.factor=lambda.factor,lambda=lambda,
                                  pf_group=pf_group,pf_sparse=pf_sparse,
                                  standardize=standardize,
                                  lower_bnd=lower_bnd,upper_bnd=upper_bnd,
                                  eps=eps,maxit=maxit)
        
        
      }
      
      
      if(cvmethod=="grid_search"){
        
        log_asp1=log(asparse1)
        log_asp2=log(asparse2)
        
        
        if(asparse1_num==1){
          asparse1_base=exp(mean(log_asp1))
        }else{
          if(asparse1_num==2){
            asparse1_base=c(exp(log_asp1[1]),exp(log_asp1[2]))
          }else{
            asparse1_base=1:asparse1_num
            asparse1_base[1]=exp(log_asp1[1])
            asparse1_base[asparse1_num]=exp(log_asp1[2])
            for(i in 1:(asparse1_num-2)){
              asparse1_base[i+1]=exp(log_asp1[1]+i*(log_asp1[2]-log_asp1[1])/(asparse1_num-1))
            }
            
          }
        }
        
        if(asparse2_num==1){
          asparse2_base=exp(mean(log_asp2))
        }else{
          if(asparse2_num==2){
            asparse2_base=c(exp(log_asp2[1]),exp(log_asp2[2]))
          }else{
            asparse2_base=1:asparse2_num
            asparse2_base[1]=exp(log_asp2[1])
            asparse2_base[asparse2_num]=exp(log_asp2[2])
            for(i in 1:(asparse2_num-2)){
              asparse2_base[i+1]=exp(log_asp2[1]+i*(log_asp2[2]-log_asp2[1])/(asparse2_num-1))
            }
          }
        }
        
        res_min=c()
        res_return=list()
        lambda_seq=c()
        alpha1_seq=c()
        alpha2_seq=c()
        minvalue_seq=c()
        n1_min=1
        n2_min=1
        minvalue=NA
        
        for(n1 in 1:asparse1_num){
          for(n2 in 1:asparse2_num){
            #print(c(asparse1_base[n1],asparse2_base[n2]))
            
            if((n1==1)&(n2==1)){
              
              first.res= cv.sparsegl(x.design.spars,y,
                                        group =group_use,family=family,
                                        cn=cn,
                                        drgix=drgix,
                                        drgiy=drgiy,
                                        cn_s=cn_s,
                                        cn_e=cn_e,
                                        intercept = intercept,
                                        pred.loss =pred.loss,
                                        subgroupnumber=subgroupnumber,partition = partition,
                                        asparse1=rep(asparse1_base[1],nlambda),
                                        asparse2=rep(asparse2_base[1],nlambda),
                                        nlambda=nlambda,
                                        lambda.factor=lambda.factor,
                                      
                                        pf_group=pf_group,pf_sparse=pf_sparse,
                                        standardize=standardize,
                                        lower_bnd=lower_bnd,upper_bnd=upper_bnd,
                                        eps=eps,maxit=maxit)
              
              alpha1_seq=asparse1_base[1]
              alpha2_seq=asparse2_base[1]
              minvalue_seq=first.res$minvalue
              minvalue=first.res$minvalue
              minobj=first.res
            }else{
              temp.res= cv.sparsegl(x.design.spars,y,
                                             group =group_use,family=family,
                                             cn=cn,
                                             drgix=drgix,
                                             drgiy=drgiy,
                                             cn_s=cn_s,
                                             cn_e=cn_e,
                                             intercept = intercept,
                                             pred.loss =pred.loss,
                                             subgroupnumber=subgroupnumber,partition = partition,
                                             asparse1=rep(asparse1_base[n1],nlambda),
                                             asparse2=rep(asparse2_base[n2],nlambda),
                                             nlambda=nlambda,
                                             lambda.factor=lambda.factor,
                                             pf_group=pf_group,pf_sparse=pf_sparse,
                                             standardize=standardize,
                                             lower_bnd=lower_bnd,upper_bnd=upper_bnd,
                                             eps=eps,maxit=maxit,foldid = first.res$foldid)
              
              alpha1_seq=c(alpha1_seq,asparse1_base[n1])
              alpha2_seq=c(alpha2_seq,asparse2_base[n2])
              minvalue_seq=c(minvalue_seq,temp.res$minvalue)
              
              if(minvalue>temp.res$minvalue){
                minvalue=temp.res$minvalue
                n1_min=n1
                n2_min=n2
                minobj=temp.res
              }
            }
          }
        }
        
        res=minobj
        res[["min_alpha1"]] = alpha1_seq[order(minvalue_seq)][1]
        res[["min_alpha2"]] = alpha2_seq[order(minvalue_seq)][1]
        res[["hier.info"]] = hier_info
        res[["X.names"]] = colnames(x)
        
        # for(n1 in 1:asparse1_num){
        #   for(n2 in 1:asparse2_num){
        #     res_trans=hierNest::cv.sparsegl(x.design.spars,y,
        #                                     group =group_use,family=family,
        #                                     cn=cn,
        #                                     drgix=drgix,
        #                                     drgiy=drgiy,
        #                                     cn_s=cn_s,
        #                                     cn_e=cn_e,
        #                                     intercept = intercept,
        #                                     pred.loss =pred.loss,
        #                                     subgroupnumber=subgroupnumber,partition = partition,
        #                                     asparse1=asparse1_base[n1],asparse2=asparse2_base[n2],
        #                                     nlambda=nlambda_single,lambda.factor=lambda.factor,lambda=lambda,
        #                                     pf_group=pf_group,pf_sparse=pf_sparse,
        #                                     standardize=standardize,
        #                                     lower_bnd=lower_bnd,upper_bnd=upper_bnd,
        #                                     eps=eps,maxit=maxit)
        #     res_min=c(res_min,res_trans$cvmin)
        #     res_return=c(res_return,list(res_trans))
        #   }
        # }
        # 
        # 
        # min_inx=order(res_min)[1]
        # 
        # res=res_return[[min_inx]]
      }
      
      if(cvmethod=="user_supply"){
        
        asparse_num=length(asparse1)
        if(length(asparse1)!=length(asparse2)){
          cli::cli_abort("For user supply sparse parameters, asparse1 must have the same length as asparse2.")
        }

        res_min=c()
        res_return=list()
        lambda_seq=c()
        alpha1_seq=c()
        alpha2_seq=c()
        minvalue_seq=c()
        n1_min=1
        n2_min=1
        minvalue=NA
        
        for(n1 in 1:asparse_num){
          
          if(n1==1){
            temp.res= cv.sparsegl(x.design.spars,y,
                                           group =group_use,family=family,
                                           cn=cn,
                                           drgix=drgix,
                                           drgiy=drgiy,
                                           cn_s=cn_s,
                                           cn_e=cn_e,
                                           intercept = intercept,
                                           pred.loss =pred.loss,
                                           subgroupnumber=subgroupnumber,partition = partition,
                                           asparse1=rep(asparse1[n1],nlambda),
                                           asparse2=rep(asparse2[n1],nlambda),
                                           nlambda=nlambda,
                                           lambda.factor=lambda.factor,
                                           
                                           pf_group=pf_group,pf_sparse=pf_sparse,
                                           standardize=standardize,
                                           lower_bnd=lower_bnd,upper_bnd=upper_bnd,
                                           eps=eps,maxit=maxit)
            
            
            alpha1_seq=asparse1[1]
            alpha2_seq=asparse2[1]
            minvalue_seq=temp.res$minvalue
            minvalue=temp.res$minvalue
            minobj=temp.res
          }else{
            temp.res= cv.sparsegl(x.design.spars,y,
                                           group =group_use,family=family,
                                           cn=cn,
                                           drgix=drgix,
                                           drgiy=drgiy,
                                           cn_s=cn_s,
                                           cn_e=cn_e,
                                           intercept = intercept,
                                           pred.loss =pred.loss,
                                           subgroupnumber=subgroupnumber,partition = partition,
                                           asparse1=rep(asparse1[n1],nlambda),
                                           asparse2=rep(asparse2[n1],nlambda),
                                           nlambda=nlambda,
                                           lambda.factor=lambda.factor,
                                           pf_group=pf_group,pf_sparse=pf_sparse,
                                           standardize=standardize,
                                           lower_bnd=lower_bnd,upper_bnd=upper_bnd,
                                           eps=eps,maxit=maxit,foldid = temp.res$foldid)
            
            alpha1_seq=c(alpha1_seq,asparse1[n1])
            alpha2_seq=c(alpha2_seq,asparse2[n1])
            minvalue_seq=c(minvalue_seq,temp.res$minvalue)
            
            if(minvalue>temp.res$minvalue){
              minvalue=temp.res$minvalue
              n1_min=n1
              n2_min=n1
              minobj=temp.res
            }
          }  
        }
        res = minobj
        res[["min_alpha1"]] = alpha1_seq[order(minvalue_seq)][1]
        res[["min_alpha2"]] = alpha2_seq[order(minvalue_seq)][1]
        res[["hier.info"]] = hier_info
        res[["X.names"]] = colnames(x)
      }
      class(res) <- "cv.hierNest"
      
      return(res)
    }
  }
  
  
  
  
}

