#' Fit Hierarchical Nested Regularization Model (hierNest)
#' @importFrom stats runif
#' @importFrom plotly plot_ly
#' @importFrom rTensor khatri_rao
#' @description
#' Fits a hierarchical nested penalized regression model for subgroup-specific effects using
#' overlapping group lasso penalties. This function encodes the hierarchical structure (e.g., MDC and DRG)
#' via a reparameterized design matrix and enables information borrowing across related subgroups.
#'
#' @param x Matrix of predictors (\eqn{n \times p}), where each row is an observation.
#' @param y Response variable (numeric for "gaussian", binary or factor for "binomial").
#' @param group Optional grouping vector (not required for "overlapping" method).
#' @param family Model family; either "gaussian" for least squares or "binomial" for logistic regression.
#' @param nlambda Number of lambda values in the regularization path (default: 100).
#' @param lambda.factor Factor for minimal value of lambda in the sequence.
#' @param lambda Optional user-supplied lambda sequence.
#' @param pf_group Penalty factor(s) for each group; defaults to sqrt(group size).
#' @param pf_sparse Penalty factors for individual predictors (L1 penalty).
#' @param intercept Logical; should an intercept be included? Default is FALSE.
#' @param asparse1 Relative weight for group-level penalty (default: 1).
#' @param asparse2 Relative weight for subgroup-level penalty (default: 0.05).
#' @param standardize Logical; standardize predictors? Default is TRUE.
#' @param lower_bnd Lower bound for coefficients (default: -Inf).
#' @param upper_bnd Upper bound for coefficients (default: Inf).
#' @param eps Convergence tolerance (default: 1e-8).
#' @param maxit Maximum number of optimization iterations (default: 3e6).
#' @param hier_info Required. Matrix encoding the hierarchical structure (see Details).
#' @param random_asparse Logical; use random sparse penalty? Default: FALSE.
#' @param method Type of hierarchical regularization ("overlapping" [default], "sparsegl", or "general").
#'
#' @details
#' This function builds a hierarchical design matrix reflecting group/subgroup structure (e.g., Major Diagnostic Categories [MDCs]
#' and Diagnosis Related Groups [DRGs]), encoding overall, group-specific, and subgroup-specific effects. It fits a penalized model
#' using overlapping group lasso, as described in Jiang et al. (2024, submitted). The main computational engine is \code{hierNest::overlapping_gl}.
#'
#' @return
#' Returns a model fit object as produced by \code{hierNest::overlapping_gl}, including selected coefficients,
#' cross-validation results, and tuning parameters.
#'
#' @references
#' Jiang, Z., Huo, L., Hou, J., Vaughan-Sarrazin, M., Smith, M. A., & Huling, J. D. (2024).
#' Heterogeneous readmission prediction with hierarchical effect decomposition and regularization. 
#'
#' @seealso [hierNest::overlapping_gl()]
#'
#' @examples
#' # Example with toy data
#'library(hierNest)
#'data("example_data")
#'fit = hierNest(example_data$X,
#'               example_data$Y,
#'               hier_info=example_data$hier_info,
#'               family="binomial",
#'               asparse1 = 1,
#'               asparse2 = 1)
#' @export
#' 
hierNest = function(x, y, 
                    group = NULL,
                    family = c("gaussian", "binomial"),
                    nlambda=100,
                    lambda.factor= NULL,
                    lambda=NULL,
                    pf_group=NULL,
                    pf_sparse=NULL,
                    intercept=FALSE,
                    asparse1=1,
                    asparse2=0.05,
                    standardize=TRUE,
                    lower_bnd = -Inf,
                    upper_bnd = Inf,
                    eps = 1e-08, 
                    maxit=3e+06,
                    hier_info=NULL, 
                    random_asparse=FALSE,
                    method="overlapping"){
 
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
      
      design=(t(rTensor::khatri_rao(t(iden_matrix),t(cbind(matrix(rep(1,NROW(x)),ncol = 1),x)))))
      
      
      
      
      
      
      trans_design.inx=1:NCOL(design)
      
      
      for(i in 1:(p+1)){
        for(j in 1:NCOL(iden_matrix)){
          trans_design.inx[(i-1)*NCOL(iden_matrix)+j]=(j-1)*(p+1)+i
        }
      }
      
      x.design=design[,trans_design.inx]
      #x.design=design[,trans_design.inx][,-1]
     
      
      # 
      # 
      # x.design=matrix(nrow = NROW(x),ncol = NCOL(iden_matrix)*(p+1)-1)
      # 
      # x.design[,1:(NCOL(iden_matrix)-1)]= (t(khatri_rao(t(iden_matrix),t(matrix(rep(1,NROW(x)),ncol = 1)))) )[,-1]
      # 
      # for(i in 1:p){
      #   x.design[,(i*NCOL(iden_matrix)):((i+1)*NCOL(iden_matrix)-1)]=(t(khatri_rao(t(iden_matrix),t(x[,i]))))
      # }
      # 
      # 
      # 
      

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
      
      
      # res=overlapping_gl(x.design,y,
      #                    group =group_use,
      #                    family=family,
      #                    cn=cn,
      #                    drgix=drgix,
      #                    drgiy=drgiy,
      #                    cn_s=cn_s,
      #                    cn_e=cn_e,
      #                    intercept = intercept,
      #                    random_asparse = random_asparse,
      #                    nlambda=nlambda,lambda.factor=lambda.factor,lambda=lambda,
      #                    pf_group=pf_group,pf_sparse=pf_sparse,
      #                    asparse1=asparse1,asparse2=asparse2,
      #                    standardize=standardize,
      #                    lower_bnd=lower_bnd,upper_bnd=upper_bnd,
      #                    eps=eps,maxit=maxit)
      
      x.design.spars=as(x.design,"sparseMatrix")
      
      res= overlapping_gl(x.design.spars,y,
                         group =group_use,
                         family=family,
                         cn=cn,
                         drgix=drgix,
                         drgiy=drgiy,
                         cn_s=cn_s,
                         cn_e=cn_e,
                         intercept = intercept,
                         random_asparse = random_asparse,
                         nlambda=nlambda,lambda.factor=lambda.factor,lambda=lambda,
                         pf_group=pf_group,pf_sparse=pf_sparse,
                         asparse1=asparse1,asparse2=asparse2,
                         standardize=standardize,
                         lower_bnd=lower_bnd,upper_bnd=upper_bnd,
                         eps=eps,maxit=maxit)
      
      res$hier.info = hier_info
      
      
      class(res) = "hierNest"
      
      return(res)
      
      
    }
    
    
    
  }
  
  
  if(method=="sparsegl"){
    
  }
  
  if(method=="general"){
    
  }
  
  
}

