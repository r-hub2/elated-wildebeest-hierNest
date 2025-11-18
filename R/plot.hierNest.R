#' Plot method for cv.hierNest objects
#'
#' @param x An object of class \code{cv.hierNest}
#' @param type Plot type, e.g. "fit" for observed vs fitted (default),
#'   others you may add later (e.g., "coef").
#' @param ... Other parameters
#' @return Invisibly returns \code{x}
#' @export
#' @method plot cv.hierNest
#' 
plot.cv.hierNest <- function(x, type = c("coefficients","Subgroup effects"),...) {
  
  type <- match.arg(type)
  
  cv.fit = x
  
  cvinx = which(cv.fit$hierNest.fit$lambda == cv.fit$lambda.min)
  hier_info = cv.fit$hier.info
  n.row = length(unique(hier_info[,1])) + length(unique(hier_info[,2])) + 1
  beta_mat = matrix(cv.fit$hierNest.fit$beta[,cvinx], nrow = n.row)
  if(is.null(cv.fit$X.names)){
    colnames(beta_mat) = c("Intercept", paste0("X",1:(NCOL(beta_mat)-1)))
  }else{
    colnames(beta_mat) = c("Intercept", cv.fit$X.names)
  }
  
  nonzero_variable_inx = which(apply(beta_mat,2, mean)!=0)
  beta_set_nonzero = beta_mat[,nonzero_variable_inx]
  

  
  if (type == "coefficients") {
    
    rowname_use = c("Overall mean effect")
    for(l1 in unique(hier_info[,1])){
      rowname_use = c(rowname_use, paste0("Group ", l1))
      for(l2 in unique(hier_info[which(hier_info[,1]==l1),2])){
        rowname_use = c(rowname_use, paste0("Subgroup ", l2))
      }
    }
    
    
    col_order = order(apply(beta_set_nonzero, 2, function(x){sum(x!=0)}),decreasing = TRUE)
    beta_set_nonzero = beta_set_nonzero[, col_order]

    mat = apply(beta_set_nonzero, 2, as.numeric)
    
    nr = nrow(mat)
    nc = ncol(mat)
    # force numeric values and explicit integer indices
    df <- data.frame(
      row   = rep(seq_len(nr), times = nc),
      col   = rep(seq_len(nc), each  = nr),
      value = as.numeric(mat)
    )
    
    eps <- 0  # dead-zone half-width around 0; set to 0 if you don't want blanking
    df$plot_val <- ifelse(abs(df$value) <= eps, NA, df$value)
    
    # symmetric clipping level from 99th percentile of |value|
    L <- stats::quantile(abs(df$plot_val), 0.95, na.rm = TRUE)
    if (!is.finite(L)) L <- max(abs(df$plot_val), na.rm = TRUE)
    
    col_labels = colnames(beta_set_nonzero)
    
    plotobj = plotly::plot_ly(
      data = df,
      x = ~col,          # same direction as ggplot
      y = ~row,
      z = ~plot_val,
      type = "heatmap",
      colorscale = list(
        c(0, "red"),     # low
        c(0.5, "white"),  # midpoint
        c(1, "blue")       # high
      ),
      zmin = -L,
      zmax = L,
      colorbar = list(title = "Value"),
      showscale = TRUE
    ) %>%
      plotly::layout(
        xaxis = list(
          title = "Predictors",
          tickmode = "array",
          tickvals = unique(df$col),     # the numeric positions
          ticktext = col_labels          # the names you want to display
        ),
        yaxis = list(
          title = "Subgroup effects",
          autorange = "reversed",      # keeps top-down order
          tickmode = "array",
          tickvals = unique(df$row),     # the numeric positions
          ticktext = rowname_use          # the names you want to display
        )
        
      )
    
  }
  
  
  
  if (type == "Subgroup effects") {
    
    rowname_use = paste0("Sub-Group ", unique(hier_info[,2]))
    
    meaneffect = beta_set_nonzero[1,]
    subgroup_effect_mat = matrix(nrow = length(unique(hier_info[,2])), ncol = NCOL(beta_set_nonzero))
    
    inxx1 = 1
    inxx2 = 1
    for(l1 in unique(hier_info[,1])){
      inxx1 = inxx1 + 1
      MDC_effect = beta_set_nonzero[inxx1,]
      for(l2 in unique(hier_info[which(hier_info[,1]==l1),2])){
        inxx1 = inxx1 + 1
        DRG_effect = beta_set_nonzero[inxx1,]
        subgroup_effect_mat[inxx2, ] = meaneffect + MDC_effect + DRG_effect
        inxx2 = inxx2 + 1
      }
    }
    
    
    mat = apply(subgroup_effect_mat, 2, as.numeric)
    nr = nrow(mat)
    nc = ncol(mat)
    # force numeric values and explicit integer indices
    df <- data.frame(
      row   = rep(seq_len(nr), times = nc),
      col   = rep(seq_len(nc), each  = nr),
      value = as.numeric(mat)
    )
    
    eps <- 0  # dead-zone half-width around 0; set to 0 if you don't want blanking
    df$plot_val <- ifelse(abs(df$value) <= eps, NA, df$value)
    
    # symmetric clipping level from 99th percentile of |value|
    L <- stats::quantile(abs(df$plot_val), 0.95, na.rm = TRUE)
    if (!is.finite(L)) L <- max(abs(df$plot_val), na.rm = TRUE)
    
    col_labels = colnames(beta_set_nonzero)
    
    plotobj = plotly::plot_ly(
      data = df,
      x = ~col,          # same direction as ggplot
      y = ~row,
      z = ~plot_val,
      type = "heatmap",
      colorscale = list(
        c(0, "red"),     # low
        c(0.5, "white"),  # midpoint
        c(1, "blue")       # high
      ),
      zmin = -L,
      zmax = L,
      colorbar = list(title = "Value"),
      showscale = TRUE
    ) %>%
      plotly::layout(
        xaxis = list(
          title = "Predictors",
          tickmode = "array",
          tickvals = unique(df$col),     # the numeric positions
          ticktext = col_labels          # the names you want to display
        ),
        yaxis = list(
          title = "Subgroup effects",
          autorange = "reversed",      # keeps top-down order
          tickmode = "array",
          tickvals = unique(df$row),     # the numeric positions
          ticktext = rowname_use          # the names you want to display
        )
        
      )
    
  }
  
  
  #plot(plotobj)

  return(plotobj)
}
