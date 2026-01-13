
##bin model logloss


log_loss <- function(y, p, eps = 1e-15, weights = NULL) {
  # y: vector of 0/1 outcomes
  # p: vector of predicted probabilities in [0,1]
  # eps: clipping for numerical stability
  p <- pmin(pmax(p, eps), 1 - eps)
  if (is.null(weights)) {
    -mean(y * log(p) + (1 - y) * log(1 - p))
  } else {
    # Weighted mean log loss
    weights <- weights / sum(weights)
    -sum(weights * (y * log(p) + (1 - y) * log(1 - p)))
  }
}




columns <- c("Model", "Formula", "Time-varying", "Family", "AIC", "cAIC", "EDF_ep","EDF_om","rho",
             "Matern range", "Spatial SD", "Hessian_positive", "Sum_loglik", 'logloss' )
mod.select <- as.data.frame( matrix(data=NA, nrow =0, ncol=length(columns), byrow=TRUE))
colnames(mod.select) <- columns

##model selection table function

mod.select.fn <- function (){
  
  c<- as.data.frame( matrix(data=NA, nrow =1, ncol=length(columns), byrow=TRUE))
  colnames(c) <- columns
  c$Model <- mod.label
  c$Formula <-m$formula [1]
  c$"Time-varying" <- ifelse(is.null (m$time_varying), NA, paste(m$time_varying[1], m$time_varying[2])  )
  c$"Family" <- paste0(m$family[1]$family, "(link = ", m$family$link[1], ")")
  c$"cAIC" = cAIC(m,what='cAIC')
 #ee = cAIC(m,what='EDF')
  #c$"EDF_ep" = ee[1]
  #c$"EDF_om" = ee[2]
  
  ##spatial model 
  c$AIC <- AIC (m)
  c$rho <- m$sd_report[[1]]["rho"]
  c$`Matern range` <- m$sd_report[[1]]["range"]
  c$`Spatial SD` <- m$sd_report[[1]]["sigma_O"]
  c$"Hessian_positive" <- m$pos_def_hessian
  
  ##model validation 
  c$"Sum_loglik" <- m_cv$sum_loglik
 m_cvTT = sdmTMBcv_tntpreds(m_cv)

  fitTT = dplyr::bind_rows(m_cvTT)
  c$logloss<-  with(fitTT[fitTT$tt=='test',],log_loss(as.numeric(pa),as.numeric(pred)))
  return(c)
}
