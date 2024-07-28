lasso1 = function(x, y, penalty = FALSE, alpha = 1, nlambda = 50,
                   lambda.min.ratio = ifelse(nrow(x) < ncol(x), 0.01, 1e-04), 
                   lambda = NULL, intercept = FALSE, standardize = TRUE, relax = FALSE, gamma = NULL) {
  
  x = as.matrix(x)
  y = as.numeric(y)
  n = nrow(x)
  p = ncol(x)
  dfmax = p
  
  if (relax == TRUE && n < (p - intercept)){
    dfmax = n - intercept
  } 
    
  if (!is.null(lambda)){
    nlambda = length(lambda) 
  } 
    
  
  w = rep(1, dfmax)
  
  # adaptive lasso
  if (penalty == TRUE){
    y1 = y - mean(y)
    x1 = scale(x)
    #fit1 = glm(y1 ~ x1)
    #coeffi <- coef(fit1)
    fit1 = cv.glmnet(x1, y1, alpha = 0, intercept = intercept) #ridge
    coeffi <- predict(fit1, type = 'coefficients', s = "lambda.min")

    if(intercept == T){
      w = 1/abs(as.numeric(coeffi))
    }
    
    else{
    w = 1/abs(as.numeric(coeffi)[-1])
    }   
  }

  if(relax == TRUE){
    obj = glmnet(x, y, alpha = 1, nlambda = 50, relax = T, intercept = F,
                       lambda.min.ratio = ifelse(nrow(x) < ncol(x), 0.01, 1e-04), 
                   lambda = NULL, standardize = TRUE)
  }


  else{
  obj = glmnet(x, y, alpha = alpha, nlambda = nlambda, relax = relax,
                       penalty.factor = w, intercept = intercept,
                       lambda.min.ratio = lambda.min.ratio,
                       lambda = lambda, standardize = standardize)
  }
  obj$relax = relax
  obj$nlambda = nlambda
  obj$intercept = intercept
  obj$x = x
  obj$y = y
  obj$gamma = gamma
  class(obj) = "lasso1"
  return(obj)
}


# NEW 
coef.lasso1 = function(object, s = NULL, gamma = NULL){
  beta.lasso1 = coef.lasso1.from.glmnet(object, s)
  if (object$intercept) 
    return(beta.lasso1)
  else return(beta.lasso1[-1, ])
}


coef.lasso1.from.glmnet = function(object, s=NULL) {
  
  if (length(object$lambda)==object$nlambda) {
    if (object$relax == F){
      class(object) = 'glmnet'
        e <- predict(object, s = s, type = 'coefficients')
        return(e)}
    
    else{
      class(object) = "relaxed"
        e <- predict(object, s = s, gamma = object$gamma, type = 'coefficients')
        return(e)
    }
  }
  
  else {
    min.lam = min(object$lambda)
    max.lam = max(object$lambda)
    svec = exp(seq(log(max.lam),log(min.lam),length=object$nlambda))
      if (object$relax == F){
        class(object) = 'glmnet'
        e <- predict(object, s = svec, type = 'coefficients')
        return(e)}
      
      else{
        class(object) = "relaxed"
        e <- predict(object, s = svec, gamma = object$gamma, type = 'coefficients')
        return(e)
      }
  }
}

predict.lasso1  = function (object, newx, s = NULL){
  if (missing(newx)) 
    newx = object$x
  if (object$intercept) 
    newx = cbind(rep(1, nrow(newx)), newx)
  return(newx %*% coef.lasso1(object, s))
}



