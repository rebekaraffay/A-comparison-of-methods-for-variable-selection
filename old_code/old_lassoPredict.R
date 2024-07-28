coef.lasso1.from.glmnet = function(object, s=NULL) {
  class(object) = "glmnet"
  if (length(object$lambda)==object$nlambda) {
    return(glmnet::coef.glmnet(object,s=s))
  }
  else {
    min.lam = min(object$lambda)
    max.lam = max(object$lambda)
    svec = exp(seq(log(max.lam),log(min.lam),length=object$nlambda))
    return(glmnet::coef.glmnet(object,s=svec))
  }
}

coef.lasso1 = function(object, s = NULL, gamma = NULL){
  beta.lasso1 = coef.lasso1.from.glmnet(object, s)
  if (object$intercept) 
    return(beta.lasso1)
  else return(beta.lasso1[-1, ])
}

predict.lasso1  = function (object, newx, s = NULL){
  if (missing(newx)) 
    newx = object$x
  if (object$intercept) 
    newx = cbind(rep(1, nrow(newx)), newx)
  return(newx %*% coef.lasso1(object, s))
}