library(MASS)
library(LaplacesDemon)

generate_beta <- function(n, p, s, type) {
  beta = rep(0,p)
  if (type==1) {
    index = round(seq(1, p, by=p/s))
    beta[index] = 1
  }
  else if (type==2) {
    beta[1:s] = 1
  }
  else if (type==3) {
    values = seq(10, 0.5, length.out = s)
    beta[1:s] = values
  }
  else {
    index = (s+1):p
    values = 0.5^(index-s)
    beta[1:s] = 1
    beta[(s+1):p] = values
  }
  beta <- as.vector(beta)
  return(beta)
}



enlist1 <- function (...) 
{
  result <- list(...)
  if ((nargs() == 1) & is.character(n <- result[[1]])) {
    result <- as.list(seq(n))
    names(result) <- n
    for (i in n) result[[i]] <- get(i)
  }
  else {
    n <- sys.call()
    n <- as.character(n)[-1]
    if (!is.null(n2 <- names(result))) {
      which <- n2 != ""
      n[which] <- n2[which]
    }
    names(result) <- n
  }
  result
}

choose.tuning.params = function(obj) {
  N = length(obj$err.test) # Number of methods
  nrep = nrow(obj$err.test[[1]]) # Number of repetitions
  tun.val = matrix(NA,nrep,N)
  tun.ora = rep(NA,N)

  # Validation tuning: based on validation error
  for (i in 1:nrep) {
    for (j in 1:N) {
      tun.val[i,j] = which.min(obj$err.val[[j]][i,])
      if (length(tun.val[i,j]) == 0) tun.val[i,j] = NA
    }
  }

  # Oracle tuning: based on average test error
  for (j in 1:N) {
    tun.ora[j] = which.min(colMeans(obj$err.test[[j]], na.rm=TRUE))
  }

  return(c(obj,enlist1(tun.val,tun.ora)))
}


generate_sigma <- function(p, rho){
  # create matrix with NA with p rows and p columns 
  sigma_ = matrix(data=NA, nrow=p, ncol=p)  
  # fill the elements with values  
  for(j in 1:p){
    for(i in 1:p) {
      sigma_[i,j] <- rho^(abs(i-j))
    }
  }
  return(sigma_)
}

# this should be a bit more computationally efficient
generate_sigma1 <- function(p, rho) {
  matrix(sapply(1:(p), function(i) sapply(1:(p), function(j) rho^abs(i-j))), nrow=p, ncol=p)
}


generate_x <- function(n, p, sigma_) {
  mu <- as.vector(rep(0, p))
  x <- MASS::mvrnorm(n, mu, sigma_)
  return(as.matrix(x))
}

sigma_square <- function(n, beta, sigma, snr){
  #diag_mat <- diag(n)
  return(diag(x = as.numeric((beta %*% sigma %*% beta) / snr), n, n))
}

generate_y <- function(n, p, s, type, rho, snr, x, sigma){
  beta <- generate_beta(n, p, s, type)
  #sigma <- generate_sigma(p, rho)
  #x <- generate_x(n, p, rho, sigma)
  mu_ <- x %*% beta
  sigma_n_by_n <- sigma_square(n, beta, sigma, snr)
  y <- as.vector(MASS::mvrnorm(1, mu_, sigma_n_by_n))
 return(y)
}

generate_student_x <- function(n, p, df, sigma){
  x <- LaplacesDemon::rmvt(n=n, mu = rep(0, p), sigma, df=df)
  #u <- rchisq(n = n, df = df)
  #normal <- generate_x(n, p, sigma)
  #df = as.matrix(rep(df, n))
  #constant = (df / as.matrix(u))^(1/2)
  #constant = as.matrix(rep(constant, p))
  #dim(constant) <- c(p, n)
  #return(as.matrix(t(constant) * normal))
  return(as.matrix(x))
}

#sigma <- generate_sigma(20, 0.5)
#generate_student_x(100, 20, 3, sigma)

generate_y_t <- function(n, p, s, beta, rho, snr, x, sigmasq){
  mu_ <- x %*% beta
  sigma_n_by_n <- diag(sigmasq, n, n)
  y <- as.numeric(MASS::mvrnorm(1, mu_, sigma_n_by_n))
 return(y)
}

sim.xy_t = function(n, p, nval, rho=0, s=5, beta.type=1, snr=1, df = 3) {
  
  # Generate predictors
  Sigma <- generate_sigma1(p, rho)
  x <- generate_student_x(n, p, df, Sigma)
  xval <- generate_student_x(nval, p, df, Sigma)
  beta <- generate_beta(n, p, s, beta.type)

  # Set snr based on sample variance on infinitely large test set
  scale = as.numeric(t(beta) %*% Sigma %*% beta) #scale of X * beta
  variance = (df/(df - 2)) * scale #variance of X * beta
  sigma = variance/snr #variance of Y is sigma

  # Generate responses
  y <- generate_y_t(n, p, s, beta, rho, snr, x, sigma)
  yval <- generate_y_t(nval, p, s, beta, rho, snr, xval, sigma)

  sigma <- sqrt(sigma) #added sqrt bc of the og function
  
  # Sigma set to the true variance of the X matrix:
  Sigma <- df/(df - 2) * Sigma
  return(list(x,y,xval,yval,Sigma,beta,sigma))
}

sim.master_t = function(n, p, nval, reg.funs, nrep=50, seed=NULL, verbose=FALSE,
                      file=NULL, file.rep=5, rho=0, s=5, beta.type=1, snr=1, df = 3, sforpredict = 0.84) {
  
  this.call = match.call()
  if (!is.null(seed)) set.seed(seed)

  N = length(reg.funs)
  reg.names = names(reg.funs)
  if (is.null(reg.names)) reg.names = paste("Method",1:N)

  err.train = err.val = err.test = prop = risk = nzs = fpos = fneg = F1 = 
    opt = runtime = vector(mode="list",length=N)
  names(err.train) = names(err.val) = names(err.test) = names(prop) =
    names(risk) = names(nzs) = names(fpos) = names(fneg) = names(F1) = 
    names(opt) = names(runtime) = reg.names
  for (j in 1:N) {
    err.train[[j]] = err.val[[j]] = err.test[[j]] = prop[[j]] = risk[[j]] =
      nzs[[j]] = fpos[[j]] = fneg[[j]] = F1[[j]] = opt[[j]] = runtime[[j]] =
      matrix(NA,nrep,1)
  }
  filled = rep(FALSE,N)
  err.null = risk.null = sigma = rep(NA,nrep)

  # Loop through the repetitions
  for (i in 1:nrep) {
    if (verbose) {
      cat(sprintf("Simulation %i (of %i) ...\n",i,nrep))
      cat("  Generating data ...\n")
    }

    # Generate x, y, xval, yval
    xy.obj = sim.xy_t(n,p,nval,rho,s,beta.type,snr,df) 
    xy.obj$x = matrix(unlist(xy.obj[1]), n, p)
    xy.obj$xval = matrix(unlist(xy.obj[3]), nval, p)
    xy.obj$y = as.numeric(unlist(xy.obj[2]))
    xy.obj$yval = as.numeric(unlist(xy.obj[4]))
    xy.obj$Sigma = matrix(unlist(xy.obj[5]), p, p)
    xy.obj$beta = as.vector(unlist(xy.obj[6]))
    xy.obj$sigma= as.numeric(unlist(xy.obj[7]))


    risk.null[i] = diag(t(xy.obj$beta) %*% xy.obj$Sigma %*% xy.obj$beta)
    err.null[i] = risk.null[i] + xy.obj$sigma^2
    sigma[i] = xy.obj$sigma

    # Loop through the regression methods
    for (j in 1:N) {
      if (verbose) {
        cat(sprintf("  Applying regression method %i (of %i) ...\n",
                    j,N))
      }

      tryCatch({
        # Apply the regression method in hand
        runtime[[j]][i] = system.time({
          reg.obj = reg.funs[[j]](xy.obj$x,xy.obj$y)
        })[1]

        # Grab the estimated coefficients, and the predicted values on the
        # training and validation sets
        betahat = as.matrix(coef(reg.obj))
        m = ncol(betahat); nc = nrow(betahat)

        # Check for intercept
        if (nc == p+1) {
          intercept = TRUE
          betahat0 = betahat[1,]
          betahat = betahat[-1,]
        }
        else intercept = FALSE

        muhat.train = as.matrix(predict(reg.obj,xy.obj$x))
        muhat.val = as.matrix(predict(reg.obj,xy.obj$xval))

        # Populate empty matrices for our metrics, of appropriate dimension
        if (!filled[j]) {
          err.train[[j]] = err.val[[j]] = err.test[[j]] = prop[[j]] =
            risk[[j]] = nzs[[j]] = fpos[[j]] = fneg[[j]] = F1[[j]] = opt[[j]] =
            matrix(NA,nrep,m)
          filled[j] = TRUE
          # N.B. Filling with NAs is important, because the filled flag could
          # be false for two reasons: i) we are at the first iteration, or ii)
          # we've failed in all previous iters to run the regression method
        }

        
        # Record all of our metrics

        err.train[[j]][i,] = colMeans((muhat.train - xy.obj$y)^2)
        
        err.val[[j]][i,] = colMeans((muhat.val - xy.obj$yval)^2)
        delta = betahat - xy.obj$beta
       
        risk[[j]][i,] = diag(t(delta) %*% xy.obj$Sigma %*% delta)
        if (intercept) risk[[j]][i,] = risk[[j]][i,] + betahat0^2
        err.test[[j]][i,] = risk[[j]][i,] + xy.obj$sigma^2
        prop[[j]][i,] = 1 - err.test[[j]][i,] / err.null[i]
        nzs[[j]][i,] = colSums(betahat!=0)
        tpos = colSums((betahat!=0)*(xy.obj$beta!=0))
        fpos[[j]][i,] = nzs[[j]][i,]-tpos
        fneg[[j]][i,] = colSums((betahat==0)*(xy.obj$beta!=0))
        F1[[j]][i,] = 2*tpos/(2*tpos+fpos[[j]][i,]+fneg[[j]][i,])
        opt[[j]][i,] = (err.test[[j]][i,] - err.train[[j]][i,]) /
          err.train[[j]][i,]
      }, error = function(err) {
        if (verbose) {
          cat(paste("    Oops! Something went wrong, see error message",
                    "below; recording all metrics here as NAs ...\n"))
          cat("    ***** Error message *****\n")
          cat(sprintf("    %s\n",err$message))
          cat("    *** End error message ***\n")
        }
        # N.B. No need to do anything, the metrics are already filled with NAs
      })
    }

    # Save intermediate results?
    if (!is.null(file) && file.rep > 0 && i %% file.rep == 0) {
      saveRDS(enlist1(err.train,err.val,err.test,err.null,prop,risk,risk.null,
                     nzs,fpos,fneg,F1,opt,sigma,runtime),file=file)
    }
  }

  # Save results now (in case of an error that might occur below)
  out = enlist1(err.train,err.val,err.test,err.null,prop,risk,risk.null,nzs,fpos,
               fneg,F1,opt,sigma,runtime)
  if (!is.null(file)) saveRDS(out, file)

  # Tune according to validation error, and according to test error
  out = choose.tuning.params(out)
  # Save final results
  out = c(out, list(rho=rho,s=s,beta.type=beta.type,snr=snr,call=this.call))
  class(out) = "sim"
  if (!is.null(file)) { saveRDS(out, file); invisible(out) }
  else return(out)
}

