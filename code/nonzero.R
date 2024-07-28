# Function to simulate from multivariate t non - adjusting for the degrees of 
# freedom 
sim.xy_t_df = function(n, p, nval, rho=0, s=5, beta.type=1, snr=1, df = 3) {
  
  # Generate predictors
  Sigma <- generate_sigma1(p, rho) # Same as normal
  x <- generate_student_x(n, p, df, Sigma)
  xval <- generate_student_x(nval, p, df, Sigma)
  beta <- generate_beta(n, p, s, beta.type)
  
  # Set snr based on sample variance on infinitely large test set
  scale = as.numeric(t(beta) %*% Sigma %*% beta) #scale of X * beta
  sigma = (df/(df - 2)) * scale #variance of X * beta
  # NOW WE DON'T ADJUST THE VARIANCE IN ORDER 
  # TO HAVE A CERTAIN SNR
  
  # Generate responses
  y <- generate_y_t(n, p, s, beta, rho, snr, x, sigma)
  yval <- generate_y_t(nval, p, s, beta, rho, snr, xval, sigma)
  
  sigma <- sqrt(sigma) #added sqrt bc of the og function
  
  # Sigma set to the true variance of the X matrix:
  Sigma <- df/(df - 2) * Sigma
  return(list(x,y,xval,yval,Sigma,beta,sigma))
}


#sim mastter for multivariate t not adjusted for the dof
sim.master_t_NoSNR = function(n, p, nval, reg.funs, nrep=50, seed=NULL, verbose=FALSE,
                              file=NULL, file.rep=5, rho=0, s=5, beta.type=1, snr = NA, df = 3, sforpredict = 0.84) {
  
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
    xy.obj = sim.xy_t_df(n,p,nval,rho,s,beta.type,snr,df) 
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



############## SIMULATION 
n = 100
p = 10
s = 5
nval = n

df_vec = seq(2.1, 10, length.out = 10)
store_results = matrix(NA, 6, 10)

for (i in 1:length(df_vec)){
  prova = sim.master_t_NoSNR(n,p,nval,reg.funs=reg.funs, rho = 0.35, nrep=20,seed=0, file = paste0("prova_T", i),
                             beta.type=2,s=5,verbose=TRUE, df = df_vec[i])
  nonz = tune.and.aggregate(prova, prova$nzs)
  store_results[,i] = nonz$z.val.ave
}

store_results = data.frame(store_results)
store_results = store_results[1:5, ]
colnames(store_results) = round(df_vec ,2)

store_results$method = c("Best Subset", "Stepwise", "Lasso", "Relaxed Lasso", "Adaptive Lasso")

store_results_log = pivot_longer(store_results, 1:10 )

ggplot(store_results_log) +
  geom_line(aes(x = as.numeric(name), y = value, group = method, 
                col = method) ,linewidth = 1) +
  geom_hline(yintercept=5, linetype = "dashed") +
  labs(y = "Number of nonzero coefficients", x = "degrees of freedom") + 
  ggtitle("Number of nonzeros for different degrees of freedom of multivariate-t") +
  scale_x_continuous(breaks = round(df_vec,2 )) + 
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.title = element_text(face = "bold")) 

