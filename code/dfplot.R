## Degrees of freedom simulations
source('newLasso.R')
source('generate_data.R')


# Set some overall simulation parameters
n = 70; p = 30 # Size of training set, and number of predictors
nval = n # Size of validation set
seed = 0 # Random number generator seed
s = 5 # Number of nonzero coefficients
beta.type = 2 # Coefficient type


## Degrees of freedom simulation
# we want to keep x fixed for the df calculations
set.seed(seed)
nrep = 500 # Number of repetitions

### RUN IF STUDENT
# set.seed(123)
# xy.obj = sim.xy_t(n,p,nval,rho=0.35,s=s,beta.type=beta.type,snr=0.7)
# x = matrix(unlist(xy.obj[1]), n, p)
# y = as.numeric(unlist(xy.obj[2]))
# mu = as.numeric(x %*% as.vector(unlist(xy.obj[6])))
# sigma = as.numeric(unlist(xy.obj[7]))

### RUN IF NORMAL
set.seed(123)
xy.obj = sim.xy(n,p,nval,rho=0.35,s=s,beta.type=beta.type,snr=0.7)
x = xy.obj$x
y = xy.obj$y
mu = as.numeric(x %*% xy.obj$beta)
sigma = xy.obj$sigma
nlam = 300

### RUN FOR BOTH 

# matrix that will contain results of simulations 
ip.las = matrix(0,nrep,(p+1))
ip.adapt =  matrix(0,nrep,(p+1))
ip.relax_las = ip.relax_las0 = matrix(0,nrep,(p+1))
ip.adap_las = matrix(0,nrep,(p+1))
p.fs = ip.bs = matrix(0,nrep,(p+1))

set.seed(123)
for (r in 1:nrep) {
  cat(r,"... ")
  eps = rnorm(n)*sigma
  y = mu + eps

  beta.las = coef(lasso1(x,y,intercept=FALSE, nlambda = nlam))
  nzs.las = colSums(beta.las != 0)
  j = nlam - rev(match(p:0, rev(nzs.las), NA)-1)
  yhat.las = (x %*% beta.las)[,j]

  #relax gamma 0.5
  beta.las = coef(lasso1(x,y,intercept=FALSE,nlambda =nlam, relax = T, gamma = 0.5))
  nzs.las = colSums(beta.las != 0)
  j = match(0:p, nzs.las, NA)
  yhat.relax_las = (x %*% beta.las)[,j]
  #
  #relax gamma 0
  beta.las = coef(lasso1(x,y,intercept=FALSE,nlambda =nlam, relax = T, gamma = 0))
  nzs.las = colSums(beta.las != 0)
  j = match(0:p, nzs.las, NA)
  yhat.relax_las0 = (x %*% beta.las)[,j]

  #adaptive
  beta.las = coef(lasso1(x,y,intercept=FALSE,nlambda =nlam, penalty = T))
  nzs.las = colSums(beta.las != 0)
  j = match(0:p, nzs.las, NA)
  yhat.adap_las = (x %*% beta.las)[,j]

  #relax gamma 0.5
  beta.las = coef(lasso1(x,y,intercept=FALSE,nlam=nlam, relax = T, gamma = 0.5))
  nzs.las = colSums(beta.las != 0)
  j = match(0:p, nzs.las, NA)
  yhat.relax_las = (x %*% beta.las)[,j]
  
  #adaptive 
  beta.las = coef(lasso1(x,y,intercept=FALSE,nlam=nlam, penalty = T))
  nzs.las = colSums(beta.las != 0)
  j = match(0:p, nzs.las, NA)
  yhat.adap_las = (x %*% beta.las)[,j]

  #elastic net
  beta.las = coef(lasso1(x,y,intercept=FALSE,nlam=nlam, alpha = 0.05))
  nzs.las = colSums(beta.las != 0)
  j = match(0:p, nzs.las, NA)
  yhat.el = (x %*% beta.las)[,j]

  
  yhat.fs = predict(fs(x,y,intercept=FALSE))[, 1:(p+1)]
  yhat.bs = predict(bs(x,y,intercept=FALSE))
  ip.las[r,] = colSums(yhat.las * eps)
  ip.reax_las[r,] = colSums(yhat.relax_las * eps)
  ip.adap_las[r,] = colSums(yhat.adap_las * eps)
  ip.fs[r,] = colSums(yhat.fs * eps)
  ip.bs[r,] = colSums(yhat.bs * eps)
  ip.el[r,] = colSums(yhat.el * eps)
}

df.las = colMeans(ip.las, na.rm=TRUE) / sigma^2
df.relax = colMeans(ip.reax_las, na.rm=TRUE) / sigma^2

df.adapt = colMeans(ip.adap_las, na.rm=TRUE) / sigma^2
df.el = colMeans(ip.el, na.rm=TRUE) / sigma^2

df.fs = colMeans(ip.fs, na.rm=TRUE) / sigma^2
df.bs = colMeans(ip.bs, na.rm=TRUE) / sigma^2

# Plot the results

dat = data.frame(x=rep(0:p, 7),
               y=c(df.las, df.relax, df.relax0, df.adapt, df.fs, df.bs, df.el),
               Method=factor(rep(c("lasso", "relaxed lasso gamma 0.5", "relaxed lasso gamma 0", "adaptive lasso", "forward", "best subset", "forward_glmnet"), rep(p+1, 7))))

dat = data.frame(x=rep(0:p, 6),
                 y=c(df.las, df.relax, df.adapt, df.fs, df.bs, df.el),
                 Method=factor(rep(c("lasso", "relaxed lasso", "adaptive lasso", "forward", "best subset", "elastic net"), rep(p+1, 6))))

ggplot(dat, aes(x=x,y=y,color=Method)) +
  xlab("Number of nonzero coefficients") +
  ylab("Degrees of freedom") +
  geom_line(lwd=0.5, color="#9b9b9b", linetype=3, aes(x,x)) +
  geom_line(lwd=1) + geom_point(pch=19) +
    theme_bw() + theme(legend.just=c(1,0), legend.pos=c(0.95,0.05))


