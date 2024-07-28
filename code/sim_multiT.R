source("generate_data.R")
source("bestsubset.R")

### testing for student

set.seed(0)
n = 100
p = 20
nval = n

# Regression functions: lasso, forward stepwise, and best subset selection
reg.funs = list()
reg.funs[["Lasso"]] = function(x,y) lasso1(x,y,intercept=FALSE)
reg.funs[["Stepwise"]] = function(x,y) fs(x,y,intercept=FALSE)
reg.funs[["Best subset"]] = function(x,y) bs(x,y,intercept=FALSE)
reg.funs[["Relaxed lasso"]] = function(x,y) lasso(x,y,intercept=FALSE,
                                                  nrelax=5)
reg.funs[["Adaptive Lasso"]] = function(x,y) lasso1(x,y,intercept=FALSE, penalty =T)

sim.obj.hisnr = sim.master_t(n,p,nval,reg.funs=reg.funs,nrep=10, seed=0,
                           beta.type=2,s=5, snr=0.7, rho = 0.35, verbose=TRUE, df = 2.1)


plot(sim.obj.hisnr, 1:3 , what ="risk")






reg.funs[["Lasso"]] = function(x,y) lasso1(x,y,intercept=FALSE)
reg.funs[["Relaxed lasso with gamma = 0.5"]] = function(x,y) lasso1(x,y, intercept=FALSE,
                                                                    relax = T, gamma = 0.5)
reg.funs[["Relaxed lasso with gamma = 1"]] = function(x,y) lasso1(x,y, intercept=FALSE,
                                                                  relax = T, gamma = 1)                                                 
reg.funs[["Adaptive Lasso"]] = function(x,y) lasso1(x,y,intercept=FALSE, penalty =T)
reg.funs[["Ridge"]] = function(x,y) lasso1(x,y, alpha = 0, intercept=FALSE)
reg.funs[["Forward stepwise"]] = function(x,y) lasso1(x,y,intercept=FALSE, gamma = 0, relax = T)
reg.funs[["Elastic net with alpha = 0.75"]] = function(x,y) lasso1(x,y,intercept=FALSE, alpha = 0.75)
reg.funs[["Elastic net with alpha = 0.25"]] = function(x,y) lasso1(x,y,intercept=FALSE, alpha = 0.25)
reg.funs[["Elastic net with alpha = 0.5"]] = function(x,y) lasso1(x,y,intercept=FALSE, alpha = 0.5)



sim_fig2 = sim.master(n,p,nval,reg.funs=reg.funs, rho = 0.35, nrep=3,seed=0, 
                      beta.type=2,s=5, snr = 0.7, verbose=TRUE)
plot(sim_fig2, what ="error")
plot(sim_fig2, main="SNR = 1")