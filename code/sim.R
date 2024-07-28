renv::init()

library(bestsubset)
library(gridExtra)
library(LaplacesDemon)
library(tidyverse)
library(glmnet)

source("newLasso.R")
source('generate_data.R')
source('bestsubset.R')

set.seed(0)
n = 100
p = 10
nval = n

# Regression functions: lasso, forward stepwise, and best subset selection
reg.funs = list()
reg.funs[["Best subset"]] = function(x,y) bs(x,y,intercept=FALSE)
reg.funs[["Stepwise"]] = function(x,y) fs(x,y,intercept=FALSE)
reg.funs[["Lasso"]] = function(x,y) lasso1(x,y,intercept=FALSE)
reg.funs[["Relaxed lasso"]] = function(x,y) lasso1(x,y,intercept=FALSE,
                                                   relax=T, gamma = 0.5)
reg.funs[["Adaptive Lasso"]] = function(x,y) lasso1(x,y,intercept=FALSE, penalty =T)
reg.funs[["Elastic Net"]] = function(x,y) lasso1(x,y,intercept=FALSE, alpha = 0.05)


# Run the master simulation function, for two different SNRs 
sim.obj.hisnr = sim.master(n,p,nval,reg.funs=reg.funs,nrep=20,seed=0, rho = 0.35,
                           beta.type=2,s=5,snr=3,verbose=TRUE)
sim.obj.losnr = sim.master(n,p,nval,reg.funs=reg.funs,nrep=20,seed=0, rho = 0.35,
                           beta.type=2,s=5,snr=0.1,verbose=TRUE)

# Print simulation results
sim.obj.hisnr
sim.obj.losnr

# Plot simulation results, excluding relaxed lasso 

fig1.1 = plot(sim.obj.hisnr, main="n = 100, p = 20, beta.type = 2, SNR = 3")
fig1.2 = plot(sim.obj.losnr, main="n = 100, p = 20, beta.type = 2, SNR = 0.1")

grid.arrange(fig1.1, fig1.2, ncol = 2 )

# Plot simulation results, including relaxed lasso 
par(mfrow=c(1,2))
plot(sim.obj.hisnr, main="SNR = 1") 
plot(sim.obj.losnr, main="SNR = 0.1") #legend.pos="topleft")

# same plot using relaxed lasso from Bestsubset Package 
reg.funs1 = list()
reg.funs1[["Best subset"]] = function(x,y) bs(x,y,intercept=FALSE)
reg.funs1[["Stepwise"]] = function(x,y) fs(x,y,intercept=FALSE)
reg.funs1[["Lasso"]] = function(x,y) lasso1(x,y,intercept=FALSE)
reg.funs1[["Relaxed lasso"]] = function(x,y) lasso1(x,y,intercept=FALSE,
                                                   nrelax = 5)
reg.funs1[["Adaptive Lasso"]] = function(x,y) lasso(x,y,intercept=FALSE, penalty =T)
reg.funs1[["Elastic Net 0.05"]] = function(x,y) lasso1(x,y,intercept=FALSE, alpha = 0.05)

sim.obj.hisnr1 = sim.master_t(n,p,nval,reg.funs=reg.funs1,nrep=20,seed=0, rho = 0.35,
                            beta.type=2,s=5,snr=3,verbose=TRUE)
sim.obj.losnr1 = sim.master(n,p,nval,reg.funs=reg.funs1,nrep=20,seed=0, rho = 0.35,
                            beta.type=2,s=5,snr=0.1,verbose=TRUE)
fig11.1 = plot(sim.obj.hisnr1, main="n = 100, p = 20, beta.type = 2, SNR = 3")
fig11.2 = plot(sim.obj.losnr1, main="n = 100, p = 20, beta.type = 2, SNR = 0.1")
fig11.1
grid.arrange(fig11.1, fig11.2, ncol = 2 )


########## FIG 5 ############
#Run the master simulation function, for two different SNRs
n = 100
p = 10
s = 5
nval = n

snr_val = exp( seq(log(0.05), log(6), length.out= 10))


setwd("C:/Users/Maria Vittoria/Desktop/EPFL/SCV/main-project-mvk/fig5")
for (i in 1:length(snr_val)){
  sim.master(n,p,nval,reg.funs=reg.funs, rho = 0.35, nrep=20,seed=0, file = paste0("prova", i),
             beta.type=2,s=5, snr = snr_val[i],verbose=TRUE)
}


file_list_fig5 = paste0("prova", 1:10)
a = plot.from.file(file_list_fig5, main=paste0("n=",n,", p=",p,", s=",s), what = "error" )
b = plot.from.file(file_list_fig5, main=paste0("n=",n,", p=",p,", s=",s), what = "prop")
c = plot.from.file(file_list_fig5, main=paste0("n=",n,", p=",p,", s=",s), what = "nonzero")
d = plot.from.file(file_list_fig5, main=paste0("n=",n,", p=",p,", s=",s), what = "F")
grid.arrange(a, b, ncol = 2)
grid.arrange(c, d, ncol = 2)


####### FIGURE 2: 
n = 100
p = 10
nval = n 

# fig 2 left 
sim_fig2 = sim.master(n,p,nval,reg.funs=reg.funs, rho = 0.35, nrep=20,seed=0, 
                  beta.type=2,s=5, snr = 0.7, verbose=TRUE)

fig2_sx = plot(sim_fig2, what ="risk", main = "SNR = 0.7, Cor=0.7")

# fig 2 right
sim_fig2.1 = sim.master(n,p,nval,reg.funs=reg.funs, rho = 0, nrep=20,seed=0, 
                      beta.type=2,s=5, snr = 2, verbose=TRUE)

fig2_dx = plot(sim_fig2.1, what ="risk", main = "SNR = 0.7, Cor=0.7")

grid.arrange(fig2_sx, fig2_dx, ncol = 2)



####### FIGURE 2 with beta type 5: 
# n = 70
# p = 30
# nval = n 
# 
# # fig 2 left 
# sim_fig2_bt5 = sim.master(n,p,nval,reg.funs=reg.funs, rho = 0.35, nrep=20,seed=0, 
#                       beta.type=5,s=5, snr = 0.7, verbose=TRUE)
# 
# fig2_bt5_sx = plot(sim_fig2_bt5, what ="risk")
# 
# # fig 2 right
# sim_fig2.1_bt5 = sim.master(n,p,nval,reg.funs=reg.funs, rho = 0, nrep=20,seed=0, 
#                         beta.type=5,s=5, snr = 2, verbose=TRUE)
# 
# fig2_bt5_dx = plot(sim_fig2.1_bt5, what ="risk")
# 
# grid.arrange(fig2_bt5_sx, fig2_bt5_dx, ncol = 2)


################################################################
######### MULTIVARIATE t part : 

########## FIG 5 - multivariate t ############
#Run the master simulation function, for two different SNRs
n = 100
p = 10
s = 5
nval = n


setwd("C:/Users/Maria Vittoria/Desktop/EPFL/SCV/main-project-mvk/fig5_t")
for (i in 1:length(snr_val)){
  sim.master_t(n,p,nval,reg.funs=reg.funs, rho = 0.35, nrep=20,seed=0, file = paste0("prova_T", i),
             beta.type=2,s=5, snr = snr_val[i],verbose=TRUE)
}


file_list_fig5_T = paste0("prova_T", 1:10)
cc = plot.from.file(file_list_fig5_T , main=paste0("n=",n,", p=",p,", s=",s), what = "nonzero")
dd = plot.from.file(file_list_fig5_T , main=paste0("n=",n,", p=",p,", s=",s), what = "F")
grid.arrange(cc, dd, ncol = 2)

#############################

