# Notes of RR

# code for adaptive lasso:
# https://stats.stackexchange.com/questions/231643/lasso-vs-adaptive-lasso 
# code for relaxed adaptive lasso:
# https://stats.stackexchange.com/questions/614896/relaxed-adaptive-lasso 
# code for elastic net (tuning alpha):
# https://rpubs.com/jmkelly91/881590

# If we know the distribution of the coordinates
# https://stats.stackexchange.com/questions/80809/linear-combination-of-multivariate-t-distribution



# for Adaptive Lasso, we just need to add a penalty term (penalty.factor)
# glmnet already has this implemented
# Options: 1) Rewrite the Lasso function of the authors
#          2) Just call glmnet with penalty.factor

# search for "penalty.factor" to find implementation of Adaptive Lasso
# https://glmnet.stanford.edu/articles/glmnet.html 

# Questions: 1) What is the difference between the relaxed lasso
#               of the authors and the glmnet one?
#               The authors use nrelax = 10

# search for "relax" to find implementation of Relaxed Lasso
# https://glmnet.stanford.edu/articles/relax.html 

# Remarks: 1) I think that we should rewrite the Lasso function
#             of the authors in lasso.R file
#             Probably the Relaxed Lasso and the Adaptive lasso
#             were implemented after the publication of the study
#             so they do not use them
#             But for us it would be easier to use them, so we can
#             just to modify some lines of their code
#          2) Elastic net is already implemented by them!
#             Check glmnet alpha and lambda parameters on how to
#             get to the Elastic net case

# Ideas: 1) We can study the Adaptive Relaxed Lasso method
#           (would be very easy with the modification of their
#            code as suggested above)
#        2) Try something else besides Gaussian samples
#        3) Use different metrics for the evaluation of the goodness of fit



##Testing 
n = 100
p = 20
nval = n
type = 2
s = 5
rho = 0.5
snr = 0.1

sigma <- generate_sigma(p, rho)
x <- generate_x(n, p, sigma)
y <- generate_y(n, p, s, type, rho, snr, x, sigma)

x_test <- generate_x(n, p, sigma)
y_test <- generate_y(n, p, s, type, rho, snr, x_test, sigma)

lasso_ = lasso1(x, y, relax = T, intercept = F, alpha = 1, , gamma = 0.8)
print('ok')
predict_ = predict(lasso_, x, s = 0.84)
print(predict_)
print('predict done')
x %*% coef.lasso1(lasso_, s = 0.84) - predict_
print('if all zeros then good')


set.seed(0)
n = 70
p = 30
nval = n

# Regression functions: lasso, forward stepwise, and best subset selection
reg.funs = list()
reg.funs[["Lasso"]] = function(x,y) lasso1(x,y,intercept=FALSE)
reg.funs[["Relaxed lasso with gamma = 0.5"]] = function(x,y) lasso1(x,y, intercept=FALSE,
                                                 relax = T, gamma = 0.5)
reg.funs[["Relaxed lasso with gamma = 0.25"]] = function(x,y) lasso1(x,y, intercept=FALSE,
                                                  relax = T, gamma = 0.25)                                                 
reg.funs[["Adaptive Lasso"]] = function(x,y) lasso1(x,y,intercept=FALSE, penalty =T)
reg.funs[["Ridge"]] = function(x,y) lasso1(x,y, alpha = 0, intercept=FALSE)
reg.funs[["Relaxed lasso with gamma = 0"]] = function(x,y) lasso1(x,y,intercept=FALSE, gamma = 0, relax = T)
reg.funs[["Elastic net with alpha = 0.75"]] = function(x,y) lasso1(x,y,intercept=FALSE, alpha = 0.75)
reg.funs[["Elastic net with alpha = 0.25"]] = function(x,y) lasso1(x,y,intercept=FALSE, alpha = 0.25)
reg.funs[["Elastic net with alpha = 0.5"]] = function(x,y) lasso1(x,y,intercept=FALSE, alpha = 0.5)
reg.funs[["Elastic net with alpha = 0.05"]] = function(x,y) lasso1(x,y,intercept=FALSE, alpha = 0.05)


 
sim_fig2 = sim.master_t(n,p,nval,reg.funs=reg.funs, rho = 0, nrep=3,seed=0, 
                  beta.type=2,s=5, snr = 2, verbose=TRUE)
plot(sim_fig2, what ="error")
plot(sim_fig2, what = 'risk')

