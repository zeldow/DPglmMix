library(MASS) # for mvnorm()
library(mvtnorm)
setwd('C:/Users/aoganisi/Box Sync/R Packages/DPglmMix/code')
source("BNP_GLM_FunctionSourceCode.r")

################################################################################
### 0) simulate data
################################################################################
set.seed(1)
n <- 500
# setting 2
x<-seq(1,10*pi, length.out = n)
y<-rnorm( n = length(x), sin(.5*x), .09*x)
plot(x,y)

# simulate data set with one outcome, y and one continuous covariate, x.
d <- data.frame(x=x, y=y)

################################################################################
### 1) Run DP mixture of GLMs
################################################################################
# should take < 10 minutes
burnin<-100
gibbs_iter<-1000
K <- 10
alpha <- 1000 # pick very high value to keep class sizes from shrinking to 0

DPglm_res<-DPglmMix(d = d, y = 'y', x = 'x', gibbs_iter = gibbs_iter, K = K, a = alpha)

################################################################################
### 2) Estimate Conditional Mean 
################################################################################
c_shell <- DPglm_res[[1]]
beta_shell <- DPglm_res[[2]]
mu_x_shell <- DPglm_res[[3]]
phi_x_shell <- DPglm_res[[4]]
psi_shell <- DPglm_res[[5]]

# x-values at which to compute mean of y
xvec<-seq(0,32,.1)

yvec<-numeric(length(xvec))

p<-table(c_shell[,gibbs_iter])/n
mupost<-rowMeans(mu_x_shell[,burnin:gibbs_iter])

M <- length(burnin:gibbs_iter)
yvec<-vector(mode = 'list', length = 3)
yvec[[1]]<-yvec[[2]]<-yvec[[3]]<- numeric(length(xvec))


# need to wrap this part into a function that can be called 
# after calling DPglmMix() to estimate parameters.

for(i in 1:length(xvec)){
  xi<-xvec[i]
  yvec_mat<-matrix(NA, nrow=M, ncol=20)
  
  cond<-numeric(length = M)
  meanv<-numeric(length = M)
  for(m in 1:M){
    cond_num<-0
    cond_denom<-0
    mean_num<-0
    mean_denom<-0
    for (c in c(1:K)){
      cond_num<-cond_num + rnorm(1, mean = c(1, xi) %*% beta_shell[[c]][ , burnin + m -1], sd=sqrt(psi_shell[c,burnin + m -1]) ) * dnorm(xi, mean=mu_x_shell[c, burnin + m -1], sqrt(phi_x_shell[c, burnin + m -1]))
      cond_denom<-cond_denom + dnorm(xi, mean=mu_x_shell[ c , burnin + m -1], sqrt(phi_x_shell[c, burnin + m -1]))
      
      mean_num<-cond_num + (c(1, xi) %*% beta_shell[[c]][ , burnin + m -1] ) * dnorm(xi, mean=mu_x_shell[c, burnin + m -1], sqrt(phi_x_shell[c, burnin + m -1]))
      mean_denom<-cond_denom + dnorm(xi, mean=mu_x_shell[ c , burnin + m -1], sqrt(phi_x_shell[c, burnin + m -1]))
      
    }
    cond[m] <- cond_num / cond_denom
    meanv[m] <- mean_num/mean_denom
  }
  yvec[[1]][i] <- mean(meanv,na.rm = T)
  yvec[[2]][i] <- quantile(cond, .05, na.rm = T)
  yvec[[3]][i] <- quantile(cond, .95, na.rm = T)
  
}

################################################################################
### 3) Plot Results
################################################################################


par(mfrow=c(3, 1))

plot(x,y, main='Original Data')

plot(x, y, col=c_shell[,gibbs_iter], pch=20, cex=.8, main = 'Classification Results')

plot(x, y, pch=20, cex=.5, main = 'Regression Line with Credible Interval')
polygon(x= c(xvec,rev(xvec)),y= c(yvec[[3]],rev(yvec[[2]])), col="grey", border=NA)
points(x, y, pch=20, cex=.5)
lines(xvec, sin(.5*xvec), col='blue', lwd=2)
lines(xvec, yvec[[1]], col='red', lwd=2)





