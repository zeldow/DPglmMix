log_cond_post_beta<-function(beta, y, x, psi, theta, s){
  lik<-sum(dnorm(y, mean = cbind(1, x) %*% beta, sd = sqrt(psi),log = T  ))
  pr <- dnorm(beta, mean = theta, sd=sqrt(s), log=T)
  post <- lik + pr
  return(post)
}

log_cond_post_alpha<-function(alpha, n, k){
  post <- (k - 3/2) * log(alpha) - 1/(2*alpha) + lgamma(alpha) - lgamma(n + alpha)
  return(post)
}

rcond_post_psi <- function(beta, y, x, g1, b1){
  n_k <- length(y)
  if(n_k==0){ 
    shape_k <- g1
    rate_k <- b1
  }else{
    shape_k <- g1 + n_k/2
    rate_k <- .5*(  sum( (y - cbind(1, x) %*% beta)^2 ) )  + b1
  }
  psi_post <- invgamma::rinvgamma(n = 1, shape = shape_k, rate = rate_k)
  return(psi_post)
}

rcond_post_phi <- function(mu_x, x, g2, b2){
  n_k <- length(x)
  if(n_k==0){ 
    shape_k <- g2
    rate_k <- b2
  }else{
    shape_k <- g2 + n_k/2
    rate_k <- .5*(  sum( (x - mu_x )^2 ) )  + b2
  }
  phi_post <- invgamma::rinvgamma(n = 1, shape = shape_k, rate = rate_k)
  return(phi_post)
}

metrop_hastings<-function(x_0, iter=1, log_post_density,
                          proposal_dist = function(x){ 
                            mvrnorm(1, mu = x, Sigma = 100*diag(length(x))  )
                          }, 
                          ... ){
  for(i in 1:iter){
    # draw from proposal distribution
    x_star <- proposal_dist(x_0)
    
    # calculate ratio of conditional posterior densities
    r_num <- do.call(log_post_density, c(list(x_star), list(...)) )
    r_denom <- do.call(log_post_density, c(list(x_0), list(...)) )
    r <- exp(r_num - r_denom)
    rmin<-min(r,1)
    if(is.na(rmin)) browser()
    # accept / reject proposal
    if(rbinom(1,1,rmin)==1){ 
      x_0<-x_star
    }
  }
  return(x_0)
}

DPglmMix<-function(d, y, x, gibbs_iter, K, a=NULL, 
                   theta=0, s=1000^2, lambda=0, tau=100^2, 
                   g1=1, b1=.1, g2=1, b2=.1){
  n<-nrow(d)
  
  y <- d[,y]
  x <- d[,x]
  
  # proposal dist for MH sampling of concentration parameter posterior
  prop_alpha<-function(x) rgamma(n = 1, shape = 1, rate = .001)
  
  ###------------------------------------------------------------------------###
  #### 0 - Create Shells for storing Gibbs Results                          ####
  ###------------------------------------------------------------------------###
  c_shell<-matrix(NA, nrow=n, ncol=gibbs_iter) # shell for indicators
  beta_shell <- vector(mode = 'list', length = K) # shell for betas
  psi_shell <- matrix(NA, nrow= K, ncol=gibbs_iter) # shell for Y's cond. var
  
  mu_x_shell <- matrix(NA, nrow=K, ncol=gibbs_iter) # shell for mean of Xs
  phi_x_shell <- matrix(NA, nrow= K, ncol=gibbs_iter) # shell for var of Xs
  
  mu_y_shell <- numeric(length = K) # shell for predicted Ys
  
  ###------------------------------------------------------------------------###
  #### 1 - Initialize Gibbs Sampler by setting starting values              ####
  ###------------------------------------------------------------------------###
  c_shell[,1] <- sample(x = 1:K, size = n, replace = T)
  
  for(k in 1:K){
    beta_shell[[k]] <- matrix(NA, ncol=gibbs_iter, nrow=2)
    beta_shell[[k]][,1]<-c(1,1)
  }
  phi_x_shell[1:K,1]<-100
  
  # if user sets concentration parameter, use it. If not, estimate it
  if(is.null(a)){
    alpha<-numeric(length = gibbs_iter)
    alpha[1] <- 2000
  }else{
    alpha<-numeric(length = gibbs_iter)
    alpha[1:gibbs_iter] <- a
  }
  
  ###------------------------------------------------------------------------###
  #### 2 - Gibbs Sampler                                                    ####
  ###------------------------------------------------------------------------###
  
  for(i in 2:gibbs_iter){
    ####### 2.0 - stratify data based on previous class indicator           ####
    class_ind<-c_shell[,i-1]
    nvec<-table(factor(class_ind,levels = 1:K) )
    
    ####### 2.1 - update concentration parameter if user did not select one ####
    if(is.null(a)){
      k<-length(unique(class_ind))
      alpha[i]<-metrop_hastings(x_0 = alpha[i-1], iter = 10,
                                log_post_density = log_cond_post_alpha,
                                proposal_dist = prop_alpha,
                                n = n, k = k )
    }
    
    ####### 2.2 - update cluster parameters                                 ####
    for(k in 1:K){
      ck_ind<-class_ind==k
      y_ck<-y[ck_ind] 
      x_ck<-x[ck_ind]
      
      # sample psi posteriors
      psi_shell[k, i] <- rcond_post_psi(beta_shell[[k]][,i-1], 
                                        y_ck, x_ck, g1, b1)
      
      # sample beta posteriors (M-H used for easy generality to other links )
      beta_shell[[k]][,i]<-metrop_hastings(x_0 = beta_shell[[k]][,i-1], 
                                           log_post_density=log_cond_post_beta, 
                                           y=y_ck, x=x_ck, psi=psi_shell[k, i], 
                                           theta=theta, s=s )
      
      
      # sample covariate parameters
      mu_mean <- (1/(1/tau + nvec[k]/phi_x_shell[k,i-1]))*
        (lambda/tau + sum(x_ck)/phi_x_shell[k,i-1])
      mu_sd <- sqrt( (1/tau + nvec[k]/phi_x_shell[k,i-1])^(-1) )
      
      mu_x_shell[k, i] <- rnorm(n = 1, mean = mu_mean, sd =  mu_sd)
      phi_x_shell[k, i] <- rcond_post_phi(mu_x_shell[k, i], x_ck, g2, b2)
    }
    
    ####### 2.4 - update class indicators                                   ####
    for(j in 1:n){
      n_min_j<-table(factor(class_ind[-j],levels = 1:K) )
      
      for(k in 1:K){  
        mu_y_shell[k] <- t(c(1, x[j])) %*% beta_shell[[k]][,i] 
      }
      
      pr <- log((n_min_j + alpha[i]/K)/(n+alpha[i]-1))
      lk_y <- dnorm(y[j], mean = mu_y_shell, 
                    sd = sqrt(psi_shell[,i]), log = T)
      
      lk_x <- dnorm(x[j], mean = mu_x_shell[,i], 
                    sd = sqrt(phi_x_shell[,i]), log = T)
      
      log_p_vec <- pr + lk_y + lk_x
      p_vec<-exp(log_p_vec - max(log_p_vec)) # rescale for numeric reasons
      
      c_shell[j,i]<-sample(x = 1:K, size = 1, replace = T, prob = p_vec)
    }
    
  }
  #### 2.5 Output Results                                                   ####
  results<-list(c_shell, beta_shell, mu_x_shell, phi_x_shell, psi_shell)
  return(results)
}

