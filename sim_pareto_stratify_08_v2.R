rm(list=ls())

#install.packages("VGAM")
#install.packages("mvtnorm")
#install.packages("LaplacesDemon")
#install.packages("cubature")
#install.packages("rmutil")

#library(VGAM)
library(mvtnorm)
library(LaplacesDemon)
library(cubature)
#library(rmutil)

m_iteration = 1000 # number of random samples to obtain the optimal mean of factor variable
n_iteration = 100000 # number of iterations
n_obligor = 1000 # number of obligors

# exposure <- sort(rpareto(n_obligor, 1.2), decreasing =  T )
# dp <- 0.02*runif(n_obligor)
# write.table(dp, "./dp.txt", sep="\t")

exposure <- read.table("./exposure_pareto_25.txt", sep="\t")$x[1:n_obligor]
dp <- read.table("./dp.txt", sep="\t")$x[1:n_obligor]

#calculate expected loss
expected_loss <- sum(exposure * dp)
sd_loss <- sqrt(sum( exposure**2 * dp *(1-dp))) 

# set the values of the parameters and functions
mu <- c(0,0,0)
sigma<-matrix(c(1,0.5,0.5,0.5,1,0.5,0.5,0.5,1),nrow=3)
b <- as.matrix(rep(1/3,3))
c <- as.numeric(sqrt(1+(t(b)%*%sigma%*%b)))*qnorm(dp)


num_take_all_obligor <- 4
threshold <- sum(exposure[1:num_take_all_obligor]) - 0.1
###########################
## crude MC              ##
###########################

tic_cmc = Sys.time()

psi <- rmvnorm(mean=mu, sigma, n = n_iteration)
loss <- rep(0, n_iteration)                     # realized loss at each iteration
for (i in 1:n_iteration){
  cond_dp <- pnorm(as.numeric(psi[i,]%*%b)+c)
  
  y <- rbern(n_obligor,cond_dp)
  loss[i] <-  sum(exposure * y)
}

indicate <- as.numeric(loss > threshold)
excess_prob_cmc <- mean(indicate)
SE_EP_cmc <- sd(indicate)/sqrt(n_iteration)  
ES_cmc <- (1/excess_prob_cmc)*mean(indicate*loss)
SE_ES_cmc <- (1/excess_prob_cmc)*sd(indicate*loss)/sqrt(n_iteration)

toc_cmc = Sys.time()
elapsed_time_cmc = as.numeric(toc_cmc - tic_cmc,units = "secs")


###########################################
## 2-step exponential twisting using CE  ##
###########################################
tic_et = Sys.time()

upper = 0.5
lower = 0.0
psi <- rmvnorm(mean=mu, sigma, n = m_iteration)
loss <- rep(0, m_iteration)                     # realized loss at each iteration
likelihood_ratio <- rep(0, m_iteration)         # likelihood ratio of given prob. measure to IS prob. measure

for(i in 1: m_iteration){
  cond_dp <- pnorm(as.numeric(psi[i,]%*%b)+c)
  expected_loss <- sum(cond_dp*exposure)
  if (expected_loss > threshold) {
    t_1 <- 0
  } else {
    expected_loss_t <-function(t){
      	q <-cond_dp/(cond_dp+ exp(-t*exposure)*(1-cond_dp))
      	sum(exposure*q) - threshold
    }
    t<-uniroot(expected_loss_t, interval=c(lower,upper),tol = 1e-10)  
    t_1 <-t$root
  }
  
  twisted_dp <- exp(t_1*exposure)*cond_dp/(exp(t_1*exposure)*cond_dp+(1-cond_dp))
  
  y <- rbern(n_obligor,twisted_dp)  
  loss[i] <- sum(y * exposure)
  likelihood_ratio[i] <- prod(cond_dp * exp(exposure*t_1) + (1-cond_dp)) * exp( - loss[i] * t_1)
}

weight <- matrix(as.numeric(loss > threshold)*likelihood_ratio, now <- 1)
mu_star <- weight %*% psi / sum(weight)

## doing the combined IS for obtaining the excess probability
## in this simulation, the factor 'psi' follows N(mu_star, sigma)

# generate N random samples
dp_ratio <- rep(0,n_iteration)
psi_g <-rmvnorm(mean=mu_star, sigma=sigma, n=n_iteration)  #sampling from new psi distribution with mu_star
psi_ratio <- dmvnorm(psi_g, mean=mu,sigma)/dmvnorm(psi_g,mean=mu_star,sigma)

for (j in 1:n_iteration){
  
  #assignment of default proability given psi_g 
  cond_dp <- pnorm(as.numeric(psi_g[j,]%*%b)+c)
  
  #assignment of twisted default probability given psi_g
  expected_loss <- cond_dp %*% exposure
  if (expected_loss > threshold) {
    t_1 <- 0
  } else {
    expected_loss_t <-function(t){
      	q <-cond_dp/(cond_dp+ exp(-t*exposure)*(1-cond_dp))
    	sum(exposure*q) - threshold 
    }
    
    # find t* given psi_g
    t<-uniroot(expected_loss_t, interval=c(lower,upper),tol = 1e-10)
    t_1 <-t$root
  }
  
  #compute the exponentially twisted default probability
  twisted_dp <-exp(t_1*exposure)*cond_dp/(exp(t_1*exposure)*cond_dp+(1-cond_dp))
  
  # generate the default events of the obligors
  y <- rbern(n_obligor,twisted_dp)
  loss[j] <- sum(exposure*y)
  
  # compute the likelihood ratio of given prob. measure to IS prob. measure
  likelihood_ratio[j] <- psi_ratio[j]*prod(cond_dp*exp(t_1*exposure)+(1-cond_dp))*exp(-loss[j]*t_1)
}

I <- as.numeric(loss> threshold) 

#estimate Pr{L>L0}
excess_prob_et = mean(likelihood_ratio*I)  
SE_EP_et <- sd(likelihood_ratio*I)/sqrt(n_iteration)  
ES_et <- (1/excess_prob_et) * mean(likelihood_ratio*I*loss)
SE_ES_et <- (1/excess_prob_et) * sd(likelihood_ratio*I*loss)/sqrt(n_iteration)

toc_et = Sys.time()
elapsed_time_et = as.numeric(toc_et - tic_et,units = "secs")


###########################
## take-all scheme       ##
###########################

tic_ta = Sys.time()

upper = 0.5
lower = 0.0
k <- num_take_all_obligor
y <- rep(0,n_obligor)

loss_s <- rep(0, m_iteration)                     # realized loss at each iteration
p_a <- rep(0, m_iteration)                     # realized loss at each iteration
q_s <- rep(0, m_iteration)                     # realized loss at each iteration
dp_ratio <- rep(0, m_iteration)         # likelihood ratio of given prob. measure to IS prob. measure


psi <- rmvnorm(mean=mu, sigma, n = m_iteration)

for(i in 1: m_iteration){
  
  	cond_dp <- pnorm(as.numeric(psi[i,]%*%b)+c)
  	expected_loss <- sum(cond_dp*exposure)
  	if (expected_loss > threshold) {
    		t_1 <- 0
  	} else {
    		expected_loss_t <-function(t){
      		q <-cond_dp/(cond_dp + exp(-t*exposure)*(1-cond_dp))
      		sum(exposure*q) - threshold
    		}

    		t_1 <-uniroot(expected_loss_t, interval=c(lower,upper),tol = 1e-10)$root  
  	}
  
  	twisted_dp <- exp(t_1*exposure)*cond_dp/(exp(t_1*exposure)*cond_dp+(1-cond_dp))
	y[1:k] <- rbern(k,twisted_dp[1:k])		
  	while (prod(y[1:k]) != 0){  
  		y[1:k] <- rbern(k,twisted_dp[1:k])  
	}
	y[(k+1):n_obligor] <- rbern(n_obligor - k,twisted_dp[(k+1):n_obligor])

	p_a[i] <- prod(cond_dp[1:k])
	q_s[i] <- 1 - prod(twisted_dp[1:k])
	loss_s[i] <- sum(y * exposure)
	dp_ratio[i] <- prod(cond_dp +  exp(-exposure*t_1)*(1-cond_dp)) * exp(sum(exposure*(1-y))*t_1)
}



I <- as.numeric(loss_s > threshold)
weight <- p_a + q_s*dp_ratio*I
mu_star <- c(sum(weight * psi[,1]), sum(weight * psi[,2]), sum(weight * psi[,3]))/sum(weight)

## doing the combined IS for obtaining the excess probability
## in this simulation, the factor 'psi' follows N(mu_star, sigma)

# generate N random samples
p_a <- rep(0,n_iteration)
q_s <- rep(0,n_iteration)
loss_a <- rep(0,n_iteration)
loss_s <- rep(0,n_iteration)
t_1 <- rep(0,n_iteration)
dp_ratio <- rep(0, n_iteration)        

psi <- rmvnorm(mean=mu_star, sigma, n = n_iteration)
for(i in 1:n_iteration){ 
  	cond_dp <- pnorm(as.numeric(psi[i,]%*%b)+c)
	
  	expected_loss <- sum(cond_dp*exposure)
  	if (expected_loss > threshold) {
    		t_1 <- 0
  	} else {
    		expected_loss_t <-function(t){
      		q <-cond_dp/(cond_dp + exp(-t*exposure)*(1-cond_dp))
      		sum(exposure*q) - threshold
    		}

		t_1 <-uniroot(expected_loss_t, interval=c(lower,upper),tol = 1e-10)$root  
  	}
 
  	twisted_dp <- cond_dp/(cond_dp + exp(-t_1*exposure)*(1-cond_dp))
	y[1:k] <- rbern(k,twisted_dp[1:k])		
  	while (prod(y[1:k]) != 0){  
  		y[1:k] <- rbern(k,twisted_dp[1:k])  
	}
	y[(k+1):n_obligor] <- rbern(n_obligor - k,twisted_dp[(k+1):n_obligor])

	p_a[i] <- prod(cond_dp[1:k])
	q_s[i] <- 1 - prod(twisted_dp[1:k])
  	loss_a[i] <- sum(exposure[1:k]) + sum(cond_dp[(k+1):n_obligor]*exposure[(k+1):n_obligor])
	loss_s[i] <- sum(y * exposure)
	dp_ratio[i] <- prod(cond_dp +  exp(-exposure*t_1)*(1-cond_dp)) * exp(sum(exposure*(1-y))*t_1)
}

psi_ratio <- dmvnorm(psi, mean=mu,sigma)/dmvnorm(psi,mean=mu_star,sigma)
w <- psi_ratio*dp_ratio
I <- as.numeric(loss_s > threshold) 

#estimate Pr{L>L0}
prob_ta <- psi_ratio*p_a + q_s*w*I 
es_ta <- psi_ratio*p_a*loss_a+q_s*w*I*loss_s
excess_prob_ta = mean(prob_ta)
SE_EP_ta <- sd(prob_ta)/sqrt(n_iteration)  
ES_ta <- (1/excess_prob_ta) * mean(es_ta)
SE_ES_ta <- (1/excess_prob_ta) * sd(es_ta)/sqrt(n_iteration)
toc_ta = Sys.time()
elapsed_time_ta = as.numeric(toc_ta - tic_ta,units = "secs")

vt_EP_cmc <- elapsed_time_cmc*SE_EP_cmc^2 
vt_EP_et <- elapsed_time_et*SE_EP_et^2 
vt_EP_ta <- elapsed_time_ta*SE_EP_ta^2

CI_EP_cmc<-c(excess_prob_cmc-2*SE_EP_cmc,excess_prob_cmc+2*SE_EP_cmc)

vt_ES_cmc <- elapsed_time_cmc*SE_ES_cmc^2 
vt_ES_et <- elapsed_time_et*SE_ES_et^2 
vt_ES_ta <- elapsed_time_ta*SE_ES_ta^2 

CI_ES_cmc<-c(ES_cmc-2*SE_ES_cmc,ES_cmc+2*SE_ES_cmc)

{
  cat("### result ###\n")
  cat("number of obligors:", n_obligor, "\n")
  cat("threshold:", threshold, "\n")
  cat("elapsed time:", elapsed_time_cmc, elapsed_time_et, elapsed_time_ta,  "\n")
  cat("number of take-all:", k, "\n")
  #cat("probability of take_all event in P measure:",mean(cond_p_all), "\n")
  #cat("probability of take_all event in Q measure:",mean(cond_q_all), "\n")
  cat("excess loss probs:", excess_prob_cmc, excess_prob_et, excess_prob_ta,  "\n")
  cat("s.e.:", SE_EP_cmc, SE_EP_et, SE_EP_ta,  "\n")
#  cat("CI_EP_cmc:",CI_EP_cmc, "\n")
#  cat("var*time:", vt_EP_cmc, vt_EP_et, vt_EP_ta, "\n")
#  cat("expected shortfall:", ES_cmc, ES_et, ES_ta, "\n")
#  cat("s.e.:", SE_ES_cmc, SE_ES_et, SE_ES_ta, "\n")
#  cat("CI_ES_cmc:",CI_ES_cmc, "\n")
#  cat("var*time:", vt_ES_cmc, vt_ES_et, vt_ES_ta, "\n")
}
