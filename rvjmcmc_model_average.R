library(tidyverse)

# Data Generation


beta_0 = 1/4
beta_1 = 0.05
beta_2 = -1

x = seq(0,5, 0.05)

sigma = 1
tau = 2

u = runif(length(x),0,1)

len = length(x)

y = (u<0.5)*(beta_0 + beta_1 * x + rnorm(len,0,sigma))+(u>=0.5)*(beta_0 + beta_1 * x + beta_2 * x^2 + rnorm(len,0,sigma))

plot(x, y, col = "blue", pch = 1, cex = 2, xlab = "x", ylab = "y")


## Reversible Jump Model Averaging

x_full = cbind(1, x, x^2)
#x_sat = cbind(1, x)

xtx = t(x_full) %*% x_full

D = svd(xtx)$d

P = svd(xtx)$u

X_orth = x_full %*% P

true_alpha = cbind(beta_0,beta_1,beta_2) %*% P

n = length(y)


miu_alpha = (y %*% X_orth) * tau^2 / (sigma^2+tau^2 * D)

tau_alpha = tau^2 * sigma^2 / (sigma^2 + tau^2 * D)


f11 = function(alpha_old)
{
  alpha_new = c(rnorm(1,miu_alpha[1],tau_alpha[1]),
                rnorm(1,miu_alpha[2],tau_alpha[2]),
                alpha_old[3])
  
  ratio = dnorm(alpha_new[1],miu_alpha[1],tau_alpha[1],log=TRUE) + 
    dnorm(alpha_new[2],miu_alpha[2],tau_alpha[2],log=TRUE)-
    dnorm(alpha_old[1],miu_alpha[1],tau_alpha[1],log = TRUE)-
    dnorm(alpha_old[2],miu_alpha[2],tau_alpha[2],log = TRUE)
  
  ratio = min(1,exp(ratio))
  
  return( list(ratio=ratio,alpha=alpha_new) ) 
}

f22 = function(alpha_old)
{
  alpha_new = c(rnorm(1,miu_alpha[1],tau_alpha[1]),
                rnorm(1,miu_alpha[2],tau_alpha[2]),
                rnorm(1,miu_alpha[3],tau_alpha[3]))
  
  ratio = dnorm(alpha_new[1],miu_alpha[1],tau_alpha[1],log=TRUE) + 
    dnorm(alpha_new[2],miu_alpha[2],tau_alpha[2],log=TRUE) + 
    dnorm(alpha_new[3],miu_alpha[3],tau_alpha[3],log=TRUE)-
    dnorm(alpha_old[1],miu_alpha[1],tau_alpha[1],log = TRUE)-
    dnorm(alpha_old[2],miu_alpha[2],tau_alpha[2],log = TRUE)-
    dnorm(alpha_old[3],miu_alpha[3],tau_alpha[3],log = TRUE)
  
  ratio = min(exp(ratio),1)
  return(list(ratio=ratio,alpha=alpha_new))
}

f12 = function(alpha_old)
{
  alpha_new = c(alpha_old[1],
                alpha_old[2],
                rnorm(1,miu_alpha[3],tau_alpha[3]))
  
  ratio = dnorm(alpha_new[3],miu_alpha[3],tau_alpha[3],log=TRUE)-
    dnorm(alpha_old[3],miu_alpha[3],tau_alpha[3],log = TRUE)
  
  ratio = min(exp(ratio),1)
  return(list(ratio=ratio,alpha=alpha_new))
}

f21 = function(alpha_old)
{
  alpha_new = c(rnorm(1,miu_alpha[1],tau_alpha[1]),
                rnorm(1,miu_alpha[2],tau_alpha[2]),
                alpha_old[3])
  
  ratio = dnorm(alpha_new[1],miu_alpha[1],tau_alpha[1],log=TRUE) + 
    dnorm(alpha_new[2],miu_alpha[2],tau_alpha[2],log=TRUE) -  
    dnorm(alpha_old[1],miu_alpha[1],tau_alpha[1],log = TRUE)-
    dnorm(alpha_old[2],miu_alpha[2],tau_alpha[2],log = TRUE)-
    dnorm(alpha_old[3],miu_alpha[3],tau_alpha[3],log = TRUE)
  
  ratio = min(exp(ratio),1)
  
  return(list(ratio=ratio,alpha=alpha_new))
}

#Burn

accept = 0
alpha_old = c(0,1,2)
I_11 = 0
I_22 = 0
I_12 = 0
I_21 = 0

for(i in 1:2000)
{
  #print(i/1000)
  u = runif(1)
  if(u<0.25){
    result = f11(alpha_old)
    v = runif(1)
    if(v<result$ratio){
      alpha_old = result$alpha
      accept = accept + 1
      I_11 = I_11 + 1
    }
  }
  else if(u<0.5){
    result = f22(alpha_old)
    v = runif(1)
    if(v<result$ratio){
      alpha_old = result$alpha
      accept = accept + 1
      I_22 = I_22 + 1
    }
  }
  else if(u<0.75){
    result = f12(alpha_old)
    v = runif(1)
    if(v<result$ratio){
      alpha_old = result$alpha
      accept = accept + 1
      I_12 = I_12 + 1
    }
  }
  else{
    result = f21(alpha_old)
    v = runif(1)
    if(v<result$ratio){
      alpha_old = result$alpha
      accept = accept + 1
      I_21 = I_21 + 1
    }
  }
}

print(accept/2000)
print(c(I_11,I_22,I_12,I_21)/accept)
alpha_old

#MCMC_stage2

accept = 0
I_11 = 0
I_22 = 0
I_12 = 0
I_21 = 0

y_frame = matrix(0,nrow=1,ncol=len) %>% as.data.frame()

for(i in 1:400)
{
  #print(i/1000)
  u = runif(1)
  if(u<0.25){
    result = f11(alpha_old)
    v = runif(1)
    if(v<result$ratio){
      alpha_old = result$alpha
      accept = accept + 1
      I_11 = I_11 + 1
    }
    y_temp = X_orth[,1:2] %*% alpha_old[1:2]
    y_temp = y_temp %>% t()
    y_frame = rbind(y_frame,y_temp)
  }
  else if(u<0.5){
    result = f22(alpha_old)
    v = runif(1)
    if(v<result$ratio){
      alpha_old = result$alpha
      accept = accept + 1
      I_22 = I_22 + 1
    }
    y_temp = X_orth %*% alpha_old
    y_temp = y_temp %>% t()
    y_frame = rbind(y_frame,y_temp)
  }
  else if(u<0.75){
    result = f12(alpha_old)
    v = runif(1)
    if(v<result$ratio){
      alpha_old = result$alpha
      accept = accept + 1
      I_12 = I_12 + 1
      y_temp = X_orth %*% alpha_old
      y_temp = y_temp %>% t()
      y_frame = rbind(y_frame,y_temp)
    }
    else{
      y_temp = X_orth[,1:2] %*% alpha_old[1:2]
      y_temp = y_temp %>% t()
      y_frame = rbind(y_frame,y_temp)
    }
  }
  else{
    result = f21(alpha_old)
    v = runif(1)
    if(v<result$ratio){
      alpha_old = result$alpha
      accept = accept + 1
      I_21 = I_21 + 1
      y_temp = X_orth[,1:2] %*% alpha_old[1:2]
      y_temp = y_temp %>% t()
      y_frame = rbind(y_frame,y_temp)
    }
    else{
      y_temp = X_orth %*% alpha_old
      y_temp = y_temp %>% t()
      y_frame = rbind(y_frame,y_temp)
    }
  }
}

print(accept/400)
print(c(I_11,I_22,I_12,I_21)/accept)
alpha_old

y_average = colMeans(y_frame[2:401,])

temp_frame = data.frame(x=x,true_y=y,estimated_y=y_average)

ggplot(temp_frame,aes(x=x,y=estimated_y)) + 
  geom_point(aes(x=x,y=true_y),color="red") + 
  geom_line(color="blue")+
  theme_minimal()
