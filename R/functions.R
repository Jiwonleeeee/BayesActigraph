## initial_position: setting initial locations to the labels of active and resting states. 
#' @param y_input T-dimensional y vector.
#' @param cutoff user-specific cutoff (threshold) value. Measurements less than the cutoff are classified to the resting states. 
#' @param end_input the last element of the vector 'rest'. We fix it as T + 0.5 throughout the entire algorithm but it varies by subjects.  
#' @return active B-dimensional vector 
#' @return rest B-dimensional vector. Each element in 'active' and 'rest' constructs the active interval by (active[b], rest[b]) and the resting interval by (rest[b], active[b+1]) for b=1,...,B.
initial_position <- function(y_input, cutoff, end_input){
  
  Signs <- sign(y_input - cutoff)
  Num  <- length(y_input)
  changePoints <- Indi <- numeric(Num)
  
  for(i in 2:Num){
    if(Signs[i]!=Signs[i-1]){
      changePoints[i] <- i
      
      Indi[i] <- ifelse(Signs[i] > Signs[i-1], "A" , "R")
      
    }
  }
  
  CPs <- changePoints[changePoints!=0]  
  Indi <- Indi[Indi!="0"]
  LthCP <- length(CPs)
  
  if(LthCP%%2 == 0){
    
    if(Indi[1]=="A"){
      active <- CPs[seq(1, LthCP-1, by =2)]
      rest <- CPs[seq(2, LthCP, by=2)]
      
      active <- active[-1]
      if(Indi[length(Indi)]=="R"){
        rest <- rest[-length(rest)]
      }
    }else{
      rest <- CPs[seq(1, LthCP-1, by =2)]
      active <- CPs[seq(2, LthCP, by=2)]
      if(Indi[length(Indi)]=="R"){
        rest <- rest[-length(rest)]
      }
    }
    
  }else{
    
    if(Indi[1]=="A"){
      active <- CPs[seq(1, LthCP, by =2)]
      rest <- CPs[seq(2, LthCP-1, by=2)]
      
      active <- active[-1]
      if(Indi[length(Indi)]=="R"){
        rest <- rest[-length(rest)]
      }
      
    }else{
      rest <- CPs[seq(1, LthCP, by =2)]
      active <- CPs[seq(2, LthCP-1, by=2)]
      if(Indi[length(Indi)]=="R"){
        rest <- rest[-length(rest)]
      }
    }
    
  }
  
  active <- c(0.5, active-0.5)
  rest <- c(rest-0.5, end_input)
  
  return(list(active=active, rest=rest))
  
}

## log_pois_lkh: return the sum of log likelihood of the mixture of Poisson distributions.
#' @param y_input scalar or vector of y of interest. 
#' @param Pi_input K-dimensional vector of mixing proportions pi. 
#' @param Lambda_input K-dimensional vector of parameters of Poisson distributions. 
#' @return sum_value calculated sum of log likelihood.
log_pois_lkh <- function(y_input, Pi_input, Lambda_input){
  
  
  if(length(y_input)==1){ # then y_input is a scalar
    sum_value <- log(sum(Pi_input * dpois(y_input, Lambda_input)))
  }else{
    
    sum_value <- 0
    for(t in 1:length(y_input)){
      sum_value <- sum_value + log(sum(Pi_input * dpois(y_input[t], Lambda_input)))
    }
    
  }
  return(sum_value)
}

## log_mvn_lkh: return the log likelihood of multivariate normal distribution (MVN).
#' @param y_input scalar or vector of y of interest. 
#' @param t_input length of y_input.
#' @param mu_input mean parameter of MVN.
#' @param sigma2_input scale parameter of MVN.
#' @param l_input range parameter of MVN.
#' @param g_input nugget term of MVN.
#' @return value calculated log likelihood
log_mvn_lkh <- function(y_input, t_input, mu_input, sigma2_input, l_input, g_input){
  d <- length(y_input)
  k_t_t_prime <- exp(-abs(outer(t_input,t_input,FUN = "-"))/l_input)+g_input*diag(length(t_input))
  value <- dmvnorm(y_input, mean = rep(mu_input, d), sigma = sigma2_input * k_t_t_prime ,log=T )
  return(value)
}



## sum_log_mvn_lkh: return the sum of log likelihood of MVN.
#' @param y_input T-dimensional y vector.
#' @param active_input B-dimensional vector
#' @param rest_input B-dimensional vector, constructing intervals like (active_input, rest_input) to filter the elements of 'y_input' in the active states.   
#' @return sum_value calculated sum of log likelihood
sum_log_mvn_lkh <- function(y_input, t_input, mu_input, sigma2_input, l_input, g_input, active_input, rest_input){
  
  sum_value <- 0
  length_a <- length(active_input)
  for(I in 1:length_a){
    
    t_set <- t_input[ t_input > active_input[I] & t_input < rest_input[I]]
    y_set <- y_input[ t_input > active_input[I] & t_input < rest_input[I]]
    
    if(length(y_set)!=0){
      sum_value <- sum_value + log_mvn_lkh(y_set, t_set, mu_input, sigma2_input, l_input, g_input)
    }
    
  }
  
  return(sum_value)
  
}

## prob_return_birth: return probabilities for choosing an active (resting) interval in which a new resting (active) interval will be proposed. 
#' @param n_active_input total number of active intervals. It is the same as the length of 'active_input' but the separate parameter 'n_active_input' is defined to accommodate more general cases.
#' @param n_rest_input total number of resting intervals. It is the same as the length of 'rest_input' - 1.
#' @return active_prob 'n_active_input'-dimensional vector with each element representing probability proportional to the length of the corresponding active interval  
#' @return rest_prob 'n_rest_input'-dimensional vector with each element representing probability proportional to the length of the corresponding resting interval  
prob_return_birth <- function(active_input, rest_input, n_active_input, n_rest_input){
  
  
  active_length <- numeric(n_active_input)
  rest_length <- numeric(n_rest_input)
  
  for(I in 1:n_active_input){
    active_length[I] <- rest_input[I]-active_input[I]
  }
  
  for(I in 1:n_rest_input){
    rest_length[I] <- active_input[I+1]-rest_input[I]
  }
  
  active_prob <- active_length/ sum(active_length)
  rest_prob <- rest_length/ sum(rest_length)
  
  return(list(active_prob=active_prob, rest_prob=rest_prob))
  
  
}


## prob_return_death: return probabilities for choosing an active (resting) interval to be removed. 
#' @return active_prob ('n_active_input'-2)-dimensional vector with each element representing probability inversely proportional to the length of the corresponding active interval. Since we fix the state of the first and last interval as active, they cannot be removed. 
#' @return rest_prob 'n_rest_input'-dimensional vector with each element representing probability inversely proportional to the length of the corresponding resting interval
prob_return_death<- function(active_input, rest_input, n_active_input, n_rest_input){
  
  
  active_length <- numeric(n_active_input-2)
  rest_length <- numeric(n_rest_input)
  
  for(I in 2:(n_active_input-1)){
    active_length[I-1] <- rest_input[I]-active_input[I]
  }
  
  for(I in 1:n_rest_input){
    rest_length[I] <- active_input[I+1]-rest_input[I]
  }
  
  
  active_length <- 1/active_length
  rest_length <- 1/rest_length
  
  active_prob <- active_length/ sum(active_length)
  rest_prob <- rest_length/ sum(rest_length)
  
  return(list(active_prob=active_prob, rest_prob=rest_prob))
  
}

## acc_birth_ar: return the updated parameters (if birth of (a,r) is accepted), return the same parameters otherwise. 
#' @param y_input T-dimensional y vector.
#' @param t_input T-dimensional time index vector (1,...,T).
#' @param mvn_inputs 4-dimensional vector of MVN parameters (mu, sigma2, l, g).
#' @param p_birth_input probability of choosing birth move.
#' @param p_death_input probability of choosing death move.
#' @param active_input B-dimensional vector
#' @param rest_input B-dimensional vector
#' @return active (B+1)-dimensional vector if the birth of (a,r) is accepted, the same B-dimensional vector otherwise
#' @return rest (B+1)-dimensional vector if the birth of (a,r) is accepted, the same B-dimensional vector otherwise
#' @return accept indicator if accept = 1 then the proposed is accepted, and accept = 0 otherwise. (for tracking the acceptance probability)
#' @return n_active number of active intervals after updating
#' @return n_rest number of resting intervals after updating

acc_birth_ar <- function(y_input, t_input, mvn_inputs, Pi_input, Lambda_input, p_birth_input, p_death_input, active_input, rest_input, n_active_input, n_rest_input){
  
  Time_T <- length(t_input)
  L <- t_input[length(t_input)]+0.5
  accept <- 0
  a_len <- length(active_input); r_len <- length(rest_input)
  
 
  ar_prob <- prob_return_birth(active_input, rest_input, n_active_input, n_rest_input)$rest_prob
  j <- sample(1:n_rest_input, 1, prob=ar_prob) 
  start <- rest_input[j]
  end <- active_input[j+1]
  
  
  if((end-start)>2){
    
    temp <- sample( seq(start+1, end-1, by=1), 2 )
    a_star <- min(temp); r_star <- max(temp) 
    
    t_star <- t_input[t_input>a_star & t_input<r_star]
    y_star <- y_input[t_input>a_star & t_input<r_star]
    
    if(length(y_star)!=0){
      k <- length(c(active_input,rest_input))
      

      mu <- mvn_inputs[1]
      sigma2 <- mvn_inputs[2]
      l <- mvn_inputs[3]
      g <- mvn_inputs[4]
      
      den <- (ar_prob[j]/sum(ar_prob)) * 2 * (1/(end-start-1))^2  * p_birth_input
      num <- (r_star - a_star)/(sum(prob_return_death(active_input, rest_input, n_active_input, n_rest_input)$active_prob) + (r_star - a_star)) * p_death_input
      
      log_lkh_ratio <- log_mvn_lkh(y_star, t_star, mu, sigma2, l, g) - log_pois_lkh(y_star, Pi_input, Lambda_input)
      log_step_ratio <- log(k) + log(k-1) - 2 * log(Time_T-1)
      log_proposal_ratio <- log(num)-log(den)
      
      log_acc <- log_lkh_ratio + log_step_ratio + log_proposal_ratio
      accept <- ifelse(log(runif(1))<log_acc,1,0) 
      
    }
    
    if(accept==1){
      if(j==1){
        
        rest <- c(rest_input[1:j],r_star,rest_input[(j+1):r_len])
        active <- c(active_input[1:j],a_star,active_input[(j+1):a_len])
        
        
      }else if(j==n_rest_input){
        
        rest <- c(rest_input[1:j],r_star,rest_input[(j+1):r_len])
        
        active <- c(active_input[1:j],a_star,active_input[(j+1):a_len])
        
      }else{
        
        rest <- c(rest_input[1:j],r_star,rest_input[(j+1):r_len])
        active <- c(active_input[1:j],a_star,active_input[(j+1):a_len])
        
      }
      a_len <- length(active); r_len <- length(rest)
      n_active <- a_len
      n_rest <- r_len-1
      
      
    }
    
  }
  return(list(active=active, rest=rest, n_active=n_active, n_rest=n_rest, accept=accept))
}

## acc_death_ar: return the updated parameters (if death of (a,r) is accepted), return the same parameters otherwise. 

acc_death_ar <- function(y_input, t_input, mvn_inputs, Pi_input, Lambda_input, p_birth_input, p_death_input, active_input, rest_input, n_active_input, n_rest_input){
  
  
  
  accept <- 0
  a_len <- length(active_input); r_len <- length(rest_input)
  L <- t_input[length(t_input)]+0.5
  Time_T <- length(t_input)
  ar_prob <- prob_return_death(active_input, rest_input, n_active_input, n_rest_input)$active_prob
  
  if(n_active_input==3){
    j <- 2
  }else{
    j <- sample(2:(n_active_input-1), 1, prob=ar_prob)
  }
  
  start <- active_input[j]
  end <- rest_input[j]
 
  
  a_star <- start; r_star <- end
  
  t_star <- t_input[t_input>a_star & t_input<r_star]
  y_star <- y_input[t_input>a_star & t_input<r_star]
  
  if(length(y_star)!=0){
    
    k <- length(c(active_input,rest_input))
    
    mu <- mvn_inputs[1]
    sigma2 <- mvn_inputs[2]
    l <- mvn_inputs[3]
    g <- mvn_inputs[4]
    
    
    new_sum <- (sum(prob_return_birth(active_input, rest_input, n_active_input, n_rest_input)$rest_prob)
                - (active_input[j] - rest_input[j-1])
                - (active_input[j+1] - rest_input[j]) 
                + (active_input[j+1]-rest_input[j-1]))
    
    den <- ar_prob[j-1] / sum(ar_prob) * p_death_input
    num <- (active_input[j+1] - rest_input[j-1])/new_sum * 2 * (1/(active_input[j+1] - rest_input[j-1]-1))^2 * p_birth_input
    
    
    log_lkh_ratio <-  log_pois_lkh(y_star, Pi_input, Lambda_input) - log_mvn_lkh(y_star, t_star, mu, sigma2, l, g) 
    log_step_ratio <- 2 * log(Time_T-1) - log(k-2) - log(k-3)
    log_proposal_ratio <- log(num)- log(den)
    
    log_acc <- log_lkh_ratio + log_step_ratio + log_proposal_ratio
    accept <- ifelse(log(runif(1))<log_acc,1,0) 
    
    if(accept==1){
      active <- active_input[-j]
      rest <- rest_input[-j]
      
      
      if(j==1){
        rest <- c(0.5,rest_input)
      }else if(j==n_active_input){
        active <- c(active_input,L)
      }
      
      a_len <- length(active); r_len <- length(rest)
      n_active <- a_len
      n_rest <- r_len-1
      
    }
    
  }
  
  
  return(list(active=active, rest=rest, n_active=n_active, n_rest=n_rest, accept=accept))
}

## acc_birth_ra: return the updated parameters (if birth of (r,a) is accepted), return the same parameters otherwise. 

acc_birth_ra <- function(y_input, t_input, mvn_inputs, Pi_input, Lambda_input, p_birth_input, p_death_input, active_input, rest_input, n_active_input, n_rest_input){
  
  accept <- 0
  L <- t_input[length(t_input)]+0.5
  Time_T <- length(t_input)
  a_len <- length(active_input); r_len <- length(rest_input)
  
  ra_prob <- prob_return_birth(active_input, rest_input, n_active_input, n_rest_input)$active_prob
  j <- sample(1:n_active_input, 1, prob=ra_prob) 
  start <- active_input[j]
  end <- rest_input[j]
  if((end-start)>2){
    
    temp <- sample( seq(start+1, end-1, by=1), 2 )
    r_star <- min(temp); a_star <- max(temp) 
    
    t_star <- t_input[t_input>r_star & t_input<a_star]
    y_star <- y_input[t_input>r_star & t_input<a_star]
    
    if(length(y_star)!=0){
      k <- length(c(active_input,rest_input))
      
      mu <- mvn_inputs[1]
      sigma2 <- mvn_inputs[2]
      l <- mvn_inputs[3]
      g <- mvn_inputs[4]
      
      
      den <- (ra_prob[j]/sum(ra_prob)) * 2 * (1/(end-start-1))^2 * p_birth_input
      num <- (a_star - r_star) / (sum(prob_return_death(active_input, rest_input, n_active_input, n_rest_input)$rest_prob) + (a_star-r_star)) * p_death_input
      
      
      log_lkh_ratio <- log_pois_lkh(y_star, Pi_input, Lambda_input) - log_mvn_lkh(y_star, t_star, mu, sigma2, l, g) 
      log_step_ratio <- log(k) + log(k-1) - 2 * log(Time_T-1)
      log_proposal_ratio <- log(num) - log(den)
      
      log_acc <- log_lkh_ratio + log_step_ratio + log_proposal_ratio
      accept <- ifelse(log(runif(1))<log_acc,1,0)
      
    }
    
    if(accept==1){
      
      
      if(j==1){
        
        rest <- c(r_star,rest_input)
        
        active <- c(active_input[1:j],a_star,active_input[(j+1):a_len])
        
      }else if(j==n_active_input){
        rest <- c(rest_input[1:(j-1)],r_star,rest_input[j:r_len])
        
        active <- c(active_input,a_star)
        
      }else{
        rest <- c(rest_input[1:(j-1)],r_star,rest_input[j:r_len])
        
        active <- c(active_input[1:j],a_star,active_input[(j+1):a_len])
      }
      
      a_len <- length(active); r_len <- length(rest)
      n_active <- a_len
      n_rest <- r_len-1
      
      
    }
    
  }
  
  return(list(active=active, rest=rest, n_active=n_active, n_rest=n_rest, accept=accept))
}

## acc_death_ra: return the updated parameters (if death of (r,a) is accepted), return the same parameters otherwise. 

acc_death_ra <- function(y_input, t_input, mvn_inputs, Pi_input, Lambda_input , p_birth_input, p_death_input, active_input, rest_input, n_active_input, n_rest_input){
  
  Time_T <- length(t_input)
  L <- t_input[length(t_input)]+0.5
  accept <- 0
  ra_prob <- prob_return_death(active_input, rest_input, n_active_input, n_rest_input)$rest_prob
  j <- sample(1:n_rest_input, 1, prob=ra_prob)
  
  start <- rest_input[j]
  end <- active_input[j+1]
  
  r_star <- start; a_star <- end
  
  t_star <- t_input[t_input>r_star & t_input<a_star]
  y_star <- y_input[t_input>r_star & t_input<a_star]
  
  if(length(y_star)!=0){
    
    k <- length(c(active_input,rest_input))
    
    mu <- mvn_inputs[1]
    sigma2 <- mvn_inputs[2]
    l <- mvn_inputs[3]
    g <- mvn_inputs[4]
    
    
    
    new_sum <- (sum(prob_return_birth(active_input, rest_input, n_active_input, n_rest_input)$active_prob)
                -(rest_input[j] - active_input[j])
                -(rest_input[j+1] - active_input[j+1])
                + rest_input[j+1]- active_input[j])
    
    den <- ra_prob[j]/sum(ra_prob) * p_birth_input
    num <- (rest_input[j+1]- active_input[j])/new_sum * 2 * (1/(rest_input[j+1]-active_input[j]))^2 * p_death_input
    
    
    log_lkh_ratio <- log_mvn_lkh(y_star, t_star, mu, sigma2, l, g) -log_pois_lkh(y_star, Pi_input, Lambda_input)
    log_step_ratio <- 2 * log(Time_T-1) - log(k-2) - log(k-3)
    log_proposal_ratio <- log(num) - log(den)
    
    log_acc <- log_lkh_ratio + log_step_ratio + log_proposal_ratio
    accept <- ifelse(log(runif(1))<log_acc,1,0)
    
    if(accept==1){
      
      rest <- rest_input[-j]
      active <- active_input[-(j+1)]
      
      a_len <- length(active); r_len <- length(rest)
      n_active <- a_len
      n_rest <- r_len-1
      
    }
    
  }
  return(list(active=active, rest=rest, n_active=n_active, n_rest=n_rest, accept=accept))
}


