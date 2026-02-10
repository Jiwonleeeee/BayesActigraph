library(tidyverse)
library(mvtnorm)

#' Initialize states based on data quantiles and duration constraints
#' @param data list containing the observed data vector 'y'.
#' @param Lmin integer, minimum duration length for a state to be valid.
#' @param q numeric, quantile probability used to determine the cutoff threshold.
#' @return A list containing the dataframe with state assignments and a vector of change points 's'.
initial_states <- function(data, Lmin, q) {
  
  y  <- data$y
  Tn <- length(y)
  margin <- 0
  
  cutoff <- quantile(y, q)
  
  cutoff_up   <- cutoff + margin
  cutoff_down <- cutoff - margin
  
  # Start with rest
  cur <- 1
  state <- integer(Tn); state[1] <- cur
  
  # Streak counters
  up_streak <- 0   # for rest -> active switch
  dn_streak <- 0   # for active -> rest switch
  
  for (t in 2:Tn) {
    yt <- y[t]
    
    if (cur == 2) {
      # Current rest: switch to active if > cutoff_up for Lmin consecutive points
      if (yt > cutoff_up) up_streak <- up_streak + 1 else up_streak <- 0
      if (up_streak >= Lmin) {
        # Prevent switch if remaining segment is too short
        if ((Tn - t + 1) >= Lmin) {
          cur <- 1
          up_streak <- 0; dn_streak <- 0
        } else {
          up_streak <- 0  # Cancel switch
        }
      }
    } else {
      # Current active: switch to rest if <= cutoff_down for Lmin consecutive points
      if (yt <= cutoff_down) dn_streak <- dn_streak + 1 else dn_streak <- 0
      if (dn_streak >= Lmin) {
        # Prevent switch if remaining segment is too short
        if ((Tn - t + 1) >= Lmin) {
          cur <- 2
          up_streak <- 0; dn_streak <- 0
        } else {
          dn_streak <- 0  # Cancel switch
        }
      }
    }
    state[t] <- cur
  }
  
  # Recalculate interval
  df <- data |>
    dplyr::mutate(
      state      = state,
      tmp_diff   = state - dplyr::lag(state),
      group_start= ifelse(tmp_diff != 0 | is.na(tmp_diff), 1L, 0L),
      interval   = cumsum(group_start)
    ) |>
    dplyr::select(y, t, state, interval)
  
  s <- c(0, which(diff(df$state) != 0), Tn)
  
  list(df = df, s = s)
}

#' Calculate Dirichlet-Multinomial log-likelihood
#' @param g vector of counts or gaps.
#' @param alpha numeric, concentration parameter.
#' @return Numeric value of the log-likelihood.
loglik_DM <- function(g, alpha) {
  K <- length(g); R <- sum(g)
  lgamma(R + 1) - sum(lgamma(g + 1)) +
    lgamma(K * alpha) - lgamma(R + K * alpha) +
    sum(lgamma(g + alpha)) - K * lgamma(alpha)
}

#' Convert change points to gap lengths
#' @param s vector of change points.
#' @param m integer, subtraction adjustment factor.
#' @return Vector 'g' representing gap lengths.
s_to_g <- function(s, m){
  
  nn <- length(s)
  g <- numeric(nn - 1)
  
  for(i in 1:(nn - 1)){
    g[i] <- s[i+1]-s[i]-m
  }
  return(g)
}

#' Calculate sum of log MVN likelihoods for active intervals
#' @param df dataframe containing 'y', 'state', and 'interval' columns.
#' @param mu numeric, mean of the MVN distribution.
#' @param sigma2 numeric, scale parameter of the MVN distribution.
#' @param l numeric, range parameter of the MVN distribution.
#' @param nugget numeric, nugget parameter of the MVN distribution.
#' @return Numeric sum of log-likelihoods for all active intervals.
sum_log_mvn_lkh <- function(df, mu, sigma2, l, nugget){
  
  active_interval <- df %>%
    filter(state==1) %>%
    distinct(interval) %>%
    pull(interval)
  
  log_sum <- 0
  for(i in seq_along(active_interval)){
    
    y_chunk <- df$y[df$interval == active_interval[i]]
    t_chunk <- 1:length(y_chunk)
    
    Cor <- exp(-abs(outer(t_chunk, t_chunk,FUN = "-"))/l) + nugget* diag(length(y_chunk))
    
    log_sum <- log_sum + dmvnorm(y_chunk, mean = rep(mu, length(y_chunk)), sigma = sigma2 * Cor, log = TRUE)
  }
  return(log_sum)
}

#' Calculate likelihood for a selected state chunk
#' @param y_chunk vector of observations for the specific interval.
#' @param state integer, current state (1 = active, other = rest).
#' @param mvn_parameters vector containing named parameters (mu, sigma2, l, nugget).
#' @param dp_parameters list containing Poisson mixture parameters (Pi, Lambda).
#' @return Numeric log-likelihood value.
selected_y_lkh <- function(y_chunk, state, mvn_parameters, dp_parameters){
  
  # mvn_parameters: mu, sigma2, l, nugget
  mu <- mvn_parameters["mu"]
  sigma2 <- mvn_parameters["sigma2"]
  l <- mvn_parameters["l"]
  nugget <- mvn_parameters["nugget"]
  
  # dp_parameters: Pi, Lambda
  Pi <- dp_parameters[["Pi"]]
  Lambda <- dp_parameters[["Lambda"]]
  
  if(state == 1){
    
    # active
    d <- length(y_chunk)
    t <- 1:d
    Cor <- exp(-abs(outer(t, t,FUN = "-"))/l) + nugget * diag(d)
    log_lkh_value <- dmvnorm(y_chunk, mean = rep(mu, d), sigma = sigma2 * Cor ,log=TRUE )
  } else {
    y_chunk <- round(y_chunk)
    # resting
    if(length(y_chunk)==1){ # then it's scalar
      log_lkh_value <- log(sum(Pi * dpois(y_chunk, Lambda)))
    } else {
      log_lkh_value <- 0
      for(t in 1:length(y_chunk)){
        log_lkh_value <-log_lkh_value + log(sum(Pi * dpois(y_chunk[t], Lambda)))
      }
    }
  }
  
  return(log_lkh_value)
}