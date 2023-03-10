library(tidyverse)


times = c(1, 2, 3, 4, 5)
p_exposed = c(0.2, 0.2, 0.2, 0.2)
inst_rate_exposed =   c(0.02, 0.02, 0.02, 0.02)/
inst_rate_unexposed = rep(0.01, 0.04)
n_start = 1000


calc_OR_RR = function(times, p_exposed, inst_rate_exposed, 
  inst_rate_unexposed, n_start, OR_matched = NA, OR_unmatched = NA){


  df = data.frame( 
    t_start = times[-length(times)],
    t_end = times[-1],
    
    p_exposed = p_exposed,
    
    inst_rate_exposed =   inst_rate_exposed,
    inst_rate_unexposed = inst_rate_unexposed,
    
    n = c(n_start, rep(NA, 3)),
    
    n_event_exposed = NA,
    n_event_unexposed = NA
  )
  
  df = df %>%
    mutate(
      delta_t = t_end - t_start,
      
      p_unexposed = 1 - p_exposed,
      
      n_no_event_exposed_begin = p_exposed * n,
      n_no_event_unexposed_begin = p_unexposed * n,
      
      n_no_event_exposed_end = NA,
      n_no_event_unexposed_end = NA,
    )
  
  simulation = simulate(df)
  
  df = simulation[1]
  results = simulation[2]
  
}




simulate = function(df){
  
  for (i in 1:nrow(df)){
    
    # The end quantities are related by exponential decay to the begin quantities
    df[i, 'n_no_event_exposed_end'] = df[i, 'n_no_event_exposed_begin'] * exp( - df[i, 'inst_rate_exposed'] * df[i, 'delta_t'] ) 
    df[i, 'n_no_event_unexposed_end'] = df[i, 'n_no_event_unexposed_begin'] * exp( - df[i, 'inst_rate_unexposed'] * df[i, 'delta_t']) 
    
    if (i == nrow(df)){}
    else {
      
      # The total number still susceptible is just the sum of who is left in 
      # exposed and unexposed groups at the end of previous period
      df[i+1, 'n'] = df[i, 'n_no_event_exposed_end'] + df[i, 'n_no_event_unexposed_end']
      
      # Then, the new number of individuals in the exposed and unexposed groups at the begining of the
      # next period is determined by p_exposed, p_unexposed in that period
      df[i+1, 'n_no_event_exposed_begin'] = df[i+1, 'n'] * df[i+1, 'p_exposed']
      df[i+1, 'n_no_event_unexposed_begin'] = df[i+1, 'n'] * df[i+1, 'p_unexposed']
    }
    
  }
  
  df = df %>%
    mutate(
      # Number of events in each exposure group is equal to number without the
      # event at the beginning, minus number without the event at the end
      n_event_exposed = n_no_event_exposed_begin - n_no_event_exposed_end,
      n_event_unexposed = n_no_event_unexposed_begin - n_no_event_unexposed_end,
      
      a1 = inst_rate_exposed * n_no_event_exposed_begin * exp(-delta_t),
      a0 = inst_rate_unexposed * n_no_event_unexposed_begin * exp(-delta_t),
      
      b0 = (a1 + a0) * (1 - p_exposed),
      b1 = (a1 + a0) * p_exposed,
      
      m10 = p_exposed * inst_rate_exposed * n_no_event_unexposed_begin * (1 - exp(-delta_t)),
      m01 = p_unexposed * inst_rate_unexposed * n_no_event_exposed_begin * (1 - exp(-delta_t)),
    )
  
  
    # The number of units of person time spent in each category
    # is given by the sum of the number of people in each category
    # at the *beginning* of each time period
    rate_exposed = sum(df$n_event_exposed) / sum(df$n_no_event_exposed_begin * df$delta_t)
    rate_unexposed = sum(df$n_event_unexposed) / sum(df$n_no_event_unexposed_begin * df$delta_t)
    
    RR = rate_exposed / rate_unexposed
    
    OR_matched = sum(df$m10) / sum(df$m01)
    OR_unmatched = (sum(df$a1) * sum(df$b0)) / (sum(df$a0) * sum(df$b1))
    
    return( list(df, c(RR, OR_matched, OR_unmatched)) )
}


