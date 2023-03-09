library(tidyverse)


# df = data.frame( 
#   scenario = c('a', 'b', 'b', 'c', 'c', 'd', 'd', 'e', 'e', 'f', 'f'),
#   t = c(10, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5),
#   
#   n_exposed   = c( rep(1000, 9), 1000, 1500),
#   n_unexposed = c( rep(1000, 9), 1000, 500),
#   
#   inst_rate_exposed =   c(2,  2, 1,  5, 1,  10, 2,  100, 10,  100, 10),
#   inst_rate_unexposed = c(1,  1, 2,  1, 1,  2,  2,  10,  10,   10, 10)
#   )
# 
# 
# df = df %>%
#   mutate(
#     t = 1:10,
#     
#     pyears_exposed = t * n_exposed,
#     pyears_unexposed = t * n_unexposed,
#     
#     # Rates are measured in number of events per thousand person years
#     n_event_exposed = pyears_exposed * inst_rate_exposed / 1000,
#     n_event_unexposed = pyears_unexposed * inst_rate_unexposed / 1000,
#     
#     p_event_exposed = n_event_exposed / n_exposed,
#     p_event_unexposed = n_event_unexposed / n_unexposed,
#     
#     # Proportion of people that are in exposed and unexposed groups
#     p_exposed   = c(0.2,  0.2, 0.5,  0.2, 0.8,  0.2, 0.8,  0.8, 0.2,  0.5, 0.7),
#     p_unexposed = 1 - p_exposed,
#     
#     inst_t_odds_exposed = n_event_exposed / (pyears_exposed - n_event_exposed),
#     inst_t_odds_unexposed = n_event_unexposed / (pyears_unexposed - n_event_unexposed),
#     
#     inst_OR = inst_t_odds_exposed / inst_t_odds_unexposed,
#     inst_RR = inst_rate_exposed / inst_rate_unexposed,
#     
#     a1 = n_exposed * inst_rate_exposed,
#     a0 = n_unexposed * inst_rate_unexposed,
#     
#     b0 = (a1 + a0) * (1 - p_exposed),
#     b1 = (a1 + a0) * p_exposed,
#     
#     m10 = n_exposed * inst_rate_exposed * p_unexposed,
#     m01 = n_unexposed * inst_rate_unexposed * p_exposed,
#     
#     #matched_weight = n_exposed * inst_rate_exposed * p_unexposed * inst_rate_unexposed
#   )
# 
# 
# df = df %>%
#   group_by(scenario) %>%
#   mutate(
#     odds_exposed = sum(n_event_exposed) / (sum(pyears_exposed) - sum(n_event_exposed)),
#     odds_unexposed = sum(n_event_unexposed) / (sum(pyears_unexposed) - sum(n_event_unexposed)),
#     
#     rate_exposed = sum(n_event_exposed) / sum(pyears_exposed),
#     rate_unexposed = sum(n_event_unexposed) / sum(pyears_unexposed),
#     
#     OR = odds_exposed / odds_unexposed,
#     RR = rate_exposed / rate_unexposed,
#     
#     OR_m = sum(m10 * t) / sum(m01 * t),
#     OR_um = (sum(a1 * t) * sum(b0 * t)) / (sum(a0 * t) * sum(b1 * t))
#   ) %>%
#   ungroup()
# 
# 
# bob = df %>% 
#   select(
#     scenario,
#     t,
#     
#     pyears_exposed,
#     pyears_unexposed,
#     
#     inst_rate_exposed,
#     inst_rate_unexposed,
#     
#     n_event_exposed,
#     n_event_unexposed,
#     
#     OR,
#     RR,
#     
#     OR_m,
#     OR_um
#   )


##########

df = data.frame( 
  t_start = 1:4,
  t_end = 2:5,
  
  p_exposed = c(0.2, 0.2, 0.2, 0.2),
  
  inst_rate_exposed =   c(20, 20, 20, 20)/1000,
  inst_rate_unexposed = rep(10, 4)/1000,
  
  n = c(1000, rep(NA, 3)),
  
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

for (i in 1:nrow(df)){
  
  # The end quantities are related by exponential decay to the begin quantities
  df[i, 'n_no_event_exposed_end'] = df[i, 'n_no_event_exposed_begin'] * exp( - df[i, 'inst_rate_exposed'] * df[i, 'delta_t'] ) 
  df[i, 'n_no_event_unexposed_end'] = df[i, 'n_no_event_unexposed_begin'] * exp( - df[i, 'inst_rate_unexposed'] * df[i, 'delta_t']) 
  
  if (i == nrow(df)){
    
  } else {
  
    # The total number still susceptible is just the sum of who is left in 
    # exposed and unexposed groups at the end of previous period
    df[i+1, 'n'] = df[i, 'n_no_event_exposed_end'] + df[i, 'n_no_event_unexposed_end']
    
    # Then, the new number of individuals in the exposed and unexposed groups at the begining of the
    # next period is determined by p_exposed, p_unexposed in that period
    df[i+1, 'n_no_event_exposed_begin'] = df[i+1, 'n'] * df[i+1, 'p_exposed']
    df[i+1, 'n_no_event_unexposed_begin'] = df[i+1, 'n'] * df[i+1, 'p_unexposed']
  }

  # df[t, 'n_event_exposed'] = df[t, 'n'] * df[t, 'p_exposed'] * df[t, 'inst_rate_exposed']
  # df[t, 'n_event_unexposed'] = df[t, 'n'] * df[t, 'p_unexposed'] * df[t, 'inst_rate_unexposed']
  # 
  # if (t!= max(df$t)){
  #   df[t+1, 'n'] = df[t, 'n'] - df[t, 'n_event_exposed'] - df[t, 'n_event_unexposed']
  # }
}



df = df %>%
  mutate(
    #n = n_no_event_exposed + n_no_event_unexposed,
    
    #n_exposed = n * p_exposed,
    #n_unexposed = n * p_unexposed,
    
    # n_event_exposed = n_exposed - n_no_event_exposed,
    # n_event_unexposed = n_unexposed - n_no_event_unexposed,
    
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
    
    # a1 = n * p_exposed * inst_rate_exposed,
    # a0 = n * p_unexposed * inst_rate_unexposed,
    # 
    # b0 = (a1 + a0) * (1 - p_exposed),
    # b1 = (a1 + a0) * p_exposed,
    # 
    # m10 = n * p_exposed* inst_rate_exposed * p_unexposed,
    # m01 = n * p_unexposed * inst_rate_unexposed * p_exposed,
    
    #matched_weight = n_exposed * inst_rate_exposed * p_unexposed * inst_rate_unexposed
  )


df = df %>%
  mutate(
    # rate_exposed = sum(n_event_exposed) / sum(n_exposed),
    # rate_unexposed = sum(n_event_unexposed) / sum(n_unexposed),
    
    # The number of units of person time spent in each category
    # is given by the sum of the number of people in each category
    # at the *beginning* of each time period
    rate_exposed = sum(n_event_exposed) / sum(n_no_event_exposed_begin * delta_t),
    rate_unexposed = sum(n_event_unexposed) / sum(n_no_event_unexposed_begin * delta_t),
    
    RR = rate_exposed / rate_unexposed,
    
    OR_m = sum(m10) / sum(m01),
    OR_um = (sum(a1) * sum(b0)) / (sum(a0) * sum(b1))
  )



df$n_event_exposed / df$n_no_event_exposed_begin * df$delta_t


