######################################################################
## Code author: Steven Kerr, steven.kerr@ed.ac.uk
######################################################################

library(tidyverse)

### Simulator function

# This function carries out a simulation of a cohort over
# a user-defined time period. There are two groups: exposed
# and unexposed. The event of interest occurs at a rate that
# is piecewise constant, i.e. it is constant for a while, then
# changes suddenly, is constant again, then changes suddenly, etc
# 
# 
# The point is to illustrate how departures from a constant
# instantaneous event ratio and constant proportion exposed
# affect the magnitude of the discrepancy between the above quantities.
# 
# Inputs:
#   - times:
#       A numeric vector that defines start and end times
#       of the time intervals where event rates are constant.
#       e.g. c(1, 4, 5) specifies a cohort covering two time periods:
#       t = 1-4, and t = 4-5
#   - num_chunks:
#       This is the discretisation parameter. The simulation involves
#       numerical integration, and num_chunks is the number of chunks
#       that each time interval is broken up into. Default is 1,000.
#   - period_event_rate_exposed:
#       The event rate in the exposed group over time intervals given
#       by times. So if times = c(1, 4, 5), then
#       period_event_rate_exposed = c(0.01, 0.02) means the event rate
#       in exposed group is 0.01 for t = 1-4, and 0.02 for t = 4-5
#   - period_event_rate_unexposed:
#       Similar to period_event_rate_exposed, but for the unexposed
#   - p_exposed_start:
#       The proportion of the cohort that is exposed at the start
#       of the simulation
#   - period_exposure_rate:
#       The rate at which unexposed individuals move to the exposed
#       grou over time intervals given by times. So if times = c(1, 4, 5),
#       then period_exposure_rate = c(0.1, 0.2) the rate of movement
#       from unexposed to exposed is 0.1 for t = 1-4, and 0.2 for t = 4-5
#
# Outputs:
# A named vector whose elements are:
#   - Event rate ratio
#   - Event odds ratio in a matched study
#   - Event odds ratio in an unmatched study
# 
# Here, the event rate is the total number of events over the
# study, divided by the total number of person years.

calc_OR_RR = function(times, num_chunks = 1000,
  period_event_rate_exposed, period_event_rate_unexposed, 
  p_exposed_start, period_exposure_rate){
  
  rate_1 = c()
  rate_0 = c()
  exposure_rate = c()
  
  for (i in 1:(length(times)-1) ){
    rate_1 = c(rate_1, rep(period_event_rate_exposed[i], num_chunks))
    rate_0 = c(rate_0, rep(period_event_rate_unexposed[i], num_chunks))
    
    exposure_rate = c(exposure_rate, rep(period_exposure_rate[i], num_chunks))
  }
  
  delta_t = 1/num_chunks
  
  df = data.frame(
    t_start = seq(times[1], times[length(times)] - delta_t, by = delta_t),

    rate_1 = rate_1,
    rate_0 = rate_0,
    
    exposure_rate = exposure_rate,
    
    n = NA,
    n_0 = NA,
    n_1 = NA,
    
    n_0_event = NA,
    n_1_event = NA,
    
    newly_exposed = NA
  )
  
  df[ , 't_end'] = df[ , 't_start'] + delta_t
  
  df[1, 'n'] = 1
  df[1, 'p_1'] = p_exposed_start
  df[1, 'n_0'] = (1 - p_exposed_start) * df[1, 'n']
  df[1, 'n_1'] = p_exposed_start * df[1, 'n']

  # Step-by-step simulation starts here
  for (i in 1:nrow(df)){
    df[i, 'n_0_event'] = df[i, 'rate_0'] * delta_t * df[i, 'n_0']
    df[i, 'n_1_event'] = df[i, 'rate_1'] * delta_t * df[i, 'n_1']

    df[i, 'newly_exposed'] = df[i, 'exposure_rate'] * delta_t * df[i, 'n_0']
    
    if (i < (nrow(df))){
    
      df[i+1, 'n_0'] = df[i, 'n_0'] - df[i, 'n_0_event'] - df[i, 'newly_exposed']
      df[i+1, 'n_1'] = df[i, 'n_1'] - df[i, 'n_1_event'] + df[i, 'newly_exposed']
    }
  }
  
  df = df %>%
    mutate(
      n = n_0 + n_1,
      
      p_0 = n_0/n,
      p_1 = 1 - p_0
    )
  
  # These are quantities that are used to estimate matched and unmatched
  # odds ratios. Taken from OR_ij and OR_m on page 549 of the following 
  # paper:
  #
  # https://doi.org/10.1093/oxfordjournals.aje.a113439
  df = df %>%
    mutate(
      a_1 = n_1 * rate_1 * delta_t,
      a_0 = n_0 * rate_0 * delta_t,
      
      b_0 = (a_0 + a_1) * p_0 * delta_t,
      b_1 = (a_0 + a_1) * p_1 * delta_t,
      
      m_10 = a_1 * p_0,
      m_01 = a_0 * p_1
    )
  
  # The number of units of person time spent in each category
  # is given by the sum of the number of people in each category
  # at the *beginning* of each time period
  rate_1 = sum(df$n_1_event) / sum(df$n_1 * delta_t)
  rate_0 = sum(df$n_0_event) / sum(df$n_0 * delta_t)
  
  RR = rate_1 / rate_0
  
  OR_matched = sum(df$m_10) / sum(df$m_01)
  OR_unmatched = (sum(df$a_1) * sum(df$b_0)) / (sum(df$a_0) * sum(df$b_1))
  
  return(c('Rate ratio' = RR, 'Odds ratio - matched' = OR_matched, 'Odds ratio - unmatched' = OR_unmatched))
}


### Simulation

# Example values
times = c(1, 2, 3, 4, 5)

delta_t = 0.001

period_event_rate_exposed = c(0.01, 0.02, 0.02, 0.01)
period_event_rate_unexposed = c(0.01, 0.04, 0.07, 0.03)

p_exposed_start = 0.3
period_exposure_rate = c(0.1, 0.2, 0.1, 0.2)

calc_OR_RR(
  times, delta_t,
  period_event_rate_exposed, period_event_rate_unexposed, 
  p_exposed_start, period_exposure_rate)
