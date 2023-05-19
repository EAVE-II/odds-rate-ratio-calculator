######################################################################
## Code author: Steven Kerr, steven.kerr@ed.ac.uk
######################################################################

library(tidyverse)

# This function carries out a simulation of a cohort over
# a user-defined time period. There are two groups: exposed
# and unexposed. The event of interest occurs at a rate that
# is piecewise constant, i.e. it is constant for a while, then
# changes suddenly, is constant again, then changes suddenly, etc.
# 
# The point is to illustrate how departures from a constant
# instantaneous event ratio and constant proportion exposed
# affect the magnitude of the discrepancy between the odds ratios,
# rate ratios and hazard ratios. If the time intervals are sufficiently 
# fine-grained, then this can be used to do calculations to any desired
# degree of accuracy
# 
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
#   - Outcome rate ratio
#   - Outcome odds ratio in a matched study
#   - Outcome odds ratio in an unmatched study
#   - Final proportion unexposed
# 
# Here, the event rate is the total number of events over the
# study, divided by the total number of person years.

calc_OR_RR = function(times, num_chunks = 1000,
  period_event_rate_exposed, period_event_rate_unexposed, 
  p_exposed_start, period_exposure_rate){
  
  # _1 refers to exposed, and _0 refers to unexposed
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
  #
  # Note that I have not multiplied by delta_t, because this factor will
  # cancel out in the OR_unmatched calculation below
  df = df %>%
    mutate(
      a_1 = n_1 * rate_1,
      a_0 = n_0 * rate_0,
      
      b_0 = (a_0 + a_1) * p_0,
      b_1 = (a_0 + a_1) * p_1,
      
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
  
  p_0 = df[nrow(df), 'p_0']
  
  return(c(
    'Rate ratio' = RR, 
    'Odds ratio - matched' = OR_matched, 
    'Odds ratio - unmatched' = OR_unmatched,
    'Proportion unexposed end' = p_0)
    )
}


### Simulation

# Example values
times = c(0, 1, 2, 3)

p_exposed_start = 0.2

# HR = (0.8, 0.4, 0.2)
period_event_rate_exposed = c(0.01, 0.01, 0.01)
period_event_rate_unexposed = c(0.0125, 0.025, 0.05)

period_exposure_rate = c(0.1, 0.1, 0.1)

calc_OR_RR(
  times = times,
  period_event_rate_exposed = period_event_rate_exposed, 
  period_event_rate_unexposed = period_event_rate_unexposed, 
  p_exposed_start = p_exposed_start, 
  period_exposure_rate = period_exposure_rate)


period_exposure_rate = c(0.1, 0.5, 0.1)

calc_OR_RR(
  times = times,
  period_event_rate_exposed =period_event_rate_exposed, 
  period_event_rate_unexposed = period_event_rate_unexposed, 
  p_exposed_start = p_exposed_start, 
  period_exposure_rate = period_exposure_rate)

period_exposure_rate = c(0, 0, 0)

calc_OR_RR(
  times = times,
  period_event_rate_exposed =period_event_rate_exposed, 
  period_event_rate_unexposed = period_event_rate_unexposed, 
  p_exposed_start = p_exposed_start, 
  period_exposure_rate = period_exposure_rate)

# HR = (0.4, 0.4, 0.4)
period_event_rate_exposed = c(0.01, 0.01, 0.01)
period_event_rate_unexposed = c(0.025, 0.025, 0.025)

period_exposure_rate = c(0.1, 0.1, 0.1)

calc_OR_RR(
  times = times,
  period_event_rate_exposed =period_event_rate_exposed, 
  period_event_rate_unexposed = period_event_rate_unexposed, 
  p_exposed_start = p_exposed_start, 
  period_exposure_rate = period_exposure_rate)

period_exposure_rate = c(0.1, 0.5, 0.1)

calc_OR_RR(
  times = times,
  period_event_rate_exposed =period_event_rate_exposed, 
  period_event_rate_unexposed = period_event_rate_unexposed, 
  p_exposed_start = p_exposed_start, 
  period_exposure_rate = period_exposure_rate)

period_exposure_rate = c(0, 0, 0)

calc_OR_RR(
  times = times,
  period_event_rate_exposed =period_event_rate_exposed, 
  period_event_rate_unexposed = period_event_rate_unexposed, 
  p_exposed_start = p_exposed_start, 
  period_exposure_rate = period_exposure_rate)



# Event rates switched for exposed/unexposed
# HR = (1.25, 2.5, 5)
period_event_rate_exposed = c(0.0125, 0.025, 0.05)
period_event_rate_unexposed = c(0.01, 0.01, 0.01)

period_exposure_rate = c(0.1, 0.1, 0.1)

calc_OR_RR(
  times = times,
  period_event_rate_exposed =period_event_rate_exposed, 
  period_event_rate_unexposed = period_event_rate_unexposed, 
  p_exposed_start = p_exposed_start, 
  period_exposure_rate = period_exposure_rate)

period_exposure_rate = c(0.1, 0.5, 0.1)

calc_OR_RR(
  times = times,
  period_event_rate_exposed =period_event_rate_exposed, 
  period_event_rate_unexposed = period_event_rate_unexposed, 
  p_exposed_start = p_exposed_start, 
  period_exposure_rate = period_exposure_rate)

period_exposure_rate = c(0, 0, 0)

calc_OR_RR(
  times = times,
  period_event_rate_exposed =period_event_rate_exposed, 
  period_event_rate_unexposed = period_event_rate_unexposed, 
  p_exposed_start = p_exposed_start, 
  period_exposure_rate = period_exposure_rate)

# HR = (2.5, 2.5, 2.5)
period_event_rate_exposed = c(0.025, 0.025, 0.025)
period_event_rate_unexposed = c(0.01, 0.01, 0.01)

period_exposure_rate = c(0.1, 0.1, 0.1)

calc_OR_RR(
  times = times,
  period_event_rate_exposed =period_event_rate_exposed, 
  period_event_rate_unexposed = period_event_rate_unexposed, 
  p_exposed_start = p_exposed_start, 
  period_exposure_rate = period_exposure_rate)

period_exposure_rate = c(0.1, 0.5, 0.1)

calc_OR_RR(
  times = times,
  period_event_rate_exposed =period_event_rate_exposed, 
  period_event_rate_unexposed = period_event_rate_unexposed, 
  p_exposed_start = p_exposed_start, 
  period_exposure_rate = period_exposure_rate)

period_exposure_rate = c(0, 0, 0)

calc_OR_RR(
  times = times,
  period_event_rate_exposed =period_event_rate_exposed, 
  period_event_rate_unexposed = period_event_rate_unexposed, 
  p_exposed_start = p_exposed_start, 
  period_exposure_rate = period_exposure_rate)

