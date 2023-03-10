# EAVE II Analysis Repository
## [Insert full title here]  

### Paper Details
> **Short title:** [Insert short title of paper here]
>
>**DOI:** [Insert DOI of paper here]
>
>**Paper authors:** Steven Kerr, Sander Greenland, Igor Rudan

### Description

This repo contains code that calculates estimators of odds and rate ratios in a simulated case-control study. These estimators needn't be equal, and this simulator allows the magnitude of the disparity to be seen in a simple simulation that uses piecewise constant event rates.

simulator.R has one function that carries out a simulation of a cohort over a user-defined time period. There are two groups: exposed and unexposed. The event of interest occurs at a rate that is piecewise constant, i.e. it is constant for a while, then changes suddenly, is constant again, then changes suddenly, etc.

#### Inputs:
  - **times**:
      A numeric vector that defines start and end times
      of the time intervals where event rates are constant.
      e.g. c(1, 4, 5) specifies a cohort covering two time periods:
      t = 1-4, and t = 4-5
  - **num_chunks**:
      This is the discretisation parameter. The simulation involves
      numerical integration, and num_chunks is the number of chunks
      that each time interval is broken up into. Default is 1,000.
  - **period_event_rate_exposed**:
      The event rate in the exposed group over time intervals given
      by times. So if **times** = c(1, 4, 5), then
      **period_event_rate_exposed** = c(0.01, 0.02) means the event rate
      in exposed group is 0.01 for t = 1-4, and 0.02 for t = 4-5
  - **period_event_rate_unexposed**:
      Similar to **period_event_rate_exposed**, but for the unexposed
  - **p_exposed_start**:
      The proportion of the cohort that is exposed at the start
      of the simulation
  - **period_exposure_rate**:
      The rate at which unexposed individuals move to the exposed
      grou over time intervals given by times. So if times = c(1, 4, 5),
      then **period_exposure_rate** = c(0.1, 0.2) the rate of movement
      from unexposed to exposed is 0.1 for t = 1-4, and 0.2 for t = 4-5

#### Outputs:
A named vector whose components are
- Event rate ratio
- Event odds ratio in a matched study
- Event odds ratio in an unmatched study

Here, the event rate is the total number of events over the study, divided by the total number of person years.

### Funding
This analysis is part of the Early Assessment of COVID-19 epidemiology and Vaccine/anti-viral Effectiveness (EAVE II) study. EAVE II is funded by the Medical Research Council (MR/R008345/1) with the support of BREATHE - The Health Data Research Hub for Respiratory Health [MC_PC_19004], which is funded through the UK Research and Innovation Industrial Strategy Challenge Fund and delivered through Health Data Research UK. Additional support has been provided through the Scottish Government DG Health and Social Care.  
[Insert any other funding information on paper here]