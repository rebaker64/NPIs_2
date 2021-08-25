##########################################################

run_rsv_model <- function(b0 = 0.015421,
                          b1 = 0.39722,
                          phi = 0.98456,
                          prop_detected_1 = 0.424, # hospitalization rate
                          prop_detected_2 = 0.088,
                          prop_detected_3 = 0.047,
                          prop_detected_4 = 0.020,
                          mixing = mixing,
                          max_t = 2000,
                          init_conds_from_file = 0, 
                          timestart = 0, 
                          timeend = 0, 
                          betared = 0, 
                          maternalimmunity = TRUE,
                          popInput = 1861923, 
                          age_threshhold = 10){
  
  # parameters from user
  omega <- 0.6 # reduced infectiousness from older age groups
  delta <- 7.632 # latency rate
  gamma <- 3.342 # infectious rate
  nu <- 0.132 # rate of waning immunity
  dt <- 0.25
  T0 <- 0

  
  # proportions of infections that are detected for age groups 0-2m, 3-5m, 6-11m, 12-23m, and >=24m
  prop_detected_5 <- 0
  # maternal protection parameters - represent scaled susceptibility to infection in first few months of life
  sigma_1 <- 0.08
  sigma_2 <- 0.45
  sigma_3 <- 0.45
  sigma_4 <- 1
  sigma_5 <- 1
  sigma_6 <- 1
  if(maternalimmunity==FALSE){
    sigma_1 <- 1
    sigma_2 <- 1
    sigma_3 <- 1
    sigma_4 <- 1
    sigma_5 <- 1
    sigma_6 <- 1
  }
  
  ##########################################################
  # ageing-related parameters
  perth_pop <- popInput # 2014 greater perth pop, 0-79 years, ABS
  age_vect_years <- c(seq(0,5,1/12), seq(10,75,5))
  nAges <- length(age_vect_years) #number of age cohorts
  final_age <- 80
  age_vect_months <- age_vect_years*12
  size_cohorts_months <- c(diff(age_vect_months), final_age*12 - age_vect_months[length(age_vect_months)])
  trans_rate <- 1/size_cohorts_months
  rel_sizes <- size_cohorts_months/sum(size_cohorts_months)
  
  ##########################################################
  #  reduced infectiousness
  omega_vect <- as.vector(rep(1, nAges))
  omega_vect[!(age_vect_years < age_threshhold)] <- omega
  
  # proportion detected
  prop_detected_vect <- as.vector(c(rep(prop_detected_1, 3),
                                    rep(prop_detected_2, 3),
                                    rep(prop_detected_3, 6),
                                    rep(prop_detected_4, 12),
                                    rep(prop_detected_5, (nAges-24))))
  
  # Reduced susceptibility in youngest cohorts to reflect natural maternally-derived immunity
  # Reduced susceptibility in youngest cohorts to reflect natural maternally-derived immunity
  sigma_vect <- matrix(1, 1, nAges)
  sigma_vect[1] <- sigma_1
  sigma_vect[2] <- sigma_2
  sigma_vect[3] <- sigma_3
  sigma_vect[4] <- sigma_4
  sigma_vect[5] <- sigma_5
  sigma_vect[6] <- sigma_6
  sigma_vect <- as.vector(c(sigma_1, sigma_2, sigma_3, sigma_4, sigma_5, sigma_6, rep(1, (nAges-6))))
  
  ##################################################################
  # initial conditions: choose whether to reset here, or read in from csv file
  if (init_conds_from_file == 1) {
    ic <- readRDS("init_conds.rds")
    S0 <- ic$S0
    E0 <- ic$E0
    I0 <- ic$I0
    R0 <- ic$R0
    Incidence0 <- rel_sizes * 0
    DetIncidence0 <- rel_sizes * 0
  } else {
    I0 <- rel_sizes * perth_pop * 0.01
    S0 <- rel_sizes * perth_pop * 0.99
    E0 <- rel_sizes * perth_pop * 0
    R0 <- rel_sizes * perth_pop * 0
    DetIncidence0 <- rel_sizes * perth_pop * 0
    Incidence0 <- R0
  }
  
  ##########################################################
  
  pars <- list(
    b0 = b0,
    b1 = b1,
    phi = phi,
    delta = delta,
    gamma = gamma,
    nu = nu,
    prop_detected_vect = prop_detected_vect,
    sigma_vect = sigma_vect,
    omega_vect = omega_vect,
    mixing = mixing,
    timestart = timestart,
    timeend = timeend, 
    betared = betared,
    S0 = S0,
    E0 = E0,
    I0 = I0,
    R0 = R0,
    Incidence0 = Incidence0,
    DetIncidence0 = DetIncidence0
  )
  
  ##########################################################
  # run model with cohort ageing
  
  mod <- x(user = pars)
  
  pop_out <- NULL
  while (T0 <= max_t){
    
    # solve the odes first
    t <- seq(from = T0, to = T0 + 1, by = dt)
    m <- mod$run(t)
    pop <- mod$transform_variables(m)
    
    if (T0 == 0){
      pop_out <- pop
    } else {
      pop_out$time <- c(pop_out$time, pop$time[5])
      pop_out$S <- rbind(pop_out$S, pop$S[5,])
      pop_out$E <- rbind(pop_out$E, pop$E[5,])
      pop_out$I <- rbind(pop_out$I, pop$I[5,])
      pop_out$R <- rbind(pop_out$R, pop$R[5,])
      pop_out$Incidence <- rbind(pop_out$Incidence, pop$Incidence[5,])
      pop_out$DetIncidence <- rbind(pop_out$DetIncidence, pop$DetIncidence[5,])
    }
    
    # cohort ageing
    
    # extract the final state from pop
    S <- as.vector(t(data.table::last(pop$S)))
    E <- as.vector(t(data.table::last(pop$E)))
    I <- as.vector(t(data.table::last(pop$I)))
    R <- as.vector(t(data.table::last(pop$R)))
    Incidence <- as.vector(t(data.table::last(pop$Incidence)))
    DetIncidence <- as.vector(t(data.table::last(pop$DetIncidence)))
    
    # initialise the new initial condition vectors
    I0 <- rel_sizes*0
    S0 <- rel_sizes*0
    E0 <- rel_sizes*0
    R0 <- rel_sizes*0
    Incidence <- rel_sizes*0
    DetIncidence <- rel_sizes*0
    
    # then fill them in
    S0[1] <- perth_pop * rel_sizes[1]
    
    for(i in c(2:nAges)){
      S0[i] = S[(i - 1)] * trans_rate[(i-1)] + S[i] - S[i] * trans_rate[i]
      E0[i] = E[(i - 1)] * trans_rate[(i-1)] + E[i] - E[i] * trans_rate[i]
      I0[i] = I[(i - 1)] * trans_rate[(i-1)] + I[i] - I[i] * trans_rate[i]
      R0[i] = R[(i - 1)] * trans_rate[(i-1)] + R[i] - R[i] * trans_rate[i]
      Incidence0[i] = Incidence[(i - 1)] * trans_rate[(i-1)] + Incidence[i] - Incidence[i] * trans_rate[i]
      DetIncidence0[i] = DetIncidence[(i - 1)] * trans_rate[(i-1)] + DetIncidence[i] - DetIncidence[i] * trans_rate[i]
    }
    
    pars <- list(
      b0 = b0,
      b1 = b1,
      phi = phi,
      delta = delta,
      gamma = gamma,
      nu = nu,
      prop_detected_vect = prop_detected_vect,
      sigma_vect = sigma_vect,
      omega_vect = omega_vect,
      mixing = mixing,
      timestart = timestart,
      timeend = timeend, 
      betared = betared,
      S0 = S0,
      E0 = E0,
      I0 = I0,
      R0 = R0,
      Incidence0 = Incidence0,
      DetIncidence0 = DetIncidence0
    )
    
    if(T0 >= timestart & T0 <= timeend){
    pars <- list(
      b0 = b0*betared,
      b1 = b1,
      phi = phi,
      delta = delta,
      gamma = gamma,
      nu = nu,
      prop_detected_vect = prop_detected_vect,
      sigma_vect = sigma_vect,
      omega_vect = omega_vect,
      mixing = mixing,
      timestart = timestart,
      timeend = timeend, 
      betared = betared,
      S0 = S0,
      E0 = E0,
      I0 = I0,
      R0 = R0,
      Incidence0 = Incidence0,
      DetIncidence0 = DetIncidence0
    )
    }
    mod <- x(user = pars) # still in while loop
    T0 <- T0 + 1
    
  }
  pop_out <- pop_out[2:8]
  return(pop_out)
}

##########################################################