#### FOR SIMULATIONS ALTERING TREATMENT START DAY ####

library(pacman)
p_load(deSolve, ggplot2, tidyverse, coda, future, future.apply)

## Initial conditions
y0_SEIVRY <- c(
  S = 1e7,
  E = 0,        # no initial infected cells in eclipse phase
  I = 1,        # small initial infected cells
  V = 10,       # initial viral load
  R = 0,        # no innate immune response
  Y = 0,        # no adaptive immune response
  GI = 0,
  C = 0,
  P = 0,
  Rtv = 0     # proxy state for Rtv to allow dosing table to work
)

## Parameters - moderate efficacy
pars_SEIVRY <- c(
  tstart = 0,     # start time of infection relative to t=0

  # Target cell limited model 
  beta = 2e-8,          # infection rate per virion per day
  k = 4,                # rate of transition from eclipse phase to productively infected phase
  delta = 1.18,         # infected cell death rate       
  p = 500,            # virion production rate per infected cell per day
  c = 6,                # virion clearance rate without immunity (per day)
  a = 0.5,                # how quickly adaptive immunity ramps up once it starts
  
  # Innate immunity
  phi = 0.000005,     # rate of conversion from susceptible to refractory cells
  rho = 0.025,      # rate of reversion from refractory to susceptible cells

  # Adaptive immunity
  tau = 7,        # onset of adaptive immunity (days since infection); 5-10 days
  m = 0.00001,     # infected cell death rate due to adaptive immunity
  sigma = 0.00002,      # virion clearance rate with adaptive immunity
  omega = 0.7,      # rate of infected cell triggering adaptive immunity
  gamma = 0.01,    # immunity cell death rate
  
  # Standard deviation of observation error
  error = 1,

  # Switch to turn treatment on or off
  treat_on = 1,   # 1 = on, 0 = off

  # PK model
  ka = 15,      # absorption rate of nirmatrelvir
  F0 = 0.5,         # bioavailability of nirmatrelvir
  alpha_F = 0.5,  # ritonavir effect on bioavailability
  alpha_cl = 0.7, # ritonavir effect on clearance
  kcp = 1.3,     # rate of distribution from central to peripheral
  kpc = 1.6,     # rate of distribution from peripheral to central
  kcl = 5,        # drug clearance rate    

  # Ritonavir proxy
  krtv = 0.3,     # ritonavir rate of decay 
  
  # PD model
  vol_c = 40000,  # volume of central compartment 
  Emax = 0.9999,  
  IC50 = 0.002, 
  Hill_coeff = 4
)

## Model function 
SEIVRY_model <- function(t,y,pars){
  with(as.list(c(y,pars)),{
  
    ## PK MODEL
    Rtv_unit <- Rtv / 100      # normalize ritonavir dose units

    F_eff <- F0 * (1 + alpha_F * Rtv_unit)    # bioavailabilty of nirmatrelvir, dependent on concentration of ritonavir
    F_eff <- min(F_eff, 1)                    # limits value to 1

    kcl_eff <- kcl / (1 + alpha_cl * Rtv_unit)  # clearance rate of nirmatrelvir dependent on concentration of ritonavir

    dGI = -ka*GI
    dC = F_eff*ka*GI + kpc*P - kcp*C - kcl_eff*C
    dP = kcp*C - kpc*P

    # Ritonavir proxy to allow dosing table to work
    dRtv = -krtv*Rtv

    ## PD MODEL
    conc = C / vol_c    # calculate concentration from 'C' compartment output 
    epsilon <- treat_on * Emax * conc^Hill_coeff / (IC50^Hill_coeff + conc^Hill_coeff) # efficacy, with treatment on/off switch
    # epsilon <- Emax * conc^Hill_coeff / (IC50^Hill_coeff + conc^Hill_coeff)

    ## WITHIN-HOST TARGET CELL LIMITED MODEL
    if(t < tstart){
      beta <- phi <- rho <- k <- delta <- m <- epsilon <- p <- c <- sigma <- omega <- gamma <- 0 
    }
    
    # Set condition where infection clears if below 1 virion
    if(V < 10) beta <- 0
    
    # g <- ifelse(t < tau, 0, 1) # time-dependent activation of adaptive immunity (Y) - hard switch
    g <- 1 / (1 + exp(-a * (t - tau))) # time-dependent activation of adaptive immunity (Y) - smooth switch

    dS = -beta*S*V - phi*S*I + rho*R
    dE = beta*S*V - k*E
    dI = k*E - delta*I - m*Y*I
    dV = (1 - epsilon)*p*I - c*V - sigma*Y*V
    dR = phi*S*I - rho*R
    dY = omega*g*I - gamma*Y
    
    return(list(c(dS, dE, dI, dV, dR, dY, dGI, dC, dP, dRtv), conc = conc, epsilon = epsilon))
  })
}

### ANTIVIRAL DRUG TREATMENT REGIME - 1 DAY AFTER INFECTION ### 
start_time <- 1     # start day of antiviral treatment 
dose_times <- start_time + seq(0, 4.5, by=0.5)  # 10 doses over 5 days (every 12 hours for 5 days)

times_out <- seq(0, 40, by = 0.01)  # simulate for 40 days, with 15-minute time-steps
times_out <- sort(unique(c(times_out, dose_times)))   # makes sure there is no repeat times between times_out and dose_times

drug_regime <- data.frame(
  var = rep(c("GI", "Rtv"), each = length(dose_times)),
  time = rep(dose_times, times = 2),
  value = c(rep(300, length(dose_times)), rep(100, length(dose_times))),
  method = rep("add", 2 * length(dose_times))
)
drug_regime <- drug_regime[order(drug_regime$time, drug_regime$var), ]  # orders rows by time to make sure 'events' works

## Solve
out_1 <- ode(y = y0_SEIVRY, times = times_out, func = SEIVRY_model, parms = pars_SEIVRY, events = list(data = drug_regime))
out_1_df <- as.data.frame(out_1)


### ANTIVIRAL DRUG TREATMENT REGIME - 3 DAYS AFTER INFECTION ### 
start_time <- 3     # start day of antiviral treatment 
dose_times <- start_time + seq(0, 4.5, by=0.5)  # 10 doses over 5 days (every 12 hours for 5 days)

times_out <- seq(0, 40, by = 0.01)  # simulate for 40 days, with 15-minute time-steps
times_out <- sort(unique(c(times_out, dose_times)))   # makes sure there is no repeat times between times_out and dose_times

drug_regime <- data.frame(
  var = rep(c("GI", "Rtv"), each = length(dose_times)),
  time = rep(dose_times, times = 2),
  value = c(rep(300, length(dose_times)), rep(100, length(dose_times))),
  method = rep("add", 2 * length(dose_times))
)
drug_regime <- drug_regime[order(drug_regime$time, drug_regime$var), ]  # orders rows by time to make sure 'events' works

## Solve
out_3 <- ode(y = y0_SEIVRY, times = times_out, func = SEIVRY_model, parms = pars_SEIVRY, events = list(data = drug_regime))
out_3_df <- as.data.frame(out_3)


### ANTIVIRAL DRUG TREATMENT REGIME - 1 DAY AFTER SYMPTOM ONSET ### 
start_time <- 5     # start day of antiviral treatment 
dose_times <- start_time + seq(0, 4.5, by=0.5)  # 10 doses over 5 days (every 12 hours for 5 days)

times_out <- seq(0, 40, by = 0.01)  # simulate for 40 days, with 15-minute time-steps
times_out <- sort(unique(c(times_out, dose_times)))   # makes sure there is no repeat times between times_out and dose_times

drug_regime <- data.frame(
  var = rep(c("GI", "Rtv"), each = length(dose_times)),
  time = rep(dose_times, times = 2),
  value = c(rep(300, length(dose_times)), rep(100, length(dose_times))),
  method = rep("add", 2 * length(dose_times))
)
drug_regime <- drug_regime[order(drug_regime$time, drug_regime$var), ]  # orders rows by time to make sure 'events' works

## Solve
out_5 <- ode(y = y0_SEIVRY, times = times_out, func = SEIVRY_model, parms = pars_SEIVRY, events = list(data = drug_regime))
out_5_df <- as.data.frame(out_5)


### ANTIVIRAL DRUG TREATMENT REGIME - 5 DAYS AFTER SYMPTOM ONSET ### 
start_time <- 9     # start day of antiviral treatment 
dose_times <- start_time + seq(0, 4.5, by=0.5)  # 10 doses over 5 days (every 12 hours for 5 days)

times_out <- seq(0, 40, by = 0.01)  # simulate for 40 days, with 15-minute time-steps
times_out <- sort(unique(c(times_out, dose_times)))   # makes sure there is no repeat times between times_out and dose_times

drug_regime <- data.frame(
  var = rep(c("GI", "Rtv"), each = length(dose_times)),
  time = rep(dose_times, times = 2),
  value = c(rep(300, length(dose_times)), rep(100, length(dose_times))),
  method = rep("add", 2 * length(dose_times))
)
drug_regime <- drug_regime[order(drug_regime$time, drug_regime$var), ]  # orders rows by time to make sure 'events' works

## Solve
out_9 <- ode(y = y0_SEIVRY, times = times_out, func = SEIVRY_model, parms = pars_SEIVRY, events = list(data = drug_regime))
out_9_df <- as.data.frame(out_9)


### ANTIVIRAL DRUG TREATMENT REGIME - 10 DAYS AFTER SYMPTOM ONSET ### 
start_time <- 14     # start day of antiviral treatment 
dose_times <- start_time + seq(0, 4.5, by=0.5)  # 10 doses over 5 days (every 12 hours for 5 days)

times_out <- seq(0, 40, by = 0.01)  # simulate for 40 days, with 15-minute time-steps
times_out <- sort(unique(c(times_out, dose_times)))   # makes sure there is no repeat times between times_out and dose_times

drug_regime <- data.frame(
  var = rep(c("GI", "Rtv"), each = length(dose_times)),
  time = rep(dose_times, times = 2),
  value = c(rep(300, length(dose_times)), rep(100, length(dose_times))),
  method = rep("add", 2 * length(dose_times))
)
drug_regime <- drug_regime[order(drug_regime$time, drug_regime$var), ]  # orders rows by time to make sure 'events' works

## Solve
out_14 <- ode(y = y0_SEIVRY, times = times_out, func = SEIVRY_model, parms = pars_SEIVRY, events = list(data = drug_regime))
out_14_df <- as.data.frame(out_14)


### NO TREATMENT ###

## Parameters
pars_notreat <- c(
  tstart = 0,     # start time of infection relative to t=0

  # Target cell limited model 
  beta = 2e-8,          # infection rate per virion per day
  k = 4,                # rate of transition from eclipse phase to productively infected phase
  delta = 1.18,         # infected cell death rate       
  p = 500,            # virion production rate per infected cell per day
  c = 6,                # virion clearance rate without immunity (per day)
  a = 0.5,                # how quickly adaptive immunity ramps up once it starts
  
  # Innate immunity
  phi = 0.000005,     # rate of conversion from susceptible to refractory cells
  rho = 0.025,      # rate of reversion from refractory to susceptible cells

  # Adaptive immunity
  tau = 7,        # onset of adaptive immunity (days since infection); 5-10 days
  m = 0.00001,     # infected cell death rate due to adaptive immunity
  sigma = 0.00002,      # virion clearance rate with adaptive immunity
  omega = 0.7,      # rate of infected cell triggering adaptive immunity
  gamma = 0.01,    # immunity cell death rate
  
  # Standard deviation of observation error
  error = 1,

  # Switch to turn treatment on or off
  treat_on = 0,   # 1 = on, 0 = off

  # PK model
  ka = 15,      # absorption rate of nirmatrelvir
  F0 = 0.5,         # bioavailability of nirmatrelvir
  alpha_F = 0.5,  # ritonavir effect on bioavailability
  alpha_cl = 0.7, # ritonavir effect on clearance
  kcp = 1.3,     # rate of distribution from central to peripheral
  kpc = 1.6,     # rate of distribution from peripheral to central
  kcl = 5,        # drug clearance rate    

  # Ritonavir proxy
  krtv = 0.3,     # ritonavir rate of decay 
  
  # PD model
  vol_c = 40000,  # volume of central compartment 
  Emax = 0.9999,  
  IC50 = 0.0005, 
  Hill_coeff = 4
)

## Solve
out_notreat <- ode(y = y0_SEIVRY, times = times_out, func = SEIVRY_model, parms = pars_notreat, events = list(data = drug_regime))
out_notreat_df <- as.data.frame(out_notreat)


### LOG TRANSFORM ###
eps <- 1e-8  # small constant to avoid log(0)

out1_log <- out_1_df %>%
  mutate(across(c(S, E, I, V, R, Y, GI, C, P, conc), 
                ~ log10(.x + eps),
                .names = "log10_{.col}")
        )

out3_log <- out_3_df %>%
  mutate(across(c(S, E, I, V, R, Y, GI, C, P, conc), 
                ~ log10(.x + eps),
                .names = "log10_{.col}")
        )

out5_log <- out_5_df %>%
  mutate(across(c(S, E, I, V, R, Y, GI, C, P, conc), 
                ~ log10(.x + eps),
                .names = "log10_{.col}")
        )

out9_log <- out_9_df %>%
  mutate(across(c(S, E, I, V, R, Y, GI, C, P, conc), 
                ~ log10(.x + eps),
                .names = "log10_{.col}")
        )

out14_log <- out_14_df %>%
  mutate(across(c(S, E, I, V, R, Y, GI, C, P, conc), 
                ~ log10(.x + eps),
                .names = "log10_{.col}")
        )

out_notreat_log <- out_notreat_df %>%
  mutate(across(c(S, E, I, V, R, Y, GI, C, P, conc), 
                ~ log10(.x + eps),
                .names = "log10_{.col}")
        )


### MAX VIRAL LOAD (PEAK) FOR ALL SCENARIOS ###
print(max(out1_log$log10_V))
print(max(out3_log$log10_V))
print(max(out5_log$log10_V))
print(max(out9_log$log10_V))
print(max(out14_log$log10_V))
print(max(out_notreat_log$log10_V))

### TIME TO CLEARANCE FORM INFECTION FOR ALL SCENARIOS ###
threshold <- 2
peak_idx <- which.max(out_notreat_log$log10_V) # index of initial peak
post_peak <- out_notreat_log[(peak_idx + 1):nrow(out_notreat_log), ] # data after the peak
clear_idx <- which(post_peak$log10_V <= threshold)[1] # first time after the peak that viral load is at or below threshold
clear_time <- post_peak$time[clear_idx]
time_from_peak <- clear_time - out_notreat_log$time[peak_idx] # time from peak to clearance threshold

clear_time
time_from_peak


### COMBINE AND PLOT ###
day1_V <- out1_log %>%
  select(time, log10_V) %>%
  mutate(group = "Day 1")

day3_V <- out3_log %>%
  select(time, log10_V) %>%
  mutate(group = "Day 3")

day5_V <- out5_log %>%
  select(time, log10_V) %>%
  mutate(group = "1 day post symptom")

day9_V <- out9_log %>%
  select(time, log10_V) %>%
  mutate(group = "5 days post symptom")

day14_V <- out14_log %>%
  select(time, log10_V) %>%
  mutate(group = "10 days post symptom")

notreat_V <- out_notreat_log %>%
  select(time, log10_V) %>%
  mutate(group = "No treatment")

# Combine and plot
combined_V <- bind_rows(
  day1_V, 
  day3_V,
  day5_V,
  day9_V,
  day14_V,
  notreat_V
)

combined_V <- combined_V %>%
  mutate(group = factor(group, levels = c(
    "Day 1",
    "Day 3",
    "1 day post symptom",
    "5 days post symptom",
    "10 days post symptom",
    "No treatment"
  )))

treatstart_plot <- ggplot(combined_V, aes(x = time, y = log10_V, color = group)) +
  geom_line(linewidth = 1) +
  labs(
    x = "Time (days)",
    y = "Log10 Viral Load (copies/mL)",
    color = "Scenario"
  ) +
  theme_bw(base_size = 20) + 
  theme(legend.position = "bottom", panel.grid = element_blank())

ggsave(filename = "treatstart_simulations.png",
       plot = treatstart_plot,
       device = "png",
       width = 15, units = "in")


### COMBINE AND PLOT FOR IMMUNE RESPONSE - DAY 1 ###
day1_V <- out1_log %>%
  select(time, log10_S, log10_I, log10_R, log10_Y) %>%
  mutate(group = "Day 1")

df_long_day1 <- day1_V %>%
  pivot_longer(
    cols = c(log10_S, log10_I, log10_R, log10_Y),
    names_to = "compartment",
    values_to = "value"
  ) %>%
  mutate(compartment = sub("^log10_", "", compartment))

day1_immune_plot <- ggplot(df_long_day1, aes(x = time, y = value, color = compartment)) +
  geom_line(linewidth = 1) +
  labs(
    x = "Time (days)",
    y = "Log10 Cell Count",
    color = "Compartment"
  ) +
  theme_bw(base_size = 18) + 
  theme(legend.position = "bottom", panel.grid = element_blank())

ggsave(filename = "efficacy_immunity_day1.png",
       plot = day1_immune_plot,
       device = "png",
       height = 5, width = 10, units = "in")


### COMBINE AND PLOT FOR IMMUNE RESPONSE - DAY 3 ###
day3_V <- out3_log %>%
  select(time, log10_S, log10_I, log10_R, log10_Y) %>%
  mutate(group = "Day 3")

df_long_day3 <- day3_V %>%
  pivot_longer(
    cols = c(log10_S, log10_I, log10_R, log10_Y),
    names_to = "compartment",
    values_to = "value"
  ) %>%
  mutate(compartment = sub("^log10_", "", compartment))

day3_immune_plot <- ggplot(df_long_day3, aes(x = time, y = value, color = compartment)) +
  geom_line(linewidth = 1) +
  labs(
    x = "Time (days)",
    y = "Log10 Cell Count",
    color = "Compartment"
  ) +
  theme_bw(base_size = 18) + 
  theme(legend.position = "bottom", panel.grid = element_blank())

ggsave(filename = "efficacy_immunity_day3.png",
       plot = day3_immune_plot,
       device = "png",
       height = 5, width = 10, units = "in")


### COMBINE AND PLOT FOR IMMUNE RESPONSE - DAY 5 ###
day5_V <- out5_log %>%
  select(time, log10_S, log10_I, log10_R, log10_Y) %>%
  mutate(group = "1 day post symptom")

df_long_day5 <- day5_V %>%
  pivot_longer(
    cols = c(log10_S, log10_I, log10_R, log10_Y),
    names_to = "compartment",
    values_to = "value"
  ) %>%
  mutate(compartment = sub("^log10_", "", compartment))

day5_immune_plot <- ggplot(df_long_day5, aes(x = time, y = value, color = compartment)) +
  geom_line(linewidth = 1) +
  labs(
    x = "Time (days)",
    y = "Log10 Cell Count",
    color = "Compartment"
  ) +
  theme_bw(base_size = 18) + 
  theme(legend.position = "bottom", panel.grid = element_blank())

ggsave(filename = "efficacy_immunity_day5.png",
       plot = day5_immune_plot,
       device = "png",
       height = 5, width = 10, units = "in")


### COMBINE AND PLOT FOR IMMUNE RESPONSE - DAY 9 ###
day9_V <- out9_log %>%
  select(time, log10_S, log10_I, log10_R, log10_Y) %>%
  mutate(group = "5 days post symptom")

df_long_day9 <- day9_V %>%
  pivot_longer(
    cols = c(log10_S, log10_I, log10_R, log10_Y),
    names_to = "compartment",
    values_to = "value"
  ) %>%
  mutate(compartment = sub("^log10_", "", compartment))

day9_immune_plot <- ggplot(df_long_day9, aes(x = time, y = value, color = compartment)) +
  geom_line(linewidth = 1) +
  labs(
    x = "Time (days)",
    y = "Log10 Cell Count",
    color = "Compartment"
  ) +
  theme_bw(base_size = 18) + 
  theme(legend.position = "bottom", panel.grid = element_blank())

ggsave(filename = "efficacy_immunity_day9.png",
       plot = day9_immune_plot,
       device = "png",
       height = 5, width = 10, units = "in")


### COMBINE AND PLOT FOR IMMUNE RESPONSE - DAY 14 ###
day14_V <- out14_log %>%
  select(time, log10_S, log10_I, log10_R, log10_Y) %>%
  mutate(group = "10 days post symptom")

df_long_day14 <- day14_V %>%
  pivot_longer(
    cols = c(log10_S, log10_I, log10_R, log10_Y),
    names_to = "compartment",
    values_to = "value"
  ) %>%
  mutate(compartment = sub("^log10_", "", compartment))

day14_immune_plot <- ggplot(df_long_day14, aes(x = time, y = value, color = compartment)) +
  geom_line(linewidth = 1) +
  labs(
    x = "Time (days)",
    y = "Log10 Cell Count",
    color = "Compartment"
  ) +
  theme_bw(base_size = 18) + 
  theme(legend.position = "bottom", panel.grid = element_blank())

ggsave(filename = "efficacy_immunity_day14.png",
       plot = day14_immune_plot,
       device = "png",
       height = 5, width = 10, units = "in")



