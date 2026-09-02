#### FOR SIMULATIONS ALTERING EFFICACY ####

library(pacman)
p_load(deSolve, ggplot2, tidyverse, coda, future)

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
pars_moderate <- c(
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

## Antiviral drug treatment regime
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
out_moderate <- ode(y = y0_SEIVRY, times = times_out, func = SEIVRY_model, parms = pars_moderate, events = list(data = drug_regime))
out_moderate_df <- as.data.frame(out_moderate)

## Mean efficacy
mean_efficacy_moderate <- mean(out_moderate_df$epsilon)
print(paste("Mean efficacy over time:", round(mean_efficacy_moderate, 4)))


### LOW EFFICACY SIMULATION ###

## Parameters - low efficacy (IC50 = 0.01)
pars_low <- c(
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
  IC50 = 0.01, 
  Hill_coeff = 4
)

## Solve
out_low <- ode(y = y0_SEIVRY, times = times_out, func = SEIVRY_model, parms = pars_low, events = list(data = drug_regime))
out_low_df <- as.data.frame(out_low)

## Mean efficacy
mean_efficacy_low <- mean(out_low_df$epsilon)
print(paste("Mean efficacy over time:", round(mean_efficacy_low, 4)))


### HIGH EFFICACY SIMULATION ###

## Parameters - high efficacy (IC50 = 0.0005)
pars_high <- c(
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
  IC50 = 0.0005, 
  Hill_coeff = 4
)

## Solve
out_high <- ode(y = y0_SEIVRY, times = times_out, func = SEIVRY_model, parms = pars_high, events = list(data = drug_regime))
out_high_df <- as.data.frame(out_high)
View(out_high_df)

## Mean efficacy
mean_efficacy_high <- mean(out_high_df$epsilon)
print(paste("Mean efficacy over time:", round(mean_efficacy_high, 4)))


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

out_moderate_log <- out_moderate_df %>%
  mutate(across(c(S, E, I, V, R, Y, GI, C, P, conc), 
                ~ log10(.x + eps),
                .names = "log10_{.col}")
        )

out_low_log <- out_low_df %>%
  mutate(across(c(S, E, I, V, R, Y, GI, C, P, conc), 
                ~ log10(.x + eps),
                .names = "log10_{.col}")
        )

out_high_log <- out_high_df %>%
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
print(max(out_moderate_log$log10_V))
print(max(out_low_log$log10_V))
print(max(out_high_log$log10_V))
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


### COMBINE AND PLOT FOR VIRAL LOAD ###
moderate_V <- out_moderate_log %>%
  select(time, log10_V) %>%
  mutate(group = "Moderate")

low_V <- out_low_log %>%
  select(time, log10_V) %>%
  mutate(group = "Low")

high_V <- out_high_log %>%
  select(time, log10_V) %>%
  mutate(group = "High")

notreat_V <- out_notreat_log %>%
  select(time, log10_V) %>%
  mutate(group = "No treatment")

# Combine and plot
combined_V <- bind_rows(
  moderate_V,
  low_V,
  high_V,
  notreat_V
)

combined_V_filtered <- combined_V %>%
  filter(log10_V > 0)   # keeps V > 1

efficacy_plot <- ggplot(combined_V_filtered, aes(x = time, y = log10_V, color = group)) +
  geom_line(linewidth = 1) +
  labs(
    x = "Time (days)",
    y = "Log10 Viral Load (copies/mL)",
    color = "Scenario"
  ) +
  theme_bw(base_size = 18) + 
  theme(legend.position = "bottom", panel.grid = element_blank())

ggsave(filename = "efficacy_simulations.png",
       plot = efficacy_plot,
       device = "png",
       height = 5, width = 10, units = "in")


### COMBINE AND PLOT FOR IMMUNE RESPONSE - MODERATE ###
moderate_immune <- out_moderate_log %>%
  select(time, log10_S, log10_I, log10_R, log10_Y) %>%
  mutate(group = "Moderate")

df_long_moderate <- moderate_immune %>%
  pivot_longer(
    cols = c(log10_S, log10_I, log10_R, log10_Y),
    names_to = "compartment",
    values_to = "value"
  ) %>%
  mutate(compartment = sub("^log10_", "", compartment))

moderate_immune_plot <-ggplot(df_long_moderate, aes(x = time, y = value, color = compartment)) +
  geom_line(linewidth = 1) +
  labs(
    x = "Time (days)",
    y = "Log10 Cell Count",
    color = "Compartment"
  ) +
  theme_bw(base_size = 18) + 
  theme(legend.position = "bottom", panel.grid = element_blank())

ggsave(filename = "efficacy_immunity_moderate.png",
       plot = moderate_immune_plot,
       device = "png",
       height = 5, width = 10, units = "in")


### COMBINE AND PLOT FOR IMMUNE RESPONSE - LOW ###
low_immune <- out_low_log %>%
  select(time, log10_S, log10_I, log10_R, log10_Y) %>%
  mutate(group = "Low")

df_long_low <- low_immune %>%
  pivot_longer(
    cols = c(log10_S, log10_I, log10_R, log10_Y),
    names_to = "compartment",
    values_to = "value"
  ) %>%
  mutate(compartment = sub("^log10_", "", compartment))

low_immune_plot <- ggplot(df_long_low, aes(x = time, y = value, color = compartment)) +
  geom_line(linewidth = 1) +
  labs(
    x = "Time (days)",
    y = "Log10 Cell Count",
    color = "Compartment"
  ) +
  theme_bw(base_size = 18) + 
  theme(legend.position = "bottom", panel.grid = element_blank())

ggsave(filename = "efficacy_immunity_low.png",
       plot = low_immune_plot,
       device = "png",
       height = 5, width = 10, units = "in")


### COMBINE AND PLOT FOR IMMUNE RESPONSE - HIGH ###
high_immune <- out_high_log %>%
  select(time, log10_S, log10_I, log10_R, log10_Y) %>%
  mutate(group = "High")

df_long_high <- high_immune %>%
  pivot_longer(
    cols = c(log10_S, log10_I, log10_R, log10_Y),
    names_to = "compartment",
    values_to = "value"
  ) %>%
  mutate(compartment = sub("^log10_", "", compartment))

high_immune_plot <- ggplot(df_long_high, aes(x = time, y = value, color = compartment)) +
  geom_line(linewidth = 1) +
  labs(
    x = "Time (days)",
    y = "Log10 Cell Count",
    color = "Compartment"
  ) +
  theme_bw(base_size = 18) + 
  theme(legend.position = "bottom", panel.grid = element_blank())

ggsave(filename = "efficacy_immunity_high.png",
       plot = high_immune_plot,
       device = "png",
       height = 5, width = 10, units = "in")

### COMBINE AND PLOT FOR IMMUNE RESPONSE - NO TREATMENT ###
notreat_immune <- out_notreat_log %>%
  select(time, log10_S, log10_I, log10_R, log10_Y) %>%
  mutate(group = "No treatment")

df_long_notreat <- notreat_immune %>%
  pivot_longer(
    cols = c(log10_S, log10_I, log10_R, log10_Y),
    names_to = "compartment",
    values_to = "value"
  ) %>%
  mutate(compartment = sub("^log10_", "", compartment))

notreat_immune_plot <- ggplot(df_long_notreat, aes(x = time, y = value, color = compartment)) +
  geom_line(linewidth = 1) +
  labs(
    x = "Time (days)",
    y = "Log10 Cell Count",
    color = "Compartment"
  ) +
  theme_bw(base_size = 18) + 
  theme(legend.position = "bottom", panel.grid = element_blank())

ggsave(filename = "efficacy_immunity_notreat.png",
       plot = notreat_immune_plot,
       device = "png",
       height = 5, width = 10, units = "in")

