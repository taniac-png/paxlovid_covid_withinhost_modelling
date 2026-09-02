#### ECLIPSE MODEL WITH ANTIVIRAL EFFECT, INNATE IMMUNITY AND ADAPTIVE IMMUNITY ####

library(pacman)
p_load(deSolve, ggplot2, tidyverse, coda, future, future.apply, envalysis)

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

## Parameters
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

## Solve and plot
out <- ode(y = y0_SEIVRY, times = times_out, func = SEIVRY_model, parms = pars_SEIVRY, events = list(data = drug_regime))

out_df <- as.data.frame(out)
View(out_df)

## Mean efficacy
mean_efficacy <- mean(out_df$epsilon)
print(paste("Mean efficacy over time:", round(mean_efficacy, 4)))

# Visualise all compartments and values over time
out_long <- as.data.frame(out) %>% pivot_longer(-time)
ggplot(out_long) + geom_line(aes(x=time, y=value, color=name)) + scale_y_log10() + coord_cartesian(ylim=c(1,1e10))

compartments <- c("S","E","I","V","R","Y","GI","C","P","Rtv", "conc", "epsilon") # to plot per compartment
out_long2 <- out_long %>%
  filter(name %in% compartments) %>%
  mutate(value = ifelse(value <= 0, NA_real_, value))  # to avoid log errors

ggplot(out_long2, aes(x = time, y = value)) +
  geom_line() +
  facet_wrap(~name, scales = "free_y", ncol = 3) +
  scale_y_log10() +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white"))

# Log transform
eps <- 1e-8  # small constant to avoid log(0)

out_log <- out_df %>%
  mutate(across(c(S, E, I, V, R, Y, GI, C, P, conc), 
                ~ log10(.x + eps),
                .names = "log10_{.col}")
        )

out_log_long <- as.data.frame(out_log) %>% pivot_longer(-time)

# plot viral load over time
viral_load <- ggplot(out_log_long %>% filter(name %in% c("log10_V")), aes(x = time, y = value, color = name)) +
  geom_line(linewidth = 1) +
  labs(
    title = "Treatment", x = "Time (days)", y = "Log10 Viral Load (copies/mL)", color = "Compartment"
  ) +
  theme_bw(base_size = 20) + 
  theme(legend.position = "none", 
        panel.grid = element_blank(), 
        plot.title = element_text(
            hjust = 0.5, 
            face = "bold"))

ggsave(filename = "viral_load_plot.png",
       plot = viral_load,
       device = "png",
       height = 5, width = 10, units = "in")

# plot drug concentration in central compartment over time
drug_conc <- ggplot(out_log_long %>% filter(name %in% c("log10_conc")), aes(x = time, y = value, color = name)) +
  geom_line(linewidth = 1) +
  labs(
    title = "Antiviral concentration", x = "Time (days)", y = "Log10 Drug Concentration (mg/L)", color = "Compartment"
  ) +
  theme_bw(base_size = 20) + 
  theme(legend.position = "none", 
        panel.grid = element_blank(), 
        plot.title = element_text(
            hjust = 0.5,  
            face = "bold"))

ggsave(filename = "drug_conc_plot.png",
       plot = drug_conc,
       device = "png",
       height = 5, width = 10, units = "in")


# plot drug efficacy over time
drug_efficacy <- ggplot(out_log_long %>% filter(name %in% c("epsilon")), aes(x = time, y = value, color = name)) +
  geom_line(linewidth = 1) +
  labs(
    title = "Antiviral efficacy", x = "Time (days)", y = "Drug Efficacy", color = "Compartment"
  ) +
  theme_bw(base_size = 20) + 
  theme(legend.position = "none", 
        panel.grid = element_blank(), 
        plot.title = element_text(
            hjust = 0.5,  
            face = "bold"))

ggsave(filename = "drug_efficacy_plot.png",
       plot = drug_efficacy, 
       device = "png",
       height = 5, width = 10, units = "in")

# plot viral load for untreated over time
viral_load_notreat <- ggplot(out_log_long %>% filter(name %in% c("log10_V")), aes(x = time, y = value, color = name)) +
  geom_line(linewidth = 1) +
  labs(
    title = "No treatment", x = "Time (days)", y = "Log10 Viral Load (copies/mL)", color = "Compartment"
  ) +
  theme_bw(base_size = 20) + 
  theme(legend.position = "none", 
        panel.grid = element_blank(), 
        plot.title = element_text(
            hjust = 0.5,  
            face = "bold"))

ggsave(filename = "viral_load_notreat_plot.png",
       plot = viral_load_notreat,
       device = "png",
       height = 5, width = 10, units = "in")


