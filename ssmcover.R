#
# State-space modeling to estimate plant cover proportion
# from cover class data
#

library(readr)
library(dplyr)
library(extraDistr)
library(ggplot2)
library(rstan)
library(cmdstanr)
library(bayesplot)

# Function to convert cover proportion to class data
conv_class <- function(x, delta = 0.05,
                      cut_points = c(0.01, 0.1, 0.25, 0.5, 0.75)) {
  alpha <- x / delta - x
  beta <- (1 - x) * (1 - delta) / delta
  p <- pbeta(cut_points, alpha, beta)
  p <- c(p, 1)
  pr <- rep(0, length(cut_points) + 1)
  pr[1] <- p[1]
  pr[2:length(pr)] <- p[2:length(pr)] - p[1:(length(pr) - 1)]
  cls <- rcat(1, pr)
  return(cls)
}



# Simulated data
## Settings
Nt <- 15                    # time span
Nq <- 10                     # number of quadrats
cut_points <- c(0.01, 0.1, 0.25, 0.5, 0.75) # cut points
theta <- rep(NA, Nt)        # temporal change in the site
r <- rep(NA, Nq)            # spatial change in the site
cover <- matrix(NA, Nt, Nq) # matrix to store cover proportion data
cls <- matrix(NA, Nt, Nq)   # matrix to store cover class data
sd <- c(0.5, 0.5)           # SDs for spatial and temporal changes, respecitively
delta <- 0.05               # uncertainty in classification

## Data generation
### temporal change
set.seed(1)
theta[1] <- -6
for (t in 2:Nt) {
  theta[t] <- rnorm(1, theta[t - 1] + 0.3, sd[1])
}

### spatial variation
set.seed(2)
r[1] <- 0
for (q in 2:Nq) {
  r[q] <- rnorm(1, r[q - 1], sd[2])
}

### convert data from propotion to class
for (t in 1:Nt) {
  for (q in 1:Nq) {
    cover[t, q] <- 1 / (1 + exp(-theta[t] - r[q]))
    cls[t, q] <- conv_class(cover[t, q], delta, cut_points)
  }
}

## View data
df <- data.frame(Cover = factor(c(cls)),
                 Time = rep(1:Nt, Nq),
                 Quadrat = rep(1:Nq, each = Nt))
ggplot(df) +
  geom_tile(aes(x = Time, y = Quadrat, fill = Cover)) +
  scale_fill_discrete(h = c(0, 180) + 15) +
  scale_x_continuous(breaks = seq(0, 15, 5), minor_breaks = NULL) +
  scale_y_continuous(breaks = c(1, 5, 10), minor_breaks = NULL) +
  coord_fixed() +
  theme_bw(base_family = "Helvetica", base_size = 10)
ggsave("sim_data.pdf", device = "pdf", width = 12, height = 8, units = "cm")

## Bind data
stan_data <- list(N_q = Nq,
                  N_y = Nt,
                  N_cls = length(cut_points) + 1,
                  N_obs = Nt,
                  Obs_y = 1:Nt,
                  CP = cut_points,
                  Y = cls)

## CmdStanR settings
cmdstanpath <- "/usr/local/cmdstan"
set_cmdstan_path(cmdstanpath)
Sys.setenv(PATH = "/usr/bin:/bin")

## Initial values
inits <- c("inits/init_1.R", "inits/init_2.R",
           "inits/init_3.R", "inits/init_4.R")

## Compile and fitting
model_file <- "ssmcover.stan"
results_file <- "sim.RData"

if (!file.exists(results_file) |
    (file.mtime(model_file) > file.mtime(results_file))) {
  model <- cmdstan_model(model_file)
  
  fit_sim <- model$sample(data = stan_data,
                          seed = 1, refresh = 200, init = inits,
                          chains = 4, parallel_chains = 4,
                          iter_sampling = 1000, iter_warmup = 1000, thin = 1,
                          adapt_delta = 0.9, max_treedepth = 20)
  ###  Diagnose
  fit_sim$cmdstan_diagnose()
  
  ### Convert to stanfit object
  stanfit_sim <- rstan::read_stan_csv(fit_sim$output_files())
  save(stanfit_sim, file = results_file)
} else {
  load(results_file)
}

## Summary
print(stanfit_sim, pars = c("delta", "sigma"))
print(stanfit_sim, pars = "phi")

## Traceplot
rstan::traceplot(stanfit_sim, pars = "sigma")

## Posterior predictive check
yrep <-  extract(stanfit_sim, pars = "yrep")[["yrep"]]
for (i in 1:15) {
  print(pp_check(cls[i, ], yrep[ , i, ], ppc_rootogram))
}
pp_check(cls[15, ], yrep[ , 15, ], ppc_rootogram) +
  ggplot2::scale_x_discrete(name = "Class", limits = as.character(1:6)) +
  theme_classic(base_family = "Helvetica")
ggsave("sim_ppc.pdf", device = "pdf", width = 12, height = 9, units = "cm")

## View simulated data and posterior median and 95% CI
class_median <- c(cut_points, 1) / 2 + c(0, cut_points) / 2
phi <- rstan::extract(stanfit_sim, pars = "phi")[[1]]
phi.median <- apply(phi, 2, median)
phi.ci <- apply(phi, 2, quantile, probs = c(0.025, 0.975))
ggplot(data.frame(Time = 1:Nt,
                  Cover = apply(cover, 1, mean))) +
  geom_line(aes(x = Time, y = Cover), size = 1) +
  geom_jitter(data = data.frame(Time = rep(1:Nt, Nq),
                                Cover = class_median[c(cls)]),
              aes(x = Time, y = Cover),
              width = 0, height = 0.01,
              color = "black", alpha = 0.6) +
  geom_line(data = data.frame(Time = 1:Nt,
                              Cover = phi.median),
            aes(x = Time, y = Cover), size = 1, colour = "red") +
  geom_ribbon(data = data.frame(Time = 1:Nt,
                                Cover_lower = phi.ci[1, ],
                                Cover_upper = phi.ci[2, ]),
              aes(x = Time, ymin = Cover_lower, ymax = Cover_upper),
              fill = "red", alpha = 0.25) +
  theme_classic()
ggsave("sim.pdf", device = "pdf", width = 12, height = 9, units = "cm")


# Real data

## CmdStanR settings
cmdstanpath <- "/usr/local/cmdstan"
set_cmdstan_path(cmdstanpath)
Sys.setenv(PATH = "/usr/bin:/bin")
model_file <- "ssmcover.stan"
model <- cmdstan_model(model_file)

## Read data
data_file <- "Table S1.txt"
data <- read_tsv(data_file) %>%
  dplyr::mutate(Cov2 = dplyr::case_when(
    Cover == "+" ~ 1,
    Cover == "1" ~ 2,
    Cover == "2" ~ 3,
    Cover == "3" ~ 4,
    Cover == "4" ~ 5,
    Cover == "5" ~ 6))

## Plots (Belt No.)
plots <- levels(factor(data$Belt))

# Species and layer
species <- "Sasa senanensis"
layer <- "shrub"

# Cut points
cut_points <- c(0.01, 0.1, 0.25, 0.5, 0.75)

pl <- plots[1] # Choose plot

# Number of quadrats
n_q <- data %>%
  dplyr::filter(Belt == as.numeric(pl)) %>%
  dplyr::select(Quadrat) %>%
  max()
# Number of years
yrs <- data %>%
  dplyr::filter(Belt == as.numeric(pl)) %>%
  dplyr::select(Year) %>%
  unlist() %>%
  factor() %>%
  levels() %>%
  as.character() %>%
  as.integer()

ss <- data %>%
  dplyr::filter(Belt == as.numeric(pl) & Layer == layer &
                  Scientific_name == species)
y <- matrix(0, nrow = length(yrs), n_q)
for (i in seq_len(nrow(ss))) {
  d <- ss[i, ]
  y[match(d$Year, yrs), d$Quadrat] <- d$Cov2
}

# View data
df <- data.frame(Cover = factor(c(y)),
                 Year = rep(yrs, n_q),
                 Quadrat = rep(1:n_q, each = length(yrs)))
ggplot(df) +
  geom_tile(aes(x = Year, y = Quadrat, fill = Cover)) +
  scale_y_continuous(breaks = c(1, 5, 10, 15, 20, 25),
                     minor_breaks = NULL) +
  scale_fill_discrete(h = c(180, 0) + 15) +
  coord_fixed() +
  theme_bw(base_family = "Helvetica", base_size = 10)
ggsave(paste0("ss", pl, "_data.pdf"), device = "pdf",
       width = 15, height = 7.5, units = "cm")

## Fitting using Stan
stan_data <- list(N_q = n_q,
                  N_y = max(yrs) - min(yrs) + 1,
                  N_cls = length(cut_points) + 1,
                  N_obs = length(yrs),
                  Obs_y = yrs - min(yrs) + 1,
                  CP = cut_points,
                  Y = y)

## Initial values
inits <- sapply(1:4, function(i) paste0("inits/init", pl, "_", i, ".R"))

## Compile and fitting
results_file <- paste0("ss", pl, ".RData")
if (!file.exists(results_file) |
    (file.mtime(model_file) > file.mtime(results_file))) {
  # Note: For chain 3 of the plot 30, setting seed=3 does not work.
  #       setting seed=5 works.
  fit <- model$sample(data = stan_data,
                      seed = 1, init = inits,
                      chains = 4, parallel_chains = 4,
                      refresh = 200, 
                      iter_sampling = 2000, iter_warmup = 2000, thin = 1,
#                      refresh = 1,  # for test-run
#                      iter_sampling = 20, iter_warmup = 20, thin = 1,
                      adapt_delta = 0.85, max_treedepth = 20)
  ###  Diagnose
  fit$cmdstan_diagnose()
  
  ### Convert to stanfit object
  stanfit_ss <- rstan::read_stan_csv(fit$output_files())
  save(stanfit_ss, file = results_file)
} else {
  load(results_file)
}

## Summary
print(stanfit_ss, pars = c("delta", "sigma"))
print(stanfit_ss, pars = "phi")

## Posterior predictive check
y[y == 0] <- 1
yrep <-  extract(stanfit_ss, pars = "yrep")[["yrep"]]
for (i in seq_along(yrs)) {
  print(pp_check(y[i, ], yrep[ , i, ], ppc_rootogram))
}
pp_check(y[21, ], yrep[ , 21, ], ppc_rootogram) +
  ggplot2::scale_x_discrete(name = "Class", limits = as.character(1:6)) +
  theme_classic(base_family = "Helvetica") +
  theme(legend.position = "none")
ggsave(paste0("ss", pl, "_ppc.pdf"), device = "pdf",
       width = 12, height = 9, units = "cm")

## Traceplot
rstan::traceplot(stanfit_ss, pars = c("delta", "sigma"))

## View simulated data and posterior median and 95% CI
class_median <- c(cut_points, 1) / 2 + c(0, cut_points) / 2
class_median <- c(0, class_median) # include zero
phi <- rstan::extract(stanfit_ss, pars = "phi")[[1]]
phi.median <- apply(phi, 2, median)
phi.ci <- apply(phi, 2, quantile, probs = c(0.025, 0.975))
ggplot(data.frame(Year = rep(yrs, n_q),
                  Cover = class_median[c(y + 1)])) +
  geom_jitter(aes(x = Year, y = Cover),
              width = 0, height = 0.01,
              color = "black", alpha = 0.6) +
  geom_line(data = data.frame(Time = min(yrs):max(yrs),
                              Cover = phi.median),
            aes(x = Time, y = Cover), size = 1, colour = "red") +
  geom_ribbon(data = data.frame(Time = min(yrs):max(yrs),
                                Cover_lower = phi.ci[1, ],
                                Cover_upper = phi.ci[2, ]),
              aes(x = Time, ymin = Cover_lower, ymax = Cover_upper),
              fill = "red", alpha = 0.25) +
  xlim(1957, 2020) +
  theme_classic()
ggsave(paste0("ss", pl, ".pdf"), device = "pdf", width = 12, height = 9, units = "cm")


