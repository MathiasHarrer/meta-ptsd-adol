# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#                                                                             #
#   META-ANALYSIS OF INTERNET INTERVENTIONS FOR PTSD IN ADOLESCENTS           #
#   2023-05-22                                                                #
#                                                                             #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

library(rjags)
library(runjags)
library(dplyr)
library(tidyr)
library(writexl)
library(readxl)
library(metafor)
library(HDInterval)
library(ggplot2)
library(ggridges)
library(sjPlot)
library(forcats)
library(tidybayes)
library(glue)
source("utils.R")


# 1. Preparation --------------------------------------------------------------

# Import data
data = read_excel("data/data.xlsx")

# Define required class
within(data,{
       n_arm1 = as.numeric(n_arm1)
       n_arm2 = as.numeric(n_arm2)
       n_bl_arm1 = as.numeric(n_bl_arm1)
       sd_arm2 = as.numeric(sd_arm2)
       m_arm2 = as.numeric(m_arm2)
       m_bl_arm2 = as.numeric(m_bl_arm2)
    }) -> data

# Calculate effect sizes
# For within-group changes, we assume r=0.8
rho = 0.8
data %>% 
  mutate(
    sd_pooled = sqrt(((n_arm1-1L)*sd_arm1^2L + 
                        (n_arm2-1)*sd_arm2^2L)/(n_arm1 + n_arm2 - 2L)),
    smd_betw = (m_arm1-m_arm2)/sd_pooled,
    se_smd_betw =sqrt(((n_arm1 + n_arm2)/(n_arm1*n_arm2)) + 
                        (smd_betw^2L/(2L*(n_arm1+n_arm2)))),
    smd_within = (m_arm1 - m_bl_arm1)/sd_bl_arm1,
    se_smd_within = sqrt(((2L*(1L-rho))/n_bl_arm1) + 
                           (smd_within^2L/(2L*n_bl_arm1)))
  ) -> data


# Save ES values
writexl::write_xlsx(data, "data/data_es.xlsx")



# 2. Between-Group Effect Sizes ------------------------------------------------

# Filter out effect sizes
data %>% 
  mutate(yi = hedgesG(smd_betw, n_arm1 + n_arm2), vi = se_smd_betw) %>% 
  drop_na(yi, vi) -> dat.between


# Define JAGS model (standard IV REM with Half-Cauchy prior)
M = "
model{

  # Likelihood
  for (i in 1:length(yi)){

    # Hierarchical model
    theta[i] ~ dnorm(mu, prec.tau)
    yi[i] ~ dnorm(theta[i], prec.sigma[i])
    
    # Plug-in estimator of sigma
    prec.sigma[i] = pow(vi[i], -1)

  }
  
  # Prior
  mu ~ dnorm(0, 0.00001)
  prec.tau = 1/pow(tau,2)
  tau ~ dt(0, pow(0.2,-2), 1)T(0,)

  # Deterministic: define I-squared
  K = length(yi)
  w = pow(vi, -1)
  v.tilde = ((K-1)*sum(w))/(pow(sum(w),2)-sum(pow(w,2)))
  i2 = pow(tau,2)/(pow(tau,2) + v.tilde)
}
"

# Run model; get samples
params = c("mu", "tau", "i2", "theta")
dat = as.list(dat.between %>% dplyr::select(yi, vi))
fit = run.jags(M, monitor = params, data = dat, n.chains = 4, 
               sample = 90000, burnin = 1000, thin = 10)
samples = combine.mcmc(fit)


# Create results data.frame
data.frame(
  analysis = "Between-Group Effects",
  k = length(dat$yi),
  g = median(samples[,"mu"] %>% as.numeric()),
  g.se = sd(samples[,"mu"] %>% as.numeric()),
  g.lo = hdi(samples[,"mu"] %>% as.numeric())[,1],
  g.hi = hdi(samples[,"mu"] %>% as.numeric())[,2],
  i2 = median(samples[,"i2"] %>% as.numeric()),
  i2.lo = hdi(samples[,"i2"] %>% as.numeric())[1,1],
  i2.hi = hdi(samples[,"i2"] %>% as.numeric())[1,2],
  tau = median(samples[,"tau"] %>% as.numeric()),
  tau.lo = hdi(samples[,"tau"] %>% as.numeric())[,1],
  tau.hi = hdi(samples[,"tau"] %>% as.numeric())[,2]
) %>% 
  mutate(
    pi.lo = g - qt(.975, k-1)*sqrt(g.se^2+tau^2),
    pi.hi = g + qt(.975, k-1)*sqrt(g.se^2+tau^2)
  ) %>% 
  write_xlsx("results/results_between_group.xlsx")

# Save for plotting
samples.between = samples



# 3. Within-Group Effect Sizes ------------------------------------------------

# Filter out effect sizes, create cluster variable
data %>% 
  mutate(yi = hedgesG(smd_within, n_arm1),
         vi = se_smd_within, cluster = as.integer(as.factor(studlab))) %>% 
  drop_na(yi, vi) -> dat.within

# Order by cluster
dat.within[order(dat.within$cluster),] -> dat.within

# Write JAGS Model: Three-Level Hierarchical with double-Half-Cauchy
M = "
model{

  # Likelihood
  for (i in 1:kcluster){

    # Hierarchical model
    # Level 2 (Between-Studies)
    kappa[i] ~ dnorm(mu, prec.tau.l2)

    # Level 1 (Within Studies)
    for (j in slots[i,1]:slots[i,2]){
      theta[j] ~ dnorm(kappa[i], prec.tau.l1)
      yi[j] ~ dnorm(theta[j], prec.sigma[j])
      
      # Plug-in estimator of sigma
      prec.sigma[j] = pow(vi[j], -1)
    }
  }
  
  # Prior
  mu ~ dnorm(0, 0.00001)
  prec.tau.l2 = 1/pow(tau.l2,2)
  tau.l2 ~ dt(0, pow(0.1,-2), 1)T(0,)
  prec.tau.l1 = 1/pow(tau.l1,2)
  tau.l1 ~ dt(0, pow(0.1,-2), 1)T(0,)

  # Deterministic: define I-squared
  K = length(yi)
  w = pow(vi, -1)
  v.tilde = ((K-1)*sum(w))/(pow(sum(w),2)-sum(pow(w,2)))
  i2.l2 = pow(tau.l2,2)/(pow(tau.l2,2) + pow(tau.l1,2) + v.tilde)
  i2.l1 = pow(tau.l1,2)/(pow(tau.l2,2) + pow(tau.l1,2) + v.tilde)
}
"

# Define parameters
params = c("mu", "kappa", "theta", "tau.l1", "tau.l2", "i2.l1", "i2.l2")
dat = as.list(dat.within %>% dplyr::select(yi, vi))
cluster = as.integer(as.factor(dat.within$studlab))
dat$kcluster = length(unique(cluster))

# Create slots matrix to loop through within clusters
dat$slots = matrix(c(1,1,2,2,3,4,5,7,8,9), ncol=2, byrow=TRUE)

# Run MCMC simulation
fit = run.jags(M, monitor = params, data = dat,
               n.chains = 4, sample = 90000,
               burnin = 1000, thin = 10)

samples = combine.mcmc(fit)

# Create results data.frame
data.frame(
  analysis = "Within-Group Effects",
  k.effects = length(dat$yi),
  k = dat$kcluster,
  g = median(samples[,"mu"]),
  g.se = sd(samples[,"mu"]),
  g.lo = hdi(samples[,"mu"] %>% as.numeric())[1,1],
  g.hi = hdi(samples[,"mu"] %>% as.numeric())[1,2],
  i2.lvl1 = median(samples[,"i2.l1"] %>% as.numeric()),
  i2.lvl1.lo = hdi(samples[,"i2.l1"] %>% as.numeric())[1,1],
  i2.lvl1.hi = hdi(samples[,"i2.l1"] %>% as.numeric())[1,2],
  i2.lvl2 = median(samples[,"i2.l2"] %>% as.numeric()),
  i2.lvl2.lo = hdi(samples[,"i2.l2"] %>% as.numeric())[1,1],
  i2.lvl2.hi = hdi(samples[,"i2.l2"] %>% as.numeric())[1,2],
  tau.lvl1 = median(samples[,"tau.l1"] %>% as.numeric()),
  tau.lvl1.lo = hdi(samples[,"tau.l1"] %>% as.numeric())[1,1],
  tau.lvl1.hi = hdi(samples[,"tau.l1"] %>% as.numeric())[1,2],
  tau.lvl2 = median(samples[,"tau.l2"] %>% as.numeric()),
  tau.lvl2.lo = hdi(samples[,"tau.l2"] %>% as.numeric())[1,1],
  tau.lvl2.hi = hdi(samples[,"tau.l2"] %>% as.numeric())[1,2]
) %>% 
  mutate(
    pi.lo = g - qt(.975, k-1)*sqrt(g.se^2+tau.lvl2^2),
    pi.hi = g + qt(.975, k-1)*sqrt(g.se^2+tau.lvl2^2)
  ) %>% 
  write_xlsx("results/results_within_group.xlsx")


# Save for plotting
samples.within = samples



# 4. Forest Plot --------------------------------------------------------------


## 4.1 Between-Group Effect Sizes ---------------------------------------------

# Create posterior draws data set
samples.between %>% 
  as.data.frame() %>% 
  select(mu, `theta[1]`:`theta[3]`) %>% 
  pivot_longer(everything()) %>% 
  mutate(name = as.factor(name),
         name = recode(name, 
                       "theta[1]" = "Cox et al. (2009)",
                       "theta[2]" = "Kassam-Adams et al. (2016)",
                       "theta[3]" = "Ruggiero et al. (2015)",
                       "mu" = "Overall Effect"),
         name = fct_relevel(name, "Overall Effect")) -> dat.plt

# Define labeling order
dat.plt$name = factor(dat.plt$name,
                      levels = c("Cox et al. (2009)", 
                                 "Kassam-Adams et al. (2016)", 
                                 "Ruggiero et al. (2015)", 
                                 "Overall Effect") %>% rev())

# Generate summary data
dat.plt.sum = samples.between %>% 
  as.data.frame() %>% 
  select(mu, `theta[1]`:`theta[3]`) %>% 
  pivot_longer(everything()) %>% 
  group_by(name) %>% 
  summarise(.est = median(value),
            .ci = hdi(value)) %>% 
  as.data.frame() %>% 
  {data.frame(name = .$name, .est = .$.est,
              .lower = .[,3][,1], .upper = .[,3][,2])} %>% 
  mutate(name = as.factor(name),
         name = recode(name, 
                       "theta[1]" = "Cox et al. (2009)",
                       "theta[2]" = "Kassam-Adams et al. (2016)",
                       "theta[3]" = "Ruggiero et al. (2015)",
                       "mu" = "Overall Effect"),
         name = fct_relevel(name, "Overall Effect"))

# Change label order for summary data
dat.plt.sum$name = factor(dat.plt.sum$name,
                      levels = c("Cox et al. (2009)", 
                                 "Kassam-Adams et al. (2016)", 
                                 "Ruggiero et al. (2015)", 
                                 "Overall Effect") %>% rev())

# Define colors
colors = c("lightskyblue", "cornflowerblue", "lightblue", "gray20") %>% rev()
colors2 = c("dodgerblue", "dodgerblue", "dodgerblue", "gray20") %>% rev()

# Start plotting
ggplot(dat.plt, aes(value, name, fill = name)) + 
  
  # Add vertical lines for pooled effect and CI
  geom_vline(xintercept = median(samples.between[,"mu"]), 
             color = "grey", linewidth = 1) +
  geom_vline(xintercept = c(hdi(as.numeric(samples.between[,"mu"]))[,1], 
                            hdi(as.numeric(samples.between[,"mu"]))[,2]), 
             color = "grey", linetype = 2) +
  geom_vline(xintercept = 0, color = "black", 
             size = 1) +
  
  # Add densities
  geom_density_ridges(color = "darkgrey", rel_min_height = 0.05, 
                      col = NA, scale = 0.9, alpha = 0.5) +
  geom_pointinterval(aes(xmin = .lower, xmax = .upper, x = .est, y = name),
                     data = dat.plt.sum, size = 1) + 
  
  # Add text and labels
  geom_label(data = mutate_if(dat.plt.sum, is.numeric, ~sprintf("%.2f", .)),
             aes(label = glue("{.est} [{.lower}; {.upper}]"), 
                 x = Inf), hjust = "inward", fill = "white", 
             label.size = NA, label.padding = unit(0.5, "lines"),
             family = "Fira Sans") +
  
  # Theming
  scale_fill_manual(values = colors2) + 
  xlab("Standardized Mean Difference (Hedges' g)") +
  ylab("") +  xlim(c(-1, 2)) + theme_minimal() +
  theme(legend.position = "none", 
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        text=element_text(family="Fira Sans")) -> plt

# Save plot
ggsave(plt, file="results/plot.between.png", bg = "white",
       height = 3.5)



## 4.2 Within-Group Effect Sizes ----------------------------------------------

# Create posterior draws data set
samples.within %>% 
  as.data.frame() %>% 
  select(mu, `theta[1]`:`theta[9]`) %>% 
  pivot_longer(everything()) %>% 
  mutate(name = as.factor(name),
         name = recode(name, 
                       "theta[1]" = "Cox et al. (2009) (IMI)",
                       "theta[2]" = "Kassam-Adams et al. (2016) (IMI)",
                       "theta[3]" = "Ruggiero et al. (2015) (IMI)",
                       "theta[4]" = "Ruggiero et al. (2015) (IMI + ASH)",
                       "theta[5]" = "Schuurmans et al. (2022) (IMI Muse)",
                       "theta[6]" = "Schuurmans et al. (2022) (IMI Daydream)",
                       "theta[7]" = "Schuurmans et al. (2022) (IMI Wild divine)",
                       "theta[8]" = "van Rosmalen-Nooijens et al. (2017) (IMI)",
                       "theta[9]" = "van Rosmalen-Nooijens et al. (2017) (IMI for controls)",
                       "mu" = "Overall Effect"),
         name = fct_relevel(name, "Overall Effect")) -> dat.plt

# Define labeling order
dat.plt$name = factor(
  dat.plt$name,
  levels = c("Cox et al. (2009) (IMI)",
             "Kassam-Adams et al. (2016) (IMI)",
             "Ruggiero et al. (2015) (IMI)",
             "Ruggiero et al. (2015) (IMI + ASH)",
             "Schuurmans et al. (2022) (IMI Muse)",
             "Schuurmans et al. (2022) (IMI Daydream)",
             "Schuurmans et al. (2022) (IMI Wild divine)",
             "van Rosmalen-Nooijens et al. (2017) (IMI)",
             "van Rosmalen-Nooijens et al. (2017) (IMI for controls)",
             "Overall Effect") %>% rev())

# Generate summary data
dat.plt.sum = samples.within %>% 
  as.data.frame() %>% 
  select(mu, `theta[1]`:`theta[9]`) %>% 
  pivot_longer(everything()) %>% 
  group_by(name) %>% 
  summarise(.est = median(value),
            .ci = hdi(value)) %>% 
  as.data.frame() %>% 
  {data.frame(name = .$name, .est = .$.est,
              .lower = .[,3][,1], .upper = .[,3][,2])} %>% 
  mutate(name = as.factor(name),
         name = recode(name, 
                       "theta[1]" = "Cox et al. (2009) (IMI)",
                       "theta[2]" = "Kassam-Adams et al. (2016) (IMI)",
                       "theta[3]" = "Ruggiero et al. (2015) (IMI)",
                       "theta[4]" = "Ruggiero et al. (2015) (IMI + ASH)",
                       "theta[5]" = "Schuurmans et al. (2022) (IMI Muse)",
                       "theta[6]" = "Schuurmans et al. (2022) (IMI Daydream)",
                       "theta[7]" = "Schuurmans et al. (2022) (IMI Wild divine)",
                       "theta[8]" = "van Rosmalen-Nooijens et al. (2017) (IMI)",
                       "theta[9]" = "van Rosmalen-Nooijens et al. (2017) (IMI for controls)",
                       "mu" = "Overall Effect"),
         name = fct_relevel(name, "Overall Effect"))
  
# Change label order for summary data
dat.plt.sum$name = factor(
  dat.plt.sum$name,
  levels = c("Cox et al. (2009) (IMI)",
             "Kassam-Adams et al. (2016) (IMI)",
             "Ruggiero et al. (2015) (IMI)",
             "Ruggiero et al. (2015) (IMI + ASH)",
             "Schuurmans et al. (2022) (IMI Muse)",
             "Schuurmans et al. (2022) (IMI Daydream)",
             "Schuurmans et al. (2022) (IMI Wild divine)",
             "van Rosmalen-Nooijens et al. (2017) (IMI)",
             "van Rosmalen-Nooijens et al. (2017) (IMI for controls)",
             "Overall Effect") %>% rev())

# Define colors
colors = c("lightskyblue", "steelblue",
           "cornflowerblue", "lightblue", "lightblue",
           "dodgerblue", "dodgerblue", "dodgerblue",
           "blue", "blue", "gray20") %>% rev()

# Start plotting
ggplot(dat.plt, aes(value, name, fill = name)) + 
  
  # Add vertical lines for pooled effect and CI
  geom_vline(xintercept = median(samples.within[,"mu"]), 
             color = "grey", linewidth = 1) +
  geom_vline(xintercept = c(hdi(as.numeric(samples.within[,"mu"]))[,1], 
                            hdi(as.numeric(samples.within[,"mu"]))[,2]), 
             color = "grey", linetype = 2) +
  geom_vline(xintercept = 0, color = "black", 
             size = 1) +
  
  # Add densities
  geom_density_ridges(color = "darkgrey", rel_min_height = 0.05, 
                      col = NA, scale = 0.9, alpha = 0.5) +
  geom_pointinterval(aes(xmin = .lower, xmax = .upper, x = .est, y = name),
                     data = dat.plt.sum, size = 1) + 
  
  # Add text and labels
  geom_label(data = mutate_if(dat.plt.sum, is.numeric, ~sprintf("%.2f", .)),
            aes(label = glue("{.est} [{.lower}; {.upper}]"), 
                x = Inf), hjust = "inward", fill = "white", 
            label.size = NA, label.padding = unit(0.5, "lines"),
            family = "Fira Sans") +
  
  # Theming
  scale_fill_manual(values = colors) + 
  xlab("Standardized Mean Difference (Hedges' g)") +
  ylab("") +  xlim(c(-1, 1)) + theme_minimal() +
  theme(legend.position = "none", 
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        text=element_text(family="Fira Sans")) -> plt

# Save plot
ggsave(plt, file="results/plot.within.png", bg = "white")



# 5. Sensitivity Analysis ------------------------------------------------------

# van Rosmalen et al. showed very strong baseline imbalance and was therefore
# excluded in the between-group effect analysis. We now conduct a sensitivity
# analysis in which this effect is included.






