# Nhat Hoang Pham MATH 7397 
# Project Proposal
library(haven)
library(ggplot2)
library(rstan)
library(rstanarm)
options(mc.cores = parallel::detectCores())
library(loo)
library(GGally)
library(bayesplot)
theme_set(bayesplot::theme_default())
library(projpred)
library(devtools)
library(BayesVarSel)
library(car)
library(coda)
library(bayesplot)
library(bayesutils)

#Data
set.seed(90)
master_data <- read_dta("Desktop/Study/CU_Denver/Spring_2023/MATH 7397/Final Project/master_data.dta")



#Data manipulation
df = master_data[,c("year","DL","VM","IN","age","pop", "unemployed","decriminalized", "MML_1", "texting_ban","VM_per_driver","fatals_total_rate")]
dim(df)
X = df[,c("year","DL","VM","IN","age","pop", "unemployed","decriminalized", "MML_1", "texting_ban","VM_per_driver")]
y_1 = df[,"fatals_total_rate"]
#Income 1k
df$IN = df$IN /1000

# Create Densityplot for Bayesian analysis
plot_title <- "Density Plot of Fatality Rate"
x_axis_label <- "Traffic Fatality Rate"
density_plot <- ggplot(df, aes(x = fatals_total_rate)) +
  geom_density(fill = "lightblue", alpha = 0.7) +
  ggtitle(plot_title) +
  xlab(x_axis_label) +
  ylab("Density")
# Graph Summary data:
library(panelView)
panelview(fatals_total_rate ~ MML_1, data = master_data[master_data$state%in%c(42,36,17,27,33,29,39), ],  index = c("state","year"), type = "outcome")
density_plot


#Data summary 

# Find mean of the data frame
mean_df <- apply(df, 2, mean)
# This applies the `mean()` function to each column (2nd dimension) of the data frame, 
# resulting in a vector of means for each column.
mean_df
# Find standard deviation of the data frame
sd_df <- apply(df, 2, sd)
# This applies the `sd()` function to each column (2nd dimension) of the data frame, 
# resulting in a vector of standard deviations for each column.
sd_df


# Collinearity 
ggpairs(df,diag = list(continuous = "barDiag"))


#Fit a general model with all variables 
stan_glm1 = stan_glm(fatals_total_rate ~ MML_1+year+DL+VM+IN+age+pop+unemployed+decriminalized+texting_ban+VM_per_driver,
                     family = gaussian(),
                     prior = normal(0.0,10),
                     prior_intercept = normal(0.0,10),
                     data = df,
                     QR = TRUE,
                     iter = 50000,chains = 4)
vif(stan_glm1)

# Assessing Collinearity
ggpairs(df,diag = list(continuous = "barDiag"))
post_1 <- as.matrix(stan_glm1)
post_1 <- post_1[,-1]
mcmc_areas(post_1,prob_outer = 0.95)
mcmc_pairs(post_1,pars = c("DL","age","pop","VM"))


# Perform variable selection  using total traffic fatalities
df.bvs_1 <- Bvs(fatals_total_rate ~ MML_1+year+IN+age+unemployed+decriminalized+texting_ban+VM_per_driver, data = df)
df.bvs_1
summary(df.bvs_1)
model_fitcv <- cv_varsel(stan_glm1,method = "L1", cv_metthod = "kfold", k_fold = 10)
# Print the selected variables



# Using (Vehtari, Gelman, and Gabry, 2017) to choose model size take too long to run
#model_fitcv <- cv_varsel(stan_glm1,method = "L1", cv_metthod = "kfold", k_fold = 10)

#plot(model_fitcv, stats = c("elpd",'rmse'))


 
stan_glm_final_1 = stan_glm(fatals_total_rate ~ MML_1+year+IN+decriminalized+texting_ban+VM_per_driver+age + unemployed,
                     family = gaussian(link = "identity"),
                     prior = normal(0.0,10),
                     prior_intercept = normal(0.0,10),
                     #prior_aux = Gamma(link = 'inverse'),
                     data = df,
                     QR = TRUE, 
                     iter = 50000,chains = 4)


sum_final_1<-summary(stan_glm_final_1,digits = 5)


waic_1 <- waic(stan_glm_final_1)
looic_1 <- loo(stan_glm_final_1)

stan_glm_final_2 = stan_glm(fatals_total_rate ~ MML_1+year+IN+decriminalized+VM_per_driver+age+texting_ban,
                            family = gaussian(link = "identity"),
                            prior = normal(0.0,10),
                            prior_intercept = normal(0.0,10),
                            #prior_aux = Gamma(link = 'inverse'),
                            data = df,
                            QR = TRUE, 
                            iter = 50000,chains = 4)
sum_final_2<-summary(stan_glm_final_2,digits = 5)


waic_2 <- waic(stan_glm_final_2)
looic_2 <- loo(stan_glm_final_2)



stan_glm_final_3 = stan_glm(fatals_total_rate ~ MML_1+year+IN+decriminalized+VM_per_driver+age+unemployed,
                            family = gaussian(link = "identity"),
                            prior = normal(0.0,10),
                            prior_intercept = normal(0.0,10),
                            #prior_aux = Gamma(link = 'inverse'),
                            data = df,
                            QR = TRUE, 
                            iter = 50000,chains = 4)

sum_final_3<-summary(stan_glm_final_3,digits = 5)


waic_3<- waic(stan_glm_final_3)
looic_3 <- loo(stan_glm_final_3)

stan_glm_final_4 = stan_glm(fatals_total_rate ~ MML_1+year+IN+decriminalized+VM_per_driver+age,
                            family = gaussian(link = "identity"),
                            prior = normal(0.0,10),
                            prior_intercept = normal(0.0,10),
                            #prior_aux = Gamma(link = 'inverse'),
                            data = df,
                            QR = TRUE, 
                            iter = 50000,chains = 4)

sum_final_4<-summary(stan_glm_final_4,digits = 5)
waic_4<- waic(stan_glm_final_4)
looic_4 <- loo(stan_glm_final_4)




# LOOIC and WAIC comparison:
loo_compare(waic_1,waic_2,waic_3,waic_4)
loo_compare(looic_1,looic_2,looic_3,looic_4)   




# Assessing convergence





df$log_fatals_total_rate <- log(df$fatals_total_rate)

stan_glm_final_5 = stan_glm(log_fatals_total_rate ~ MML_1+year+IN+decriminalized+VM_per_driver+age,
                            family = gaussian(link = "identity"),
                            prior = normal(0.0,10),
                            prior_intercept = normal(0.0,10),
                            #prior_aux = Gamma(link = 'inverse'),
                            data = df,
                            QR = TRUE, 
                            iter = 50000,chains = 4)


sum_final_5<-summary(stan_glm_final_5,digits = 5)
summary_final <- as.data.frame(summary(stan_glm_final_5, probs = c(0.025, 0.975)))
color_scheme_set("brightblue") # see help("color_scheme_set")
summary_final[,c(1, 3, 4,5)]
#Assesing convergence 
mcmc_rhat(rhat(stan_glm_final_5))

mcmc_neff(neff_ratio(stan_glm_final_5), size = 2) # worry if less than 0.1

post_final<- as.matrix(stan_glm_final_5)

heidel.diag(post_final)

raftery.diag(post_final)

geweke.diag(post_final)

post_final <- post_final[,-1]

mcmc_intervals(post_final,prob_outer = 0.99)


mcmc_trace(post_final)



#Outlier:

# Combine variable selection and assessing collinearity check, we have : 
#fatals_total_rate ~ MML_1+year+IN+decriminalized+texting_ban+VM_per_driver


# Model checking
library(car)

# compute leave-one-out predictive checks
B = 100000 # number of simulations

# length, mean, and sd of data
n = length(df$log_fatals_total_rate)
m = mean(df$log_fatals_total_rate)
s = sd(df$log_fatals_total_rate)

### posterior interval for mu (simulation)
# sample from posterior of sigmasq
sigmasq = rinvchisq(B, df = n - 1, scale = s^2)
# sample from conditional posterior of mu
mu = rnorm(B, mean = m, sd = sqrt(sigmasq/n))

# store nyrep samples of yrep
# from from posterior predictive distribution
nyrep = 200
yrep = matrix(0, nrow = nyrep, ncol = n)

# rename for convenience
y = as.vector(df$log_fatals_total_rate)
# sample 510 observations from posterior predictive distribution
for (i in seq_len(nyrep)) {
  yrep[i, ] = rnorm(n, mean = mu[i], sd = sqrt(sigmasq[i]))
}

# compare minimum of observed data to minimum from
# replicated samples

# minimum of replicated sample
mins = apply(yrep, 1, min)
# estimated p-value
(sum(mins <= min(y)) + 1)/(length(mins) + 1)
# histogram comparing T(y) and T(yrep)
ppc_stat(y, yrep, stat = "min")


# compare histograms of y and yrep
ppc_hist(y, yrep[sample(nyrep, 8), ])

ppc_boxplot(y, yrep[sample(nyrep, 8), ])

# compare densities of y and yrep
ppc_dens_overlay(y, yrep[sample(nyrep, 20), ])
# compare ecdfs of y and yrep
ppc_ecdf_overlay(y, yrep[sample(nyrep, 20) , ])
# compare histograms of y - yrep
ppc_error_hist(y, yrep[sample(nyrep, 9) , ])
# compare scatterplots of y vs yrep
ppc_scatter(y, yrep[sample(nyrep, 9) , ])
# compare scatterplots of y vs y - yrep
ppc_error_scatter(y, yrep[sample(nyrep, 9) , ])

# marginal predictive checks
# comparison of observed data and 90% predictive cpis
ppc_intervals(y, yrep)


# Outlier check
outlier_check <-lm(log_fatals_total_rate ~ MML_1+year+IN+decriminalized+VM_per_driver+age,data = df)
influencePlot(outlier_check)
library(BAS)
outliers <- Bayes.outlier(outlier_check, k = 3)
plot(outliers$prob.outlier, type = "h", ylab = "Posterior Probability")





nyrep = 10000
yrep = matrix(0, nrow = nyrep, ncol = n)
loglik_yrep = matrix(0, nrow = nyrep, ncol = n)
for (i in seq_len(nyrep)) {
  yrep[i, ] = rnorm(n, mean = mu[i], sd = sqrt(sigmasq[i]))
  loglik_yrep[i,] = dnorm(yrep[i, ],
                          mean = mu[i],
                          sd = sqrt(sigmasq[i]),
                          log = TRUE)
}

# compute relative effective MCMC sample size
# divided by the total sample size for likelihood
r_eff = relative_eff(exp(loglik_yrep),
                     chain_id = rep(1, nyrep))
# compute leave-one-out information
loo_info = loo(loglik_yrep, r_eff = r_eff, save_psis = TRUE)

# construct leave-one-out prediction intervals
ppc_loo_intervals(y, yrep,
                  psis_object = loo_info$psis_object)

# construct leave-one-out quantiles
ppc_loo_pit_qq(y, yrep,
               lw = weights(loo_info$psis_object))

# ppo vs y
PPO = colMeans(exp(loglik_yrep))
plot(PPO ~ y, ylab = "PPO")
dy = density(y)
# scale density to match scale of PPO
dy$y = dy$y/max(dy$y)*(max(PPO) - min(PPO)) + min(PPO)
lines(dy)



pp_check(stan_glm_final_5, plotfun = "stat", stat = "mean")
pp_check(stan_glm_final_5, plotfun = "stat_2d", stat = c("mean", "sd"))
# Create graph

ppc_dens_overlay(y = as.vector(stan_glm_final_5$y),
                 yrep = posterior_predict(stan_glm_final_5, draws = 200))


