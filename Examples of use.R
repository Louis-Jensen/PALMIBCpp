# Packages ----
# --------------------------------------------------
# Please install the following packages.
library(tidyverse)
library(pracma)
library(spatstat)
library(markovchain)
library(GET)
#library(magrittr)
#library(wrapr)
#library(zeallot)
#library(janitor)
#library(purrr)

# Loading in the dependencies ----
# --------------------------------------------------
# let 'root' be the full path to the IBCpp folder - e.g. root = "C:/Users/Me/Documents/PALMIBCpp/"
#root = "C:/Users/louis/Documents/PhD/Produced Work/Articles/2019 - A spatio-temporal point-process model for correction of photoblinking artifacts/PALMIBCpp/"
root = ""
ffs = c("Estimation.R", "Gamma1.R", "Misc.R", "Space-Time Simulation.R", "Temporal Simulation.R")
for(f in ffs) {
  fn = paste(root, "Depend/", f, sep = "")
  source(fn)
}
rm(ffs, fn, f)


# INTRODUCTION ----
# --------------------------------------------------
# We present here some ways the PALM-IBCpp model can be used to model PALM data.
# EXAMPLE 1 shows basic usage for estimation of kinetic rates, and obtaining various derived statistics.
# EXAMPLE 2 shows how to perform the blinking corrected CSR test.
# EXAMPLE 3 shows some ways to simulate blinking data.




# EXAMPLE 1: Estimating kinetic rates, derived blinking statistics ----
# --------------------------------------------------
framerate = 25
delt = 1/framerate

# Reading in clustered example dataset
fn = paste(root, "Example data/ClusterPar1.csv", sep = "")
dt = read.csv(fn)

# Converting to a point process object
pp = dt_to_pp(dt)
plot(pp, use.marks = F, pch = 19)

# Estimating kinectic rates. Eta = 1 since we assume no background noise here.
fit = estimate.IBCpp(pp, framerate, dt$sd, eta = 1)

# Parameters stored as (r_F, r_B, r_D, r_R).
params = fit$par
params # True values (0.004, 3, 6, 1).

# Mean of G, the total number of reappearances per PA-FP.
meanG = moment.approx(params[2:4], delt)[1]
meanG # True value of 11.3.

# Expected number of proteins in image.
Nprotein = pp$n/meanG
Nprotein # True value is 500. I obtained 489 when I ran this.

# Mean number of blinks-1 (total number of F state visits)
nblink = pars.to.terms(params[2:4])$mninf
nblink

# Bleaching probability
p = params[2]/sum(params[2:3])
p

# Approximate distribution of G, showing the long tail.
G = approx.gdist(1e4, params[2:4], delt)
hist(G, breaks = 100);abline(v = meanG, lty = 2, col = 2)

# CDF of time from activation to permanent photobleaching.
curve(bleach_cdf(params[2:4])(x), 0, 30, xlab = "Time (s)", ylab = "CDF of PA-FP lifetime")

# Quantile function of the bleaching time CDF is also available.
bleach_qf(params[2:4])(c(0.25,0.5,0.75,0.99))




# EXAMPLE 2: Blinking corrected CSR testing ----
# --------------------------------------------------
framerate = 25
delt = 1/framerate

# Reading in CSR and clustered data
fn_csr = paste(root, "Example data/CSRPar1.csv", sep = "")
dt_csr = read.csv(fn_csr)
fn_clus = paste(root, "Example data/ClusterPar1.csv", sep = "")
dt_clus = read.csv(fn_clus)

# Converting to a point process object
pp_csr = dt_to_pp(dt_csr)
pp_clus = dt_to_pp(dt_clus)
plot(pp_csr, use.marks = F, pch = 19)
plot(pp_clus, use.marks = F, pch = 19)

# Fitting the model to each dataset.
fit_csr  = estimate.IBCpp(pp_csr , framerate, dt_csr$sd , eta = 1)
fit_clus = estimate.IBCpp(pp_clus, framerate, dt_clus$sd, eta = 1)

params_csr  = fit_csr$par
params_clus = fit_clus$par

# Testing for CSR on CSR data (null model = true). Might take a minute or two per test.
# Might want to use larger values of nsim if time is not an issue (say, 500).
CSR_test_null = BlinkingCSRTest(pp_csr , dt_csr$sd , params_csr , 25, nsim = 100)
CSR_test_clus = BlinkingCSRTest(pp_clus, dt_clus$sd, params_clus, 25, nsim = 100)

# Test result when ground truth is CSR.
plot(CSR_test_null) # The observed L(r)-r function is inside the envelope - I get a pvalue of 12%. Not significant.

# Test result when ground truth is clustered.
plot(CSR_test_clus) # The observed L(r)-r function is very clearly outside the envelope - I get a pvalue of < 1%. Highly significantly clustered.

# Thus, we did not reject CSR when proteins were CSR, but we did reject CSR when proteins were clustered.




# EXAMPLE 3: Simulating blinking data ----
# --------------------------------------------------
# Simulation of blinking data is done by adding blinking clusters to a given protein sample.

# We simulate first proteins. Here we make them CSR, but could be anything.
win = square(1e4)
lambdax = 1e-5 # intensity corresponds to roughly 1000 proteins.
proteins = rpoispp(lambdax, win = win)
plot(proteins)

# We set up the blinking dynamics.
pars = c(0.004,2,8,0.5) # Blinking rates as (r_F, r_B, r_D, r_R).
framerate = 25 # Desired framerate
eta = 0.95 # 1 minutes the desired fraction of noise points (eta = 0.95 corresponds to 5% noise).
sd = abs(rnorm(1e3, 25, 20)) # We simulate some localization uncertainties from a desired distribution. These will be used in simulation.

a = sim_blinking_pp(function() proteins, 
                pars.to.simulator(pars[2:4], 0, framerate, T), 
                function() rexp(1, pars[1]), framerate, noise = 1-eta, s = sd)

# The resulting IBCpp can be accessed in the O output
plot(a$O, use.marks = F, pch = 19, cex = 0.25, main = "IBCpp O")

# Can also color points according to blinking cluster membership
points(a$O$y~a$O$x, pch = 19, cex = 0.25, col = a$Labels, asp = 1)

# We can plot the frames against space, to get an idea of the blinking dynamics
par(mfrow = c(1,2))
plot(a$O$marks~a$O$x, use.marks = F, pch = 19, cex = 0.25, main = "IBCpp O", xlab = "x", ylab = "Time (s)")
plot(a$O$marks~a$O$x, use.marks = F, pch = 19, cex = 0.25, main = "IBCpp O", xlab = "x", ylab = "Time (s)", ylim = c(0,500))
par(mfrow = c(1,1))

# Blinking clusters without noise are in Z
plot(a$Z, use.marks = F, pch = 19, cex = 0.25, main = "Z")

# The nosie points are in E
plot(a$E, use.marks = F, pch = 19, cex = 0.25, main = "E")


