library(tidyverse)
library(magrittr)
library(ggplot2)
library(wrapr)
library(zeallot)
library(janitor)
library(pracma)
library(spatstat)
library(purrr)
library(markovchain)

dir = "C:/Users/louis/Documents/PhD/Produced Work/Articles/2019 - A spatio-temporal point-process model for correction of photoblinking artifacts/STPP4BC/R"
ffs = c("Estimation.R", "Gamma1.R", "Misc.R", "Space-Time Simulation.R", "Temporal Simulation.R")

for(f in ffs) {
  fn = paste(dir, "/", f, sep = "")
  source(fn)
}

rm(ffs, fn, dir, f)
