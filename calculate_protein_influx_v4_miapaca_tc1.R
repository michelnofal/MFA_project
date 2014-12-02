library(data.table)
library(reshape2)
library(ggplot2)
library(dplyr)
library(broom)
library(tidyr)
options(stringsAsFactors = FALSE)

# setwd("/Users/Michel/Desktop/Research/CRISPR screen/")
# source("code/useful R/libraries_and_themes_1_1025.R")
theme_set(figure_theme)
setwd("/Users/Michel/Desktop/Research/MFA/code/")
# source("protein_influx_functions_v2.R")

file_handle <- "miapaca_timecourse1"

# import growth data
growth.dat <- read.table("first_timecourses/miapaca2_tc1_pcv.csv", header=TRUE, sep=",")
# import intracellular data
intra.dat <- read.table("first_timecourses/miapaca2_tc1_intracellullar.csv", header=TRUE, sep=",")
# import extracellular data
extra.dat <- read.table("first_timecourses/miapaca2_tc1_extracellullar.csv", header=TRUE, sep=",")

# 1 - calculate growth predictions
# 2 - calculate fitted growth curve parameters
# 3 - calculate area under growth curve
# 4 - generate growth plot
growth.list <- growth.func(growth.dat)

# 1 - return organized intra- and extracellular data
# 2 - calculate fitted amino acid labeling curve parameters and predictions
# 3 - generate amino acid labeling plots
aa.lab.list <- aa.lab.func(intra.dat, extra.dat)

# 1 - calculate amino acid uptake rates
# 2 - generate plot of predicted pcv vs amino acid levels in medium (slope = uptake rate / k_growth)
aa_table <- data.table(compound = c("lysine","threonine","tyrosine"), DMEM.conc = c(800,800,500))
uptake.list <- calc.uptake.rates(aa_table, aa.lab.list$data, growth.list$preds, growth.list$curve.params)

# calculate umol unlab amino acid at final time point
final.umol.unlab <- get.final.umol.unlab(aa_table, aa.lab.list$data)

# calculate integral of product of exponentials describing growth rate and amino acid labeling
integral.prod.curves <- integrate.prod.exponentials(growth.list$curve.params, growth.list$integrals$max, aa.lab.list$curve.params) 

# 1 - calculate flux using MFA equation
# 2 - plot flux estimates from different aas
# 3 - plot influx of amino acids through monomeric uptake vs through protein catabolism 
alpha.df <- data.frame(compound = c("lysine","threonine","tyrosine"), alpha = c(59, 32, 21))
fluxes.list <- compute.flux(alpha.df, final.umol.unlab, uptake.list$uptake.rates, growth.list$integrals, integral.prod.curves)

save.plots(growth.list$plot, aa.lab.list$plots, uptake.list$plot, fluxes.list$plot.fluxes, fluxes.list$plot.compare)

growth.list$plot
aa.lab.list$plots
uptake.list$plot
fluxes.list$plot.fluxes
fluxes.list$plot.compare

miapaca_timecourse1_fluxes <- fluxes.list$fluxes
# save(miapaca_timecourse1_fluxes, file="miapaca_timecourse1_fluxes.Rda")
