library(data.table)
library(reshape2)
library(ggplot2)
library(dplyr)

# setwd("/Users/Michel/Desktop/Research/CRISPR screen/")
# source("code/useful R/libraries_and_themes_1_1025.R")
setwd("/Users/Michel/Desktop/Research/MFA/code/")
theme_set(figure_theme)
options(stringsAsFactors = FALSE)

library(broom)
library(tidyr)

# import growth data
growth.dat <- read.table("first_timecourses/PCV_timecourse_2.csv", header=TRUE, sep=",")
# import intracellular data
intra.dat <- read.table("../MFA code/first_timecourses/intracellular_data_timecourse2.csv", header=TRUE, sep=",")
# import extracellular data
extra.dat <- read.table("../MFA code/first_timecourses/extracellular_data_timecourse2.csv", header=TRUE, sep=",")

###########################
#### GENERAL FUNCTIONS ####
###########################

# generate sequence of xvalues to plot and predict y values using fit
generate.predictions <- function (fit, dat) { 
  max.x <- max(dat$time.pt)
  xvals <- seq(0, max.x, max.x/200)
  data.frame(time.pt=xvals, predicted=predict(fit, data.frame(time.pt=xvals)))
}

# integrate exponential curve using model parameters
do.integral <- function (intercept, slope, max) {
  intercept <- as.numeric(intercept)
  slope <- as.numeric(slope)
  max <- as.numeric(max)
  as.numeric(integrate(function (t) exp(intercept + slope*t), 0, max)[1])
}

# integrate product of 2 exponentials
do.integral.product.exponentials <- function (intercept1, slope1, intercept2, slope2, max) {
  intercept1 <- as.numeric(intercept1)
  slope1 <- as.numeric(slope1)
  intercept2 <- as.numeric(intercept2)
  slope2 <- as.numeric(slope2)
  max <- as.numeric(max)
  as.numeric(integrate(function (t) exp((intercept1 + slope1*t) + (intercept2 + slope2*t)), 0, max)[1])
}

# plug calculated concentrations / rates / integrals into final MFA equation
do.equation <- function (alpha, final.unlab.med, uptake.rate, growth.int, prod.int) {
  alpha <- as.numeric(alpha)
  final.unlab.med <- as.numeric(final.unlab.med)
  uptake.rate <- as.numeric(uptake.rate)
  growth.int <- as.numeric(growth.int)
  prod.int <- as.numeric(prod.int)
  numerator <- final.unlab.med + uptake.rate*prod.int
  denominator <- growth.int - prod.int
  (1/alpha)*numerator/denominator
}

# convert extracellular amino acid counts to umol in medium
get.extra.umol <- function(aa_table, aa.dt) {
  init.extra.aa <- aa.dt %>% filter(time.pt == 0, loc == "extracellular") %>% group_by(compound) %>% summarize(mean=mean(sum))
  aa.conv.dt <- merge(init.extra.aa, aa_table, by="compound") 
  aa.conv.dt <- aa.conv.dt %>% mutate(conv.factor = DMEM.conc*.003/mean) %>% select(-mean, -DMEM.conc)
  
  extra.sums.dt <- aa.dt %>% filter(loc == "extracellular") %>% select(compound, time.pt, sum, frac)
  extra.umol.dt <- merge(extra.sums.dt, aa.conv.dt)
  extra.umol.dt <- extra.umol.dt %>% mutate(umol = sum*conv.factor) %>% select(-sum, -conv.factor) %>% mutate(umol.unlab = umol*frac)
  
  return(extra.umol.dt)
}

###########################
##### OTHER FUNCTIONS #####
###########################

growth.func <- function(data) {
  output.list <- list() # any object to be returned will be saved in this list
  
  growth.dt <- tbl_df(data)
  growth.fits <- growth.dt %>% do(fit=nls(pcv ~ exp(intercept + rate*time.pt), data=., start=list(intercept=0, rate=.1)), data=(.))
  
  growth.augment <- growth.fits %>% do(augment(.$fit[[1]], .$data[[1]]))
  growth.preds <- growth.augment %>% mutate(predicted = .fitted) %>% select(time.pt, pcv, predicted)
  output.list$preds <- growth.preds
  
  growth.tidy <- growth.fits %>% do(tidy(.$fit[[1]]))
  curve.params <- data.frame(intercept = growth.tidy[growth.tidy$term == "intercept",]$estimate, 
                             rate = growth.tidy[growth.tidy$term == "rate",]$estimate)
  output.list$curve.params <- curve.params
  
  growth.max <- growth.dt %>% do(data.frame(max = max(.$time.pt)))
  growth.integrals <- merge(curve.params, growth.max)
  growth.integrals$growth.integral <- apply(growth.integrals, 1, function (x) { do.integral(x['intercept'], x['rate'], x['max']) } )
  output.list$integrals <- growth.integrals
  
  # predict pcv using fit
  growth.preds.curve <- growth.fits %>% do(generate.predictions(.$fit[[1]], .$dat[[1]]))
  # plotting growth data and fitted curves
  growth.plot <- ggplot() + 
    geom_point(data=growth.dt, aes(x=time.pt, y=pcv)) + 
    geom_line(data=growth.preds.curve, aes(x=time.pt, y=predicted))
  output.list$plot <- growth.plot
  
  return(output.list)
}

# calculate growth predictions
# calculate fitted growth curve parameters
# calculate area under growth curve
# generate growth plot
growth.list <- growth.func(growth.dat)

########
########

aa.lab.func <- function (intra.data, extra.data) {
  output.list <- list()
  
  intra.dt <- tbl_df(intra.data)
  intra.dt <- intra.dt %>% mutate(sum=unlabeled+labeled, frac=unlabeled/sum, loc="intracellular")
  
  extra.dt <- tbl_df(extra.data)
  extra.dt <- extra.dt %>% mutate(sum=unlabeled+labeled, frac=unlabeled/sum, loc="extracellular")
  
  aa.dt <- rbind(intra.dt, extra.dt)
  output.list$data <- aa.dt

  aa.fits <- aa.dt %>% group_by(compound, loc) %>% 
    do(fit=nls(frac ~ exp(intercept + rate*time.pt), data=., start=list(intercept=0, rate=.1)), data=(.))
  
  aa.aug <- aa.fits %>% group_by(compound, loc) %>% do(augment(.$fit[[1]], .$data[[1]]))
  aa.tidy <- aa.fits %>% group_by(compound, loc) %>% do(tidy(.$fit[[1]]))
  aa.tidy$estimate <- as.numeric(aa.tidy$estimate)
  aa.tidy.spread <- aa.tidy %>% select(-(stderror:p.value)) %>% spread(term, estimate)
  output.list$curve.params <- aa.tidy.spread
  
  aa.preds <- aa.fits %>% group_by(compound, loc) %>% do(generate.predictions(.$fit[[1]], .$dat[[1]]))
  aa.dt$loc <- factor(aa.dt$loc, levels = c("intracellular","extracellular"))
  aa.preds$loc <- factor(aa.preds$loc, levels = c("intracellular","extracellular"))
  output.list$preds <- aa.preds
  
  aa.plots <- ggplot () + 
    geom_point(data=aa.dt, aes(x=time.pt, y=frac)) + 
    geom_line(data=aa.preds, aes(x=time.pt, y=predicted)) + 
    facet_grid(compound ~ loc, scales="free_y")
  output.list$plots <- aa.plots
  
  return(output.list)
}

# 1 - return organized intra- and extracellular data
# 2 - calculate fitted amino acid labeling curve parameters and predictions
# 3 - generate amino acid labeling plots
aa.lab.list <- aa.lab.func(intra.dat, extra.dat)

########
########

calc.uptake.rates <- function (aa_table, aa.dt, growth.preds, growth.curve.params) {
  output.list <- list()
  
  # convert amino acid 
  extra.umol.dt <- get.extra.umol(aa_table, aa.dt)
  
  # make data frame containing pcv predictions and extracellular 
  pcv.v.extra <- merge(unique(growth.preds %>% select(-pcv)), extra.umol.dt)
  
  # make model to predict extracellular aa abundance using predicted pcv
  pve.fits <- pcv.v.extra %>% group_by(compound) %>% do(fit=lm(umol ~ predicted, data=.), data=(.))
  pve.tidy <- pve.fits %>% group_by(compound) %>% do(tidy(.$fit[[1]]))
  pve.tidy$term[pve.tidy$term == "(Intercept)"] = "intercept"
  pve.tidy$term[pve.tidy$term == "predicted"] = "slope"
  pve.spread <- pve.tidy %>% select(-(stderror:p.value)) %>% spread(term, estimate)
  
  growth.k <- growth.curve.params$rate
  uptake.rates <- pve.spread %>% mutate(uptake.rate = -growth.k*slope) %>% select(compound, uptake.rate)
  output.list$uptake.rates <- uptake.rates
  
  pve.aug <- pve.fits %>% group_by(compound) %>% do(augment(.$fit[[1]], .$dat[[1]]))
  pve.aug <- pve.aug %>% mutate(pred.umol = .fitted) %>% select(predicted, compound, umol, pred.umol)
  pve.plot <- ggplot(pve.aug) + geom_point(aes(x=predicted,y=umol)) + geom_line(aes(x=predicted, y=pred.umol)) +
    facet_grid( ~ compound) + scale_y_continuous(limits=c(0, ceiling(max(pve.aug$umol))))
  output.list$plot <- pve.plot
  
  return(output.list)
}

# 1 - calculate amino acid uptake rates
# 2 - generate plot of predicted pcv vs amino acid levels in medium (slope = uptake rate / k_growth)
aa_table <- data.table(compound = c("lysine","threonine","tyrosine"), DMEM.conc = c(800,800,500))
uptake.list <- calc.uptake.rates(aa_table, aa.lab.list$data, growth.list$preds, growth.list$curve.params)

########
########

get.final.umol.unlab <- function(aa_table, aa.dt) {
  extra.umol.dt <- get.extra.umol(aa_table, aa.dt)
  
  # calculate umol unlab amino acid at final time point
  final.umol.unlab.dt <- merge(extra.umol.dt, all.max) %>% filter(loc == "extracellular", time.pt == max) %>% 
    group_by(compound) %>% summarize(mean.umol.unlab=mean(umol.unlab))
  
  final.umol.unlab.dt <- merge(extra.umol.dt, all.max) %>% filter(loc == "extracellular", time.pt %in% c(0,max)) %>% 
    group_by(compound, time.pt) %>% summarize(mean.umol.unlab=mean(umol.unlab))
  final.umol.unlab.dt$time.pt[final.umol.unlab.dt$time.pt == 0] = "time0"
  final.umol.unlab.dt$time.pt[final.umol.unlab.dt$time.pt != "time0"] = "endtime"
  final.umol.unlab.dt <- final.umol.unlab.dt %>% spread(time.pt, mean.umol.unlab)
  final.umol.unlab.dt <- final.umol.unlab.dt %>% mutate(net.umol.unlab = endtime - time0) %>% select(-time0, -endtime)
  
  return(final.umol.unlab.dt)
}

# calculate umol unlab amino acid at final time point
final.umol.unlab <- get.final.umol.unlab(aa_table, aa.lab.list$data)
  
########
########

integrate.prod.exponentials <- function (growth.curve.params, growth.max, aa.params) {
  growth.dt <- growth.curve.params %>% mutate(max = growth.max)
  intra.params <- aa.params %>% filter(loc == "intracellular")
  intra.params <- intra.params %>% mutate(aa.intercept = intercept, aa.slope = rate) %>% 
    select(compound, aa.intercept, aa.slope)
  integrals <- merge(intra.params, growth.dt)
  integrals$integral.prods <- apply(integrals, 1, function (x) {
    do.integral.product.exponentials(x['aa.intercept'], x['aa.slope'], x['intercept'], x['rate'], x['max']) } )
  integrals <- integrals %>% select(compound, integral.prods)
  
  return(integrals)
}

# calculate integral of product of exponentials describing 1 - growth rate and 2 - amino acid labeling
integral.prod.curves <- integrate.prod.exponentials(growth.list$curve.params, growth.integrals$max, aa.lab.list$curve.params) 

########
########

alpha.df <- data.frame(compound = c("lysine","threonine","tyrosine"), alpha = c(59, 32, 21))

compute.flux <- function (alpha.df, final.umol.unlab, uptake.rates, growth.integral, integral.prods) {
  output.list <- list()
  
  eqn.vars <- merge(alpha.df, merge(final.umol.unlab, uptake.rates))
  eqn.vars$growth.integral <- growth.integral$growth.integral
  eqn.vars <- merge(eqn.vars, integral.prods)
  eqn.vars$umol.protein.flux <- apply(eqn.vars, 1, function (x) {
    do.equation(x['alpha'], x['net.umol.unlab'], x['uptake.rate'], x['growth.integral'], x['integral.prods']) } )
  output.list$fluxes <- eqn.vars
  
  compare.uptakes <- eqn.vars %>% select(-net.umol.unlab, -growth.integral, -integral.prods) %>%
    gather(uptake.route, flux, c(uptake.rate, umol.protein.flux))
  compare.uptake.plot <- ggplot(compare.uptake.dt, aes(x=compound, y=flux, fill=uptake.route)) + 
    geom_bar(stat="identity", position="dodge", color="black") + 
    scale_fill_manual(values=c("grey50","tomato2")) + 
    geom_abline(height=0, slope=0) +
    xlab("") + ylab("flux (umol / uL cell / hr)")
  output.list$plot <- compare.uptake.plot
  
  return(output.list)
}

fluxes.list <- compute.flux(alpha.df, final.umol.unlab, uptake.list$uptake.rates, growth.list$integrals, integral.prod.curves)