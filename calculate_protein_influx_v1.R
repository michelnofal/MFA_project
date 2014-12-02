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
growth.dat <- read.table("first_timecourses/PCV_timecourses_1_and_2.csv", header=TRUE, sep=",")
growth.dt <- tbl_df(growth.dat)
growth.dt <- growth.dt %>% mutate(log.pcv = log(pcv))

# import intracellular data
intra.dat <- read.table("../MFA code/first_timecourses/intracellular_data_timecourse2.csv", header=TRUE, sep=",")
intra.dt <- tbl_df(intra.dat2)
intra.dt <- intra.dt %>% mutate(sum=unlabeled+labeled, frac=unlabeled/sum, loc="intracellular")

# import extracellular data
extra.dat <- read.table("../MFA code/first_timecourses/extracellular_data_timecourse2.csv", header=TRUE, sep=",")
extra.dt <- tbl_df(extra.dat2)
extra.dt <- extra.dt %>% mutate(sum=unlabeled+labeled, frac=unlabeled/sum, loc="extracellular")

# combine intracellular and extracellular data
all.aa <- rbind(intra.dt, extra.dt)

###########################
######## FUNCTIONS ########
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
integrate.product.exponentials <- function (intercept1, slope1, intercept2, slope2, max) {
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

###########################
###########################

growth.fits <- growth.dt %>% group_by(expt) %>% do(fit=lm(log.pcv ~ time.pt, data=.), data=(.))
growth.preds.curve <- growth.fits %>% group_by(expt) %>% do(generate.predictions(.$fit[[1]], .$dat[[1]]))
growth.preds.curve <- growth.preds.curv %>% mutate(predicted = exp(predicted))

growth.augment <- growth.fits %>% group_by(expt) %>% do(augment(.$fit[[1]], .$data[[1]])) %>% select(-.rownames)
growth.preds <- growth.augment %>% mutate(predicted = exp(.fitted)) %>% select(time.pt, pcv, expt, predicted)
growth.tidy <- growth.fits %>% group_by(expt) %>% do(tidy(.$fit[[1]]))
growth.tidy$term[growth.tidy$term == "(Intercept)"] <- "intercept" 
growth.tidy$term[growth.tidy$term == "time.pt"] <- "slope" 

growth.tidy.spread <- growth.tidy %>% select(-(stderror:p.value)) %>% spread(term, estimate)
growth.max <- growth.dt %>% group_by(expt) %>% do(data.frame(max = max(.$time.pt)))
growth.integrals <- merge(growth.tidy.spread, growth.max)
growth.integrals$growth.integral <- apply(growth.integrals, 1, function (x) { do.integral(x['intercept'], x['slope'], x['max']) } )

# plotting growth data and fitted curves
growth.plots <- ggplot() + 
  geom_point(data=growth.dt, aes(x=time.pt, y=pcv)) + 
  geom_line(data=growth.preds.curve, aes(x=time.pt, y=predicted)) +
  facet_wrap(~ expt) #, scales="free_x")
growth.plots

########
########

all.fits <- all.aa %>% group_by(compound, loc) %>% 
  do(fit=nls(frac ~ exp(intercept + rate*time.pt), data=., start=list(intercept=0, rate=.1)), data=(.))
all.preds <- all.fits %>% group_by(compound, loc) %>% do(generate.predictions(.$fit[[1]], .$dat[[1]]))

all.aug <- all.fits %>% group_by(compound, loc) %>% do(augment(.$fit[[1]], .$data[[1]]))
all.tidy <- all.fits %>% group_by(compound, loc) %>% do(tidy(.$fit[[1]]))
all.tidy$estimate <- as.numeric(all.tidy$estimate)

all.tidy.spread <- all.tidy %>% select(-(stderror:p.value)) %>% spread(term, estimate)
all.max <- all.aa %>% group_by(compound, loc) %>% do(data.frame(max = max(.$time.pt)))
#all.integrals <- merge(all.tidy.spread, all.max)
#all.integrals$all.integral <- apply(all.integrals, 1, function (x) { do.integral(x['intercept'], x['rate'], x['max']) } )

# plotting amino acid labeling data and fitted curves
all.aa$loc <- factor(all.aa$loc, levels = c("intracellular","extracellular"))
all.preds$loc <- factor(all.preds$loc, levels = c("intracellular","extracellular"))
all.aa.plots <- ggplot () + 
  geom_point(data=all.aa, aes(x=time.pt, y=frac)) + 
  geom_line(data=all.preds ,aes(x=time.pt, y=predicted)) + 
  facet_grid(compound ~ loc, scales="free_y")
all.aa.plots
# how do i add geom_abline???

output.list <- list()
output.list$spread <- all.tidy.spread
output.list$aa.plot <- all.aa.plots

########
########

# convert counts to umol
#aa_table <- data.table(compound = c("lysine","threonine","tyrosine"), DMEM.conc = c(800,800,500))
#save(aa_table, file="aa_DMEM_conc_table.Rda")
load(file="aa_DMEM_conc_table.Rda") 
aa.conv.dt <- all.aa %>% filter(time.pt == 0, loc == "extracellular") %>% group_by(compound) %>% summarize(mean=mean(sum))
aa.conv.dt <- merge(aa.conv.dt, aa_table, by="compound") 
aa.conv.dt <- aa.conv.dt %>% mutate(conv.factor = DMEM.conc*.003/mean) %>% select(-mean, -DMEM.conc)

extra.sums.dt <- all.aa %>% filter(loc == "extracellular") %>% select(compound, time.pt, sum, frac)
extra.umol.dt <- merge(extra.sums.dt, aa.conv.dt)
extra.umol.dt <- extra.umol.dt %>% mutate(umol = sum*conv.factor) %>% select(-sum, -conv.factor) %>% mutate(umol.unlab = umol*frac)

# calculate uptake rates
pcv.v.extra.dt <- merge(growth.preds %>% filter(expt == 2), extra.umol.dt)

pve.fits <- pcv.v.extra.dt %>% group_by(compound) %>% do(fit=lm(umol ~ predicted, data=.))
pve.tidy <- pve.fits %>% group_by(compound) %>% do(tidy(.$fit[[1]]))
pve.tidy$term[pve.tidy$term == "(Intercept)"] = "intercept"
pve.tidy$term[pve.tidy$term == "predicted"] = "slope"
pve.spread <- pve.tidy %>% select(-(stderror:p.value)) %>% spread(term, estimate)

growth.k.dt <- (growth.tidy %>% filter(expt == 2))$estimate[2]
uptake.rates.dt <- pve.spread %>% mutate(uptake.rate = -slope*growth.k.dt) %>% select(compound, uptake.rate)

# calculate umol unlab amino acid at final time point
final.umol.unlab.dt <- merge(extra.umol.dt, all.max) %>% filter(loc == "extracellular", time.pt == max) %>% 
  group_by(compound) %>% summarize(mean.umol.unlab=mean(umol.unlab))

final.umol.unlab.dt <- merge(extra.umol.dt, all.max) %>% filter(loc == "extracellular", time.pt %in% c(0,max)) %>% 
  group_by(compound, time.pt) %>% summarize(mean.umol.unlab=mean(umol.unlab))
final.umol.unlab.dt$time.pt[final.umol.unlab.dt$time.pt == 0] = "time0"
final.umol.unlab.dt$time.pt[final.umol.unlab.dt$time.pt != "time0"] = "endtime"
final.umol.unlab.dt <- final.umol.unlab.dt %>% spread(time.pt, mean.umol.unlab)
final.umol.unlab.dt <- final.umol.unlab.dt %>% mutate(net.umol.unlab = endtime - time0) %>% select(-time0, -endtime)

# integrating two curves
growth.tidy.spread.dt <- merge(growth.tidy.spread, growth.max) %>% filter(expt == 2) %>% select(-expt)
intra.tidy.spread <- all.tidy.spread %>% filter(loc == "intracellular")
intra.tidy.spread <- intra.tidy.spread %>% mutate(aa.intercept = intercept, aa.slope = rate) %>% 
  select(compound, aa.intercept, aa.slope)
product.integrals <- merge(intra.tidy.spread, growth.tidy.spread.dt)
product.integrals$product.integral <- apply(product.integrals, 1, 
                                    function (x) { integrate.product.exponentials(
                                      x['aa.intercept'], x['aa.slope'], x['intercept'], x['slope'], x['max']) } )

########
########

alpha.df <- data.frame(compound = c("lysine","threonine","tyrosine"), alpha = c(59, 32, 21))

# gather equation variables
eqn.vars.dt <- merge(alpha.df, merge(final.umol.unlab.dt, uptake.rates.dt))
eqn.vars.dt$growth.integral <- (growth.integrals %>% filter(expt == 2))$growth.integral
eqn.vars.dt <- merge(eqn.vars.dt, product.integrals %>% select(compound, product.integral))
# plug variables into equation
eqn.vars.dt$umol.prot.flux <- apply(eqn.vars.dt, 1, 
                                    function(x) { 
                                      do.equation(x['alpha'], x['net.umol.unlab'], x['uptake.rate'], x['growth.integral'], x['product.integral']) } )

eqn.vars.dt <- eqn.vars.dt %>% mutate(resulting.aa.flux = alpha*umol.prot.flux)

compare.uptake.dt <- eqn.vars.dt %>% select(-net.umol.unlab, -(growth.integral:umol.prot.flux)) %>%
  gather(uptake.route, flux, c(uptake.rate, resulting.aa.flux))

compare.uptake.plot <- ggplot(compare.uptake.dt, aes(x=compound, y=flux, fill=uptake.route)) + 
  geom_bar(stat="identity", position="dodge", color="black") + 
  scale_fill_manual(values=c("grey50","tomato2")) + 
  geom_abline(height=0, slope=0) +
  xlab("") + ylab("flux (umol / uL cell / hr)") + theme(legend.position = "right")
compare.uptake.plot

########################################
########################################
