library(data.table)
library(reshape2)
library(ggplot2)
library(dplyr)

# source("code/useful R/libraries_and_themes_1_1025.R")
setwd("/Users/Michel/Desktop/Research/MFA/code/")
theme_set(figure_theme)
options(stringsAsFactors = FALSE)

library(broom)
library(tidyr)

# import growth data
growth.dat <- read.table("../MFA code/first_timecourses/PCV_timecourses_1_and_2.csv", header=TRUE, sep=",")
growth.dt <- tbl_df(growth.dat)
growth.dt <- growth.dt %>% mutate(log.pcv = log(pcv))
growth.dt1 <- growth.dt %>% filter(expt == 1)
growth.dt2 <- growth.dt %>% filter(expt == 2)

# generate sequence of xvalues to plot and predict y values using fit
generate.predictions <- function (fit, dat) { 
  min.x <- min(dat$time.pt)
  max.x <- max(dat$time.pt)
  xvals <- seq(min.x, max.x, (max.x-min.x)/200)
  data.frame(time.pt=xvals, predicted=exp(predict(fit, data.frame(time.pt=xvals))))
}
# find extreme x values and return them in data frame
get.min.max <- function (dat) {
  min = min(dat$time.pt)
  max = max(dat$time.pt)
  data.frame(min=min(dat$time.pt), max=max(dat$time.pt))
}
# integrate exponential curve using model parameters
do.integral <- function (intercept, slope, min, max) {
  as.numeric(integrate(function (t) exp(intercept + slope*t), min, max)[1])
}

growth.fits <- growth.dt %>% group_by(expt) %>% do(fit=lm(log.pcv ~ time.pt, data=.), data=(.))
growth.preds <- growth.fits %>% group_by(expt) %>% do(generate.predictions(.$fit[[1]], .$dat[[1]]))

growth.augment <- growth.fits %>% group_by(expt) %>% do(augment(.$fit[[1]], .$data[[1]])) %>% select(-.rownames)
growth.tidy <- growth.fits %>% group_by(expt) %>% do(tidy(.$fit[[1]]))
growth.tidy$term[growth.tidy$term == "(Intercept)"] <- "intercept" 
growth.tidy$term[growth.tidy$term == "time.pt"] <- "slope" 

tidy.spread <- growth.tidy %>% select(-(stderror:p.value)) %>% spread(term, estimate)
growth.min.max <- growth.dt %>% group_by(expt) %>% do(get.min.max(.))
growth.integrals <- merge(tidy.spread, growth.min.max)
growth.integrals$integral <- apply(growth.with.fits, 1, function (x) { do.integral(x['intercept'], x['slope'], x['min'], x['max']) } )

growth.plots <- ggplot() + 
  geom_point(data=growth.dt, aes(x=time.pt, y=pcv)) + 
  geom_line(data=growth.preds, aes(x=time.pt, y=predicted)) +
  facet_grid(~ expt, scales="free")
growth.plots

########

# import extracellular data
extra.dat2 <- read.table("../MFA code/first_timecourses/extracellular_data_timecourse2.csv", header=TRUE, sep=",")
extra.dt2 <- tbl_df(extra.dat2)

### LINEAR MODEL of log-transformed growth data
###### fit log_PCV to linear model -- return intercept and slope
fit.to.exponential <- function (time.pts, measurements) {
  lin.mod <- lm(log(measurements) ~ time.pts)
  intercept <- as.numeric(lin.mod$coef[1])
  slope <- as.numeric(lin.mod$coef[2])
  c(intercept, slope)
}
###### return predictive function using intercept and slope
get.func <- function (time.pts, measurements) {
  params <- fit.to.exponential(time.pts, measurements)
  intercept <- params[1]
  slope <- params[2]
  function (x) exp(intercept + slope*x)
}

growth.k1 <- fit.to.exponential(growth.dt1$time.pt, growth.dt1$pcv)[2]
func1 <- get.func(growth.dt1$time.pt, growth.dt1$pcv)
growth.dt1 <- growth.dt1 %>% mutate(pred.pcv = func1(time.pt))
plot1 <- ggplot(growth.dt1) + geom_point(aes(x=time.pt, y=log(pcv))) + geom_line(aes(x=time.pt, y=log(pred.pcv)))
pcv.auc1 <- as.numeric(integrate(func1,0,24)[1])

growth.k2 <- fit.to.exponential(growth.dt2$time.pt, growth.dt2$pcv)[2]
func2 <- get.func(growth.dt2$time.pt, growth.dt2$pcv)
growth.dt2 <- growth.dt2 %>% mutate(pred.pcv = func2(time.pt))
plot2 <- ggplot(growth.dt2) + geom_point(aes(x=time.pt, y=pcv)) + geom_line(aes(x=time.pt, y=pred.pcv))
pcv.auc2 <- as.numeric(integrate(func2,0,12)[1])

### calculate uptake rate by plotting uptake vs pcv
extra.dt2 <- extra.dt2 %>% mutate(sum = unlabeled + labeled)
#aa_table <- data.table(compound = c("lysine","threonine","tyrosine"), DMEM.conc = c(800,800,500))
#save(aa_table, file="aa_DMEM_conc_table.Rda")
load(file="aa_DMEM_conc_table.Rda") # aa_conc defined here
aa.conv.dt <- extra.dt2 %>% filter(time.pt == 0) %>% group_by(compound) %>% summarize(mean=mean(sum))
aa.conv.dt <- merge(aa.conv.dt, aa_table, by="compound") 
aa.conv.dt <- aa.conv.dt %>% mutate(conv.factor = DMEM.conc*.003/mean)

full.dt2 <- merge(extra.dt2, aa.conv.dt, by="compound")
select.dt2 <- full.dt2 %>% mutate(mol=sum*conv.factor) %>% dplyr::select(compound,time.pt,mol)

merged.dt2 <- merge(select.dt2, unique(growth.dt2 %>% dplyr::select(time.pt, pred.pcv))) 

aas <- unique(merged.dt2$compound)
means <- c()
stderrs <- c()
ints <- c()
for (aa in aas) {
  temp.dt <- merged.dt2 %>% filter(compound == aa)
  temp.summ <- summary(lm(temp.dt$mol ~ temp.dt$pred.pcv))
  means <- c(means, temp.summ$coefficients[2,1])
  stderrs <- c(stderrs, temp.summ$coefficients[2,2])
  ints <- c(ints, temp.summ$coefficients[1,1])
}
summs.dt2 <- data.table(compound=aas, model.mean=means, model.stderr=stderrs, model.int=ints)
summs.dt2 <- summs.dt2 %>% mutate(uptake.mean=model.mean*growth.k2, uptake.se=model.stderr*growth.k2)

big.dt2 <- merge(merged.dt2, summs.dt2)
tidy.dt2 <- big.dt2 %>% mutate(pred.mol = pred.pcv*model.mean + model.int) %>% 
  dplyr::select(-c(model.mean, model.stderr, model.int, uptake.mean, uptake.se))

uptake.plots <- ggplot(tidy.dt2) + geom_point(aes(x=pred.pcv, y=mol)) + geom_line(aes(x=pred.pcv, y=pred.mol)) + 
  facet_wrap( ~ compound, scale="free") + xlab("Packed cell volume (uL)") + ylab("mmol AA in medium")
uptake.plots

### need to fit full equation using metabolomics data.....
frac.dt2 <- extra.dt2 %>% mutate(frac = unlabeled/sum)
frac.quickplot <- ggplot(frac.dt2) + geom_point(aes(x=time.pt, y=frac)) + 
  facet_wrap( ~ compound, scale="free_y")
frac.quickplot



# foo <- frac.dt2 %>% group_by(compound) %>% do(lm=lm(frac) ~ time.pt), data=.)
# bar <- foo %>% group_by(compound) %>% summarize(slope=summary(lm[[1]])$coefficients[2,1]) #, intercept=lm[1]$coefficients[1])

library(broom)

performfit <- function(dat, type="lm") {
  # returns an lm or nls fit
  if (type == "lm") {
    lm(log(frac) ~ time.pt, data=dat)
  } else if (type == "nls") {
    nls(frac ~ exp(intercept + rate * time.pt), data=dat, start=list(intercept=0, rate=.1))
  }
}

# double it up for two types of fits
doubled <- frac.dt2 %>% data.frame(type=c("lm", "nls")) %>% group_by(type) %>% do(frac.dt2)
fits <- doubled %>% group_by(type, compound) %>% do(fit=performfit(., type=.$type[[1]]), data=(.))

augmented <- fits %>% group_by(type, compound) %>% do(augment(.$fit[[1]], .$data[[1]]))
preds <- fits %>% group_by(type, compound) %>%
  do(data.frame(time.pt=xvals, predicted=predict(.$fit[[1]], newdata=data.frame(time.pt=xvals))))
preds <- preds %>% mutate(predicted=ifelse(type == "lm", exp(predicted), predicted))

ggplot(doubled, aes(time.pt, frac)) + geom_point() + facet_wrap(type ~ compound) +
  geom_line(aes(time.pt, predicted), data=preds)

allcoefs <- fits %>% group_by(type, compound) %>% do(tidy(.$fit[[1]]))
allcoefs$term[allcoefs$term == "(Intercept)"] = "intercept"
allcoefs$term[allcoefs$term == "time.pt"] = "rate"

####
glances <- fits %>% group_by(type, compound) %>% do(glance(.$fit[[1]]))

# compare nls and lm values
library(tidyr)

estimates <- allcoefs %>% dplyr::select(type, compound, term, estimate) %>% spread(type, estimate)

ggplot(estimates, aes(lm, nls, color=compound)) + geom_point() + theme(legend.position="right") +
  geom_abline() + facet_wrap(~ term, scales="free")

augmented <- frac.dt2 %>% group_by(compound) %>% do(augment(lm(log(frac) ~ time.pt, data=.), .))

ggplot(augmented, aes(time.pt, frac)) + geom_point() + facet_wrap(~ compound) +
  geom_line(aes(y=exp(.fitted))) + scale_y_log10() 

allcoefs_nls <- frac.dt2 %>% group_by(compound) %>% do(tidy(nls(frac ~ exp(b + m * time.pt), data=., start=list(b=0, m=.1))))
augmented_nls <- frac.dt2 %>% group_by(compound) %>% do(augment(nls(frac ~ exp(b + m * time.pt), data=., start=list(b=0, m=.1)), .))

ggplot(augmented_nls, aes(time.pt, frac)) + geom_point() + facet_wrap(~ compound) +
  geom_line(aes(y=.fitted))

#, formula=I(log(y)) ~ x)

foo %>% gather(metric, value, estimate:p.value)
