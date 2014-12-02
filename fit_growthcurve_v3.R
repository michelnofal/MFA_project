library(data.table)
library(reshape2)
library(ggplot2)
library(dplyr)

# source("code/useful R/libraries_and_themes_1_1025.R")
setwd("/Users/Michel/Desktop/Research/MFA/code/")
theme_set(figure_theme)
options(stringsAsFactors = FALSE)

# import growth data
growth.dat <- read.table("../MFA code/first_timecourses/PCV_timecourses_1_and_2.csv", header=TRUE, sep=",")
growth.dt <- tbl_df(growth.dat)
growth.dt <- growth.dt %>% mutate(log.pcv = log(pcv))
growth.dt1 <- growth.dt %>% filter(expt == 1)
growth.dt2 <- growth.dt %>% filter(expt == 2)

# import extracellular data
extra.dat2 <- read.table("../MFA code/first_timecourses/extracellular_data_timecourse2.csv", header=TRUE, sep=",")
extra.dt2 <- tbl_df(extra.dat2)

### prediction using smooth.spline on untransformed data
#spl1 <- smooth.spline(growth.dt1$time.pt, growth.dt1$pcv)
#plot(growth.dt1$time.pt, growth.dt1$pcv, col = "gray"); lines(predict(spl.1, c(0,4,8,12)), col = 2)
#integrate(,lower=0,upper=24)

### prediction using lm to fit log-transformed data

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
plot1 <- ggplot(growth.dt1) + geom_point(aes(x=time.pt, y=pcv)) + geom_line(aes(x=time.pt, y=pred.pcv))
pcv.auc1 <- as.numeric(integrate(func1,0,24)[1])

growth.k2 <- fit.to.exponential(growth.dt2$time.pt, growth.dt2$pcv)[2]
func2 <- get.func(growth.dt2$time.pt, growth.dt2$pcv)
growth.dt2 <- growth.dt2 %>% mutate(pred.pcv = func2(time.pt))
plot2 <- ggplot(growth.dt2) + geom_point(aes(x=time.pt, y=pcv)) + geom_line(aes(x=time.pt, y=pred.pcv))
pcv.auc2 <- as.numeric(integrate(func2,0,12)[1])

### calculate uptake rate by plotting uptake vs pcv
extra.dt2 <- extra.dt2 %>% mutate(sum = unlabeled + labeled)
merged.dt2 <- merge(extra.dt2, unique(growth.dt2 %>% dplyr::select(time.pt, pred.pcv))) 


lys.dt2 <- merged.dt2 %>% filter(compound == "lysine")
thr.dt2 <- merged.dt2 %>% filter(compound == "threonine")

lys.init.mol <- 800*0.003 # (umol)
lys.conv2 <- lys.init.mol/mean((lys.dt2 %>% filter(time.pt == 0))$sum)
lys.dt2 <- lys.dt2 %>% mutate(mol = sum*lys.conv2)

lys.lm2 <- lm(lys.dt2$mol ~ lys.dt2$pred.pcv)
lys.summ2 <- summary(lys.lm2)
lys.est2 <- lys.summ2$coefficients[2,1]*growth.k2
lys.stderr2 <- lys.summ2$coefficients[2,2]*growth.k2
lys.int2 <- lys.summ2$coefficients[1,1]
lys.dt2 <- lys.dt2 %>% mutate(pred.mol = pred.pcv*lys.est2/growth.k2 + lys.int2)

lys.plot2 <- ggplot(lys.dt2) + geom_point(aes(x=pred.pcv, y=mol)) + geom_line(aes(x=pred.pcv, y=pred.mol))
lys.plot2

thr.init.mol <- 800*0.003 # (umol)
thr.conv2 <- thr.init.mol/mean((thr.dt2 %>% filter(time.pt == 0))$sum)
thr.dt2 <- thr.dt2 %>% mutate(mol = sum*thr.conv2)

thr.lm2 <- lm(thr.dt2$mol ~ thr.dt2$pred.pcv)
thr.summ2 <- summary(thr.lm2)
thr.est2 <- thr.summ2$coefficients[2,1]*growth.k2
thr.stderr2 <- thr.summ2$coefficients[2,2]*growth.k2
thr.int2 <- thr.summ2$coefficients[1,1]
thr.dt2 <- thr.dt2 %>% mutate(pred.mol = pred.pcv*thr.est2/growth.k2 + thr.int2)


thr.plot2 <- ggplot(thr.dt2) + geom_point(aes(x=pred.pcv, y=mol)) + geom_line(aes(x=pred.pcv, y=pred.mol))
thr.plot2




### need to fit full equation using metabolomics data.....
