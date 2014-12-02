library(data.table)
library(reshape2)
library(ggplot2)
library(dplyr)

# source("code/useful R/libraries_and_themes_1_1025.R")
setwd("/Users/Michel/Desktop/Research/MFA/code/")
theme_set(figure_theme)
options(stringsAsFactors = FALSE)

growth.dat <- read.table("../MFA code/first_timecourses/PCV_timecourses_1_and_2.csv", header=TRUE, sep=",")
growth.dt <- tbl_df(growth.dat)
growth.dt <- growth.dt %>% mutate(log.pcv = log(pcv))
growth.dt1 <- growth.dt %>% filter(expt == 1)
growth.dt2 <- growth.dt %>% filter(expt == 2)

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


