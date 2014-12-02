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
prelim.dat <- read.table("../data/D3_Ras_data.csv", header=TRUE, sep=",")
prelim.df <- tbl_df(prelim.dat)

gathered <- prelim.df %>% gather(sample, ion.count, -compound)
separated <- gathered %>% separate(sample, into = c("cell.line", "condition", "replicate"), sep = "_")

common.data <- separated %>% filter(!(cell.line %in% c("D3","Ras"))) %>% filter(condition != "blank") %>% 
  mutate(name = paste(cell.line, condition, sep=".")) %>% arrange(compound, name, replicate)
aa40 <- c("Lysine","Histidine","Phenylalanine","Threonine","Valine","Glutamine", "Serine", "Tyrosine")
common.data40 <- common.data %>% filter(compound %in% aa40)

common.means <- common.data40 %>% group_by(compound, name) %>% 
  summarize(length=n(), mean = mean(ion.count), stderr = sqrt(var(ion.count))/length, se.low = mean-stderr, se.high = mean+stderr)

common.plot <- ggplot(common.means, aes(x=compound, y=mean, fill=name)) + geom_bar(stat="identity", position="dodge", color="black") +
  geom_abline(slope=0, lwd=1) + geom_errorbar(aes(ymin=se.low, ymax=se.high), position=position_dodge(width=.9), width=.5)
common.plot

common.dat.tomerge <- common.data40 %>% filter(cell.line == "fresh40B")
common.dat.tomerge$condition = "fresh"
common.dat.tomerge <- common.dat.tomerge %>% data.frame(type=c("D3","Ras")) %>% group_by(type) %>% do(common.dat.tomerge)
common.dat.tomerge <- common.dat.tomerge %>% mutate(cell.line = type) %>% ungroup %>% select(-type, -name)

d3.ras.data <- separated %>% filter(cell.line %in% c("D3","Ras")) %>% filter(compound %in% aa40)

compiled.dat <- rbind(common.dat.tomerge, d3.ras.data)
temp.dat <- compiled.dat %>% filter(condition=="fresh") %>% group_by(compound, cell.line) %>% summarize(divide.by = mean(ion.count))
compiled.dat <- merge(compiled.dat, temp.dat)
compiled.dat <- compiled.dat %>% mutate(ion.count = ion.count/divide.by) %>% select(-divide.by)

helperfunc <- function (x) {
  bla <- 
}
compiled.dat <- compiled.dat %>% mutate(blah = ifelse(condition == fresh, 1, 0)) %>% group_by(compound, cell.line) %>%  norm=())

compiled.dat %>% group_by(compound, cell.line) %>% mutate(mean = mean(ion.count[condition == "fresh"])) %>% arrange(condition)


compiled.dat <- compiled.dat %>% group_by(compound, cell.line, condition) %>%
  summarize(length=n(), mean = mean(ion.count), stderr = sqrt(var(ion.count))/length, se.low = mean-stderr, se.high = mean+stderr)

compiled.dat$condition <- factor(compiled.dat$condition, level=c("fresh","1.25","1.5","2","2.5","3"))

plot1 <- ggplot(compiled.dat %>% filter(cell.line=="D3"), aes(x=compound, y=mean, fill=condition)) + geom_bar(stat="identity", position="dodge", color="black") +
  geom_abline(slope=0, lwd=1) + geom_errorbar(aes(ymin=se.low, ymax=se.high), position=position_dodge(width=.9), width=.5)

plot1 + theme(legend.position="right")





