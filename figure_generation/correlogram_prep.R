setwd("")
library(tidyverse)
library(viridis)
library(reshape2)
library(corrplot)
theme_set(theme_classic())

misincorporation_data <- read_csv("collated_misincorporation.csv")
#manually collated csv containing average misincorporation
#at every position of each tRNA across all experiments
misincorporation_data <- na.omit(misincorporation_data)
misincorporation_correlation <- cor(misincorporation_data[a, b])
#select only data containing columns

plot.new()
pdf("misincorporation_correlogram.pdf")
corrplot(
    misincorporation_correlation,
    type = "lower",
    order = "hclust",
    addCoef.col = "black",
    col = viridis(100, option = "mako")
    )
dev.off()

termination_data <- read_csv("collated_termination.csv")
#manually collated csv containing average termination
#at every position of each tRNA across all experiments
termination_data <- na.omit(termination_data)
termination_correlation <- cor(termination_data[a, b])
#select only data containing columns
plot.new()
pdf("termination_correlogram.pdf")
corrplot(
    termination_correlation,
    type = "lower",
    order = "hclust",
    addCoef.col = "black",
    col = viridis(100, option = "mako")
    )
dev.off()