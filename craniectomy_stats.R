library(tidyr)
library(RPostgreSQL)
library(plyr)
library(dplyr)
library(stringr)
library(lubridate)
library(ggplot2)
library(ggsci)
library(gridExtra)
library(cowplot)
library(here)

dataFile <- "C:/Users/david/UW/Ryan Kellogg - Kempe files/volume_overview_threshold.csv"
handCurated <- read.csv(file=dataFile, header=TRUE, sep=",",na.strings=c(""))


plot1 <- ggplot(handCurated,aes(x=Absolute.difference)) +
  geom_histogram(fill="blue") + 
  labs(x="Absolute difference in mm^3",y="Count",title = "Distribution of Volume Differences") +
  xlim(0,max(handCurated$Absolute.difference))
print(plot1)


plot2 <- ggplot(handCurated,aes(x=X..Dif)) +
  geom_histogram(fill="blue") + 
  labs(x="% Difference from Pre-op",y="Count",title = "Distribution of % Differences from Pre-op") +
  xlim(0,max(handCurated$X..Dif))
print(plot2)