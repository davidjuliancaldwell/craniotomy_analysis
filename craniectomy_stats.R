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
library(emmeans)
library(here)

here()

data_file_excel <- "Kempe-CTScans_DATA_2020-04-23_1712.csv"
crani_e_data <- read.csv(file=data_file_excel, header=TRUE, sep=",",na.strings=c(""))


crani_e_data$abs_diff <- crani_e_data$crani_size_inner - crani_e_data$crani_size_outer
cols <- c("crani_type2","crani_type","pathology","mechanism")
crani_e_data[cols] <- lapply(crani_e_data[cols], factor)

summary_data <- crani_e_data %>% group_by(crani_type2) %>% summarize(min_diff=min(abs_diff,na.rm = TRUE),
                                                                     max_diff = max(abs_diff,na.rm = TRUE),
                                                                     min_perc=min(crani_dif,na.rm = TRUE),
                                                                     max_perc=max(crani_dif,na.rm = TRUE),
                                                                     med_diff=median(abs_diff,na.rm = TRUE),
                                                                     med_perc=median(crani_dif,na.rm = TRUE))

data_model_ancova <- crani_e_data %>% anova_test(crani_size_inner)

data_model <- lm(crani_size_outer ~ crani_size_inner + crani_type2 + age,data=crani_e_data)
summary(data_model)

pwc <- emmeans(data_model,specs=pairwise ~crani_type2:crani_size_inner)

plot1 <- ggplot(crani_e_data,aes(x=abs_diff,group=crani_type2,fill=crani_type2)) +
  geom_histogram(position="dodge") + 
  labs(x="Absolute difference in mm^3",y="Count",title = "Distribution of Volume Differences") +
  xlim(0,max(crani_e_data$abs_diff)) + 
  theme_bw() + 
  scale_fill_discrete(name="Type of\nSurgery",breaks=c("0","1"),labels=c("Reverse ? Mark","Kempe"))
print(plot1)


plot2 <- ggplot(crani_e_data,aes(x=crani_dif,group=crani_type2,fill=crani_type2)) +
  geom_histogram(position="dodge") + 
  labs(x="Percent Difference",y="Count",title = "Distribution of Percent Differences") +
  xlim(0,max(crani_e_data$crani_dif)) + 
  theme_bw()  + 
  scale_fill_discrete(name="Type of\nSurgery",breaks=c("0","1"),labels=c("Reverse ? Mark","Kempe"))
print(plot2)

plot3 <- ggplot(crani_e_data,aes(x=crani_size_inner,group=crani_type2,fill=crani_type2)) +
  geom_histogram(position="dodge") + 
  labs(x="Initial volume",y="Count",title = "Distribution of Starting Volumes") +
  xlim(0,max(crani_e_data$crani_dif)) + 
  theme_bw()  + 
  scale_fill_discrete(name="Type of\nSurgery",breaks=c("0","1"),labels=c("Reverse ? Mark","Kempe"))
print(plot3)
