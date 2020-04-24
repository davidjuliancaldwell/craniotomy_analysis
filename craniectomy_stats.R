library(tidyr)
library(RPostgreSQL)
library(plyr)
library(dplyr)
library(stringr)
library(ggplot2)
library(ggsci)
library(gridExtra)
library(cowplot)
library(emmeans)
library(here)
library(CBCgrps)

here()

# if wanting to select only good scans 
select_subset = FALSE
save_plot = TRUE

data_file_excel <- "Kempe-CTScansWithDemograph_DATA_2020-04-24_1621.csv"
crani_e_data <- read.csv(file=data_file_excel, header=TRUE, sep=",",na.strings=c(""))


crani_e_data$abs_diff <- crani_e_data$crani_size_inner - crani_e_data$crani_size_outer
cols <- c("crani_type2","crani_type","pathology","mechanism")
crani_e_data[cols] <- lapply(crani_e_data[cols], factor)

crani_e_data <- crani_e_data %>% mutate(ebl_cplasty = as.numeric(ebl_cplasty),
                                        or_time = as.numeric(or_time),
                        ebl = as.numeric(ebl),
                        surgery_duration = as.numeric(surgery_duration),
                        daysto_copasty = as.numeric((daysto_copasty)))

if (select_subset){
crani_e_data <- crani_e_data %>% filter(thin_ct == 1)
} else{
 # crani_e_data <- crani_e_data %>% filter(!is.na(thin_ct))
}

twogrps(crani_e_data,"crani_type2",workspace=2e7,ShowStatistic = TRUE)

# force kruskal on cplasty ebl 
kruskal.test(ebl_cplasty ~ crani_type2,data = crani_e_data)

summary_data <- crani_e_data %>% group_by(crani_type2) %>% summarize(min_diff=min(abs_diff,na.rm = TRUE),
                                                                     max_diff = max(abs_diff,na.rm = TRUE),
                                                                     min_perc=min(crani_dif,na.rm = TRUE),
                                                                     max_perc=max(crani_dif,na.rm = TRUE),
                                                                     med_diff=median(abs_diff,na.rm = TRUE),
                                                                     med_perc=median(crani_dif,na.rm = TRUE),
                                                                     std_diff=sd(abs_diff,na.rm=TRUE),
                                                                     std_perc=sd(crani_dif,na.rm=TRUE),
                                                                     mean_diff=mean(abs_diff,na.rm=TRUE),
                                                                     mean_perc=mean(crani_dif,na.rm=TRUE))

summary_data_demo <- crani_e_data %>% group_by(crani_type2) %>% summarise(med =sprintf("%0.1f",mean(ebl_cplasty,na.rm=TRUE)),std = sprintf("%0.1f",sd(ebl_cplasty,na.rm=TRUE)))


data_model <- lm(crani_size_outer ~ crani_size_inner + crani_type2 + age,data=crani_e_data)
summary(data_model)

data_model_abs_diff <- lm(abs_diff ~ crani_type2,data=crani_e_data)
summary(data_model_abs_diff)

pwc <- emmeans(data_model,specs=pairwise ~crani_type2:crani_size_inner)

plot1 <- ggplot(crani_e_data,aes(x=abs_diff,group=crani_type2,fill=crani_type2)) +
  geom_histogram(position="dodge") + 
  labs(x="Absolute difference in mm^3",y="Count",title = "Distribution of Volume Differences") +
  xlim(0,max(crani_e_data$abs_diff)) + 
  theme_bw() + 
  scale_fill_discrete(name="Type of\nSurgery",breaks=c("0","1"),labels=c("Reverse ? Mark","Kempe"))
print(plot1)

if (save_plot){
  if (select_subset){
  ggsave("abs_diffs_subset.png", units="in", width=5, height=4, dpi=600)
  }
  else{
    ggsave("abs_diffs.png", units="in", width=5, height=4, dpi=600)
  }
}



plot2 <- ggplot(crani_e_data,aes(x=crani_dif,group=crani_type2,fill=crani_type2)) +
  geom_histogram(position="dodge") + 
  labs(x="Percent Difference",y="Count",title = "Distribution of Percent Differences") +
  xlim(0,max(crani_e_data$crani_dif)) + 
  theme_bw()  + 
  scale_fill_discrete(name="Type of\nSurgery",breaks=c("0","1"),labels=c("Reverse ? Mark","Kempe"))
print(plot2)

if (save_plot){
  if (select_subset){
    ggsave("percents_subset.png", units="in", width=5, height=4, dpi=600)
  }
  else{
    ggsave("percents.png", units="in", width=5, height=4, dpi=600)
  }
}


plot3 <- ggplot(crani_e_data,aes(x=crani_size_inner,group=crani_type2,fill=crani_type2)) +
  geom_histogram(position="dodge") + 
  labs(x="Initial volume",y="Count",title = "Distribution of Starting Volumes") +
  xlim(0,max(crani_e_data$crani_dif)) + 
  theme_bw()  + 
  scale_fill_discrete(name="Type of\nSurgery",breaks=c("0","1"),labels=c("Reverse ? Mark","Kempe"))
print(plot3)

if (save_plot){
  if (select_subset){
    ggsave("distrib_starts_subset.png", units="in", width=5, height=4, dpi=600)
  }
  else{
    ggsave("distrib_starts.png", units="in", width=5, height=4, dpi=600)
  }
}


plot4 <- ggplot(crani_e_data,aes(x=crani_size_inner,y=crani_size_outer,color=factor(crani_type2,labels=c("Reverse ? Mark","Kempe")))) +
  geom_point() + 
  geom_smooth(method=lm) +
  labs(x="Initial volume (mm^3)",y="Final Volume (mm^3)",title = "Final vs. Starting Volume by Surgery",color="Type of\nSurgery") +
  xlim(min(crani_e_data$crani_size_inner,crani_e_data$crani_size_outer,na.rm = TRUE),max(crani_e_data$crani_size_inner,crani_e_data$crani_size_outer,na.rm = TRUE)) +
  ylim(min(crani_e_data$crani_size_inner,crani_e_data$crani_size_outer,na.rm = TRUE),max(crani_e_data$crani_size_inner,crani_e_data$crani_size_outer,na.rm = TRUE)) +
  theme_bw()
print(plot4)

if (save_plot){
  if (select_subset){
    ggsave("scatter_vols_subset.png", units="in", width=5, height=4, dpi=600)
  }
  else{
    ggsave("scatter_vols.png", units="in", width=5, height=4, dpi=600)
  }
}


