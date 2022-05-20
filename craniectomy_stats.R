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
select_subset = TRUE
save_plot = FALSE

data_file_excel <- "Kempe-CTScansWithDemograph_DATA_2020-05-18_1359.csv"
crani_e_data <- read.csv(file=data_file_excel, header=TRUE, sep=",",na.strings=c(""))

#crani_e_data <- crani_e_data %>% filter(record_id != 125 & record_id != 119)
crani_e_data <- crani_e_data %>% filter(record_id != 125)


crani_e_data$abs_diff <- crani_e_data$crani_size_inner - crani_e_data$crani_size_outer
cols <- c("crani_type2","crani_type","pathology","mechanism")
crani_e_data[cols] <- lapply(crani_e_data[cols], factor)

crani_e_data <- crani_e_data %>% mutate(ebl_cplasty = as.numeric(ebl_cplasty),
                                        or_time = as.numeric(or_time),
                        ebl = as.numeric(ebl),
                        surgery_duration = as.numeric(surgery_duration),
                        daysto_copasty = as.numeric((daysto_copasty)),
                        sex = as.factor(sex))

if (select_subset){
crani_e_data <- crani_e_data %>% filter(thin_ct == 1)
} else{
  crani_e_data <- crani_e_data %>% filter(!is.na(thin_ct))
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


data_model <- lm(crani_size_outer ~ crani_size_inner + crani_type2 + sex + age ,data=crani_e_data)
summary(data_model)

data_model_abs_diff <- lm(abs_diff ~ crani_type2,data=crani_e_data)
summary(data_model_abs_diff)

data_final_only <- lm(crani_size_outer ~ crani_type2,data=crani_e_data)
summary(data_final_only)

data_initial_only <-  lm(crani_size_inner ~ crani_type2,data=crani_e_data)
summary(data_initial_only)

a<- kruskal.test(abs_diff ~ crani_type2,data=crani_e_data)

data_model.emm <- emmeans(data_model,"crani_type2")
pwpm(data_model.emm)
pwpp(data_model.emm)

plot1 <- ggplot(crani_e_data,aes(x=abs_diff,group=crani_type2,fill=crani_type2)) +
  geom_histogram(position="stack",bins=20) + 
  labs(x="Absolute difference in mm^3",y="Count",title = "Distribution of Volume Differences") +
  coord_cartesian(xlim=c(0,max(crani_e_data$abs_diff)+1000))  + 
  theme_bw() + 
  scale_fill_discrete(name="Type of\nSurgery",breaks=c("0","1"),labels=c("Reverse ? Mark","Kempe"))
print(plot1)

if (save_plot){
  if (select_subset){
  ggsave("abs_diffs_subset_v2.png", units="in", width=5, height=4, dpi=600)
  }
  else{
    ggsave("abs_diffs_v2.png", units="in", width=5, height=4, dpi=600)
  }
}

plot2 <- ggplot(crani_e_data,aes(x=crani_dif,group=crani_type2,fill=crani_type2)) +
  geom_histogram(position="stack",bins=20) + 
  labs(x="Percent Difference",y="Count",title = "Distribution of Percent Differences") +
  coord_cartesian(xlim=c(0,max(crani_e_data$crani_dif)+5))  + 
  theme_bw()  + 
  scale_fill_discrete(name="Type of\nSurgery",breaks=c("0","1"),labels=c("Reverse ? Mark","Kempe"))
print(plot2)

if (save_plot){
  if (select_subset){
    ggsave("percents_subset_v2.png", units="in", width=5, height=4, dpi=600)
  }
  else{
    ggsave("percents_v2.png", units="in", width=5, height=4, dpi=600)
  }
}


plot3 <- ggplot(crani_e_data,aes(x=crani_size_inner,group=crani_type2,fill=crani_type2)) +
  geom_histogram(position="stack",bins=10) + 
  labs(x="Initial volume",y="Count",title = "Distribution of Starting Volumes") +
  coord_cartesian(xlim=c(0,max(crani_e_data$crani_size_inner)+1000))  + 
  theme_bw()  + 
  scale_fill_discrete(name="Type of\nSurgery",breaks=c("0","1"),labels=c("Reverse ? Mark","Kempe"))
print(plot3)

if (save_plot){
  if (select_subset){
    ggsave("distrib_starts_subset_v2.png", units="in", width=5, height=4, dpi=600)
  }
  else{
    ggsave("distrib_starts_v2.png", units="in", width=5, height=4, dpi=600)
  }
}


plot4 <- ggplot(crani_e_data,aes(x=crani_size_inner,y=crani_size_outer,color=factor(crani_type2,labels=c("Reverse ? Mark","Kempe")))) +
  geom_point() + 
  geom_smooth(method=lm) +
  labs(x="Initial volume (mm^3)",y="Final Volume (mm^3)",title = "Final vs. Starting Volume by Surgery",color="Type of\nSurgery") +
  coord_cartesian(xlim=c(min(crani_e_data$crani_size_inner,crani_e_data$crani_size_outer,na.rm = TRUE),max(crani_e_data$crani_size_inner,crani_e_data$crani_size_outer,na.rm = TRUE)),ylim=c(min(crani_e_data$crani_size_inner,crani_e_data$crani_size_outer,na.rm = TRUE),max(crani_e_data$crani_size_inner,crani_e_data$crani_size_outer,na.rm = TRUE)))  + 
  theme_bw()
print(plot4)

if (save_plot){
  if (select_subset){
    ggsave("scatter_vols_subset_v2.png", units="in", width=5, height=4, dpi=600)
  }
  else{
    ggsave("scatter_vols_v2.png", units="in", width=5, height=4, dpi=600)
  }
}

plot5 <- ggplot(crani_e_data,aes(x=crani_size_outer,group=crani_type2,fill=crani_type2)) +
  geom_histogram(position="stack",bins=10) + 
  labs(x="Final volume",y="Count",title = "Distribution of Final Volumes") +
  coord_cartesian(xlim=c(0,max(crani_e_data$crani_size_outer)+1000))  + 
  theme_bw()  + 
  scale_fill_discrete(name="Type of\nSurgery",breaks=c("0","1"),labels=c("Reverse ? Mark","Kempe"))
print(plot5)

if (save_plot){
  if (select_subset){
    ggsave("distrib_final_subset_v2.png", units="in", width=5, height=4, dpi=600)
  }
  else{
    ggsave("distrib_final_v2.png", units="in", width=5, height=4, dpi=600)
  }
}

plot6 <- ggplot(crani_e_data,aes(x=crani_size_inner,group=crani_type2,fill=crani_type2)) +
  geom_histogram(position="stack",bins=20) + 
  labs(x="Absolute difference in mm^3",y="Count",title = "Distribution of Volume Differences") +
#  coord_cartesian(xlim=c(0,max(crani_e_data$abs_diff)+1000))  + 
  theme_bw() + 
  scale_fill_discrete(name="Type of\nSurgery",breaks=c("0","1"),labels=c("Reverse ? Mark","Kempe")) +
  facet_grid(vars(sex))
print(plot6)

if (save_plot){
  if (select_subset){
    ggsave("distrib_gender_subset_v2.png", units="in", width=5, height=4, dpi=600)
  }
  else{
    ggsave("distrib_gender_v2.png", units="in", width=5, height=4, dpi=600)
  }
}



