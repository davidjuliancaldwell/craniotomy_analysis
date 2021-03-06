library(tidyr)
library(RPostgreSQL)
library(plyr)
library(dplyr)
library(stringr)
library(ggplot2)
library(ggExtra)
library(ggsci)
library(gridExtra)
library(cowplot)
library(emmeans)
library(here)
library(CBCgrps)
library(tidyverse)

here()

save_plot = TRUE

data_file_excel <- "kempe_data_josh.csv"
crani_e_data <- read.csv(file=data_file_excel, header=TRUE, sep=",",na.strings=c(""))

crani_e_data <- crani_e_data %>% distinct(crani_type,crani_type2) %>% rowid_to_column("unique_id") %>% left_join(crani_e_data)
crani_e_data$unique_id <- as.factor(crani_e_data$unique_id)
crani_e_data$dos <- as.Date(crani_e_data$dos)

crani_e_data <- crani_e_data %>% filter(dos >= "2015-01-01")

plot1 <- ggplot(crani_e_data,aes(x=dos,y=unique_id,shape=unique_id)) +
  geom_jitter(width = 0.2, height = 0.2) + 
  labs(x="Date of Surgery",y="Surgical Approach",title = "Penetrance of Kempe Incision") +
  scale_x_date(date_breaks="1 year",date_labels="%b-%Y") + 
  scale_shape_discrete(name="Type of\nSurgery",labels=c("Kempe-CraniE","Kempe-Crani","RQM-Crani","RQM-CraniE")) +
  theme_bw() + 
  removeGridY() + 
  theme(axis.ticks.y = element_blank(),
                     axis.text.y = element_blank(),
                     axis.text.x = element_text(angle=45,vjust=0.5),
        text = element_text(family="sans")) + 
  theme(legend.title=element_blank())
print(plot1)

if (save_plot){
    ggsave("time_series_jitter_revision.png", units="in", width=7, height=4, dpi=600)
  ggsave("time_series_jitter_revision.svg", units="in", width=7, height=4, dpi=600)
  
}
