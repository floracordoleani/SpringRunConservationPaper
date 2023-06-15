# Code created by Flora Cordoleani to:
# 1. Compare Mill, Deer and Butte creek spring-run adult abundances pre- and post- floodplain restoration
# 2. Compare Mill, Deer and Butte creek spring-run coefficient of variation pre- and post- floodplain restoration
# 3. Assess stock complex stability by estimating the combined coefficient of variation pre- and post- floodplain restoration
# 4. Compare Mill, Deer and Butte creek spring-run juvenile migratory strategies when leaving the natal reaches

library(dplyr)
library(ggplot2)
library(lubridate)
library(tidyverse)
library(viridis)
library(cowplot)
library(gridExtra)
library(ggpubr)
library(ggalt)
library(reshape2)

# Load and prepare spring-run escapement data ------------------------------------------------------------
### Prepare escapement data for figures and analysis
# Select years for analysis
year_vec <- data.frame(Year=c(1970:2019))

ButteEscap <-read.csv("Data/ButteAdultAbundance.csv")
MillEscap <- read.csv("Data/MillAdultAbundance.csv")
DeerEscap <- read.csv("Data/DeerAdultAbundance.csv")

ButteEscap_ally_withNA  <- merge(ButteEscap,year_vec,by="Year",all.y=TRUE) %>% 
                           mutate(PopID = replace_na(PopID,'Butte Creek'))

MillEscap_ally_withNA <- merge(MillEscap,year_vec,by="Year",all.y=TRUE) %>% 
  mutate(PopID = replace_na(PopID,'Mill Creek'))

DeerEscap_ally_withNA <- merge(DeerEscap,year_vec,by="Year",all.y=TRUE) %>% 
  mutate(PopID = replace_na(PopID,'Deer Creek'))

# Replace NA with zero values for Figure 2
ButteEscap_ally_withzero <- ButteEscap_ally_withNA %>% 
                            mutate(Total = replace_na(Total,0))

MillEscap_ally_withzero <- MillEscap_ally_withNA %>% 
  mutate(Total = replace_na(Total,0))

DeerEscap_ally_withzero <- DeerEscap_ally_withNA %>% 
  mutate(Total = replace_na(Total,0))

# Remove years with no estimates
ButteEscap_noNA <- ButteEscap %>%  
                  filter(Year != 1991)

MillEscap_noNA <- MillEscap %>%  
  filter(Year %in% c(1970:1975,1977,1978,1980,1982,1984:2019))

DeerEscap_noNA <- DeerEscap %>%  
  filter(Year %in% c(1970:1975,1977,1978,1980,1982,1983,1985:2019))

# Select years where Butte, Mill and Deer Creek populations have escapement estimates for the CV analysis
ButteEscap_noNA_CV <- ButteEscap %>%  
                      filter(Year %in% c(1970:1975,1977,1978,1980,1982,1985:1990,1992:2019))

MillEscap_noNA_CV <- MillEscap %>%  
                     filter(Year %in% c(1970:1975,1977,1978,1980,1982,1985:1990,1992:2019))

DeerEscap_noNA_CV <- DeerEscap %>%  
                     filter(Year %in% c(1970:1975,1977,1978,1980,1982,1985:1990,1992:2019))

# Combine population escapement estimates for figure 2 and stats estimation
MDBEscap_ally_withzero <- data.frame(rbind(DeerEscap_ally_withzero,MillEscap_ally_withzero,
                                           ButteEscap_ally_withzero))

MDBEscap_7094 <- data.frame(rbind(DeerEscap_noNA,MillEscap_noNA,ButteEscap_noNA)) %>% 
                 filter(Year %in% c(1970:1994))

MDBEscap_9519 <- data.frame(rbind(DeerEscap_noNA,MillEscap_noNA,ButteEscap_noNA)) %>% 
                 filter(Year %in% c(1995:2019))

# Combine population escapement estimates for pre and post-restoration CV calculation
BDEscap_CV <- merge(ButteEscap_noNA_CV,DeerEscap_noNA_CV,by="Year") %>% 
  mutate(Total = Total.x + Total.y,
         PopID = "Butte + Deer Creek") %>% 
  select(PopID, Year,Total) 

BMEscap_CV <- merge(ButteEscap_noNA_CV,MillEscap_noNA_CV,by="Year") %>% 
  mutate(Total = Total.x + Total.y,
         PopID = "Butte + Mill Creek") %>% 
  select(PopID, Year,Total) 

MDEscap_CV <- merge(MillEscap_noNA_CV,DeerEscap_noNA_CV,by="Year") %>% 
  mutate(Total = Total.x + Total.y,
         PopID = "Mill + Deer Creek") %>% 
  select(PopID, Year,Total) 

MDBEscap_CV <- merge(MDEscap_CV[,c('PopID','Year','Total')],
                          ButteEscap_noNA_CV,by="Year") %>% 
  mutate(PopID = "Mill + Deer + Butte Creek") %>% 
  mutate(Total = Total.x + Total.y) %>% 
  select(c(PopID,Year,Total)) 

MDBEscap_7094_CV <- data.frame(rbind(ButteEscap_noNA_CV,MillEscap_noNA_CV,DeerEscap_noNA_CV,
                                     BDEscap_CV,BMEscap_CV,MDEscap_CV,MDBEscap_CV)) %>% 
                               filter(Year %in% c(1970:1994))

MDBEscap_9519_CV <- data.frame(rbind(ButteEscap_noNA_CV,MillEscap_noNA_CV,DeerEscap_noNA_CV,
                                     BDEscap_CV,BMEscap_CV,MDEscap_CV,MDBEscap_CV)) %>% 
                                filter(Year %in% c(1995:2019))


# Spring-run population escapement statistics --------------------------------------------------------------
MDBEscap_7094_stats <- MDBEscap_7094 %>% 
                       group_by(PopID) %>% 
                       dplyr::summarise(Mean = mean(Total),
                                        Var= var(Total),
                                        Sd = sd(Total))

MDBEscap_9519_stats <- MDBEscap_9519 %>% 
                       group_by(PopID) %>% 
                       dplyr::summarise(Mean = mean(Total),
                                        Var= var(Total),
                                        Sd = sd(Total))

MDBEscap_7094_CV_stats <- MDBEscap_7094_CV %>% 
                          group_by(PopID) %>% 
                          dplyr::summarise(Mean = mean(Total),
                                            Var= var(Total),
                                            Sd = sd(Total),
                                            CV = Sd/Mean)

MDBEscap_9519_CV_stats <- MDBEscap_9519_CV %>% 
                          group_by(PopID) %>% 
                          dplyr::summarise(Mean = mean(Total),
                                           Var= var(Total),
                                           Sd = sd(Total),
                                           CV = Sd/Mean)

# Estimate Correlation coefficients between pops
Mill_7094 <- MDBEscap_7094_CV %>% filter(PopID=="Mill Creek")
Deer_7094 <- MDBEscap_7094_CV %>% filter(PopID=="Deer Creek")
Butte_7094 <- MDBEscap_7094_CV %>% filter(PopID=="Butte Creek")

(MD_cor_7094 <- cor(Deer_7094$Total,Mill_7094$Total,method="pearson"))
(MB_cor_7094 <- cor(Butte_7094$Total,Mill_7094$Total,method="pearson"))
(BD_cor_7094 <- cor(Deer_7094$Total,Butte_7094$Total,method="pearson"))


Mill_9519 <- MDBEscap_9519_CV %>% filter(PopID=="Mill Creek")
Deer_9519 <- MDBEscap_9519_CV %>% filter(PopID=="Deer Creek")
Butte_9519 <- MDBEscap_9519_CV %>% filter(PopID=="Butte Creek")

(MD_cor_9519 <- cor(Deer_9519$Total,Mill_9519$Total,method="pearson"))
(MB_cor_9519 <- cor(Butte_9519$Total,Mill_9519$Total,method="pearson"))
(BD_cor_9519 <- cor(Deer_9519$Total,Butte_9519$Total,method="pearson"))

# Figures 2 and 6 ----------------------------------------------------------------
# Figure 2
(EscapStackPlot <- ggplot(data=MDBEscap_ally_withzero,aes(x=Year,y=Total,fill=PopID),size=2)+
   geom_area(color = "black")+
   annotate("text", x=1976.3, y=1000, label= "*") + 
   annotate("text", x=1979, y=1000, label= "*") +
   annotate("text", x=1981, y=1000, label= "*") +
   annotate("text", x=1983.1, y=1000, label= "*") +    
   annotate("text", x=1984, y=1000, label= "*") +
   annotate("text", x=1991, y=1000, label= "*") +
   theme_bw()+
   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
         panel.background = element_blank(), axis.line = element_line(colour = "black"),
         axis.text.x = element_text(angle=45, hjust = 1),
         legend.position = 'top',
         text = element_text(size=18),legend.title = element_blank())+
   ylab('Escapement')+ xlab("")+
   scale_x_continuous(breaks = seq(1970,2019,by = 4),limits=c(1970,2019))+
   scale_fill_viridis_d()
)

png("Figures/Fig2.png", 
    family = "serif", width = 10, height= 6, units = "in", res =300)

EscapStackPlot

dev.off()

# Figure 6
Escap_CV_stats <- merge(MDBEscap_7094_CV_stats,MDBEscap_9519_CV_stats,by="PopID")
colnames(Escap_CV_stats) <- c("PopID","Mean","Var","Sd","CV","Mean_end","Var_end","Sd_end","CV_end")

Escap_CV_stats$PopID <- factor(Escap_CV_stats$PopID ,
                            levels = c("Deer Creek", "Mill Creek","Butte Creek",
                                       "Mill + Deer Creek","Butte + Deer Creek",
                                       "Butte + Mill Creek","Mill + Deer + Butte Creek"))

Escap_CV_stats_melt <- melt(Escap_CV_stats[,c("PopID","CV","CV_end")],id="PopID")

(EscapCVPlot  <- ggplot(Escap_CV_stats_melt,aes(x = PopID, y = value)) +
    geom_line()+
    geom_point(aes(x = PopID, y = value,shape=variable,color=variable), size = 3)+
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          text = element_text(size=16),legend.title = element_blank(),
          axis.title.y = element_text(margin = unit(c(0, 4, 0, 0), "mm")),
          legend.position="top")+
    scale_shape_manual(values = c(16,15),labels=c('Pre-Butte Creek restoration (1970-1994)', 
                                                  'Post-Butte Creek restoration (1995-2019)'))+
    scale_color_manual(values = c("grey60", "black"),
                       labels=c('Pre-Butte Creek restoration (1970-1994)', 
                        'Post-Butte Creek restoration (1995-2019)'))+
    ylab('CV')+ xlab("")+
    
    scale_x_discrete(labels=c('Deer', 'Mill', 'Butte',
                                  'M + D','B + D','B + M', 'M + D + B'))
)

Escap_Abund_stats_melt <- melt(Escap_CV_stats[,c("PopID","Mean","Mean_end")],id="PopID")

(EscapAbundPlot  <- ggplot(Escap_Abund_stats_melt,aes(x = PopID, y = value)) +
    geom_line()+
    geom_point(aes(x = PopID, y = value,shape=variable,color=variable), size = 3)+
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.title.y = element_text(margin = unit(c(0, 4, 0, 0), "mm")),
          # axis.text.x = element_text(angle=45, hjust = 1),
          text = element_text(size=16),legend.title = element_blank(),
          legend.position="top")+
    scale_shape_manual(values = c(16,15),labels=c('Pre-Butte Creek restoration (1970-1994)', 
                                                  'Post-Butte Creek restoration (1995-2019)'))+
    scale_color_manual(values = c("grey60", "black"),
                       labels=c('Pre-Butte Creek restoration (1970-1994)', 
                                'Post-Butte Creek restoration (1995-2019)'))+
    ylab('Adult Abundance')+ xlab("")+
    scale_x_discrete(labels=c('Deer', 'Mill', 'Butte',
                              'M + D','B + D','B + M', 'M + D + B'))
)

png("figures/Fig6.png", 
    family = "serif", width = 8, height= 8, units = "in", res =300)

ggarrange(EscapCVPlot,EscapAbundPlot,
          nrow=2, common.legend=TRUE,
          labels=c('(a)','(b)'),font.label=list(size=18),
          align="v")

dev.off()

# Load RST data  ---------------------------------------------------
### Mill/Deer Creek RST
RSTdata_MDC <- read.csv('Data/RSTChinMillDeer.csv',header=T)

RSTdata_MDC$MonthDay <- format(as.Date(RSTdata_MDC$Date), format="%m-%d")

RSTdata_MDC  <- RSTdata_MDC %>%
  mutate(Type = case_when(Month > 9 & Length > 50 |
                            Month <= 2 & Length > 60 |
                            Month == 3 & Length > 76 |
                            Month == 4 & Day >= 1 & Day <15 & Length > 85 |
                            Month == 4 & Day >= 15 & Day <30 & Length > 95 |
                            Month >= 4 & Month < 7 & Length > 100 ~ 'Yearling',
                          TRUE ~ 'YoY'))

# Sample size estimate for each defined life history strategy
RST_MDC.summary <- RSTdata_MDC %>% group_by(Type) %>% 
  dplyr::summarise(count=n())%>% 
  ungroup() %>% 
  mutate(prop=count/sum(count)*100)

### Butte Creek RST
RSTdata_BC  <- read.csv('Data/RSTChinButteCreek.csv',header=T)

## Butte Creek report RST figure 
RSTdata_BC$MonthDay <- format(as.Date(RSTdata_BC$Date), format="%m-%d")

RSTdata_BC  <- RSTdata_BC %>%
  mutate(Type = case_when(Month > 8 & ForkLength > 60 |
                            Month <= 2 & ForkLength > 79 |
                            Month == 3 & ForkLength > 95 |
                            # Month == 4 & Day >= 1 & Day <15 & Length > 85 |
                            Month == 4  & ForkLength > 110 |
                            Month >= 5 & Month < 8 & ForkLength > 110 ~ 'Yearling',
                          TRUE ~ 'YoY'))

RSTdata_BC <- RSTdata_BC %>% 
  filter(!is.na(Type))

# Remove fish that do not belong to the juvenile class
RSTdata_BC <- RSTdata_BC %>%
  filter(ForkLength > 0 & ForkLength < 250)

# Sample size estimate for each defined life history strategy
RST_BC.summary <- RSTdata_BC %>% group_by(Type) %>% 
  dplyr::summarise(count=n()) %>% 
  ungroup() %>% 
  mutate(prop=count/sum(count)*100)

# Figure S4 --------------------------------------------------------------------------------
# Mill and Deer Creek figure
date_ord <-  seq(as.Date("2001-09-01"), as.Date("2002-07-30"), by="days")
date_ord <- format(date_ord, format="%m-%d")

(pRST_MDC <- ggplot(RSTdata_MDC,aes(x=factor(MonthDay,levels = date_ord),
                                    y=Length, colour = Type))+
    geom_point(alpha=0.2)+xlab("")+ylab("Length (mm)")+theme_bw()+
    ggtitle("Mill/Deer Creek")+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))+
    scale_x_discrete(breaks=c("09-01","10-01","11-01" ,"12-01","01-01","02-01","03-01","04-01","05-01","06-01"),
                     labels=c("September","October","November","December","January","February","March", "April","May","June"))+
    scale_y_continuous(breaks=seq(20,160,20),limits = c(20, 160))+
    theme(axis.text.x = element_text(angle=45, hjust = 1),
          text = element_text(size=18),
          legend.position = "none")+
    scale_color_manual(values = c("mediumpurple","goldenrod2"))
)

# Butte Creek figure
date_ord <-  seq(as.Date("2001-09-01"), as.Date("2002-07-30"), by="days")
date_ord <- format(date_ord, format="%m-%d")

(pRST_BC <- ggplot(RSTdata_BC,aes(x=factor(MonthDay,levels = date_ord),
                                  y=ForkLength, colour = Type))+
    geom_point(alpha=0.2)+xlab("")+ylab("")+theme_bw()+
    ggtitle("Butte Creek")+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))+
    scale_x_discrete(breaks=c("09-05","10-03","11-01" ,"12-01","01-01","02-01","03-01","04-01","05-01","06-01","07-01"),
                     labels=c("September","October","November","December","January","February","March", "April","May","June","July"))+
    scale_y_continuous(breaks=seq(20,160,20),limits = c(20, 160))+
    theme(axis.text.x = element_text(angle=45, hjust = 1),
          text = element_text(size=18),
          legend.position = "none")+
    scale_color_manual(values = c("mediumpurple","goldenrod2"))
)

# Mill/Deer and Butte Creek combined
(pRST_MDB <- ggarrange(pRST_MDC,pRST_BC,labels=c("(a)","(b)"),nrow=1,
                       font.label=list(size =18,color="black"))
)

png("Figures/FigS4.png", 
    family = "serif", width = 10, height= 5, units = "in", res =300)

pRST_MDB

dev.off()

