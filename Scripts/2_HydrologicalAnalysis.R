# Code created by Flora Cordoleani to: 
# 1. Estimate the relationship between floodplain inundation proportion/Sacramento River flow/Delta flow 
# and Butte, Mill and Deer Creek spring-run adult abundances
# 2. Compare Butte, Mill, and Deer Creek stream temperatures in both spawning and rearing grounds 
# 3. Estimate the number of days water from the Sacramento River was flowing into Butte Creek floodplain

library(dplyr)
library(ggplot2)
library(lubridate)
library(tidyverse)
library(ggpubr)
library(reshape2)

# Flow/Inundation and adult abundance relationship ------------------------------------------------------------------
### Load escapement data and estimate lagged abundances
ButteEscap <- read.csv("Data/ButteAdultAbundance.csv")

ButteEscap_withlag <- ButteEscap %>% 
                      mutate(LagTotal=lead(Total,2),
                             PopID="Butte Creek") %>% 
                      dplyr::select(PopID,Year,Total,LagTotal) %>% 
                      filter(Year %in% c(2001:2020))

MillEscap <- read.csv("Data/MillAdultAbundance.csv")

MillEscap_withlag <- MillEscap %>% 
                     mutate(LagTotal=lead(Total,2),
                            PopID="Mill Creek") %>% 
                     dplyr::select(PopID,Year,Total,LagTotal) %>% 
                     filter(Year %in% c(2001:2020))


DeerEscap <- read.csv("Data/DeerAdultAbundance.csv")

DeerEscap_withlag <- DeerEscap %>% 
                     mutate(LagTotal=lead(Total,2),
                            PopID="Deer Creek") %>% 
                     dplyr::select(PopID,Year,Total,LagTotal) %>% 
                     filter(Year %in% c(2001:2020))

# Combine escapement data
Escap_withlag <- data.frame(rbind(ButteEscap_withlag,MillEscap_withlag,DeerEscap_withlag))

### Load and clean Butte Creek floodplain stage/inundation data
FloodplainStage <- read.csv('Data/BSLStage.csv')  

FloodplainStage_0218 <- FloodplainStage %>% 
                        group_by(Year,Month,Day) %>% 
                        dplyr::summarise(MeanStage=mean(VALUE,na.rm=TRUE)) %>% 
                        ungroup() %>% 
                        group_by(Year) %>% 
                        mutate(WY = ifelse(Month %in% c(10,11,12),Year + 1,Year),
                        Date = ymd(paste(Year,'-',Month,'-',Day,sep=""))) %>% 
                        ungroup() %>% 
                        filter(Month %in% c(10,11,12,1,2,3,4) & WY %in% c(2002:2018)) %>% 
                        # Add floodplain inundation proportion
                        mutate(inund_prop = 1*exp(-exp(-0.4514547*(MeanStage-46.7045641)))) 

FloodplainInund_index <- FloodplainStage_0218 %>%
                         group_by(WY) %>% 
                         dplyr::summarise(Index=mean(inund_prop,na.rm=TRUE)) %>% 
                         ungroup()

### Load and clean Sacramento flow data
SacFlow <- read.csv('Data/BTCFlows.csv') 

SacFlow_0218 <- SacFlow %>% 
  group_by(Year,Month,Day) %>% 
  dplyr::summarise(MeanFlow=mean(Flow,na.rm=TRUE)) %>% 
  ungroup() %>% 
  group_by(Year) %>% 
  mutate(WY = ifelse(Month %in% c(10,11,12),Year + 1,Year),
         Date = ymd(paste(Year,'-',Month,'-',Day,sep=""))) %>% 
  ungroup() %>% 
  filter(Month %in% c(10,11,12,1,2,3,4) & WY %in% c(2002:2018))

SacFlow_index <- SacFlow_0218 %>%
                 group_by(WY) %>% 
                 dplyr::summarise(Index=mean(MeanFlow,na.rm=TRUE)) %>% 
                 ungroup()

### Load and clean Delta flow data
DeltaFlow <- read.csv('Data/DeltaFlows.csv')  

DeltaFlow_0218 <- DeltaFlow %>% 
                  mutate(Date=as.Date(Date,format="%m/%d/%Y"),
                  Day=day(Date)) %>% 
                  group_by(Year) %>% 
                  mutate(WY = ifelse(Month %in% c(10,11,12),Year + 1,Year)) %>% 
                  ungroup() %>% 
                  filter(Month %in% c(10,11,12,1,2,3,4) & WY %in% c(2002:2018)) %>% 
                  select(Date,Year,Month,Day,WY,Flow) 

DeltaFlow_index <- DeltaFlow_0218 %>%
                   group_by(WY) %>% 
                   dplyr::summarise(Index=mean(Flow,na.rm=TRUE)) %>% 
                   ungroup()

### Combine hydrologic and escapement data
# Assign DWR water year type to each WY
WY_type <- c(rep('D',3),rep('AN',3),rep('BN',3),rep('AN',3),rep('W',3),rep('D',3),rep('C',3),rep('D',3),
             rep('BN',3),rep('W',3),rep('BN',3),rep('D',3),rep('C',3),rep('C',3),rep('BN',3),rep('W',3),
             rep('BN',3)) 

# Combine floodplain inundation with escapement 
EscapInund <- merge(FloodplainInund_index,Escap_withlag,by.x="WY", by.y="Year",all.x=TRUE) %>% 
              mutate(WY_type=WY_type,
              EscapYear=WY+2) %>% 
              na.omit() %>% 
              select(WY,WY_type,Index,EscapYear,PopID,LagTotal)

# Combine Sacramento River Flow with escapement 
EscapSacFlow <- merge(SacFlow_index ,Escap_withlag,by.x="WY",by.y="Year",all.x=TRUE) %>%
                mutate(WY_type=WY_type,              
                       EscapYear=WY+2) %>% 
                na.omit() %>% 
                select(WY,WY_type,Index,EscapYear,PopID,LagTotal)

# Combine Delta Flow with escapement 
EscapDeltaFlow <- merge(DeltaFlow_index ,Escap_withlag,by.x="WY",by.y="Year",all.x=TRUE) %>%
                  mutate(WY_type=WY_type,
                         EscapYear=WY+2) %>% 
                  na.omit() %>% 
                  select(WY,WY_type,Index,EscapYear,PopID,LagTotal)


# Figures S1, S2, and 3 ----------------------------------------------------------------------------------
# Figure S1
(FloodplainStageFig <- ggplot(FloodplainStage_0218,aes(x=Date,y=MeanStage))+
                       geom_line()+
                       theme_bw()+
                       theme(panel.grid.major = element_blank(),
                       panel.background = element_blank(),
                       axis.line = element_line(colour = "black"),
                       axis.title.y = element_text(margin = unit(c(0, 4, 0, 0), "mm")),
                       text = element_text(size=18))+
                       facet_wrap(WY ~ .,scales="free_x")+
                       labs(x="",y='Stage at Meridian (ft)')
)

(InundStageFig <- ggplot(FloodplainStage_0218,aes(x=MeanStage,y=inund_prop))+
                  geom_line(size=1)+ 
                  theme(panel.grid.major = element_blank(),
                  panel.background = element_blank(),
                  axis.line = element_line(colour = "black"),
                  text = element_text(size=18),
                  axis.title.y = element_text(margin = unit(c(0, 4, 0, 0), "mm")))+
                  labs(x="Stage at Meridian (ft)",y='Inundation proportion')
)

png("Figures/FigS1.png", 
    family = "serif", width = 8, height= 11, units = "in", res =300)

ggarrange(FloodplainStageFig,InundStageFig,
          nrow=2,labels = c("(a)","(b)"), 
          heights=c(4,2),
          font.label=list(size =18,color="black"))

dev.off()

# Figure S2
(EscapButteSacFlow_Fig <- ggplot(EscapSacFlow %>% filter(PopID=="Butte Creek"),
                                 aes(x=Index,y=LagTotal,label=WY,col=WY_type))+
    ggtitle('Butte Creek')+
    xlab('')+ ylab('Escapement')+
    geom_smooth(color = 'black',method="lm") +
    geom_point(size=3)+ 
    geom_text(vjust = -0.5,  size=5,show.legend = F)+
    scale_x_continuous(limits=c(3000,35000), breaks=seq(5000,34000,5000),expand=c(0,0))+
    scale_y_continuous(limits=c(0,24000), breaks=seq(0,24000,5000),expand=c(0,0))+
    stat_cor(p.accuracy = 0.001, r.accuracy = 0.01, 
             aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
             label.x = 5500,label.y=22500,size=6,col="black")+
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          text = element_text(size=20))+
    scale_color_manual(values = c("red3","orange","seagreen","turquoise3", "royalblue1"),
                       breaks=c('C', 'D', 'BN','AN','W'))+
    guides(color = guide_legend(title = "Water year"))
)

(EscapMillSacFlow_Fig  <- ggplot(EscapSacFlow %>% filter(PopID=="Mill Creek"),
                                 aes(x=Index,y=LagTotal,label=WY,col=WY_type))+
    xlab('Sacramento River flow')+ ylab('')+ ggtitle('Mill Creek')+
    geom_smooth(color = 'black',method="lm") +
    geom_point(size=3)+
    geom_text(vjust = -0.5,  size=5,show.legend = F)+
    scale_x_continuous(limits=c(3000,35000), breaks=seq(5000,34000,5000),expand=c(0,0))+
    scale_y_continuous(limits=c(0,1399), breaks=seq(0,1399,200),expand=c(0,0))+
    stat_cor(p.accuracy = 0.001, r.accuracy = 0.01, 
             aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
             label.x = 5500,label.y=1320,size=6,col="black")+
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          text = element_text(size=20))+
    scale_color_manual(values = c("red3","orange","seagreen","turquoise3", "royalblue1"),
                       breaks=c('C', 'D', 'BN','AN','W'))+
    guides(color = guide_legend(title = "Water year"))
)

(EscapDeerSacFlow_Fig <- ggplot(EscapSacFlow %>% filter(PopID=="Deer Creek"),
                                aes(x=Index,y=LagTotal,label=WY,col=WY_type))+
    xlab('')+ylab('')+
    ggtitle('Deer Creek')+
    geom_smooth(color = 'black',method="lm") +
    geom_point(size=3)+ 
    geom_text(vjust = -0.5,  size=5,show.legend = F)+
    scale_x_continuous(limits=c(3000,35000), breaks=seq(5000,34000,5000),expand=c(0,0))+
    scale_y_continuous(limits=c(-100,2850), breaks=seq(-100,2850,400),expand=c(0,0))+
    stat_cor(p.accuracy = 0.001, r.accuracy = 0.01, 
             aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
             label.x = 5500,label.y=2700,size=6,col="black")+
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          text = element_text(size=20))+
    scale_color_manual(values = c("red3","orange","seagreen","turquoise3", "royalblue1"),
                       breaks=c('C', 'D', 'BN','AN','W'))+
    guides(color = guide_legend(title = "Water year"))
)

(EscapButteDeltaFlow_Fig <- ggplot(EscapDeltaFlow %>% filter(PopID=="Butte Creek"),
                               aes(x=Index,y=LagTotal,label=WY,col=WY_type))+
    ggtitle('Butte Creek')+
    xlab('')+ ylab('Escapement')+
    geom_smooth(color = 'black',method="lm") +
    geom_point(size=3)+ 
    geom_text(vjust = -0.5, size=5,show.legend = F)+
    scale_x_continuous(limits=c(1000,105000), breaks=seq(5000,103000,20000),expand=c(0,0))+
    stat_cor(p.accuracy = 0.001, r.accuracy = 0.01, 
             aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
             label.x = 7000,label.y=25000,size=6,col="black")+
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          text = element_text(size=20))+
    scale_color_manual(values = c("red3","orange","seagreen","turquoise3", "royalblue1"),
                       breaks=c('C', 'D', 'BN','AN','W'))+
    guides(color = guide_legend(title = "Water year"))
)

(EscapMillDeltaFlow_Fig <- ggplot(EscapDeltaFlow %>% filter(PopID=="Mill Creek"),
                                  aes(x=Index,y=LagTotal,label=WY,col=WY_type))+
    ggtitle('Mill Creek')+
    xlab('Delta Flow')+ ylab('')+
    geom_smooth(color = 'black',method="lm") +
    geom_point(size=3)+
    geom_text(vjust = -0.5, size=5,show.legend = F)+
    scale_x_continuous(limits=c(1000,105000), breaks=seq(5000,103000,20000),expand=c(0,0))+
    stat_cor(p.accuracy = 0.001, r.accuracy = 0.01, 
             aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
             label.x = 7000,label.y=1300,size=6,col="black")+
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          text = element_text(size=20))+
    scale_color_manual(values = c("red3","orange","seagreen","turquoise3", "royalblue1"),
                       breaks=c('C', 'D', 'BN','AN','W'))+
    guides(color = guide_legend(title = "Water year"))
)

(EscapDeerDeltaFlow_Fig <- ggplot(EscapDeltaFlow %>% filter(PopID=="Deer Creek"),
                                  aes(x=Index,y=LagTotal,label=WY,col=WY_type))+
    ggtitle('Deer Creek')+
    xlab('')+ylab('')+
    geom_smooth(color = 'black',method="lm") +
    geom_point(size=3)+ 
    geom_text(vjust = -0.5, size=5,show.legend = F)+
    scale_x_continuous(limits=c(1000,105000), breaks=seq(5000,103000,20000),expand=c(0,0))+
    scale_y_continuous(limits=c(-500,2800), breaks=seq(0,2800,400))+
    stat_cor(p.accuracy = 0.001, r.accuracy = 0.01, 
             aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
             label.x = 7000,label.y=2800,size=6,col="black")+
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          text = element_text(size=20))+
    scale_color_manual(values = c("red3","orange","seagreen","turquoise3", "royalblue1"),
                       breaks=c('C', 'D', 'BN','AN','W'))+
    guides(color = guide_legend(title = "Water year"))
)


png("Figures/FigS2.png", 
    family = "serif", width = 16, height= 9 , units = "in", res =300)

ggarrange(ggarrange(EscapButteSacFlow_Fig,EscapMillSacFlow_Fig,EscapDeerSacFlow_Fig,
                    nrow=1,common.legend = TRUE, legend="top",
                    labels = "(a)",font.label=list(size =20,color="black")),
          ggarrange(EscapButteDeltaFlow_Fig,EscapMillDeltaFlow_Fig,EscapDeerDeltaFlow_Fig,
                    nrow=1,legend="none",
                    labels = "(b)",font.label=list(size=20,color="black")),
          nrow=2)

dev.off()

# Figure 3
(EscapButteInund_Fig <- ggplot(EscapInund %>% filter(PopID=="Butte Creek"),
                               aes(x=Index,y=LagTotal,label=WY,col=WY_type))+
    ggtitle('Butte Creek')+
    xlab('')+ylab('Escapement')+
    geom_smooth(color = 'black',method="lm") +
    geom_point(size=3)+ 
    geom_text(vjust = -0.5,size=5,show.legend = F)+
    scale_x_continuous(limits=c(0,0.52), breaks=seq(0.05,0.52,0.1),expand=c(0,0))+
    scale_y_continuous(expand=c(0,0))+
    stat_cor(p.accuracy = 0.001, r.accuracy = 0.01, 
             aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
             label.x = 0.01,label.y=21000,size=6,col="black")+
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          text = element_text(size=18),
          plot.title =element_text(size=18, face='bold'),
          legend.position="none")+
    scale_color_manual(values = c("red3","orange","seagreen","turquoise3", "royalblue1"),
                       breaks=c('C', 'D', 'BN','AN','W'))+
    guides(color = guide_legend(title = "Water year"))
)

(EscapMillInund_Fig <- ggplot(EscapInund %>% filter(PopID=="Mill Creek"),
                              aes(x=Index,y=LagTotal,label=WY,col=WY_type))+
    ggtitle('Mill Creek')+
    xlab('Butte Creek floodplain inundated proportion')+ ylab('')+
    geom_smooth(color = 'black',method="lm") +
    stat_cor(p.accuracy = 0.001, r.accuracy = 0.01, 
             aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
             label.x = 0.01,label.y=1220,size=6,col="black")+
    geom_point(size=3)+ 
    geom_text(vjust = -0.5,size=5,show.legend = F)+
    scale_x_continuous(limits=c(0,0.52), breaks=seq(0.05,0.52,0.1),expand=c(0,0))+
    scale_y_continuous(limits=c(0,1300), breaks=seq(0,1200,200),expand=c(0,0))+
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          text = element_text(size=18),
          plot.title =element_text(size=18, face='bold'),
          legend.position="none")+
    scale_color_manual(values = c("red3","orange","seagreen","turquoise3", "royalblue1"),
                       breaks=c('C', 'D', 'BN','AN','W'))+
    guides(color = guide_legend(title = "Water year"))
)

(EscapDeerInund_Fig <- ggplot(EscapInund %>% filter(PopID=="Deer Creek"),
                              aes(x=Index,y=LagTotal,label=WY,col=WY_type))+
    ggtitle('Deer Creek')+
    xlab('')+ylab('')+
    geom_smooth(color = 'black',method="lm") +
    stat_cor(p.accuracy = 0.001, r.accuracy = 0.01, 
             aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
             label.x = 0.01,label.y=2650,size=6,col="black")+
    geom_point(size=3)+ 
    scale_x_continuous(limits=c(0,0.52), breaks=seq(0.05,0.52,0.1),expand=c(0,0))+
    scale_y_continuous(limits=c(-100,2790), breaks=seq(0,2790,400),expand=c(0,0))+
    geom_text(vjust = -0.5,size=5,show.legend = F)+
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          text = element_text(size=18),
          plot.title =element_text(size=18, face='bold'),
          legend.position="right")+
    scale_color_manual(values = c("red3","orange","seagreen","turquoise3", "royalblue1"),
                       breaks=c('C', 'D', 'BN','AN','W'))+
    guides(color = guide_legend(title = "Water year"))
  
)


png("Figures/Fig3.png", 
    family = "serif", width = 18, height= 5, units = "in", res =300)

ggarrange(EscapButteInund_Fig,EscapMillInund_Fig,EscapDeerInund_Fig,
          nrow=1,common.legend = TRUE, legend="top")

dev.off()


# Butte, Mill and Deer Creek spawning and rearing grounds temperatures comparison ------------------------------------------------------------------
# Load temperature data
TempUpperButte_0414 <- read.csv("Data/TempUpperButte_0414.csv")
TempLowerButte_0208 <- read.csv("Data/TempLowerButte_0208.csv")
TempUpperDeer_0414 <- read.csv("Data/TempUpperDeer_0414.csv")
TempLowerDeer_0414 <- read.csv("Data/TemplowerDeer_0414.csv")
TempUpperMill_0414 <- read.csv("Data/TempUpperMill_0414.csv")
TempLowerMill_0414 <- read.csv("Data/TemplowerMill_0414.csv")

# Calculate monthly temperature metrics for each tributary and reach
TempUpperButte_monthly <- TempUpperButte_0414 %>% 
                          group_by(Tributary,Area,Month) %>%
                          summarize(MeanTemp = mean(Temp,na.rm=TRUE),
                                    SdTemp = sd(Temp,na.rm=TRUE),
                                    MinTemp = min(Temp,na.rm=TRUE),
                                    MaxTemp = max(Temp,na.rm=TRUE),
                                    Q5thTemp = quantile(Temp, probs = 0.05,na.rm=TRUE),
                                    Q95thTemp = quantile(Temp, probs = 0.95,na.rm=TRUE)) 

TempLowerButte_monthly <- TempLowerButte_0208 %>% 
                          group_by(Tributary,Area,Month) %>%
                          summarize(MeanTemp = mean(Temp,na.rm=TRUE),
                                    SdTemp = sd(Temp,na.rm=TRUE),
                                    MinTemp = min(Temp,na.rm=TRUE),
                                    MaxTemp = max(Temp,na.rm=TRUE),
                                    Q5thTemp = quantile(Temp, probs = 0.05,na.rm=TRUE),
                                    Q95thTemp = quantile(Temp, probs = 0.95,na.rm=TRUE)) 

TempLowerDeer_monthly <- TempLowerDeer_0414 %>% 
                         group_by(Tributary,Area,Month) %>%
                         summarize(MeanTemp = mean(Temp,na.rm=TRUE),
                                   SdTemp = sd(Temp,na.rm=TRUE),
                                   MinTemp = min(Temp,na.rm=TRUE),
                                   MaxTemp = max(Temp,na.rm=TRUE),
                                   Q5thTemp = quantile(Temp, probs = 0.05,na.rm=TRUE),
                                   Q95thTemp = quantile(Temp, probs = 0.95,na.rm=TRUE)) 

TempUpperDeer_monthly <- TempUpperDeer_0414 %>% 
                         group_by(Tributary,Area,Month) %>%
                         summarize(MeanTemp = mean(Temp,na.rm=TRUE),
                                   SdTemp = sd(Temp,na.rm=TRUE),
                                   MinTemp = min(Temp,na.rm=TRUE),
                                   MaxTemp = max(Temp,na.rm=TRUE),
                                   Q5thTemp = quantile(Temp, probs = 0.05,na.rm=TRUE),
                                   Q95thTemp = quantile(Temp, probs = 0.95,na.rm=TRUE)) 

TempLowerMill_monthly <- TempLowerMill_0414 %>% 
                         group_by(Tributary,Area,Month) %>%
                         summarize(MeanTemp = mean(Temp,na.rm=TRUE),
                                   SdTemp = sd(Temp,na.rm=TRUE),
                                   MinTemp = min(Temp,na.rm=TRUE),
                                   MaxTemp = max(Temp,na.rm=TRUE),
                                   Q5thTemp = quantile(Temp, probs = 0.05,na.rm=TRUE),
                                   Q95thTemp = quantile(Temp, probs = 0.95,na.rm=TRUE)) 

TempUpperMill_monthly <- TempUpperMill_0414 %>% 
                         group_by(Tributary,Area,Month) %>%
                         summarize(MeanTemp = mean(Temp,na.rm=TRUE),
                                   SdTemp = sd(Temp,na.rm=TRUE),
                                   MinTemp = min(Temp,na.rm=TRUE),
                                   MaxTemp = max(Temp,na.rm=TRUE),
                                   Q5thTemp = quantile(Temp, probs = 0.05,na.rm=TRUE),
                                   Q95thTemp = quantile(Temp, probs = 0.95,na.rm=TRUE)) 

# Mill/Deer/Butte monthly temperatures combined
TempMDB_monthly <- data.frame(rbind(TempUpperButte_monthly,TempLowerButte_monthly,
                                    TempUpperDeer_monthly,TempLowerDeer_monthly,
                                    TempUpperMill_monthly,TempLowerMill_monthly))

# Figure S5 -----------------------------------------------
TempMDB_monthly$Area <- factor(TempMDB_monthly$Area,levels=c('Upper','Lower'))

TempMDB_monthly$Tributary <- factor(TempMDB_monthly$Tributary,
                                levels=c("Deer Creek","Mill Creek","Butte Creek")) 

(UpperTempFig <- ggplot(data=TempMDB_monthly %>% filter(Area=="Upper"))+
    geom_line(aes(x=Month,y=MeanTemp,colour=Tributary),size=2)+
    geom_ribbon(aes(x=Month,ymax = Q95thTemp, ymin = Q5thTemp,
                    fill=Tributary,alpha=Tributary))+ 
    scale_alpha_manual(values=c(0.6,0.3,0.3)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.text.x = element_text(angle=45, hjust = 1),
          axis.title.y = element_text(margin = unit(c(0, 4, 0, 0), "mm")),
          text = element_text(size=20),legend.title = element_blank())+
    ylab('Temperature (C)')+ xlab("") + 
    scale_y_continuous(limits=c(1,21), breaks=seq(1,21,4))+
    scale_x_continuous(breaks = seq(1,12,by =2),
                       labels=c('January','March','May','July','September','November'))+
    scale_fill_manual(values = c('darkcyan','goldenrod','purple3'))+  
    scale_color_manual(values = c('darkcyan','goldenrod','purple3'))
  
)

(LowerTempFig <- ggplot(data=TempMDB_monthly %>%  filter(Area=="Lower"))+
  geom_line(aes(x=Month,y=MeanTemp,colour=Tributary),size=2)+
  geom_ribbon(aes(x=Month,ymax = Q95thTemp, ymin = Q5thTemp,
                  fill=Tributary,alpha=Tributary))+ 
    scale_alpha_manual(values=c(0.6,0.3,0.3)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle=45, hjust = 1),
        axis.title.y = element_text(margin = unit(c(0, 4, 0, 0), "mm")),
        text = element_text(size=20),legend.title = element_blank())+
  ylab('Temperature (C)')+ xlab("") + 
  scale_y_continuous(limits=c(3,27), breaks=seq(3,27,4))+
  scale_x_continuous(breaks = seq(1,12,by =2),
                     labels=c('January','March','May','July','September','November'))+
    scale_fill_manual(values = c('darkcyan','goldenrod','purple3'))+  
    scale_color_manual(values = c('darkcyan','goldenrod','purple3'))
)

png("figures/FigS5.png", 
    family = "serif", width = 9, height= 9, units = "in", res =300)

ggarrange(UpperTempFig,LowerTempFig,
          nrow=2,
          common.legend = TRUE,legend="right",
          labels=c('(a)','(b)'),font.label=list(size=18))

dev.off()


# Sutter Bypass weir overtopping analysis ------------------------------------------------------------------
# Load weirs data
filenames <- list.files("Data/SutterWeir/",pattern = ".csv")

weir_list <- vector(mode = "list", length = 17)

j=1
for(i in filenames){
  filepath <- file.path("Data/SutterWeir/",paste(i,"",sep=""))
  weir_list[[j]] <- assign(i, read.csv(filepath))
  j=j+1
}

# Estimate number of days with Colusa weir overtopping during Dec-May period
weir.df <- data.frame(matrix(ncol=2,nrow=0))

for (i in 1:17){
  weir_temp <- cbind(weir_list[[i]][[1]],weir_list[[i]][[6]])
  weir.df <- rbind(weir.df,weir_temp)
}
colnames(weir.df) <- c('Date','Value')

weir.df$Year <- year(as.Date(weir.df$Date,format="%m/%d/%Y"))
weir.df$Month <- month(as.Date(weir.df$Date,format="%m/%d/%Y"))

weir.df <- weir.df %>% 
  filter(Month %in% c(12,1:5)) %>% 
  arrange(-desc(Year))

ColWeirdata <- data.frame(ncol=2,nrow=0)
for (i in 1:(length(unique(weir.df$Year))-1)){
datasubset <- weir.df %>% 
  filter(Year == unique(weir.df$Year)[i] & Month == 12 |
         Year == unique(weir.df$Year)[i+1] & Month  %in% c(1:5))
ColWeirdata[i,1] <- unique(weir.df$Year)[i+1]
ColWeirdata[i,2] <- length(which(!is.na(datasubset$Value)))

}
colnames(ColWeirdata) <- c('Year','NDays')

ColWeirdata
