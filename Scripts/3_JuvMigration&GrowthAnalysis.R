# Code created by Flora Cordoleani to:
# 1. compare Sr profiles and movements of Mill/Deer vs Butte fish
# 2. Compare Butte vs Mill/Deer Creek early life growth in freshwater 

library(dplyr)
library(ggplot2)
library(ggpubr)
library(reshape2)
library(tidyr)
library(NbClust)
library(car)
library(FSA)
library(rstatix)

# Load datasets needed in this script ---------------------------------------------------
otoSrBC <- read.csv('Data/OtoSrBC.csv',stringsAsFactors=FALSE)
otoSrMDC <- read.csv('Data/MillDeerOtoliths.csv',stringsAsFactors = FALSE)
NatalFWExitBC <- read.csv('Data/NatalFWExitBC.csv',stringsAsFactors=FALSE)
NatalFWExitMDC <- read.csv('Data/NatalFWExitMDC.csv',stringsAsFactors=FALSE) 
incMDC <-  read.csv('Data/IncMDC.csv',stringsAsFactors=FALSE)
incBC <- read.csv('Data/IncBC.csv',stringsAsFactors=FALSE)

# Run hierarchical clustering analysis based on otolith radius at natal exit --------------------------------------------
# Prepare data 
clust_BC <- NatalFWExitBC  %>%
  dplyr::select(sample, year,tributary,NatalExit_OR)%>%
  distinct()

oto_forcluster <- clust_BC %>%
  dplyr::select(sample, NatalExit_OR)

oto_forcluster_matrix <- oto_forcluster %>%
  dplyr::select(-sample) %>%
  scale()

# Dissimilarity matrix
dist_mat <- dist(oto_forcluster_matrix, method = "euclidean")

#Find optimal number of clusters with Nbclust
res.nbclust <- NbClust(data = oto_forcluster_matrix, diss = NULL, 
                       distance = "euclidean",
                       min.nc = 2, max.nc = 15, method ="ward.D2")

# Using Ward's method with optimal cluster number = 3
oto_k=3
hc_oto <- hclust(dist_mat, method = "ward.D2" )
sub_grp_oto <- cutree(hc_oto, k = oto_k)

# Number of members in each cluster
table(sub_grp_oto)

# Add clusters back to data
oto_clust_results <- oto_forcluster %>%
                    cbind(sub_grp_oto)

otoSrBC_clust <- merge(otoSrBC,oto_clust_results, by ="sample") %>% 
                 mutate(reartype = case_when(sub_grp_oto == '1' ~ 'Late',
                                             sub_grp_oto == '2' ~ 'Early',
                                             sub_grp_oto == '3' ~'Intermediate')) %>% 
                  dplyr::select(-sub_grp_oto) 

NatalFWExitBC <- merge(NatalFWExitBC,unique(otoSrBC_clust[,c('sample','reartype')]),by='sample') 

# Select columns of interest for both Butte Creek and Mill/Deer Creek data
NatalFWExitBC_sub <- NatalFWExitBC %>% 
  dplyr::select(sample,NatalExit_IncNum,FWExit_IncNum,
                NatalExit_OR,FWExit_OR,reartype,year) %>% 
  mutate(tributary="Butte Creek")

NatalFWExitMDC_sub <- NatalFWExitMDC %>% 
  dplyr::select(sample,NatalExit_IncNum, NatalExit_OR,
                FWExit_IncNum,FWExit_OR,reartype,year) %>% 
  mutate(tributary="Mill/Deer Creek")

# Calculation of FL at natal and FW exit based on OR ------------------------------------------------
### Broken stick FL calibration model
calib = read.csv("Data/OR_FL_FINALforR.csv")
calib = subset(calib, select=c("Sample_ID","OR","FL"))
FL <- calib$FL
OR <- calib$OR
forced.intercept <- 30
res.lm <- lm(I(FL-forced.intercept) ~ 0 + OR)
res.bs <- segmented::segmented(res.lm, seg.Z = ~ 0 + OR)
FL.bs <- predict(res.bs) + forced.intercept
calib$Predict <- FL.bs 

### Application on Butte Creek population
ORNatalBC <- data.frame(OR=as.double(NatalFWExitBC_sub$NatalExit_OR))
NatalFWExitBC_sub$NatalExitFL <- segmented::predict.segmented(res.bs,newdata = ORNatalBC) + 
  forced.intercept

ORFWBC <- data.frame(OR=as.double(NatalFWExitBC_sub$FWExit_OR))
NatalFWExitBC_sub$FWExitFL <- segmented::predict.segmented(res.bs,newdata = ORFWBC) + 
  forced.intercept

# Figure S3
(c1_BC <- ggplot() + geom_point(data = calib,aes(OR,FL),col="black")+
    geom_point(data=NatalFWExitBC_sub ,aes(NatalExit_OR,NatalExitFL),col='red',size=2)+
    geom_line(data=calib,aes(OR,Predict),col='blue',size=1)+
    labs(y=(expression(paste('Fork Length (mm)'))), fill="")+
    labs(x=(expression(paste('Otolith radius (',mu,'m)',sep = ''))))+
     scale_x_continuous(limits=c(150,1000),breaks=seq(150,1000,200))+
    scale_y_continuous(limits=c(30,160),breaks=seq(30,160,20))+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          text = element_text(size=20),legend.title = element_blank(),
          legend.position = "none", plot.title = element_text(face = "bold"))
)

(c2_BC  <- ggplot() + geom_point(data = calib,aes(OR,FL),col="black")+
    geom_point(data=NatalFWExitBC_sub,aes(FWExit_OR,FWExitFL),col='red',size=2)+
    geom_line(data=calib,aes(OR,Predict),col='blue',size=1)+
    labs(x=(expression(paste('Otolith radius (',mu,'m)',sep = ''))),y="")+
    scale_x_continuous(limits=c(150,1000),breaks=seq(150,1000,200))+
    scale_y_continuous(limits=c(30,160),breaks=seq(30,160,20))+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          text = element_text(size=20),legend.title = element_blank(),
          legend.position = "none", plot.title = element_text(face = "bold"))
)

png("figures/FigS3.png", 
    family = "serif", width = 9, height= 4, units = "in", res =300)

ggarrange(c1_BC,c2_BC,nrow=1,labels = c("(a)","(b)"), font.label=list(size =18,color="black"))

dev.off()

### Application on Mill/Deer Creek population
ORNatalMDC <- data.frame(OR=as.double(NatalFWExitMDC_sub$NatalExit_OR))
NatalFWExitMDC_sub$NatalExitFL <- segmented::predict.segmented(res.bs,newdata = ORNatalMDC) + 
  forced.intercept

ORFWMDC <- data.frame(OR=as.double(NatalFWExitMDC_sub$FWExit_OR))
NatalFWExitMDC_sub$FWExitFL <- segmented::predict.segmented(res.bs,newdata = ORFWMDC) + 
  forced.intercept

# Size and Age at natal and FW exit summary statistics -----------------------------------------------------------
# Combine Natal and FW exit data for Mill, Deer and Butte Creek
NatalFWExitMDB <- data.frame(rbind(NatalFWExitMDC_sub,NatalFWExitBC_sub))

NatalExit.summary <- NatalFWExitMDB %>% 
  group_by(tributary) %>% 
  dplyr::summarise(MeanOR = round(mean(NatalExit_OR, na.rm=TRUE)),
                   MaxOR = round(max(NatalExit_OR, na.rm=TRUE)),
                   MinOR = round(min(NatalExit_OR, na.rm=TRUE)),
                   SdOR = round(sd(NatalExit_OR, na.rm=TRUE)),
                   MedOR = round(median(NatalExit_OR, na.rm=TRUE)),
                   Q1OR = round(quantile(NatalExit_OR, probs = 0.25,na.rm=TRUE)),
                   Q3OR = round(quantile(NatalExit_OR, probs = 0.75,na.rm=TRUE)),
                   MeanFL = round(mean(NatalExitFL, na.rm=TRUE)),
                   MaxFL = round(max(NatalExitFL, na.rm=TRUE)),
                   MinFL = round(min(NatalExitFL, na.rm=TRUE)),
                   SdFL = round(sd(NatalExitFL, na.rm=TRUE)),
                   MedFL = round(median(NatalExitFL, na.rm=TRUE)),
                   Q1FL = round(quantile(NatalExitFL, probs = 0.25,na.rm=TRUE)),
                   Q3FL = round(quantile(NatalExitFL, probs = 0.75,na.rm=TRUE)),
                   MeanNatalInc = round(mean(NatalExit_IncNum, na.rm=TRUE)),
                   MaxNatalInc = round(max(NatalExit_IncNum, na.rm=TRUE)),
                   MinNatalInc = round(min(NatalExit_IncNum, na.rm=TRUE)),
                   SdNatalInc = round(sd(NatalExit_IncNum, na.rm=TRUE)),
                   MedNatalInc = round(median(NatalExit_IncNum, na.rm=TRUE)),
                   Q1NatalInc = round(quantile(NatalExit_IncNum, probs = 0.25,na.rm=TRUE)),
                   Q3NatalInc = round(quantile(NatalExit_IncNum, probs = 0.75,na.rm=TRUE))
  )

NatalExit.summary_reartype <- NatalFWExitMDB %>% 
  group_by(tributary,reartype) %>% 
  dplyr::summarise(MeanOR = round(mean(NatalExit_OR, na.rm=TRUE)),
                   MaxOR = round(max(NatalExit_OR, na.rm=TRUE)),
                   MinOR = round(min(NatalExit_OR, na.rm=TRUE)),
                   SdOR = round(sd(NatalExit_OR, na.rm=TRUE)),
                   MedOR = round(median(NatalExit_OR, na.rm=TRUE)),
                   Q1OR = round(quantile(NatalExit_OR, probs = 0.25,na.rm=TRUE)),
                   Q3OR = round(quantile(NatalExit_OR, probs = 0.75,na.rm=TRUE)),
                   MeanFL = round(mean(NatalExitFL, na.rm=TRUE)),
                   MaxFL = round(max(NatalExitFL, na.rm=TRUE)),
                   MinFL = round(min(NatalExitFL, na.rm=TRUE)),
                   SdFL = round(sd(NatalExitFL, na.rm=TRUE)),
                   MedFL = round(median(NatalExitFL, na.rm=TRUE)),
                   Q1FL = round(quantile(NatalExitFL, probs = 0.25,na.rm=TRUE)),
                   Q3FL = round(quantile(NatalExitFL, probs = 0.75,na.rm=TRUE)),
                   MeanNatalInc = round(mean(NatalExit_IncNum, na.rm=TRUE)),
                   MaxNatalInc = round(max(NatalExit_IncNum, na.rm=TRUE)),
                   MinNatalInc = round(min(NatalExit_IncNum, na.rm=TRUE)),
                   SdNatalInc = round(sd(NatalExit_IncNum, na.rm=TRUE)),
                   MedNatalInc = round(median(NatalExit_IncNum, na.rm=TRUE)),
                   Q1NatalInc = round(quantile(NatalExit_IncNum, probs = 0.25,na.rm=TRUE)),
                   Q3NatalInc = round(quantile(NatalExit_IncNum, probs = 0.75,na.rm=TRUE))
  )

FWExit.summary <-NatalFWExitMDB %>% 
  group_by(tributary) %>%
  dplyr::summarise(MeanOR = round(mean(FWExit_OR, na.rm=TRUE)),
                   MaxOR= round(max(FWExit_OR, na.rm=TRUE)),
                   MinOR = round(min(FWExit_OR, na.rm=TRUE)),
                   SdOR = round(sd(FWExit_OR, na.rm=TRUE)),
                   MedOR = round(median(FWExit_OR, na.rm=TRUE)),
                   Q1OR = round(quantile(FWExit_OR, probs = 0.25,na.rm=TRUE)),
                   Q3OR = round(quantile(FWExit_OR, probs = 0.75,na.rm=TRUE)),
                   MeanFL = round(mean(FWExitFL, na.rm=TRUE)),
                   MaxFL = round(max(FWExitFL, na.rm=TRUE)),
                   MinFL = round(min(FWExitFL, na.rm=TRUE)),
                   SdFL = round(sd(FWExitFL, na.rm=TRUE)),
                   MedFL = round(median(FWExitFL, na.rm=TRUE)),
                   Q1FL = round(quantile(FWExitFL, probs = 0.25,na.rm=TRUE)),
                   Q3FL = round(quantile(FWExitFL, probs = 0.75,na.rm=TRUE)),
                   MeanFWInc = round(mean(FWExit_IncNum, na.rm=TRUE)),
                   MaxFWInc = round(max(FWExit_IncNum, na.rm=TRUE)),
                   MinFWInc = round(min(FWExit_IncNum, na.rm=TRUE)),
                   SdFWInc = round(sd(FWExit_IncNum, na.rm=TRUE)),
                   MedFWInc = round(median(FWExit_IncNum, na.rm=TRUE)),
                   Q1FWInc = round(quantile(FWExit_IncNum, probs = 0.25,na.rm=TRUE)),
                   Q3FWInc = round(quantile(FWExit_IncNum, probs = 0.75,na.rm=TRUE))
                   )

FWExit.summary_reartype <-NatalFWExitMDB %>% 
  group_by(tributary,reartype) %>%
  dplyr::summarise(MeanOR = round(mean(FWExit_OR, na.rm=TRUE)),
                   MaxOR= round(max(FWExit_OR, na.rm=TRUE)),
                   MinOR = round(min(FWExit_OR, na.rm=TRUE)),
                   SdOR = round(sd(FWExit_OR, na.rm=TRUE)),
                   MedOR = round(median(FWExit_OR, na.rm=TRUE)),
                   Q1OR = round(quantile(FWExit_OR, probs = 0.25,na.rm=TRUE)),
                   Q3OR = round(quantile(FWExit_OR, probs = 0.75,na.rm=TRUE)),
                   MeanFL = round(mean(FWExitFL, na.rm=TRUE)),
                   MaxFL = round(max(FWExitFL, na.rm=TRUE)),
                   MinFL = round(min(FWExitFL, na.rm=TRUE)),
                   SdFL = round(sd(FWExitFL, na.rm=TRUE)),
                   MedFL = round(median(FWExitFL, na.rm=TRUE)),
                   Q1FL = round(quantile(FWExitFL, probs = 0.25,na.rm=TRUE)),
                   Q3FL = round(quantile(FWExitFL, probs = 0.75,na.rm=TRUE)),
                   MeanFWInc = round(mean(FWExit_IncNum, na.rm=TRUE)),
                   MaxFWInc = round(max(FWExit_IncNum, na.rm=TRUE)),
                   MinFWInc = round(min(FWExit_IncNum, na.rm=TRUE)),
                   SdFWInc = round(sd(FWExit_IncNum, na.rm=TRUE)),
                   MedFWInc = round(median(FWExit_IncNum, na.rm=TRUE)),
                   Q1FWInc = round(quantile(FWExit_IncNum, probs = 0.25,na.rm=TRUE)),
                   Q3FWInc = round(quantile(FWExit_IncNum, probs = 0.75,na.rm=TRUE))
  )

### Stats test 
# Levene test for testing homogeneity of data of fish size and age at natal and FW exit
leveneTest(NatalExitFL~ reartype*as.factor(year), data = NatalFWExitBC)
leveneTest(FWExitFL~ reartype*as.factor(year), data = NatalFWExitBC)

leveneTest(NatalExit_IncNum ~ reartype*as.factor(year), data = NatalFWExitBC)
leveneTest(FWExit_IncNum ~ reartype*as.factor(year), data = NatalFWExitBC)

# Anova and Kruskal-Wallis test
(res.kw.NatalFL <- kruskal.test(NatalExitFL ~ reartype , data =NatalFWExitBC))
(res.kw.FWFL <- kruskal.test(FWExitFL ~ reartype , data =NatalFWExitBC))

NatalIncdata <-   NatalFWExitBC %>% filter(NatalExit_IncNum!='NA')
res.aov.Natalinc <- aov(NatalExit_IncNum~ reartype , data = NatalIncdata)
summary(res.aov.Natalinc)

FWIncdata <-   NatalFWExitBC %>% filter(FWExit_IncNum!='NA')
res.aov.FWinc <- aov(FWExit_IncNum~ reartype , data = FWIncdata)
summary(res.aov.FWinc)

# Tukey and non-parametric Dunn tests
(dunn.res.NatalFL <- dunnTest(NatalExitFL ~ reartype,data =NatalFWExitBC,
                              method="bonferroni"))

(dunn.res.FWFL <- dunnTest(FWExitFL ~ reartype,data =NatalFWExitBC,
                           method="bonferroni"))

(tukey.res.Natalinc <- TukeyHSD(res.aov.Natalinc))
plot(tukey.res.Natalinc)

(tukey.res.FWinc <- TukeyHSD(res.aov.FWinc))
plot(tukey.res.FWinc)

# Figure 4 ------------------------------------------------------------------------------
##### Figure 4a
isoscape <- read.csv("Data/Isoscape.csv")
x = c(0,1200)
position_trib <- data.frame(x=x,y= c(isoscape$Srvalue[1],isoscape$Srvalue[2])) #range of Butte Creek watershed Sr values
position_sac <- data.frame(x=x,y=c(isoscape$Srvalue[2],isoscape$Srvalue[3])) #range of lower Sacramento River Sr values
position_delta <- data.frame(x=x,y=c(isoscape$Srvalue[3],isoscape$Srvalue[4])) #range of Delta Sr values
position_bay<- data.frame(x=x,y=c(isoscape$Srvalue[4],isoscape$Srvalue[5])) #range of Bay & Ocean Sr values
position_incub <- data.frame(x=c(0,200),y=c(isoscape$Srvalue[1],isoscape$Srvalue[5])) # show incubation period

(SrProfilesBCFig <- ggplot() +
    geom_rect(data= position_trib, inherit.aes = FALSE,
              aes(xmin=-Inf, xmax=+Inf, ymin=y[1], ymax=y[2]),
              fill="#001B87", alpha=0.6)+
    geom_rect(data= position_sac, inherit.aes = FALSE,
              aes(xmin=-Inf, xmax=+Inf, ymin=y[1], ymax=y[2]),
              fill= "#0064D6",alpha=0.5)+
    geom_rect(data= position_delta, inherit.aes = FALSE,
              aes(xmin=-Inf, xmax=+Inf, ymin=y[1], ymax=y[2]),
              fill="#00A6D7", alpha=0.2)+
    geom_rect(data= position_bay, inherit.aes = FALSE,
              aes(xmin=-Inf, xmax=+Inf, ymin=y[1], ymax=y[2]),
              fill="#C1E6E0", alpha=1)+
    geom_line(data=otoSrBC_clust,aes(distance,otoSr,group=sample,color=reartype,
                                     alpha=reartype)) + 
    geom_vline(xintercept = 200,linetype="dashed",size=1)+
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          text = element_text(size=18),
          axis.title.y = element_text(margin = unit(c(0, 3, 0, 0), "mm")),
          legend.position="none")+
    ylab(expression(paste({}^"87","Sr/",{}^"86","Sr")))+ 
    xlab(expression(paste('Otolith radius (',mu,'m)',sep = '')))+
    scale_x_continuous(limits=c(0,1200), breaks=seq(0,1000,200),expand=c(0,0))+
    scale_y_continuous(limits=c(0.7042,0.7095), breaks=seq(0.7042,0.7095,0.001),
                       expand=c(0,0))+
    scale_alpha_manual(values=c(0.8,0.3,0.6)) +
    scale_color_manual(values = c("mediumvioletred","grey48","goldenrod2"))+
    
    annotate(geom="text",x=1080, y=0.7050,label="Butte Cr. & \n Sutter Bypass",col="grey", 
             size =6,fontface="bold") +
    annotate(geom="text",x=1050, y=0.7059,label="Sacramento River", col="black",
             size = 6,fontface="bold") +
    annotate(geom="text",x=1100, y=0.7070,label="Delta",  col="black",
             size = 6,fontface="bold")+
    annotate(geom="text",x=1080, y=0.7082,label="Bay & Ocean",  col="black", 
             size = 6,fontface="bold")+
    annotate(geom="text",x=90, y=0.7092,label="Incubation", 
             col="black", size = 6,fontface="bold")
  
)

##### Figure 4b-e 
(FLNatalDens <- ggplot()+
    geom_density(data=NatalFWExitMDB,
                 aes(x=NatalExitFL,
                     fill = reartype,color=reartype),
                 alpha=0.3,adjust=1, size=1) +
    facet_grid(~as.factor(tributary))+ 
    labs(x="",y="Density")+
    ggtitle('Natal Exit ')+
    scale_x_continuous(limits=c(25,135), breaks=seq(25,135,50))+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          text = element_text(size=18),
          axis.title.y = element_text(margin = unit(c(0, 3, 0, 0), "mm")),
          legend.title = element_blank(),
          legend.position='none',
          plot.title = element_text(face = "bold"))+
    scale_color_manual(values = c("mediumvioletred","grey48","goldenrod2"))+
    scale_fill_manual(values = c("mediumvioletred","grey48","goldenrod2"))
)

(FLFWDens <- ggplot(NatalFWExitMDB)+
    geom_density(aes(x=as.double(FWExitFL),
                     fill = reartype,color= reartype),
                 alpha=0.3,adjust=1, size=1) +
    facet_grid(~as.factor(tributary))+ 
    labs(x="Fork Length (mm)",y="Density")+
    ggtitle('Freshwater Exit ')+
    scale_x_continuous(limits=c(40,160), breaks=seq(40,160,40))+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          text = element_text(size=18),
          axis.title.y = element_text(margin = unit(c(0, 3, 0, 0), "mm")),
          legend.title = element_blank(),
          legend.position='none',
          plot.title = element_text(face = "bold"))+
    scale_color_manual(values = c("mediumvioletred","grey48","goldenrod2"))+
    scale_fill_manual(values = c("mediumvioletred","grey48","goldenrod2"))
)

(IncNatalDens <- ggplot(NatalFWExitMDB)+
    geom_density(aes(as.double(NatalExit_IncNum), 
                     fill = reartype,color= reartype),
                 alpha=0.3,adjust=1.2, size=1) +
    facet_grid(~as.factor(tributary))+ 
    ggtitle('Natal Exit ')+
    labs(x="",y="")+
    scale_x_continuous(limits=c(0,300), breaks=seq(0,300,100))+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          text = element_text(size=18),
          axis.title.y = element_text(margin = unit(c(0, 3, 0, 0), "mm")),
          legend.title = element_blank(),
          legend.position='none',
          plot.title = element_text(face = "bold")) +
    scale_color_manual(values = c("mediumvioletred","grey48","goldenrod2"))+
    scale_fill_manual(values = c("mediumvioletred","grey48","goldenrod2"))
)

(IncFWDens <- ggplot(NatalFWExitMDB)+
    geom_density(aes(x=as.double(FWExit_IncNum),
                     fill = reartype,color= reartype),
                 alpha=0.3,adjust=1, size=1) +
    facet_grid(~as.factor(tributary))+ 
    ggtitle('Freshwater Exit ')+
    labs(x=expression('Otolith increment number'),y="")+
    scale_x_continuous(limits=c(0,300), breaks=seq(0,300,100))+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          text = element_text(size=18),
          axis.title.y = element_text(margin = unit(c(0, 3, 0, 0), "mm")),
          legend.title = element_blank(),
          legend.position="none",
          plot.title = element_text(face = "bold"))+
    scale_color_manual(values = c("mediumvioletred","grey48","goldenrod2"))+
    scale_fill_manual(values = c("mediumvioletred","grey48","goldenrod2"))
)

#### Figure 4 all panels combined 
(p1 <- ggarrange(SrProfilesBCFig,labels=c("(a)"),
                font.label=list(size =18,color="black"))
)

(p2 <- ggarrange(ggarrange(FLNatalDens, IncNatalDens,
                           FLFWDens, IncFWDens,
                           labels = c("(b)", "(c)","(d)","(e)"),
                           nrow=2, ncol=2,
                           font.label=list(size =18,color="black"),
                           common.legend = TRUE, legend = "bottom"))
)

(ProfDensPlot <- ggarrange(p1, p2, ncol = 1,heights=c(2,3),align = "hv")
  
)

png("figures/Fig4.png",family = "serif", width = 9, height=11, units = "in", res =300)

ProfDensPlot

dev.off()

# Rearing location comparison across watershed and years ----------------------------------------------
#### Mill/Deer Creek
RearLocMDC <- data.frame()

# Generate unique fish IDs
FishID <-unique(NatalFWExitMDC$sample)

for (i in 1:length(unique(incMDC$sample))){ 
  dataSubset_inc <-subset(incMDC,sample== unique(incMDC$sample)[i])
  dataSubset_OR <- subset(NatalFWExitMDC,sample==unique(incMDC$sample)[i])
  
  NatalExit_OR <- NatalFWExitMDC$NatalExit_OR[which(NatalFWExitMDC$sample==FishID[i])]
  FWExit_OR <- NatalFWExitMDC$FWExit_OR[which(NatalFWExitMDC$sample==FishID[i])] 
  
  endpoint <- length(dataSubset_inc$inc_distance)
  
  # Add rearing habitat type 
  if(dim(dataSubset_OR)[1]==0){ # No OR distance measured 
    dataSubset_inc$habitat <- rep(NA,endpoint)
  } else {
    tribhabitat <- which.min(abs(dataSubset_inc$inc_distance - dataSubset_OR$NatalExit_OR))
    sacdeltahabitat <- which.min(abs(dataSubset_inc$inc_distance - dataSubset_OR$FWExit_OR)) - 
      which.min(abs(dataSubset_inc$inc_distance - dataSubset_OR$NatalExit_OR))
    bayoceanhabitat <- endpoint - which.min(abs(dataSubset_inc$inc_distance - dataSubset_OR$FWExit_OR))
    dataSubset_inc$habitat <- c(rep("Mill/Deer Cr.",tribhabitat),rep("Sac. River/Delta",sacdeltahabitat),
                                rep("Bay/Ocean",bayoceanhabitat))
  }
  RearLocMDC <- rbind.data.frame(RearLocMDC, dataSubset_inc)
}

rearlocMDC.summary <-  RearLocMDC %>% 
  group_by(year,sample) %>% 
  count(habitat) %>%
  mutate(Freq = n/sum(n)) %>%
  drop_na()

# plot rearing location proportions per year
rearlocMDC.summary_cast <- dcast(rearlocMDC.summary, year + sample ~ habitat,
                                      value.var=c("n"),fun.aggregate = sum)
rearlocMDC.summary_cast$MillDeervalues <- rearlocMDC.summary_cast$`Mill/Deer Cr.`
rearlocMDC.summary_cast$Totvalue <- (rearlocMDC.summary_cast$`Mill/Deer Cr.`+
                                     rearlocMDC.summary_cast$`Sac. River/Delta` +
                                     rearlocMDC.summary_cast$`Bay/Ocean`)

rearlocMDC.summary_melt <- melt(rearlocMDC.summary_cast, 
                                     id.vars = c("year","sample","MillDeervalues","Totvalue"),
                                     measure.vars = c("Mill/Deer Cr.",
                                                      "Sac. River/Delta",
                                                      "Bay/Ocean"))

rearlocMDC.summary_melt$variable <-factor(rearlocMDC.summary_melt$variable,
                                               levels=c("Mill/Deer Cr.",
                                                        "Sac. River/Delta",
                                                        "Bay/Ocean")) 
rearlocMDC_stats <- rearlocMDC.summary_melt %>% group_by(variable) %>% 
  dplyr::summarise(Mean = round(mean(value, na.rm=TRUE)),
                   Sd=round(sd(value, na.rm=TRUE)),
                   Min=round(min(value, na.rm=TRUE)),
                   Max=round(max(value, na.rm=TRUE)),
                   Med = round(median(value, na.rm=TRUE)),
                   Q1 = round(quantile(value, probs = 0.25,na.rm=TRUE)),
                   Q3 = round(quantile(value, probs = 0.75,na.rm=TRUE)))


(RearlocMDCFig <- rearlocMDC.summary_melt %>% group_by(year)%>%
    filter(!variable=="Bay/Ocean") %>% 
    arrange(-desc(MillDeervalues)) %>% 
    ungroup() %>%
    mutate(sample= factor(sample, levels=unique(sample))) %>% 
    ggplot(aes(x = sample,y = value,fill = variable)) +
    geom_col(position = position_stack(reverse = TRUE),width=3)+
    labs(y ="No. Days") +
    facet_grid(~as.factor(year),scale="free_x")+ 
    theme(panel.grid.major = element_blank(), 
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.x = element_blank(),
          text = element_text(size=18),
          axis.title.y = element_text(margin = unit(c(0, 4, 0, 0), "mm")),
          legend.position = "bottom",
          legend.title = element_blank())+
    scale_fill_manual(values = c("#001B87","#0064D6","#C1E6E0"))+
    scale_y_continuous(expand=c(0,0),limits=c(0,300), breaks=seq(0,300,50))+
    scale_x_discrete(expand=c(0,0))
)

### Butte Creek
RearLocBC <- data.frame()

# Generate unique fish IDs
FishID <-unique(NatalFWExitBC$sample)

for (i in 1:length(unique(incBC$sample))){ 
  dataSubset_inc <-subset(incBC,sample== unique(incBC$sample)[i])
  dataSubset_OR <- subset(NatalFWExitBC,sample==unique(incBC$sample)[i])
  
  SutterExit_OR <- NatalFWExitBC$NatalExit_OR[which(NatalFWExitBC$sample==FishID[i])]
  FWExit_OR <- NatalFWExitBC$FWExit_OR[which(NatalFWExitBC$sample==FishID[i])] 
  
  endpoint <- length(dataSubset_inc$inc_distance)
  
  # Add rearing habitat type 
  if(length(dataSubset_OR)==0){ # No OR distance measured by George
    dataSubset_inc$habitat <- rep(NA,endpoint)
  } else {
    tribhabitat <- which.min(abs(dataSubset_inc$inc_distance - dataSubset_OR$NatalExit_OR))
    sacdeltahabitat <- which.min(abs(dataSubset_inc$inc_distance - dataSubset_OR$FWExit_OR)) - 
      which.min(abs(dataSubset_inc$inc_distance - dataSubset_OR$NatalExit_OR))
    bayoceanhabitat <- endpoint - which.min(abs(dataSubset_inc$inc_distance - dataSubset_OR$FWExit_OR))
    dataSubset_inc$habitat <- c(rep("Butte Cr. & Sutter Bypass",tribhabitat),rep("Sac. River/Delta",sacdeltahabitat),
                                rep("Bay/Ocean",bayoceanhabitat))
  }
  RearLocBC <- rbind.data.frame(RearLocBC, dataSubset_inc)
}

rearlocBC.summary <-  RearLocBC  %>% 
  group_by(year,sample) %>% 
  count(habitat) %>%
  mutate(Freq = n/sum(n))

### plot rearing location proportions per year
rearlocBC.summary_cast <- dcast(rearlocBC.summary, year + sample ~ habitat,
                                    value.var=c("n"),fun.aggregate = sum)

rearlocBC.summary_cast$Buttevalues <- rearlocBC.summary_cast$`Butte Cr. & Sutter Bypass`

rearlocBC.summary_cast$Totvalue <- (rearlocBC.summary_cast$`Butte Cr. & Sutter Bypass` +
                                    rearlocBC.summary_cast$`Sac. River/Delta`+
                                    rearlocBC.summary_cast$`Bay/Ocean`)

rearlocBC.summary_melt <- melt(rearlocBC.summary_cast, 
                                   id.vars = c("year","sample","Buttevalues","Totvalue"),
                                   measure.vars = c("Butte Cr. & Sutter Bypass",
                                                    "Sac. River/Delta",
                                                    "Bay/Ocean"))

rearlocBC.summary_melt$variable <-factor(rearlocBC.summary_melt$variable,
                                             levels=c("Butte Cr. & Sutter Bypass",
                                                      "Sac. River/Delta",
                                                      "Bay/Ocean")) 

rearlocBC_stats <- rearlocBC.summary_melt %>% group_by(variable) %>% 
  dplyr::summarise(Mean = round(mean(value, na.rm=TRUE)),
                   Sd=round(sd(value, na.rm=TRUE)),
                   Min=round(min(value, na.rm=TRUE)),
                   Max=round(max(value, na.rm=TRUE)),
                   Med = round(median(value, na.rm=TRUE)),
                   Q1 = round(quantile(value, probs = 0.25,na.rm=TRUE)),
                   Q3 = round(quantile(value, probs = 0.75,na.rm=TRUE)))

(RearlocBCFig <- rearlocBC.summary_melt %>% group_by(year)%>% 
    filter(!variable=="Bay/Ocean") %>% 
    arrange(-desc(Buttevalues)) %>% 
    ungroup() %>%
    mutate(sample= factor(sample, levels=unique(sample))) %>% 
    ggplot(aes(x = sample,y = value,fill = variable)) +
    geom_col(position = position_stack(reverse = TRUE),width=3)+
    labs(y ="") +
    facet_grid(~as.factor(year),scale="free_x")+ 
    theme(panel.grid.major = element_blank(), 
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.x = element_blank(),
          text = element_text(size=18),
          axis.title.y = element_text(margin = unit(c(0, 4, 0, 0), "mm")),
          legend.position = "bottom",
          legend.title = element_blank())+
    scale_fill_manual(values=c("#001B87","#0064D6","#C1E6E0"))+
    scale_y_continuous(expand=c(0,0),limits=c(0,300), breaks=seq(0,300,50))+
    scale_x_discrete(expand=c(0,0))
)

# Increment growth comparison across watershed ------------------------------------------------------------------------
# Combine Mill, Deer and Butte data
growthBC <- merge(NatalFWExitBC,incBC,by=c("sample","year"))
growthMDC <- merge(NatalFWExitMDC,incMDC,by=c("sample","year")) %>% 
             mutate(watershed="Mill/Deer Creek")

growthMDB <- data.frame(rbind(growthBC[,c('sample', 'watershed','year','reartype',
                                          'inc_num','inc_distance','inc_width')],
                               growthMDC[,c('sample', 'watershed','year','reartype',
                                            'inc_num','inc_distance','inc_width')]))

# Make sure the increment numbers are ordered
growthMDB <- growthMDB %>% 
             group_by(sample) %>% 
             arrange(inc_num) %>% 
             ungroup()
  
# Estimate FL at each increment (i.e for each day)
ORdata <- data.frame(OR=as.double(growthMDB$inc_distance)) 
growthMDB$FL <- segmented::predict.segmented(res.bs,newdata = ORdata) + 
  forced.intercept

growthMDB <- growthMDB %>% 
             group_by(sample) %>% 
             mutate(FLdiff = FL -lag(FL)) %>% 
             ungroup()

# Estimate FL growth rate statistics for each watershed and life history type
growthMDB.summary <- growthMDB %>% 
                     group_by(watershed) %>% 
                     summarise(meanFLdiff= mean(FLdiff, na.rm=TRUE),
                               minFLdiff= min(FLdiff, na.rm=TRUE),
                               maxFLdiff= max(FLdiff, na.rm=TRUE),
                               sdFLdiff= sd(FLdiff, na.rm=TRUE),
                               medFLdiff= median(FLdiff, na.rm=TRUE))

growthMDB.summary_reartype <- growthMDB %>% 
                              group_by(watershed,reartype) %>% 
                              summarise(meanFLdiff= mean(FLdiff, na.rm=TRUE),
                                        minFLdiff= min(FLdiff, na.rm=TRUE),
                                        maxFLdiff= max(FLdiff, na.rm=TRUE),
                                        sdFLdiff= sd(FLdiff, na.rm=TRUE),
                                        medFLdiff= median(FLdiff, na.rm=TRUE))

# Levene test for testing homogeneity of data
# for first 81 days (median age at freshwater exit for Butte Creek fish)
growthMDB_81 <- growthMDB %>% 
                filter(inc_num < 82) 

growthMDB_81$watershed <- as.factor(growthMDB_81$watershed)

growthMDB_81.summary <- growthMDB_81 %>% 
                        group_by(watershed) %>% 
                        summarise(meanFLdiff= mean(FLdiff, na.rm=TRUE),
                        minFLdiff= min(FLdiff, na.rm=TRUE),
                        maxFLdiff= max(FLdiff, na.rm=TRUE),
                        sdFLdiff= sd(FLdiff, na.rm=TRUE),
                        medianFLdiff= median(FLdiff, na.rm=TRUE))

leveneTest(FLdiff ~ watershed *as.factor(year), data = growthMDB_81)

# Kruskal-Wallis test
(res.kw.growth <- kruskal.test(FLdiff ~ watershed , data = growthMDB_81))

# Figures 5b,c
(GrowthMDBFig <- ggplot(growthMDB)+ 
    stat_smooth(aes(x=as.double(inc_num),y = as.double(FLdiff),group=sample,
                    color=watershed),geom='line', alpha=0.08, se=TRUE)+
    geom_smooth(aes(x=as.double(inc_num),y = as.double(FLdiff),
                    color=watershed))+ 
    scale_color_manual(values=c("purple3","springgreen3"))+
    geom_vline(data=FWExit.summary,aes(xintercept=MedFWInc,colour=as.factor(tributary)),lty=2, size=1)+
    theme(panel.grid.major = element_blank(), 
          axis.text.x = element_text(colour = "black",angle = 45, hjust = 1),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          text = element_text(size=18),
          axis.title.y = element_text(margin = unit(c(0, 4, 0, 0), "mm")),
          legend.position = "bottom",
          legend.title = element_blank())+
    labs(x="No. Days ",y='Growth rate (mm/day)')+
    scale_x_continuous(limits=c(0,350), breaks=seq(0,350,50))
)

# Add stat test to boxplot
res.wx.growth <- growthMDB_81 %>% 
  pairwise_wilcox_test(FLdiff ~ watershed)

(GrowthMDBBoxplot <- ggplot(growthMDB_81,aes(x = watershed, y =  as.double(FLdiff))) +
    geom_boxplot(outlier.shape = NA,alpha=0.5,aes(fill=watershed)) +
    scale_fill_manual(values=c("purple3","springgreen3")) +
    stat_pvalue_manual(res.wx.growth,y.position =1.35)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          text = element_text(size=18),legend.title = element_blank(),
          axis.title.y = element_text(margin = unit(c(0,4, 0, 0), "mm")),
          axis.text.x=element_blank(),axis.ticks.x=element_blank(),
          legend.position = "bottom") +
    scale_y_continuous(limits=c(0,1.35), breaks=seq(0,1.35,0.2))+
    ylab('Growth rate in first 81 days')+ xlab("")
)

# Figure 5 ------------------------------------------------------------------------
(RearlocGrowth_MDB <- ggarrange(RearlocMDCFig, RearlocBCFig, GrowthMDBFig,GrowthMDBBoxplot,
                                  nrow=2,ncol=2,
                                  widths = c(2,3), heights=c(2,3),
                                  labels = c("(a)", "(b)", "(c)","(d)"),
                                  font.label=list(size =18,color="black"))
)


png("Figures/Fig5.png", 
    family = "sans serif", width =12, height=8, units = "in", res =200)

RearlocGrowth_MDB 

dev.off()

