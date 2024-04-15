#Load packages 
library(tidyverse)
library(ggforce)
library(R2jags)
library(tidybayes)
library(DHARMa)
library(gridExtra)
library(ggridges)
library(viridisLite)
library(cowplot)

#Source the functions for models, summaries, etc 
source("Code/SH_Functions_Final.r")

#Run models? If not just read in .RDS files or only read in the summary files 

#Models <- "Run" 
#Models <- "Read"
Models <- "Neither"

#Summaries <- "Run"
Summaries <- "Read"

#Read in Test fishery catch data
Catch <- Load_Catch_Data()

#Read in Annual return to Albion including data up to 2022
annual_return<-read.csv(file="Data/annual_SH_return_to_Albion_1983-2022 DFO version Aug 1 to Nov 30_updated.csv",header=T)

#These are useful to have on hand
minDay <- min(Catch$Day)
minYear <- min(Catch$Year)
maxYear <- max(Catch$Year)

#Data inputs - tried to make more flexible
DataList <- list ("n_years"=length(unique(Catch$Year)), 
                  "n_fisheries" = length(unique(Catch$Fishery)),
                  "n_obs" = dim(Catch)[[1]],
                  "years" = (Catch$Year-minYear+1),
                  "days" = (Catch$Day-minDay+1),
                  "fishery" = ifelse(Catch$Fishery == "CN_test", 1, 2), 
                  "obs_catch"=Catch$SH_Catch,
                  "annual_return"=annual_return$N,
                  # also include prior specification
                  "rt_m_mu" = 280-minDay, "rt_m_sig" = 40, # mean date of Oct 7
                  "halft_s" = 15  , "halft_df" = 2, 
                  "logit_q_mean" = -3, "logit_q_tau" = 0.8)


# Final Set of JAGS models
#-------------------------------------------------------------------------------------------------------
#   1) Similar to Bison (2021) - Non Hier, Poisson, q hier (this might be a little different than his)
#      Non_Hier_Poisson_Q

if(Models == "Run"){
  JagsFit_Indep_Poisson <- jags.parallel(DataList, model.file =Non_Hier_Poisson_Q, 
                                     n.chains =3, n.iter = 20000, n.burnin = 5000, n.thin = 5,
                                     parameters.to.save = c("rt_m", "rt_sd", "q", "pred_abund"))
  # Check convergence
  All_Ests <-   data.frame(JagsFit_Indep_Poisson$BUGSoutput$summary)
  All_Ests$Param <- row.names(All_Ests) 
  All_Ests %>% filter(Rhat>1.01)
  # Save JAGS output
  saveRDS(JagsFit_Indep_Poisson, "Outputs/1_ModelRuns/JagsFit_Indep_Poisson.RDS")
  # Look if catchabilities look similar to Bison's
  
  ggplot(Summary_Indep_Poisson$Pred_Catch_Ints) +
    geom_pointrange(aes(x=Year, y=q, ymin = q.lower, ymax=q.upper), size = 0.25) +
    ylab("Catchability")+
    theme_light()+
    facet_wrap(~Fishery)+
    ggtitle("Catchability Across Years, Like-Bison Model")
  ggsave("Outputs/3_Plots/Bison_Catchability.png", width = 11, height = 8, units = "in")
  
  Summary_Indep_Poisson$Pred_Catch_Ints %>% group_by(Fishery) %>% summarise(mean(q))
  # <chr>            <dbl>
  #  CM_test         0.0288
  #  CN_test         0.0127
  # Bison estimates are 0.032, 0.011, so fairly similar
  
  # what about mean dates?
  Summary_Indep_Poisson$Date_Dists_Ints %>%  summarise(mean(rt_m)+minDay-1)
  # average is 282, vs Bison's 284
  Summary_Indep_Poisson$Date_Dists_Ints %>% summarise(mean(rt_sd))
  # mean is 21.5, similar to Bison's "corrected" value of 21
  
} else if(Models == "Read") {
  # read in JAGs output and summary from RDS
  JagsFit_Indep_Poisson <- readRDS("Outputs/1_ModelRuns/JagsFit_Indep_Poisson.RDS")
}
if(Summaries == "Run"){
  # summarize Data
  Summary_Indep_Poisson <- Summarize_Data(JagsFit_Indep_Poisson, Catch, Hier=F, ANorm=F,
                                          Timing_Params = "rt_m[Year], rt_sd[Year]")
  # Save Data Summary
  saveRDS(Summary_Indep_Poisson, "Outputs/1_ModelRuns/Summary_Indep_Poisson.RDS")
  # Create plots and print
  Create_Plots(Summary_Data=Summary_Indep_Poisson)
} else if(Summaries == "Read") {
  Summary_Indep_Poisson <- readRDS("Outputs/1_ModelRuns/Summary_Indep_Poisson.RDS")
}



#----------------------------------------------------------------------------------------------

#2  Neg Binom normal, q hier, mean, sd hierarchical
if(Models == "Run"){
  JagsFit_Hier_NB_HQ <- jags.parallel(DataList, model.file = Hier_NB_HQ, 
                             n.chains =3, n.iter=20000, n.burnin = 5000, n.thin = 5, 
                             parameters.to.save = c("rt_m", "rt_sd", "r", "q", "pred_abund", "rt_m_m", "rt_sd_m"))
  
  All_Ests_Hier_NB_HQ <-   data.frame(JagsFit_Hier_NB_HQ$BUGSoutput$summary)
  All_Ests_Hier_NB_HQ$Param <- row.names(All_Ests_Hier_NB_HQ) 
  All_Ests_Hier_NB_HQ %>% filter(Rhat>1.01)
  #none flagged!
  
  # Save JAGS output
  saveRDS(JagsFit_Hier_NB_HQ, "Outputs/1_ModelRuns/Hier_NB_HQ.RDS")
  
} else if(Models == "Read"){ 
  JagsFit_Hier_NB_HQ <- readRDS("Outputs/1_ModelRuns/Hier_NB_HQ.RDS")
}
if(Summaries == "Run"){
  # Summarize JAGS outputs
  Summary_Hier_NB_HQ <- Summarize_Data(JagsFit_Hier_NB_HQ, Catch, Timing_Params = "rt_m[Year], rt_sd[Year], rt_m_m, rt_sd_m", Hier=T, ANorm=F)
  # Save Summary
  saveRDS(Summary_Hier_NB_HQ, "Outputs/1_ModelRuns/Summary_Hier_NB_HQ.RDS")
  # Create Plots
  Create_Plots(Summary_Data=Summary_Hier_NB_HQ)
} else if(Summaries == "Read"){ 
  Summary_Hier_NB_HQ <- readRDS("Outputs/1_ModelRuns/Summary_Hier_NB_HQ.RDS")
}


#-----------------------------------------------------------------------
# 3 Neg Binom Asym. norm, q hier, mean, sd hierarchical 
if(Models == "Run"){
  JagsFit_Hier_NB_ANorm <- jags.parallel(DataList, model.file = Hier_NB_ANorm, 
                                          n.chains =3, n.iter=20000, n.burnin = 5000, n.thin = 5, 
                                          parameters.to.save = c("rt_m", "rt_sd", "rt_m_m", "rt_sd_m", "r", "q", "pred_abund"))
  All_Ests_Hier_NB_ANorm <-   data.frame(JagsFit_Hier_NB_ANorm$BUGSoutput$summary)
  All_Ests_Hier_NB_ANorm$Param <- row.names(All_Ests_Hier_NB_ANorm) 
  All_Ests_Hier_NB_ANorm %>% filter(Rhat>1.01) %>% filter(!str_detect(Param, "^pred_abund"))
  
  # Save JAGS output
  saveRDS(JagsFit_Hier_NB_ANorm, "Outputs/1_ModelRuns/Hier_NB_ANorm.RDS")
 
} else if(Models == "Read"){ 
  JagsFit_Hier_NB_ANorm <- readRDS("Outputs/1_ModelRuns/Hier_NB_ANorm.RDS")
}
if(Summaries == "Run"){
  # Summarize JAGS outputs
  Summary_Hier_NB_ANorm <- Summarize_Data(JagsFit_Hier_NB_ANorm, Catch, Hier=T, ANorm=T,
                                          Timing_Params = "rt_m[Year], rt_sd[Side,Year], rt_m_m, rt_sd_m[Side]")
  # Save Summary
  saveRDS(Summary_Hier_NB_ANorm, "Outputs/1_ModelRuns/Summary_Hier_NB_ANorm.RDS")
  # Create Plots
  Create_Plots(Summary_Data=Summary_Hier_NB_ANorm)
} else if(Summaries == "Read") { 
  Summary_Hier_NB_ANorm <- readRDS("Outputs/1_ModelRuns/Summary_Hier_NB_ANorm.RDS")
}


#--------------------------------------------------------------------------------------------
#Model comparison
#--------------------------------------------------------------------------------------------
# Can we compare these Models with DIC?
if(Models == "Run" | Models == "Read"){
  JagsMods <- list(
    "Hierarchical Normal"= JagsFit_Hier_NB_HQ$BUGSoutput$DIC,
    "Hierarchical Asymmetric Normal" = JagsFit_Hier_NB_ANorm$BUGSoutput$DIC,
    "Independent Normal" = JagsFit_Indep_Poisson$BUGSoutput$DIC
  )
  
  DIC_DF <- data.frame(Mod = names(JagsMods), DIC = unlist(JagsMods))
  
  JagsMods_pD <- list(
    "Hierarchical Normal"= JagsFit_Hier_NB_HQ$BUGSoutput$pD,
    "Hierarchical Asymmetric Normal" = JagsFit_Hier_NB_ANorm$BUGSoutput$pD,
    "Independent Normal"  = JagsFit_Indep_Poisson$BUGSoutput$pD
  )
  
  DIC_DF$pD <- unlist(JagsMods_pD)
  write.csv(DIC_DF, "Outputs/2_DataOut/DIC_DF.csv", row.names=F)
}




#------------------------------------------------------------------------------------------
#  Summary Plots
#----------------------------------------------------------------------------------------
# Plot showing differences in timing across years for each model

#Get the windows for each year from the summary objects 
Yearly1 <- Summary_Indep_Poisson$Yearly_Windows %>% mutate(Model = "Indep_Poisson")
Yearly2 <- Summary_Hier_NB_HQ$Yearly_Windows %>% mutate(Model = "Hier_NB_HQ") 
Yearly3 <- Summary_Hier_NB_ANorm$Yearly_Windows %>% mutate(Model = "Hier_NB_ANorm") 
Yearly_Winds <- rbind(Yearly3, Yearly2, Yearly1) %>% mutate_if(is.numeric, round) %>% mutate(Year = Year+minYear-1)
#write to a csv 
write.csv(Yearly_Winds, "Outputs/2_DataOut/Yearly_Windows.csv")

#Get the global windows from the summary objects 
Global1 <- Summary_Indep_Poisson$Global_Windows %>% mutate(Model = "Indep_Poisson") %>% select(-starts_with("."))
Global2 <- Summary_Hier_NB_HQ$Global_Windows %>% mutate(Model = "Hier_NB_HQ") %>% select(-Year) %>% select(-starts_with("."))
Global3 <- Summary_Hier_NB_ANorm$Global_Windows %>% mutate(Model = "Hier_NB_ANorm") %>% select(-Year) %>% select(-starts_with("."))
Global_Winds <- rbind(Global3, Global2, Global1) %>% mutate_if(is.numeric, round)
#Convert day to date for excel formatting
toDate <- function(x){as.Date(x-1, origin = as.Date("2023-01-01"))}
Global_Winds_Dates <- Global_Winds %>% mutate_if(is.numeric, toDate)
write.csv(Global_Winds, "Outputs/2_DataOut/Global_Windows.csv")
write.csv(Global_Winds_Dates, "Outputs/2_DataOut/Global_Windows_Dates.csv")

# Simplify catch to get rid of fishery, other_catch
Total_Catch <- Catch %>% group_by(Year, Day) %>% summarise(Catch = sum(SH_Catch)) %>%
  filter(Catch > 0)

#Get list of models 
Mods <- unique(Global_Winds$Model)
#Create empty data frame to get catch values outside "global" window
outs <- data.frame()
#Loop through each model to get catch values outside "global" window
for(mm in 1:length(Mods)){
  Dvals <- Global_Winds %>% filter(Model==Mods[mm]) %>%
    select(Model, hD_025, hD_975)
  outliers2 <- Total_Catch %>%
    filter(Day < round(Dvals$hD_025) | Day > round(Dvals$hD_975))
  outliers2$Model <- Mods[mm]
  outs <- rbind(outs, outliers2)
}

# Change Model names to make more readable and for plotting
outs <- ungroup(outs) %>%
  mutate(Model = fct_recode(as.factor(Model),
                            "Independent Normal" = "Indep_Poisson",
                            "Hierarchical Normal" = "Hier_NB_HQ",
                            "Hierarchical Asymmetric Normal" = "Hier_NB_ANorm"))

#Change the column names for yearly windows for merging and plotting
names(Yearly_Winds) <- gsub(x = names(Yearly_Winds), pattern = "h", replacement = "")  

#Join all the windows together 
All_Windows <- Yearly_Winds %>% 
  full_join(Global_Winds, join_by(Model), relationship = "many-to-one") %>%
  mutate(Model = fct_recode(as.factor(Model),
                            "Independent Normal" = "Indep_Poisson",
                            "Hierarchical Normal" = "Hier_NB_HQ",
                            "Hierarchical Asymmetric Normal" = "Hier_NB_ANorm"))
ANorm_025 <- Global_Winds$hD_025[Global_Winds$Model == "Hier_NB_ANorm"]
ANorm_975 <- Global_Winds$hD_975[Global_Winds$Model == "Hier_NB_ANorm"]

#plot windows on top of eachother
p0 <- ggplot(All_Windows) +
  geom_rect(aes(ymin=min(Year), ymax=max(Year), 
                xmin = hD_025, xmax=hD_975), fill='grey', alpha=0.1) +
  geom_point(aes(y=Year, x=D_50), size = 0.5) +
  geom_segment(aes(y=Year, x=D_025, xend=D_975, yend=Year), size = 0.25) +
  #geom_line(aes(y=Year, x=hD_50), size = 0.25, col="red") +
  geom_vline(xintercept=ANorm_025, linetype="dashed", color = "red")+
  geom_vline(xintercept=ANorm_975, linetype="dashed", color = "red")+
  geom_point(data = outs, aes(y=Year, x=Day, size=Catch), col="blue", alpha=0.5) +
  facet_wrap(~factor(Model, levels = c('Independent Normal', 'Hierarchical Normal', 'Hierarchical Asymmetric Normal')), ncol=1)  +
  xlab("Date")+
  theme_bw()+
  theme(legend.position = c(.9,.88), 
        legend.background = element_blank())+
  scale_x_continuous(breaks=c(214, 244, 274, 305, 335),
                     labels=c("Aug 1", "Sept 1", "Oct 1", "Nov 1", "Dec 1")) +
  ggtitle("95% Migration Window")

p0
#ggsave("Outputs/3_Plots/Window_comparison_All.png", width = 8, height = 11, units = "in")


#Create a single panel plot for the asym model only for the SAR
p0_asym <- ggplot(All_Windows[All_Windows$Model == "Hierarchical Asymmetric Normal",]) +
  geom_rect(aes(ymin=min(Year), ymax=max(Year), 
                xmin = hD_025, xmax=hD_975), fill='grey', alpha=0.1) +
  geom_point(aes(y=Year, x=D_50), size = 0.5) +
  geom_segment(aes(y=Year, x=D_025, xend=D_975, yend=Year), size = 0.25) +
  geom_point(data = outs[outs$Model == "Hierarchical Asymmetric Normal",], aes(y=Year, x=Day, size=Catch), col="blue", alpha=0.5) +
  xlab("Date")+
  theme_bw()+
  scale_x_reverse(breaks=c(214, 244, 274, 305, 335),
                  labels=c("Aug 1", "Sept 1", "Oct 1", "Nov 1", "Dec 1")) +
  coord_flip()
p0_asym
ggsave("Outputs/3_Plots/Run_Timing_asym.png", width = 5, height = 3, units = "in")


#names in french
#Normale asym. hiér., binomiale nég.
#Normale hiér., binomiale nég.
#Normale indép., poisson

All_Windows$Model_french <- recode_factor(All_Windows$Model, "Hierarchical Asymmetric Normal" = "Normale asymmétrique hiérarchique",
                                       "Hierarchical Normal" = "Normale hiérarchique",
                                       "Independent Normal" = "Normale indépendente")

outs$Model_french <- recode_factor(outs$Model, "Hierarchical Asymmetric Normal" = "Normale asymmétrique hiérarchique",
                                  "Hierarchical Normal" = "Normale hiérarchique",
                                  "Independent Normal" = "Normale indépendente")

p0_french <- ggplot(All_Windows) +
  geom_rect(aes(ymin=min(Year), ymax=max(Year), 
                xmin = hD_025, xmax=hD_975), fill='grey', alpha=0.1) +
  geom_point(aes(y=Year, x=D_50), size = 0.5) +
  geom_segment(aes(y=Year, x=D_025, xend=D_975, yend=Year), size = 0.25) +
  #geom_line(aes(y=Year, x=hD_50), size = 0.25, col="red") +
  geom_vline(xintercept=ANorm_025, linetype="dashed", color = "red")+
  geom_vline(xintercept=ANorm_975, linetype="dashed", color = "red")+
  geom_point(data = outs, aes(y=Year, x=Day, size=Catch), col="blue", alpha=0.5) +
  facet_wrap(~factor(Model_french, levels = c("Normale indépendente", "Normale hiérarchique", "Normale asymmétrique hiérarchique")), ncol=1)  +
  xlab("Date")+
  ylab("Année")+
  labs(size  = "Prise")+
  theme_bw()+
  theme(legend.position = c(.9,.88), 
        legend.background = element_blank())+
  scale_x_continuous(breaks=c(214, 244, 274, 305, 335),
                     labels=c("1 août", "1 sept", "1 oct", "1 nov", "1 déc")) +
  ggtitle("Période de 95 % de la montaison")

#periode de 95% de la migration?
p0_french

#Plot just the asymmetric in french
p0_asym_french <- ggplot(All_Windows[All_Windows$Model == "Hierarchical Asymmetric Normal",]) +
  geom_rect(aes(ymin=min(Year), ymax=max(Year), 
                xmin = hD_025, xmax=hD_975), fill='grey', alpha=0.1) +
  geom_point(aes(y=Year, x=D_50), size = 0.5) +
  geom_segment(aes(y=Year, x=D_025, xend=D_975, yend=Year), size = 0.25) +
  geom_point(data = outs[outs$Model == "Hierarchical Asymmetric Normal",], aes(y=Year, x=Day, size=Catch), col="blue", alpha=0.5) +
  xlab("Date")+
  ylab("Année")+
  labs(size  = "Prise")+
  theme_bw()+
  scale_x_reverse(breaks=c(214, 244, 274, 305, 335),
                  labels=c("1 août", "1 sept", "1 oct", "1 nov", "1 déc")) +
  coord_flip()
p0_asym_french
ggsave("Outputs/3_Plots/Run_Timing_asym_french.png", width = 5, height = 3, units = "in")


#Calculate proportion of catch outside window
Yearly_Catch <- Catch %>% group_by(Year) %>% summarise(Yearly_Catch = sum(SH_Catch))
Outs_Summary_Yearly <- outs %>% group_by(Model, Year) %>% 
                         summarise(Catch_Out = sum(Catch)) %>%
                         left_join(Yearly_Catch) %>%
                         mutate(Prop_Out = Catch_Out/Yearly_Catch)
Outs_Summary_Yearly %>% filter(Prop_Out>0.05)
Outs_Summary_Yearly %>% group_by(Model) %>% summarise(median(Prop_Out))
write.csv(Outs_Summary_Yearly, "Outputs/2_DataOut/Catch_Outside_Windows_Yearly.csv")
Outs_Summary <- outs %>% group_by(Model) %>% 
                  summarise(Catch_Out = sum(Catch))
Outs_Summary$Prop_Catch <- Outs_Summary$Catch_Out/sum(Catch$SH_Catch)
write.csv(Outs_Summary, "Outputs/2_DataOut/Catch_Outside_Windows.csv")

# I think the two hier Models are much more sensible, independent one gives crazy windows
# for some years, going into July

# create same plot with only two hier Models
# doesn't seem to like separate, join all together
Mods_In_Order <- c("Hierarchical Asymmetric Normal", "Hierarchical Normal")

# Plots curves on toNorm_Hparams$rt_m_m + minDay-1p
# grab global curve params from Anorm_Params and Norm_Params - or draws from posterior?
Norm_Params <- Summary_Hier_NB_HQ$Date_Dists_Ints %>% filter(Year == 1983)
ANorm_Params <- Summary_Hier_NB_ANorm$Date_Dists_Ints %>% filter(Year == 1983)
     
NormDist <- data.frame(
    x = seq(220, 350, by = 0.1),
    y = dnorm(seq(220, 350, by = 0.1), mean = Norm_Params$rt_m_m+ minDay-1,  sd = Norm_Params$rt_sd_m),
    Model = "Hierarchical Normal",
    Model_french = "Normale hiérarchique"
)
ANormDist <- data.frame(
  x = seq(220, 350, by = 0.1),
  y = ddnorm(seq(220, 350, by = 0.1), mean = unique(ANorm_Params$rt_m_m) + minDay-1, sigma = as.numeric(ANorm_Params$rt_sd)),
  Model = "Hierarchical Asymmetric Normal",
  Model_french = "Normale asymmétrique hiérarchique"
)

#Change model names for plotting
Global_Winds_Plot <- Global_Winds %>%
  mutate(Model = fct_recode(as.factor(Model),
                            "Independent Normal" = "Indep_Poisson",
                            "Hierarchical Normal" = "Hier_NB_HQ",
                            "Hierarchical Asymmetric Normal" = "Hier_NB_ANorm"))

# Subset for diff wins
Dists <- rbind(NormDist, ANormDist) %>%
  left_join(Global_Winds_Plot, join_by(Model), relationship = "many-to-many") %>%
  mutate(Win95 = y*ifelse(x > hD_025 & x < hD_975, 1, 0)) %>%
  mutate(Win90 = y*ifelse(x > hD_05 & x < hD_95, 1, 0)) %>%
  mutate(Win80 = y*ifelse(x > hD_10 & x < hD_90, 1, 0))

#Plot the global windows for the two hierarchical models
p2 <- ggplot(data = Dists) +
  geom_polygon(aes(x = x, y = Win95), fill = "black", alpha=0.2) +
  geom_polygon(aes(x = x, y = Win90), fill = "black", alpha=0.3) +
  geom_polygon(aes(x = x, y = Win80), fill = "black", alpha=0.4) +
  geom_line(aes(x,y)) +
  facet_wrap(~factor(Model, levels = c('Hierarchical Normal', 'Hierarchical Asymmetric Normal')), ncol=1)+
  scale_x_continuous(breaks=c(244, 274, 305, 335),expand = c(0,0))+
  scale_y_continuous(breaks=c(0),
                     labels=c(""))+
  #ggtitle("Average Run-Timing with 95%, 90%, and 80% Windows") +
  theme(plot.title = element_text(hjust = 1))+
  xlab("Date") +
  theme_bw()+
  theme(plot.margin=unit(c(2.5,1.9,1.7,5), "mm"),
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank())  
p2


#French translation
p2_french <- ggplot(data = Dists) +
  geom_polygon(aes(x = x, y = Win95), fill = "black", alpha=0.2) +
  geom_polygon(aes(x = x, y = Win90), fill = "black", alpha=0.3) +
  geom_polygon(aes(x = x, y = Win80), fill = "black", alpha=0.4) +
  geom_line(aes(x,y)) +
  facet_wrap(~factor(Model_french, levels = c("Normale hiérarchique", "Normale asymmétrique hiérarchique")), ncol=1)+
  scale_x_continuous(breaks=c(244, 274, 305, 335),expand = c(0,0))+
  scale_y_continuous(breaks=c(0),
                     labels=c(""))+
  #ggtitle("Average Run-Timing with 95%, 90%, and 80% Windows") +
  theme(plot.title = element_text(hjust = 1))+
  xlab("Date") +
  theme_bw()+
  theme(plot.margin=unit(c(2.5,1.9,1.7,5), "mm"),
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank())  
p2_french

#Montaison moyenne, périodes de 95 %, 90 %, et 80 %


#Need to extend the dates for the plots 
ANormDist_extended <- data.frame(
  x = seq(213, 350, by = 0.1),
  y = ddnorm(seq(213, 350, by = 0.1), mean = unique(ANorm_Params$rt_m_m) + minDay-1, sigma = as.numeric(ANorm_Params$rt_sd)),
  Model = "Hierarchical Asymmetric Normal"
)

ANormDist_extended <- ANormDist_extended %>%
  left_join(Global_Winds_Plot, join_by(Model), relationship = "many-to-many") %>%
  mutate(Win95 = y*ifelse(x > hD_025 & x < hD_975, 1, 0)) %>%
  mutate(Win90 = y*ifelse(x > hD_05 & x < hD_95, 1, 0)) %>%
  mutate(Win80 = y*ifelse(x > hD_10 & x < hD_90, 1, 0))

#Create a version of the above plot that only has the top panel for the Asym
p2_asym <- ggplot(data = ANormDist_extended) +
  geom_polygon(aes(x = x, y = Win95), fill = "black", alpha=0.2) +
  geom_polygon(aes(x = x, y = Win90), fill = "black", alpha=0.3) +
  geom_polygon(aes(x = x, y = Win80), fill = "black", alpha=0.4) +
  geom_line(aes(x,y)) +
  scale_x_continuous(breaks=c(213,244, 274, 305, 335),
                     labels=c("1 août", "1 sept", "1 oct", "1 nov", "1 déc"), limits = c(213,335))+
  #scale_y_continuous(breaks=c(0),labels=c(""))+
  theme(plot.title = element_text(hjust = 1))+
  xlab("Date")+
  ylab("Density")+
  theme_bw()+
  theme(plot.margin=unit(c(2.1,1.9,4.5,1.5), "mm"))

p2_asym

#Create a version of the above plot that only has the top panel for the Asym
p2_asym_french <- ggplot(data = ANormDist_extended) +
  geom_polygon(aes(x = x, y = Win95), fill = "black", alpha=0.2) +
  geom_polygon(aes(x = x, y = Win90), fill = "black", alpha=0.3) +
  geom_polygon(aes(x = x, y = Win80), fill = "black", alpha=0.4) +
  geom_line(aes(x,y)) +
  scale_x_continuous(breaks=c(213,244, 274, 305, 335),
                     labels=c("Aug 1", "Sept 1", "Oct 1", "Nov 1", "Dec 1"), limits = c(213,335))+
  scale_y_continuous(labels = scales::label_comma(big.mark = ".",
                                                  decimal.mark = ","))+
  theme(plot.title = element_text(hjust = 1))+
  xlab("Date")+
  ylab("Densité")+
  theme_bw()+
  theme(plot.margin=unit(c(2.1,1.9,4.5,1.5), "mm"))

p2_asym_french

#Get run timing curves for the independent poisson model
IndepParams <- Summary_Indep_Poisson$Date_Dists_Ints
IndepNormDist <- list()
for(i in 1:length(unique(IndepParams$Year))){
  year <- unique(IndepParams$Year)[i]
   probs <- data.frame(
    x = seq(220, 350, by = 0.1),
    y = dnorm(seq(220, 350, by = 0.1), mean = IndepParams$rt_m[IndepParams$Year == year] + minDay-1,  sd = IndepParams$rt_sd[IndepParams$Year == year]),
    Year = year)
  IndepNormDist[[i]] <- probs
}
IndepNormDist<-do.call(rbind, IndepNormDist)

#Get global window for independent poisson model
Indep_Mean_Window <- Global_Winds %>% filter(Model == "Indep_Poisson")

#Plot the run timing curves for the independent poisson model 
p3 <- ggplot(IndepNormDist) +
  geom_line(aes(x,y, group = Year), alpha = 0.5)+
  geom_segment(data = Indep_Mean_Window, aes(x = hD_025, y = -0.002 , yend = -0.002, xend = hD_975), alpha = 0.2, linewidth = 2)+
  geom_segment(data = Indep_Mean_Window, aes(x = hD_05, y = -0.002 , yend = -0.002, xend = hD_95), alpha = 0.3, linewidth = 2)+
  geom_segment(data = Indep_Mean_Window, aes(x = hD_10, y = -0.002 , yend = -0.002, xend = hD_90), alpha = 0.4, linewidth = 2)+
  scale_x_continuous(limits=c(220,350), expand = c(0,0))+
  scale_y_continuous(breaks=c(0),labels=c(""), limits = c(-0.0035, 0.052), expand = c(NA, 0))+
  #xlab("Date")+
  ylab("")+
  ggtitle("Average Run-Timing with 95%, 90%, and 80% Windows") +
  theme_bw()+
  theme(plot.margin=unit(c(2,1.9,-1,4.2), "mm"), axis.title = element_blank(),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank(), axis.text.x = element_blank())+
  annotate("rect", xmin = 220, xmax = 350, ymin = 0.047, ymax = 0.052,
           color="black", fill = "#D9D9D9")+ 
  annotate(label = "Independent Normal", x= 283, y=0.0495, geom="text", size = 3)
p3
#ggsave(file= "Outputs/3_Plots/indep_poisson.png")

p3_french <- ggplot(IndepNormDist) +
  geom_line(aes(x,y, group = Year), alpha = 0.5)+
  geom_segment(data = Indep_Mean_Window, aes(x = hD_025, y = -0.002 , yend = -0.002, xend = hD_975), alpha = 0.2, linewidth = 2)+
  geom_segment(data = Indep_Mean_Window, aes(x = hD_05, y = -0.002 , yend = -0.002, xend = hD_95), alpha = 0.3, linewidth = 2)+
  geom_segment(data = Indep_Mean_Window, aes(x = hD_10, y = -0.002 , yend = -0.002, xend = hD_90), alpha = 0.4, linewidth = 2)+
  scale_x_continuous(limits=c(220,350), expand = c(0,0))+
  scale_y_continuous(breaks=c(0),labels=c(""), limits = c(-0.0035, 0.052), expand = c(NA, 0))+
  #xlab("Date")+
  ylab("")+
  ggtitle("Montaison moyenne, périodes de 95 %, 90 %, et 80 %") +
  theme_bw()+
  theme(plot.margin=unit(c(2,1.9,-1,4.2), "mm"), axis.title = element_blank(),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank(), axis.text.x = element_blank())+
  annotate("rect", xmin = 220, xmax = 350, ymin = 0.047, ymax = 0.052,
           color="black", fill = "#D9D9D9")+ 
  annotate(label = "Normale indépendente", x= 283, y=0.0495, geom="text", size = 3)
p3_french


#Create 6-panel comparison plot with p1, p2, p3
p4 <- grid.arrange(p3, p2, ncol = 1, heights = c(1,2))
p4_french <- grid.arrange(p3_french, p2_french, ncol = 1, heights = c(1,2))

png("Outputs/3_Plots/Run_Timing_6panel_labelled.png", width = 10, height = 9, units = "in", res = 300)
plot_grid(p0,p4,ncol = 2, labels = c('A', 'B'))
dev.off()

png("Outputs/3_Plots/Run_Timing_6panel_labelled_french.png", width = 10, height = 9, units = "in", res = 300)
plot_grid(p0_french,p4_french,ncol = 2, labels = c('A', 'B'))
dev.off()

#=====================================================================
#Multipanel plot with all the priors 
# normal 280, 40
x1<- seq(280-120,280+120, by = 1)
y1 <- dnorm(x1, 280, sd=40)
Priors_DF <- data.frame(x=x1, y=y1, Dist = "Normal(280, 40)")
# half-t 15,2
x2 <- seq(0, 50, by=0.1)
y2 <- extraDistr::dht(x2, nu=2, sigma=15)
Priors_DF <- rbind(Priors_DF, data.frame(x=x2, y=y2, Dist="Half-t(15, 2)"))
# half-t 7.5, 2
x3 <- seq(0, 50, by=0.1)
y3 <- extraDistr::dht(x3, nu=2, sigma=7.5)
Priors_DF <- rbind(Priors_DF, data.frame(x=x3, y=y3, Dist="Half-t(7.5, 2)"))
x4 <- seq(0.0001, 0.2, by = 0.0001) 
y4 <- dilogit(x4, mu = -3, tau = 0.8) 
Priors_DF <- rbind(Priors_DF, data.frame(x=x4, y=y4, Dist="logit(Normal(-3, 1/sqrt(0.8)))")) 

png("Outputs/3_Plots/Priors_Multipanel.png", width=10, height=8, units = "in", res = 300)
ggplot(Priors_DF) +
  geom_line(aes(x, y)) +
  facet_wrap(~Dist, scales="free")+
  scale_y_continuous(breaks=c(0),
                     labels=c("")) +
  ylab("Density") +
  xlab("Value") +
  theme_bw()
dev.off()

#French
#Use scales::label_comma for french decimal display
Priors_DF$Dist_french <- recode_factor(Priors_DF$Dist, "Normal(280, 40)" = "Normale(280, 40)",
                                       "Half-t(15, 2)" = "Demi-Student(15, 2)",
                                       "Half-t(7.5, 2)" = "Demi-Student(7.5, 2)",               
                                       "logit(Normal(-3, 1/sqrt(0.8)))" = "logit(Normale(-3, 1/racine(0.8)))")

# "Normal(280, 40)" = Normale(280, 40)
# "Half-t(15, 2)" = "Demi-Student(15, 2)
# "Half-t(7.5, 2)" = "Demi-Student(7.5, 2)               
# "logit(Normal(-3, 1/sqrt(0.8)))" = "logit(Normale(-3, 1/racine(0.8)))"

#Change these to the right order 

png("Outputs/3_Plots/Priors_Multipanel_french.png", width=10, height=8, units = "in", res = 300)
ggplot(Priors_DF) +
  geom_line(aes(x, y)) +
  facet_wrap(~ factor(Dist_french,levels = c("Demi-Student(15, 2)", "Demi-Student(7.5, 2)",
                                               "logit(Normale(-3, 1/racine(0.8)))", "Normale(280, 40)")), scales="free")+
  scale_y_continuous(breaks=c(0),
                     labels=c("")) +
  scale_x_continuous(labels = scales::label_comma(big.mark = ".",
                                                  decimal.mark = ","))+
  ylab("Densité") +
  xlab("Valeur") +
  theme_bw()
dev.off()


#=========================================================
#Make multipanel figure for catch
#Plot all non-zero catch values across all years 
a <- ggplot(Catch[Catch$SH_Catch>0,])+
  geom_point(aes(x=Day, y=SH_Catch, fill = Fishery, color = Fishery,shape = Fishery),size = 2, alpha = 0.5)+
  theme_bw()+
  xlab("Date")+
  ylab("Steelhead Catch")+
  scale_x_continuous(breaks=c(213,244, 274, 305, 335),
                     labels=c("Aug 1", "Sept 1", "Oct 1", "Nov 1", "Dec 1"), limits = c(213,335))+
  theme(legend.position = c(.15, .8))+
  labs(fill  = "Fishery", color = "Fishery", shape = "Fishery")+
  scale_shape_discrete(labels = c("Chum Test", "Chinook Test"))+
  scale_fill_discrete(labels = c("Chum Test", "Chinook Test"))+
  scale_color_discrete(labels = c("Chum Test", "Chinook Test"))
a

a_french <- ggplot(Catch[Catch$SH_Catch>0,])+
  geom_point(aes(x=Day, y=SH_Catch, fill = Fishery, color = Fishery,shape = Fishery),size = 2, alpha = 0.5)+
  theme_bw()+
  xlab("Date")+
  ylab("Prise de truite arc-en-ciel anadrome")+
  scale_x_continuous(breaks=c(213,244, 274, 305, 335),
                     labels=c("1 août", "1 sept", "1 oct", "1 nov", "1 déc"), limits = c(213,335))+
  theme(legend.position = c(.15, .8))+
  labs(fill  = "Pêche d'essai", color = "Pêche d'essai", shape = "Pêche d'essai")+
  scale_shape_discrete(labels = c("Saumon kéta", "Saumon quinnat"))+
  scale_fill_discrete(labels = c("Saumon kéta", "Saumon quinnat"))+
  scale_color_discrete(labels = c("Saumon kéta", "Saumon quinnat"))

a_french 

#Ridgeline plots of catch for every year 
b <- ggplot(Catch, aes(x=Day, y=Year, height = SH_Catch, group = Year))+ 
  geom_ridgeline(aes(fill = Year, color = Year), alpha = 0.5)+ 
  scale_x_continuous(breaks=c(213,244, 274, 305, 335),
                     labels=c("Aug 1", "Sept 1", "Oct 1", "Nov 1", "Dec 1"), limits = c(213,335))+
  scale_y_reverse()+
  xlab("Date")+
  scale_color_viridis_c(direction = -1,name="Year")+
  scale_fill_viridis_c(direction = -1, name = "Year")+
  theme_bw()+
  theme(legend.position = c(.2, .85), legend.direction = "horizontal")
b

b_french <- ggplot(Catch, aes(x=Day, y=Year, height = SH_Catch, group = Year))+ 
  geom_ridgeline(aes(fill = Year, color = Year), alpha = 0.5)+ 
  scale_x_continuous(breaks=c(213,244, 274, 305, 335),
                     labels=c("1 août", "1 sept", "1 oct", "1 nov", "1 déc"), limits = c(213,335))+
  scale_y_reverse()+
  xlab("Date")+
  ylab("Année")+
  scale_color_viridis_c(direction = -1, name = "Année")+
  scale_fill_viridis_c(direction = -1, name = "Année")+
  theme_bw()+
  theme(legend.position = c(.2, .85), legend.direction = "horizontal")
b_french

#Convert year to numeric
annual_return$Year <- as.numeric(annual_return$Year)
#Plot annual return time series 
c <- ggplot(annual_return, aes(x=Year, y = N))+
  geom_point()+
  geom_line()+
  theme_bw()+
  scale_x_continuous(breaks=c(1983, 1993, 2003, 2013,2022))+
  ylab("Annual Return")+
  theme(plot.margin = unit(c(0.25,0.25,0.25,0.65), "cm"))
c

c_french <- ggplot(annual_return, aes(x=Year, y = N))+
  geom_point()+
  geom_line()+
  theme_bw()+
  scale_x_continuous(breaks=c(1983, 1993, 2003, 2013,2022))+
  ylab("Remontes annuelles")+
  xlab("Année")+
  theme(plot.margin = unit(c(0.25,0.25,0.25,0.65), "cm"))
c_french

#Create multipanel plot
plot_grid(a,b,c,ncol = 1, rel_heights = c(2,2,1), labels = "AUTO")
ggsave(filename = "Outputs/3_Plots/Catch_multipanel.png", height = 8, width = 5.5, units = "in")

#Create french figure 
plot_grid(a_french,b_french,c_french,ncol = 1, rel_heights = c(2,2,1), labels = "AUTO")
ggsave(filename = "Outputs/3_Plots/Catch_multipanel_french.png", height = 8, width = 5.5, units = "in")

#Create a figure for the SAR that combines catch and run timing curve panels
plot_grid(c, b, p2_asym, ncol = 1, rel_heights = c(1.3,2.5,2), labels = "AUTO")
ggsave(filename = "Outputs/3_Plots/Catch_Run_Timing_multipanel.png", height = 8, width = 5.5, units = "in")

plot_grid(c_french, b_french, p2_asym_french, ncol = 1, rel_heights = c(1.3,2.5,2), labels = "AUTO")
ggsave(filename = "Outputs/3_Plots/Catch_Run_Timing_multipanel_french.png", height = 8, width = 5.5, units = "in")




#===========================================
#Plot residuals 
pdf(file= "Outputs/3_Plots/Residuals_Hier_NB_ANorm.pdf", width = 8.5, height = 4)
plot(Summary_Hier_NB_ANorm$sims, title = "DHARMa residual diagnostics: Hierarchical Asymmetric Normal")
dev.off()

pdf(file= "Outputs/3_Plots/Residuals_Hier_NB_HQ.pdf", width = 8.5, height = 4)
plot(Summary_Hier_NB_HQ$sims, title = "DHARMa residual diagnostics: Hierarchical Normal")
dev.off()

pdf(file= "Outputs/3_Plots/Residuals_Indep_Poisson.pdf", width = 8.5, height = 4)
plot(Summary_Indep_Poisson$sims, title = "DHARMa residual diagnostics: Independent Normal")
dev.off()

#These plots don't seem to work with the multipanelling packages, so need to paste them together manually


#=========================================
#Write table with date_dist_ints
write.csv(Summary_Hier_NB_ANorm$Date_Dists_Ints, "Outputs/2_DataOut/ANorm_Date_Dists.csv")
write.csv(Summary_Hier_NB_HQ$Date_Dists_Ints, "Outputs/2_DataOut/Norm_Date_Dists.csv")
write.csv(Summary_Indep_Poisson$Date_Dists_Ints, "Outputs/2_DataOut/Indep_Poisson_Dists.csv")

#=========================================
#Additional plots 

#Read in spawner data from Thompson and Chilcotin
spawners <- read.csv("Data/Spawners from RBison 230630.csv")

#Subset to 1983 -
spawners$BroodYear <- as.numeric(spawners$BroodYear)
spawners <- spawners[spawners$BroodYear >= 1983,]
spawners <- spawners[,-c(4,5)]
spawners$Thompson <- as.numeric(spawners$Thompson)
spawners$Chilcotin <- as.numeric(spawners$Chilcotin)
spawners$Total <- spawners$Thompson + spawners$Chilcotin
spawners$prop_chilcotin <- spawners$Chilcotin/spawners$Total

spawners_bar <- data.frame("BroodYear" = rep(spawners$BroodYear, 2),
              "Stock" = c(rep("Thompson", nrow(spawners)), rep("Chilcotin", nrow(spawners))),
              "Count" = c(spawners$Thompson, spawners$Chilcotin))


spawners1 <- ggplot(spawners_bar)+
  geom_bar(aes(x=BroodYear, y = Count, fill = Stock), position = "fill", stat = "identity")+
  theme_bw()+
  ylab("Proportion")+
  scale_fill_manual(values = c("deepskyblue4", "darkorange1"))
spawners1

spawners2 <- ggplot(spawners)+
  geom_point(aes(x=BroodYear, y = Thompson), color = "deepskyblue4")+
  geom_line(aes(x=BroodYear, y = Thompson), color = "deepskyblue4")+
  geom_point(aes(x=BroodYear, y = Chilcotin, color = "darkorange1"))+
  geom_line(aes(x=BroodYear, y = Chilcotin, color = "darkorange1"))+
  theme_bw()+
  ylab("Count")+
  theme(plot.margin = unit(c(0,3.35,0.25,0.25), "cm"), legend.position = "none")
spawners2

plot_grid(spawners1, spawners2, nrow = 2)
ggsave(file = "Outputs/3_Plots/Thompson_Chilcotin_spawners.png", width = 8, height = 7, units = "in")

#Plot catchability values 
q_ests_Hier_NB_ANorm <- Summary_Hier_NB_ANorm$Pred_Catch_Ints %>% 
  group_by(Fishery, Year) %>% 
  summarise(q=mean(q), q.lower=mean(q.lower), q.upper=mean(q.upper)) %>% 
  mutate(Model = "Hierarchical Asymmetric Normal")

#Read in the object with q_mean
q_mean <- readRDS("Outputs/1_ModelRuns/Q_Summary_Hier_NB_Anorm.RDS")
qmedchin <- median(InvLogit(q_mean$q_mean[q_mean$Fishery == "1"]))
qmedchum <- median(InvLogit(q_mean$q_mean[q_mean$Fishery == "2"]))

q <- ggplot()+
  geom_point(data = q_ests_Hier_NB_ANorm, aes(x=q, y = Year, color = Fishery), size = 2.5, alpha = 0.6)+
  geom_segment(data = q_ests_Hier_NB_ANorm, aes(x= q.lower, xend = q.upper, y = Year, yend = Year, color = Fishery), linewidth = 1.2, alpha = 0.6)+
  geom_vline(xintercept = qmedchin, color = "#00BFC4")+
  geom_vline(xintercept = qmedchum, color = "#F8766D")+
  theme_bw()+
  theme(plot.title = element_text(size=10))+
  xlab("Catchability (q)")+
  ggtitle("Hierarchical Asymmetric Normal")+
  scale_color_discrete(labels = c("Chum Test", "Chinook Test"))
q
ggsave(file = "Outputs/3_Plots/Catchability_Hier_Asym.png", width = 5, height = 4)

#Catchability in french
q_french <- ggplot()+
  geom_point(data = q_ests_Hier_NB_ANorm, aes(x=q, y = Year, color = Fishery), size = 2.5, alpha = 0.6)+
  geom_segment(data = q_ests_Hier_NB_ANorm, aes(x= q.lower, xend = q.upper, y = Year, yend = Year, color = Fishery), linewidth = 1.2, alpha = 0.6)+
  geom_vline(xintercept = qmedchin, color = "#00BFC4")+
  geom_vline(xintercept = qmedchum, color = "#F8766D")+
  theme_bw()+
  theme(plot.title = element_text(size=10))+
  xlab("Capturabilité (q)")+
  ylab("Année")+
  labs(color = "Pêche d'essai")+
  ggtitle("Normale asymmétrique hiérarchique")+
  scale_color_discrete(labels = c("Saumon kéta", "Saumon quinnat"))+
  scale_x_continuous(labels = scales::label_comma(big.mark = ".",
                                                  decimal.mark = ","))
q_french
ggsave(file = "Outputs/3_Plots/Catchability_Hier_Asym_french.png", width = 5, height = 4)

#Plotting catchability for all the other models 

# q_ests_Hier_NB_HQ <- Summary_Hier_NB_HQ$Pred_Catch_Ints %>% 
#   group_by(Fishery, Year) %>% 
#   summarise(q=mean(q), q.lower=mean(q.lower), q.upper=mean(q.upper)) %>%
#   mutate(Model = "Hierarchical Normal")
# 
# q_ests_Indep_Poisson <- Summary_Indep_Poisson$Pred_Catch_Ints %>% 
#   group_by(Fishery, Year) %>% 
#   summarise(q=mean(q), q.lower=mean(q.lower), q.upper=mean(q.upper)) %>%
#   mutate(Model = "Independent Normal")

# q_all <- bind_rows(q_ests_Hier_NB_ANorm, q_ests_Hier_NB_HQ, q_ests_Indep_Poisson)
# 
# ggplot(q_all)+
#   geom_point(aes(x=q, y = Year, color = Fishery), size = 1.5, alpha = 0.6)+
#   geom_segment(aes(x= q.lower, xend = q.upper, y = Year, yend = Year, color = Fishery), linewidth = 0.75, alpha = 0.6)+
#   theme_bw()+
#   theme(plot.title = element_text(size=10))+
#   xlab("Catchability (q)")+
#   facet_wrap(~Model, ncol = 1)+
#   scale_color_discrete(labels = c("Chum Test", "Chinook Test"))
# ggsave(file = "Outputs/3_Plots/Catchability.png", width = 5, height = 7.5) 
# 
#   


# 
# q1 <- ggplot()+
#   geom_point(data = q_ests_Indep_Poisson, aes(x=q, y = Year, color = Fishery), size = 2.5, alpha = 0.6)+
#   geom_segment(data = q_ests_Indep_Poisson, aes(x= q.lower, xend = q.upper, y = Year, yend = Year, color = Fishery), linewidth = 1.2, alpha = 0.6)+
#   theme_bw()+
#   theme(legend.position = "none", plot.title = element_text(size=10))+
#   xlab("")+
#   ggtitle("Independent Normal")
# q1  
# 
# 
# q2 <- ggplot()+
#   geom_point(data = q_ests_Hier_NB_HQ, aes(x=q, y = Year, color = Fishery), size = 2.5, alpha = 0.6)+
#   geom_segment(data = q_ests_Hier_NB_HQ, aes(x= q.lower, xend = q.upper, y = Year, yend = Year, color = Fishery), linewidth = 1.2, alpha = 0.6)+
#   theme_bw()+
#   theme(legend.position = "none", axis.title.y = element_blank(), plot.title = element_text(size=10))+
#   xlab("Catchability")+
#   ggtitle("Hierarchical Normal")
# q2

# q3 <- ggplot()+
#   geom_point(data = q_ests_Hier_NB_ANorm, aes(x=q, y = Year, color = Fishery), size = 2.5, alpha = 0.6)+
#   geom_segment(data = q_ests_Hier_NB_ANorm, aes(x= q.lower, xend = q.upper, y = Year, yend = Year, color = Fishery), linewidth = 1.2, alpha = 0.6)+
#   theme_bw()+
#   theme(axis.title.y = element_blank(), plot.title = element_text(size=10))+
#   xlab("")+
#   ggtitle("Hierarchical Asymmetric Normal")+
#   scale_color_discrete(labels = c("Chum Test", "Chinook Test"))
# q3

# plot_grid(q1, q2, q3, ncol = 3, rel_widths = c(1,1,1.4))
# ggsave(file = "Outputs/3_Plots/Catchability.png", width = 10, height = 6) 

#Make vertical 
plot_grid(q1, q2, q3, ncol = 1, align = "v", axis = "b")
ggsave(file = "Outputs/3_Plots/Catchability_vert.png", width = 3, height = 10) 


# ggplot()+
#   geom_point(data = q_ests_Hier_NB_ANorm, aes(x=Year, y = q, color = Fishery), size = 2.5, alpha = 0.7, )+
#   geom_segment(data = q_ests_Hier_NB_ANorm, aes(y= q.lower, yend = q.upper, x = Year, xend = Year, color = Fishery), linewidth = 1.2, alpha = 0.7, position = position_dodge(width = 2))+
#   theme_bw()+
#   xlab("Catchability")
