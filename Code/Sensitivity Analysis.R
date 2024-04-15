# Sensitivity analysis

library(tidyverse)
library(ggforce)
library(R2jags)
library(ggdist)
library(tidybayes)
library(DHARMa)
library(gridExtra)
library(extras)
library(ggpubr)

source("Code/SH_Functions_Final.r")

# Read in Test fishery catch data
Catch <- Load_Catch_Data()

#Read in Annual return to Albion including data up to 2022
annual_return<-read.csv(file="Data/annual_SH_return_to_Albion_1983-2022 DFO version Aug 1 to Nov 30_updated.csv",header=T)

# These are useful to have on hand
minDay <- min(Catch$Day)
minYear <- min(Catch$Year)

# Test affect of return data
# Fit models by increasing/Reducing return data
# maybe try at extremes first, doubling and cutting in half?

# Test affect of low return years/ low observation years
# Bison found that low data years biased the window shorter
# I think hierarchical structure (on q and timing) model should ameliorate this by sharing info across years

# Interestingly got estimates similar to his adjusted values for our Indep. Poiss model

# Compare return and number of non-zero steelhead obs

Yearly_Catch <- Catch %>% group_by(Year) %>% summarise(Catch = sum(SH_Catch), N_Obs= sum(SH_Catch != 0))
Return <- annual_return %>% left_join(Yearly_Catch)
plot(Return$N, Return$Catch)
plot(Return$N, Return$N_Obs)

# 13 years with les than 20 fish caught
# 10 years less than 1000 abundance -- what Bison talked about as low abund
# all of these years have less than 20 caught -- 3 would spill over if used this cutoff
# maybe num of non-zero obs?
# less than 20 catch or less than 15 obs is same 
# means we remove 2008 and 1997 (abund > 1000) which are wonky anyways (1997 because removed multipanel)

# Bison removed 1993, 2008 and 2019 -- 1993 not removed here, but 2008 and 2019 are

# call this "RLO" remove low obs
Years_RLO <- Return %>% filter(Catch > 20) %>% select(Year)

# only going to do the reccommended model - Hier NB ANorm
# fewer runs to try and save some time

# Data inputs - tried to make more flexible

Catch_RLO <- Catch
Catch_RLO$SH_Catch[!(Catch_RLO$Year %in% Years_RLO$Year)] <- NA

annual_return_RLO <- annual_return 
annual_return_RLO$N[!(annual_return_RLO$Year %in% Years_RLO$Year)] <- NA

DataList <- list ("n_years"=length(unique(Catch_RLO$Year)), 
                  "n_fisheries" = length(unique(Catch_RLO$Fishery)),
                  "n_obs" = dim(Catch_RLO)[[1]],
                  "years" = (Catch_RLO$Year-minYear+1),
                  "days" = (Catch_RLO$Day-minDay+1),
                  "fishery" = ifelse(Catch_RLO$Fishery == "CN_test", 1, 2), 
                  "obs_catch"=Catch_RLO$SH_Catch,
                  "annual_return"=annual_return_RLO$N,
                  # also include prior specification
                  "rt_m_mu" = 280-minDay, "rt_m_sig" = 40, # mean date of Oct 7
                  "halft_s" = 15  , "halft_df" = 2, 
                  "logit_q_mean" = -3, "logit_q_tau" = 0.8)

JagsFit_Hier_NB_ANorm_RLO <- jags.parallel(DataList, model.file = Hier_NB_ANorm, 
                                           n.chains =3, n.iter=10000, n.burnin = 2000, n.thin = 5, 
                                           parameters.to.save = c("rt_m", "rt_sd", "rt_m_m", "rt_sd_m", "r", "q", "pred_abund"))
#This modifies the annual_return into a vector??

All_Ests <-   data.frame(JagsFit_Hier_NB_ANorm_RLO$BUGSoutput$summary)
All_Ests$Param <- row.names(All_Ests) 
All_Ests %>% filter(Rhat>1.01) %>% filter(!str_detect(Param, "^pred_abund"))

# Save JAGS output
saveRDS(JagsFit_Hier_NB_ANorm_RLO, "Outputs/1_ModelRuns/Hier_NB_ANorm_RLO.RDS")
JagsFit_Hier_NB_ANorm_RLO <- readRDS( "Outputs/1_ModelRuns/Hier_NB_ANorm_RLO.RDS")

# Summarize JAGS outputs
Summary_Hier_NB_ANorm_RLO <- Summarize_Data(JagsFit_Hier_NB_ANorm_RLO, Catch_RLO, Hier=T, ANorm=T,
                                            Timing_Params = "rt_m[Year], rt_sd[Side,Year], rt_m_m, rt_sd_m[Side]")
# Save Summary
saveRDS(Summary_Hier_NB_ANorm_RLO, "Outputs/1_ModelRuns/Summary_Hier_NB_ANorm_RLO.RDS")
# Create Plots
Create_Plots(Summary_Data=Summary_Hier_NB_ANorm_RLO)

#-----------------------------------------------------------------
# What about changing returns?
annual_return<-read.csv(file="Data/annual_SH_return_to_Albion_1983-2022 DFO version Aug 1 to Nov 30_updated.csv",header=T)

# Increase returns------------------
# put catch back to full data set
DataList$obs_catch = Catch$SH_Catch
DataList$annual_return = round(annual_return$N*1.5)
# IR50 = increase returns 50%
JagsFit_Hier_NB_ANorm_IR50 <- jags.parallel(DataList, model.file = Hier_NB_ANorm, 
                                            n.chains =3, n.iter=15000, n.burnin = 3000, n.thin = 5, 
                                            parameters.to.save = c("rt_m", "rt_sd", "rt_m_m", "rt_sd_m", "r", "q", "pred_abund"))
All_Ests <-   data.frame(JagsFit_Hier_NB_ANorm_IR50$BUGSoutput$summary)
All_Ests$Param <- row.names(All_Ests) 
All_Ests %>% filter(Rhat>1.01) %>% filter(!str_detect(Param, "^pred_abund"))

# Save JAGS output
saveRDS(JagsFit_Hier_NB_ANorm_IR50, "Outputs/1_ModelRuns/Hier_NB_ANorm_IR50.RDS")
JagsFit_Hier_NB_ANorm_IR50 <- readRDS("Outputs/1_ModelRuns/Hier_NB_ANorm_IR50.RDS")

# Summarize JAGS outputs
Summary_Hier_NB_ANorm_IR50 <- Summarize_Data(JagsFit_Hier_NB_ANorm_IR50, Catch, Hier=T, ANorm=T,
                                             Timing_Params = "rt_m[Year], rt_sd[Side,Year], rt_m_m, rt_sd_m[Side]")
# Save Summary
saveRDS(Summary_Hier_NB_ANorm_IR50, "Outputs/1_ModelRuns/Summary_Hier_NB_ANorm_IR50.RDS")
# Create Plots
Create_Plots(Summary_Data=Summary_Hier_NB_ANorm_IR50)                

# Decrease returns ---------------
# I don't know why annual_return keeps changing to a vector?
annual_return<-read.csv(file="Data/annual_SH_return_to_Albion_1983-2022 DFO version Aug 1 to Nov 30_updated.csv",header=T)
DataList$annual_return = round(annual_return$N*0.5)
# DR50 = increase returns 50%
JagsFit_Hier_NB_ANorm_DR50 <- jags.parallel(DataList, model.file = Hier_NB_ANorm, 
                                            n.chains =3, n.iter=15000, n.burnin = 3000, n.thin = 5, 
                                            parameters.to.save = c("rt_m", "rt_sd", "rt_m_m", "rt_sd_m", "r", "q", "pred_abund"))
All_Ests <-   data.frame(JagsFit_Hier_NB_ANorm_DR50$BUGSoutput$summary)
All_Ests$Param <- row.names(All_Ests) 
All_Ests %>% filter(Rhat>1.01) %>% filter(!str_detect(Param, "^pred_abund"))

# Save JAGS output
saveRDS(JagsFit_Hier_NB_ANorm_DR50, "Outputs/1_ModelRuns/Hier_NB_ANorm_DR50.RDS")
JagsFit_Hier_NB_ANorm_DR50 <-  readRDS( "Outputs/1_ModelRuns/Hier_NB_ANorm_DR50.RDS")

# Summarize JAGS outputs
Summary_Hier_NB_ANorm_DR50 <- Summarize_Data(JagsFit_Hier_NB_ANorm_DR50, Catch, Hier=T, ANorm=T,
                                             Timing_Params = "rt_m[Year], rt_sd[Side,Year], rt_m_m, rt_sd_m[Side]")
# Save Summary
saveRDS(Summary_Hier_NB_ANorm_DR50, "Outputs/1_ModelRuns/Summary_Hier_NB_ANorm_DR50.RDS")
# Create Plots
Create_Plots(Summary_Data=Summary_Hier_NB_ANorm_DR50)    
#-----------------------------------------------------------------------------
# Test for linear trend in mean date
annual_return<-read.csv(file="Data/annual_SH_return_to_Albion_1983-2022 DFO version Aug 1 to Nov 30_updated.csv",header=T)
DataList$annual_return = round(annual_return$N)

JagsFit_Hier_TestMu <- jags.parallel(DataList, model.file = Hier_NB_ANorm_TestMu, 
                                     n.chains =3, n.iter=20000, n.burnin = 2000, n.thin = 5, 
                                     parameters.to.save = c("rt_m", "rt_sd", "rt_sd_m", "r", "q", "rt_m_Int", "rt_m_Slope"))
All_Ests <-   data.frame(JagsFit_Hier_TestMu$BUGSoutput$summary)
All_Ests$Param <- row.names(All_Ests) 
All_Ests %>% filter(Rhat>1.01) %>% filter(!str_detect(Param, "^pred_abund"))

All_Ests %>% filter(str_detect(Param, "^rt_m_Slope"))

JagsFit_Hier_TestMu$BUGSoutput$DIC

# Save JAGS output
saveRDS(JagsFit_Hier_TestMu, "Outputs/1_ModelRuns/JagsFit_Hier_TestMu.RDS")
JagsFit_Hier_TestMu <-  readRDS( "Outputs/1_ModelRuns/JagsFit_Hier_TestMu.RDS")

# plot rt_m over time, with fit line, also posterior of rt_m_Slope

rt_m_s_sims <- JagsFit_Hier_TestMu$BUGSoutput$sims.matrix[, "rt_m_Slope"]
low = All_Ests$X2.5.[All_Ests$Param =="rt_m_Slope" ]
high = All_Ests$X97.5.[All_Ests$Param =="rt_m_Slope" ]
pdf("Outputs/3_Plots/RBeta_Slope.pdf", height=6, width=6) 
hist(rt_m_s_sims, xlab = bquote(beta[1]))
abline(v=c(low, high), col = "red", lty="dashed")
dev.off()


#French version
#Change decimal to comma for base r plot
options(OutDec = ",")
pdf("Outputs/3_Plots/RBeta_Slope_french.pdf", height=6, width=6) 
hist(rt_m_s_sims, xlab = bquote(beta[1]), ylab = "Fréquence", main =NULL)
abline(v=c(low, high), col = "red", lty="dashed")
dev.off()
#Change back to decimal
options(OutDec = ".")


JagsFit_MCMC <- as.mcmc(JagsFit_Hier_TestMu)

Timing_Params <- "rt_m[Year], rt_sd[Side, Year], rt_sd_m[Side]"

minYear <- min(Catch$Year)
minDay <- min(Catch$Day)

# don't want to have to complicate summarize_data function for this, will just manually pull out what I want
Date_Dist_Draws <- eval(parse(text = paste("tidybayes::spread_draws(JagsFit_MCMC,", Timing_Params, ")")))
Date_Dist_Ints <- point_interval(Date_Dist_Draws) %>%
  # Put years and dates back to normal
  mutate(Year = Year + minYear - 1) 
Date_Dist_Ints$Model <- "Timing Trend"

# compare to best model
Summary_Hier_NB_ANorm <- readRDS("Outputs/1_ModelRuns/Summary_Hier_NB_ANorm.RDS")
Date_Dists_Best <- Summary_Hier_NB_ANorm$Date_Dists_Ints %>% mutate(Model = "Hier. Timing")  

# Add indep too
Summary_Indep_Poisson <- readRDS("Outputs/1_ModelRuns/Summary_Indep_Poisson.RDS")
Date_Dists_Indep <- Summary_Indep_Poisson$Date_Dists_Ints %>% mutate(Model = "Indep Poiss.")
Date_Dists_Comp <- rbind(Date_Dist_Ints, Date_Dists_Best, Date_Dists_Indep) 
# join with return
# also just join all together
RT_All <- Return %>% left_join(Date_Dists_Comp, join_by(Year))

pdf("Outputs/3_Plots/RT_TrendvsBest.pdf", height=4, width=8)  
ggplot(RT_All%>% filter(Model %in% c("Hier. Timing" ,"Timing Trend"))) +
  geom_pointrange(aes(x=Year, y=rt_m, ymin=rt_m.lower, ymax = rt_m.upper, col=Model, alpha=N_Obs, shape=Model))+
  ylab("Mean Return Date") +
  theme_bw()
dev.off()

#French
pdf("Outputs/3_Plots/RT_TrendvsBest_french.pdf", height=4, width=8)  
ggplot(RT_All%>% filter(Model %in% c("Hier. Timing" ,"Timing Trend"))) +
  geom_pointrange(aes(x=Year, y=rt_m, ymin=rt_m.lower, ymax = rt_m.upper, col=Model, alpha=N_Obs, shape=Model))+
  ylab("Date de montaison moyenne") +
  xlab("Année")+
  theme_bw()+
  labs(fill  = "Modèle", color = "Modèle", shape = "Modèle", alpha = "Nombre d'obs.")+
  scale_shape_discrete(labels = c("Hiérarchique", "Tendence linéaire"))+
  scale_fill_discrete(labels = c("Hiérarchique", "Tendence linéaire"))+
  scale_color_discrete(labels = c("Hiérarchique", "Tendence linéaire"))
dev.off()

# Get window for most recent year, just to see what looks like
percs <- c(0.025, 0.05, 0.1, 0.5, 0.9, 0.95, 0.975)
perc_names <- c("D_025", "D_05", "D_10", "D_50", "D_90", "D_95", "D_975")
# If Anorm spread so that have sd+ and sd- all in same row

Param_Draws <- Date_Dist_Draws %>%
  filter(Year == 40) %>%
  select(!c(rt_sd_m)) %>%
  pivot_wider(values_from = rt_sd, names_from = Side, names_prefix = "rt_sd") 

N_Draws <- dim(Param_Draws)[1]
Win <- matrix(nrow=N_Draws, ncol = length(percs))
for(i in 1:N_Draws){
  Win[i, ] <- qdnorm(p = percs, mean = Param_Draws$rt_m[i], sigma= as.numeric(Param_Draws[i,c("rt_sd1", "rt_sd2")])) + minDay-1
}
Win_DF <- as.data.frame(Win)
names(Win_DF) <- perc_names
Window_Draws <- cbind(Param_Draws, Win_DF) %>% select(!c("rt_m", "rt_sd1", "rt_sd2", "Year"))
Global_Windows <- point_interval(Window_Draws)

Pred_Draws <- tidybayes::spread_draws(JagsFit_MCMC, q[Fishery,Year], pred_abund[Day, Year], r, ndraws = Dharma_draws)  %>%
  # Now get to get Pred_Catch 
  mutate(Pred_Catch = pred_abund*q) %>% 
  mutate(p_NB = r/(Pred_Catch + r)) %>%
  mutate(Year = Year + MinYear - 1) %>%
  mutate(Day = Day + MinDay -1) 

# Now simulate to Post_Pred_Catch
Pred_Draws$Post_Pred_Catch <- rnbinom(dim(Pred_Draws)[1], prob=Pred_Draws$p_NB, size = Pred_Draws$r)
Pred_Draws$Fishery = ifelse(Pred_Draws$Fishery==1, "CN_test", "CM_test")

Pred_Ints <- point_interval( Pred_Draws) %>%
  left_join(Catch[, c("Year", "Day", "Fishery", "SH_Catch")])

# Now pull in actual catches and merge together
All_Sims <- Pred_Draws %>% 
  group_by(Fishery, Year, Day) %>%
  select(Fishery, Year, Day, Post_Pred_Catch, .draw) %>% 
  pivot_wider(values_from = Post_Pred_Catch, names_from= .draw) %>%
  # now add true catches, and filter for those, and also add pred_catch
  left_join(Pred_Ints[c("Year", "Day", "Fishery", "Pred_Catch")]) %>%
  left_join(Catch[, c("Year", "Day", "Fishery", "SH_Catch")]) %>%
  relocate(Pred_Catch, .after = Day) %>%
  relocate(SH_Catch, .after = Pred_Catch) %>%
  filter(is.na(SH_Catch) == F)

Sim_Resp_Matrix <- as.matrix(All_Sims[, 6:(Dharma_draws+5)])

sims <- DHARMa::createDHARMa(simulatedResponse = Sim_Resp_Matrix, 
                             observedResponse = All_Sims$SH_Catch, 
                             fittedPredictedResponse = All_Sims$Pred_Catch, 
                             integerResponse = T)


# Summarize these three scenarios --------------------------------------------------------
# get windows - just do global for now, don't really use yearly ones

Summary_Hier_NB_ANorm_RLO <- readRDS( "Outputs/1_ModelRuns/Summary_Hier_NB_ANorm_RLO.RDS")
Summary_Hier_NB_ANorm_IR50 <- readRDS( "Outputs/1_ModelRuns/Summary_Hier_NB_ANorm_IR50.RDS")
Summary_Hier_NB_ANorm_DR50 <- readRDS( "Outputs/1_ModelRuns/Summary_Hier_NB_ANorm_DR50.RDS")


#Pull out windows from summaries 
Yearly1 <- Summary_Hier_NB_ANorm_RLO$Yearly_Windows %>% mutate(Model = "ANorm_RLO")
Yearly2 <- Summary_Hier_NB_ANorm_IR50$Yearly_Windows %>% mutate(Model = "ANorm_IR50") 
Yearly3 <- Summary_Hier_NB_ANorm_DR50$Yearly_Windows %>% mutate(Model = "ANorm_DR50") 
Yearly_Winds <- rbind(Yearly3, Yearly2, Yearly1) %>% mutate_if(is.numeric, round) %>% mutate(Year = Year+MinYear-1)
write.csv(Yearly_Winds, "Outputs/2_DataOut/Yearly_Windows_Sens.csv")

Global1 <- Summary_Hier_NB_ANorm_RLO$Global_Windows %>% mutate(Model = "ANorm_RLO") %>% select(-Year) %>% select(-starts_with("."))
Global2 <- Summary_Hier_NB_ANorm_IR50$Global_Windows %>% mutate(Model = "ANorm_IR50") %>% select(-Year) %>% select(-starts_with("."))
Global3 <- Summary_Hier_NB_ANorm_DR50$Global_Windows %>% mutate(Model = "ANorm_DR50") %>% select(-Year) %>% select(-starts_with("."))

Global_Winds <- rbind(Global3, Global2, Global1) %>% mutate_if(is.numeric, round)
toDate <- function(x){as.Date(x-1, origin = as.Date("2023-01-01"))}
Global_Winds_Dates <- Global_Winds %>% mutate_if(is.numeric, toDate)
write.csv(Global_Winds, "Outputs/2_DataOut/Global_Windows_Sens.csv")
write.csv(Global_Winds_Dates, "Outputs/2_DataOut/Global_Windows_Dates_Sens.csv")

# also print date dists for each
write.csv(Summary_Hier_NB_HQ$Date_Dists_Ints, "Outputs/2_DataOut/Date_Dist_Ints_Norm.csv")
write.csv(Summary_Hier_NB_ANorm$Date_Dists_Ints, "Outputs/2_DataOut/Date_Dist_Ints_ANorm.csv")
write.csv(Summary_Indep_Poisson$Date_Dists_Ints, "Outputs/2_DataOut/Date_Dist_Ints_Indep.csv")

#-------------------------------------------------------------------------------------
# want to look at trends in sigma, mean over time
# start with 2-piece normal

RT_ANorm <- Return %>% left_join(Date_Dists_Best, join_by(Year))

ggplot(RT_ANorm, aes(x=N, y=rt_sd)) +
  geom_pointrange(aes(ymax=rt_sd.upper, ymin = rt_sd.lower)) +
  facet_wrap(~Side) +
  geom_smooth(method = "lm", se = TRUE) +
  stat_cor(label.x = max(RT_ANorm$N)*0.6, label.y = max(RT_ANorm$rt_sd.upper)) +
  stat_regline_equation(label.x = max(RT_ANorm$N)*0.6, label.y = max(RT_ANorm$rt_sd.upper)-5)

# left side quite constant
# right side more variable
pdf("Outputs/3_Plots/ANorm_sd_N_Obs.pdf", height=4, width=8) 
ggplot(RT_ANorm, aes(x=N_Obs, y=rt_sd)) +
  geom_pointrange(aes( ymax=rt_sd.upper, ymin = rt_sd.lower)) +
  facet_wrap(~Side) +
  geom_smooth(method = "lm", se = TRUE) +
  stat_cor(label.x = max(RT_ANorm$N_Obs)*0.6, label.y = max(RT_ANorm$rt_sd.upper)) +
  stat_regline_equation(label.x = max(RT_ANorm$N_Obs)*0.6, label.y = max(RT_ANorm$rt_sd.upper)-5)+
  theme_bw() +
  ylab(bquote(sigma[y])) + xlab("Number of Non-zero Obs.")
# similar pattern
dev.off()

# what if remove "bad" years
ggplot(RT_ANorm %>% filter(Year %in% Years_RLO$Year),aes(x=N, y=rt_sd) ) +
  geom_pointrange(aes( ymax=rt_sd.upper, ymin = rt_sd.lower)) +
  facet_wrap(~Side)+
  geom_smooth(method = "lm", se = TRUE) +
  stat_cor(label.x = max(RT_ANorm$N_Obs)*0.6, label.y = max(RT_ANorm$rt_sd.upper)) +
  stat_regline_equation(label.x = max(RT_ANorm$N_Obs)*0.6, label.y = max(RT_ANorm$rt_sd.upper)-5)

# do we see trend over time in mean dates or sigmas?
ggplot(RT_ANorm, aes(x=Year, y=rt_sd)) +
  geom_pointrange(aes(ymax=rt_sd.upper, ymin = rt_sd.lower)) +
  facet_wrap(~Side)+
  geom_smooth(method = "lm", se = TRUE) +
  stat_cor(label.x = max(RT_ANorm$Year)-20, label.y = max(RT_ANorm$rt_sd.upper)) +
  stat_regline_equation(label.x = max(RT_ANorm$Year)-20, label.y = max(RT_ANorm$rt_sd.upper)-5)

ggplot(RT_ANorm, aes(x=Year, y=rt_m)) +
  geom_pointrange(aes(ymax=rt_m.upper, ymin = rt_m.lower)) +
  geom_smooth(method = "lm", se = TRUE) +
  stat_cor(label.x = max(RT_ANorm$Year)-20, label.y = max(RT_ANorm$rt_m.upper)) +
  stat_regline_equation(label.x = max(RT_ANorm$Year)-20, label.y = max(RT_ANorm$rt_m.upper)-5)


# remove "bad" years?
ggplot(RT_ANorm%>% filter(Year %in% Years_RLO$Year), aes(x=Year, y=rt_m)) +
  geom_pointrange(aes( ymax=rt_m.upper, ymin = rt_m.lower)) +
  geom_smooth(method = "lm", se = TRUE) +
  stat_cor(label.x = max(RT_ANorm$Year)-20, label.y = max(RT_ANorm$rt_m.upper)) +
  stat_regline_equation(label.x = max(RT_ANorm$Year)-20, label.y = max(RT_ANorm$rt_m.upper)-5)

ggplot(RT_ANorm%>% filter(Year %in% Years_RLO$Year), aes(x=Year, y=rt_m)) +
  geom_pointrange(aes( ymax=rt_m.upper, ymin = rt_m.lower)) +
  geom_smooth(method = "lm", se = TRUE) +
  stat_cor(label.x = max(RT_ANorm$Year)-20, label.y = max(RT_ANorm$rt_m.upper)) +
  stat_regline_equation(label.x = max(RT_ANorm$Year)-20, label.y = max(RT_ANorm$rt_m.upper)-5)

# if remove years of low obs, it goes away

ggplot(RT_ANorm %>% filter(!(Year %in% c(2008, 2022))), aes(x=Year, y=rt_m)) +
  geom_pointrange(aes( ymax=rt_m.upper, ymin = rt_m.lower)) +
  geom_smooth(method = "lm", se = TRUE) +
  stat_cor(label.x = max(RT_ANorm$Year)-20, label.y = max(RT_ANorm$rt_m.upper)) +
  stat_regline_equation(label.x = max(RT_ANorm$Year)-20, label.y = max(RT_ANorm$rt_m.upper)-5)


# In paper we talk about this right after the indep Poiss. model, maybe should look at that along?

RT_Indep <- RT_All %>% filter(Model == "Indep Poiss.")

# do we see trend over time in mean dates or sigmas?
pdf("Outputs/3_Plots/Indep_Poiss_sd_over_time.pdf", height=4, width=8)  
ggplot(RT_Indep, aes(x=Year, y=rt_sd)) +
  geom_pointrange(aes(ymax=rt_sd.upper, ymin = rt_sd.lower)) +
  geom_smooth(method = "lm", se = TRUE) +
  stat_cor(label.x = max(RT_Indep$Year)-20, label.y = max(RT_Indep$rt_sd.upper)) +
  stat_regline_equation(label.x = max(RT_Indep$Year)-20, label.y = max(RT_Indep$rt_sd.upper)-5)+
  ylab(bquote(sigma[y])) 
theme_bw()
# minimal positive trend,not significant
dev.off()

pdf("Outputs/3_Plots/Indep_Poiss_sd_N_Obs.pdf", height=4, width=8)  
ggplot(RT_Indep, aes(x=N_Obs, y=rt_sd)) +
  geom_pointrange(aes(ymax=rt_sd.upper, ymin = rt_sd.lower)) +
  geom_smooth(method = "lm", se = TRUE) +
  stat_cor(label.x = max(RT_Indep$N_Obs)*0.6, label.y = max(RT_Indep$rt_sd.upper)) +
  stat_regline_equation(label.x = max(RT_Indep$N_Obs)*0.6, label.y = max(RT_Indep$rt_sd.upper)-5)+
  ylab(bquote(sigma[y])) +  xlab("Number of Non-zero Obs.")+
  theme_bw()
dev.off()

pdf("Outputs/3_Plots/Indep_Poiss_m_over_time.pdf", height=4, width=8)
ggplot(RT_Indep, aes(x=Year, y=rt_m)) +
  geom_pointrange(aes(ymax=rt_m.upper, ymin = rt_m.lower)) +
  geom_smooth(method = "lm", se = TRUE) +
  stat_cor(label.x = max(RT_Indep$Year)-20, label.y = max(RT_Indep$rt_m.upper)) +
  theme_light() +
  stat_regline_equation(label.x = max(RT_Indep$Year)-20, label.y = max(RT_Indep$rt_m.upper)-5) +
  ylab(bquote(M[y])) +
  theme_bw()
#  negative trend, but not significant
dev.off()

# do we see trend over time in mean dates or sigmas

# curious if remove 2008
ggplot(RT_Indep %>% filter(Year != 2008), aes(x=Year, y=rt_m)) +
  geom_pointrange(aes(ymax=rt_m.upper, ymin = rt_m.lower)) +
  geom_smooth(method = "lm", se = TRUE) +
  stat_cor(label.x = max(RT_Indep$Year)-20, label.y = max(RT_Indep$rt_m.upper)) +
  stat_regline_equation(label.x = max(RT_Indep$Year)-20, label.y = max(RT_Indep$rt_m.upper)-5)

# remove years Bison did
ggplot(RT_Indep %>% filter(!(Year %in% c(1993, 2008, 2019))), aes(x=Year, y=rt_m)) +
  geom_pointrange(aes(ymax=rt_m.upper, ymin = rt_m.lower)) +
  geom_smooth(method = "lm", se = TRUE) +
  stat_cor(label.x = max(RT_Indep$Year)-20, label.y = max(RT_Indep$rt_m.upper)) +
  stat_regline_equation(label.x = max(RT_Indep$Year)-20, label.y = max(RT_Indep$rt_m.upper)-5)
# significant decline here

# remove years with less than 20 obs
ggplot(RT_Indep %>% filter(Year %in% Years_RLO$Year), aes(x=Year, y=rt_m)) +
  geom_pointrange(aes(ymax=rt_m.upper, ymin = rt_m.lower)) +
  geom_smooth(method = "lm", se = TRUE) +
  stat_cor(label.x = max(RT_Indep$Year)-20, label.y = max(RT_Indep$rt_m.upper)) +
  stat_regline_equation(label.x = max(RT_Indep$Year)-20, label.y = max(RT_Indep$rt_m.upper)-5)
# not sig





#####Extras sensitivity analyses for CSAS meeting: 
#By Anna

#Restarting the code here so original functions and datalist are read in
#Load packages 
library(cowplot)

#Source the functions for models, summaries, etc 
source("Code/SH_Functions_Final.r")

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

#Read in and get the windows for the original models to compare to
Summary_Indep_Poisson <- readRDS("Outputs/1_ModelRuns/Summary_Indep_Poisson.RDS")
#Get windows
windows_indep <- Summary_Indep_Poisson$Yearly_Windows
#Recalc the year 
windows_indep$Year <- windows_indep$Year+minYear-1 

Summary_Hier_NB_HQ <- readRDS("Outputs/1_ModelRuns/Summary_Hier_NB_HQ.RDS")
windows_nb <- Summary_Hier_NB_HQ$Yearly_Windows
windows_nb$Year <- windows_nb$Year+minYear-1 

Summary_Hier_NB_ANorm <- readRDS("Outputs/1_ModelRuns/Summary_Hier_NB_ANorm.RDS")
windows_asym <- Summary_Hier_NB_ANorm$Yearly_Windows
windows_asym$Year <- windows_asym$Year+minYear-1 




#Sensitivity of prior for mean date 
#Non-hierarchical poisson with uniform prior on rt_m
Non_Hier_Poisson_Q_unif <- function(){
  
  # global mean catchability priors
  for(i in 1:n_fisheries){
    q_mean[i] ~ dnorm(logit_q_mean, logit_q_tau)
    # precision in catchabilities across years
    q_tau[i] ~ dscaled.gamma(halft_s, halft_df)
  }
  
  #Likelihood
  for (y in 1:n_years){
    #Priors
    # Independent average and sd for each year
    rt_m[y] ~ dunif(1,365)
    # scaled gamma on tau gives half-t on sd
    rt_tau[y] ~ dscaled.gamma(halft_s, halft_df)
    rt_sd[y] <- 1/sqrt(rt_tau[y])
    # Prior for Catchability
    logit_q_1[y] ~ dnorm(q_mean[1], q_tau[1])
    logit_q_2[y] ~ dnorm(q_mean[2], q_tau[2])
    q[1,y] <- ilogit(logit_q_1[y])
    q[2,y] <- ilogit(logit_q_2[y])
  }
  
  
  for (k in 1:n_obs) {
    # Observation process
    obs_catch[k] ~ dpois(q[fishery[k],years[k]]*pred_abund[days[k],years[k]])
  }
  
  #Predicted abundance (state process)
  for(y in 1:n_years){
    NC[y] <-  (1/(rt_sd[y]*sqrt(44/7)))
    for(d in 1:max(days)){
      pred_abund[d,y] <- annual_return[y]* NC[y] * exp(-pow((d-rt_m[y])/rt_sd[y], 2)/2)
    }}
  
} 


JagsFit_Indep_Poisson <- jags.parallel(DataList, model.file =Non_Hier_Poisson_Q_unif, 
                                       n.chains =3, n.iter = 20000, n.burnin = 5000, n.thin = 5,
                                       parameters.to.save = c("rt_m", "rt_sd", "q", "pred_abund"))
# Check convergence
All_Ests <-   data.frame(JagsFit_Indep_Poisson$BUGSoutput$summary)
All_Ests$Param <- row.names(All_Ests) 
All_Ests %>% filter(Rhat>1.01)

Summary_Indep_Poisson_unif <- Summarize_Data(JagsFit_Indep_Poisson, Catch, Hier=F, ANorm=F,
                                             Timing_Params = "rt_m[Year], rt_sd[Year]")
saveRDS(Summary_Indep_Poisson_unif, "Outputs/1_ModelRuns/Summary_Indep_Poisson_unif.RDS")

#Get windows
windows_indep_unif <- Summary_Indep_Poisson_unif$Yearly_Windows
windows_indep_unif$Year <- windows_indep_unif$Year+minYear-0.8
#Plot comparison
ggplot()+
  geom_pointrange(data = windows_indep, aes(x = Year, y = hD_50, ymin = hD_025, ymax = hD_975, color = "Normal Prior"))+
  geom_pointrange(data = windows_indep_unif, aes(x = Year, y = D_50, ymin = D_025, ymax = D_975, color = "Uniform Prior"))+
  ggtitle("Indep. Normal, Poisson")+
  ylab("Date")+
  scale_y_continuous(breaks=c(213,244, 274, 305, 335),
                     labels=c("Aug 1", "Sept 1", "Oct 1", "Nov 1", "Dec 1"))+
  theme_bw()+
  scale_colour_manual("", values = c("Normal Prior"="dodgerblue", "Uniform Prior"="orange"))
ggsave(file = "Outputs/3_Plots/uniform_prior_Indep_Poisson.png", width = 8, height = 5, units = "in")



#---------------------------------------------------------------------------

# Hier timing, NB error model

Hier_NB_HQ_unif <- function(){
  
  #Likelihood
  # Global dists
  # Mean and sd around mean date global dist
  rt_m_m ~ dunif(1,365)
  # spread in means between years lower than 
  # spread in dates
  rt_m_tau ~ dscaled.gamma(halft_s/2, halft_df)
  rt_m_sd <- 1/sqrt(rt_m_tau)
  # Mean and sd around sdrt_m_tau global dist
  rt_tau_m ~ dscaled.gamma(halft_s, halft_df)
  rt_sd_m <- 1/sqrt(rt_tau_m)
  # also going to be smaller than spread, so reduce scale
  # most of weight will be between 0 and 7.5
  rt_sd_tau ~ dscaled.gamma(halft_s/2, halft_df)
  rt_sd_sd <- 1/sqrt(rt_sd_tau)
  
  # global mean catchability priors
  for(i in 1:n_fisheries){
    q_mean[i] ~ dnorm(logit_q_mean, logit_q_tau)
    # precision in catchabilities across years
    q_tau[i] ~ dscaled.gamma(halft_s, halft_df)
  }
  
  r ~ dunif(0.01,10)
  
  for (y in 1:n_years){
    #2a. Priors
    # m and sd from "global" shared dist
    rt_m[y] ~ dnorm(rt_m_m, rt_m_tau) # had full range for wider priors
    rt_sd[y] ~ dnorm(rt_sd_m, rt_sd_tau) # had up to 100 in wider prior
    # Priors for Catchability
    logit_q_1[y] ~ dnorm(q_mean[1], q_tau[1])
    logit_q_2[y] ~ dnorm(q_mean[2], q_tau[2])
    q[1,y] <- ilogit(logit_q_1[y])
    q[2,y] <- ilogit(logit_q_2[y])
  }
  
  for (k in 1:n_obs) {
    # Observation process
    obs_catch[k] ~ dnegbin(p[k],r)
    p[k] <- r/(r + q[fishery[k], years[k]]*pred_abund[days[k],years[k]]) 
  }
  
  #Predicted abundance (state process)
  for(y in 1:n_years){
    NC[y] <- (1/(rt_sd[y]*sqrt(44/7)))
    for(d in 1:max(days)){
      pred_abund[d,y] <- annual_return[y]* NC[y] *exp(-pow((d-rt_m[y])/rt_sd[y], 2)/2)
    }}
  
} 

JagsFit_Hier_NB_HQ <- jags.parallel(DataList, model.file = Hier_NB_HQ_unif, 
                                    n.chains =3, n.iter=20000, n.burnin = 5000, n.thin = 5, 
                                    parameters.to.save = c("rt_m", "rt_sd", "r", "q", "pred_abund", "rt_m_m", "rt_sd_m"))

All_Ests_Hier_NB_HQ <-   data.frame(JagsFit_Hier_NB_HQ$BUGSoutput$summary)
All_Ests_Hier_NB_HQ$Param <- row.names(All_Ests_Hier_NB_HQ) 
All_Ests_Hier_NB_HQ %>% filter(Rhat>1.01)

Summary_Hier_NB_HQ_unif <- Summarize_Data(JagsFit_Hier_NB_HQ, Catch, Timing_Params = "rt_m[Year], rt_sd[Year], rt_m_m, rt_sd_m", Hier=T, ANorm=F)
saveRDS(Summary_Hier_NB_HQ_unif, "Outputs/1_ModelRuns/Summary_Hier_NB_HQ_unif.RDS")
rm(JagsFit_Hier_NB_HQ)

#Get windows
windows_nb_unif <- Summary_Hier_NB_HQ_unif$Yearly_Windows
#Recalc the year 
windows_nb_unif$Year <- windows_nb_unif$Year+minYear-0.8
#Plot comparison
ggplot()+
  geom_pointrange(data = windows_nb, aes(x = Year, y = hD_50, ymin = hD_025, ymax = hD_975, color = "Normal Prior"))+
  geom_pointrange(data = windows_nb_unif, aes(x = Year, y = D_50, ymin = D_025, ymax = D_975, color = "Uniform Prior"))+
  ggtitle("Hier. Normal, Neg. Binomial")+
  ylab("Date")+
  scale_y_continuous(breaks=c(213,244, 274, 305, 335),
                     labels=c("Aug 1", "Sept 1", "Oct 1", "Nov 1", "Dec 1"))+
  theme_bw()+
  scale_colour_manual("", values = c("Normal Prior"="dodgerblue", "Uniform Prior"="orange"))
ggsave(file = "Outputs/3_Plots/uniform_prior_Hier_NB_HQ.png", width = 8, height = 5, units = "in")




#----------------------------------------------------------
#Asymmetric normal, Hier timing, NB Obs model
Hier_NB_ANorm_unif <- function(){
  
  # Global dists
  # Mean and sd around mean date global dist
  rt_m_m ~ dunif(1, 365)
  # spread in means between years lower than 
  # spread in dates
  rt_m_tau ~ dscaled.gamma(halft_s/2, halft_df)
  rt_m_sd <- 1/sqrt(rt_m_tau)
  
  r ~ dunif(0.01,10)
  
  #Priors
  for (y in 1:n_years){
    # q's for each fishery from "global" dist
    logit_q_1[y] ~ dnorm(q_mean[1], q_tau[1])
    logit_q_2[y] ~ dnorm(q_mean[2], q_tau[2])
    q[1,y] <- ilogit(logit_q_1[y])
    q[2,y] <- ilogit(logit_q_2[y])
    # For each year select SD's and means from global dist
    #2a. Priors
    # m and sd from "global" shared dist
    rt_m[y] ~ dnorm(rt_m_m, rt_m_tau) 
    rt_sd[1, y] ~ dnorm(rt_sd_m[1], rt_sd_tau[1]) 
    rt_sd[2, y] ~ dnorm(rt_sd_m[2], rt_sd_tau[2])
  }
  
  # now we have two sd's -- one for either side of normal
  # also have two fisheries -- not super flexible, but efficient for this case to have one loop
  for(i in 1:2){
    # Mean and sd around sdrt_m_tau global dist
    rt_tau_m[i] ~ dscaled.gamma(halft_s, halft_df)
    rt_sd_m[i] <- 1/sqrt(rt_tau_m[i])
    # also going to be smaller than spread, so reduce scale
    rt_sd_tau[i] ~ dscaled.gamma(halft_s/2, halft_df)
    rt_sd_sd[i] <- 1/sqrt(rt_sd_tau[i])
    # might as well just put qs in here too
    q_mean[i] ~ dnorm(logit_q_mean, logit_q_tau)
    # precision in catchabilities across years
    q_tau[i] ~ dscaled.gamma(halft_s, halft_df)
  }
  
  for (k in 1:n_obs) {
    # Observation process
    obs_catch[k] ~ dnegbin(p[k],r)
    p[k] <- r/(r + q[fishery[k], years[k]]*pred_abund[days[k],years[k]])
  }
  
  
  #Predicted abundance (state process)
  for(y in 1:n_years){
    iNC[y] <- (sqrt(2*22/7)*(rt_sd[1,y] + rt_sd[2,y])/2)^-1 
    for(d in 1:max(days)){
      varObs[d,y] <- ifelse(d <= rt_m[y], rt_sd[1,y]^2, rt_sd[2,y]^2)
      daily_prop[d,y] <- iNC[y] * exp(-0.5*(d-rt_m[y])^2/varObs[d,y])
      pred_abund[d,y] <- annual_return[y]*daily_prop[d,y]
    }
  } # end year loop
  
} 


JagsFit_Hier_NB_ANorm <- jags.parallel(DataList, model.file = Hier_NB_ANorm_unif, 
                                       n.chains =3, n.iter=30000, n.burnin = 5000, n.thin = 5, 
                                       parameters.to.save = c("rt_m", "rt_sd", "rt_m_m", "rt_sd_m", "r", "q", "pred_abund"))
All_Ests_Hier_NB_ANorm <-   data.frame(JagsFit_Hier_NB_ANorm$BUGSoutput$summary)
All_Ests_Hier_NB_ANorm$Param <- row.names(All_Ests_Hier_NB_ANorm) 
All_Ests_Hier_NB_ANorm %>% filter(Rhat>1.01) %>% filter(!str_detect(Param, "^pred_abund"))


Summary_Hier_NB_ANorm_unif <- Summarize_Data(JagsFit_Hier_NB_ANorm, Catch, Hier=T, ANorm=T,
                                             Timing_Params = "rt_m[Year], rt_sd[Side,Year], rt_m_m, rt_sd_m[Side]")
saveRDS(Summary_Hier_NB_ANorm_unif, "Outputs/1_ModelRuns/Summary_Hier_NB_ANorm_unif.RDS")

#Get windows
windows_asym_unif <- Summary_Hier_NB_ANorm_unif$Yearly_Windows
#Recalc the year 
windows_asym_unif$Year <- windows_asym_unif$Year+minYear-0.8 #trick to get the offset in the plot
#Plot comparison
ggplot()+
  geom_pointrange(data = windows_asym, aes(x = Year, y = hD_50, ymin = hD_025, ymax = hD_975, color = "Normal Prior"))+
  geom_pointrange(data = windows_asym_unif, aes(x = Year, y = D_50, ymin = D_025, ymax = D_975, color = "Uniform Prior"))+
  ggtitle("Hier. Asym. Normal, Neg. Binomial")+
  ylab("Date")+
  scale_y_continuous(breaks=c(213,244, 274, 305, 335),
                     labels=c("Aug 1", "Sept 1", "Oct 1", "Nov 1", "Dec 1"))+
  theme_bw()+
  scale_colour_manual("", values = c("Normal Prior"="dodgerblue", "Uniform Prior"="orange"))

ggsave(file = "Outputs/3_Plots/uniform_prior_Hier_NB_ANorm.png", width = 8, height = 5, units = "in")



#Get the windows from the summaries to compare the unif to the normal priors 

# #Get the windows for each year from the summary objects 
# Yearly1_unif <- Summary_Indep_Poisson_unif$Yearly_Windows %>% mutate(Model = "Indep_Poisson")
# Yearly2_unif <- Summary_Hier_NB_HQ_unif$Yearly_Windows %>% mutate(Model = "Hier_NB_HQ") 
# Yearly3_unif <- Summary_Hier_NB_ANorm_unif$Yearly_Windows %>% mutate(Model = "Hier_NB_ANorm") 
# Yearly_Winds_unif <- rbind(Yearly3_unif, Yearly2_unif, Yearly1_unif) %>% mutate_if(is.numeric, round) %>% mutate(Year = Year+minYear-1)

#Get the global windows from the summary objects 
Global1 <- Summary_Indep_Poisson$Global_Windows %>% mutate(Model = "Indep_Poisson") %>% select(-starts_with("."))
Global2 <- Summary_Hier_NB_HQ$Global_Windows %>% mutate(Model = "Hier_NB_HQ") %>% select(-Year) %>% select(-starts_with("."))
Global3 <- Summary_Hier_NB_ANorm$Global_Windows %>% mutate(Model = "Hier_NB_ANorm") %>% select(-Year) %>% select(-starts_with("."))
Global_Winds <- rbind(Global3, Global2, Global1) %>% mutate_if(is.numeric, round)

Global1_unif <- Summary_Indep_Poisson_unif$Global_Windows %>% mutate(Model = "Indep_Poisson") %>% select(-starts_with("."))
Global2_unif <- Summary_Hier_NB_HQ_unif$Global_Windows %>% mutate(Model = "Hier_NB_HQ") %>% select(-Year) %>% select(-starts_with("."))
Global3_unif <- Summary_Hier_NB_ANorm_unif$Global_Windows %>% mutate(Model = "Hier_NB_ANorm") %>% select(-Year) %>% select(-starts_with("."))
Global_Winds_unif <- rbind(Global3_unif, Global2_unif, Global1_unif) %>% mutate_if(is.numeric, round)
#Convert day to date for excel formatting
toDate <- function(x){as.Date(x-1, origin = as.Date("2023-01-01"))}
Global_Winds_Dates_unif <- Global_Winds_unif %>% mutate_if(is.numeric, toDate)
write.csv(Global_Winds_Dates_unif, file = "Outputs/2_DataOut/Global_Windows_Dates_uniform_prior.csv", row.names = F)

#Plot comparison of windows 
colnames(Global_Winds_unif) <- paste0(colnames(Global_Winds_unif),"unif")
Global_Winds_unif_joined <- cbind(Global_Winds,Global_Winds_unif)


ggplot(Global_Winds_unif_joined)+
  geom_point(aes(x=hD_50, y = hD_50unif, shape = Model, color = "Median"), alpha = 0.7, size = 5)+
  geom_point(aes(x=hD_025, y = hD_025unif, shape = Model, color = "2.5%"), alpha = 0.7, size = 5)+
  geom_point(aes(x=hD_975, y = hD_975unif, shape = Model, color = "97.5%"), alpha = 0.7, size = 5)+
  geom_abline(intercept = 0,slope = 1, color = "darkgrey")+
  theme_bw()+
  scale_y_continuous(breaks=c(213,244, 274, 305, 335),
                     labels=c("Aug 1", "Sept 1", "Oct 1", "Nov 1", "Dec 1"))+
  scale_x_continuous(breaks=c(213,244, 274, 305, 335),
                     labels=c("Aug 1", "Sept 1", "Oct 1", "Nov 1", "Dec 1"))+
  ylab("Date from model with uniform prior")+
  xlab("Date from model with normal prior")+
  ggtitle("Comparison of Global Windows for each model")+
  scale_colour_manual("Quantiles", values = c("2.5%" = "darkred", "Median"="darkorange", "97.5%"="darkblue"))
ggsave(file = "Outputs/3_Plots/uniform_prior_global_windows_comparison.png", width = 6, height = 5, units = "in")



# ggplot()+
#   geom_point(aes(x=Global_Winds$hD_50, y = Global_Winds_unif$hD_50, color = "Median"),alpha = 0.7,size = 5)+
#   geom_point(aes(x=Global_Winds$hD_025, y = Global_Winds_unif$hD_025, color = "2.5%"),alpha = 0.7, size =5)+
#   geom_point(aes(x=Global_Winds$hD_975, y = Global_Winds_unif$hD_975, color = "97.5%"),alpha = 0.7, size = 5)+
#   geom_abline(intercept = 0,slope = 1, color = "darkgrey")+
#   theme_bw()+
#   scale_y_continuous(breaks=c(213,244, 274, 305, 335),
#                      labels=c("Aug 1", "Sept 1", "Oct 1", "Nov 1", "Dec 1"))+
#   scale_x_continuous(breaks=c(213,244, 274, 305, 335),
#                      labels=c("Aug 1", "Sept 1", "Oct 1", "Nov 1", "Dec 1"))+
#   ylab("Date from model with uniform prior")+
#   xlab("Date from model with normal prior")+
#   ggtitle("Comparison of Global Windows for each model")+
#   scale_colour_manual("Quantiles", values = c("2.5%" = "darkred", "Median"="darkorange", "97.5%"="darkblue"))
#ggsave(file = "Outputs/3_Plots/uniform_prior_global_windows_comparison.png", width = 6, height = 5, units = "in")
#Maybe add which model
#Add 1:1 label




#### run the model with data cut off 
#Cut the data off at Aug 20th Nov 20th 
Catch_sub <- Catch[Catch$Day >= 233 & Catch$Day <= 325,]

#These are useful to have on hand
minDay <- min(Catch_sub$Day)
minYear <- min(Catch$Year)
maxYear <- max(Catch$Year)

#Data inputs - tried to make more flexible
DataList_sub <- list ("n_years"=length(unique(Catch_sub$Year)), 
                      "n_fisheries" = length(unique(Catch_sub$Fishery)),
                      "n_obs" = dim(Catch_sub)[[1]],
                      "years" = (Catch_sub$Year-minYear+1),
                      "days" = (Catch_sub$Day-minDay+1),
                      "fishery" = ifelse(Catch_sub$Fishery == "CN_test", 1, 2), 
                      "obs_catch"=Catch_sub$SH_Catch,
                      "annual_return"=annual_return$N,
                      # also include prior specification
                      "rt_m_mu" = 280-minDay, "rt_m_sig" = 40, # mean date of Oct 7
                      "halft_s" = 15  , "halft_df" = 2, 
                      "logit_q_mean" = -3, "logit_q_tau" = 0.8)


#Indep 
JagsFit_Indep_Poisson <- jags.parallel(DataList_sub, model.file =Non_Hier_Poisson_Q, 
                                       n.chains =3, n.iter = 20000, n.burnin = 5000, n.thin = 5,
                                       parameters.to.save = c("rt_m", "rt_sd", "q", "pred_abund"))
# Check convergence
All_Ests <-   data.frame(JagsFit_Indep_Poisson$BUGSoutput$summary)
All_Ests$Param <- row.names(All_Ests) 
All_Ests %>% filter(Rhat>1.01)

Summary_Indep_Poisson_sub <- Summarize_Data(JagsFit_Indep_Poisson, Catch_sub, Hier=F, ANorm=F,
                                            Timing_Params = "rt_m[Year], rt_sd[Year]")
saveRDS(Summary_Indep_Poisson_sub, "Outputs/1_ModelRuns/Summary_Indep_Poisson_sub.RDS")
rm(JagsFit_Indep_Poisson)

#Get windows
windows_indep_sub <- Summary_Indep_Poisson_sub$Yearly_Windows
#Recalc the year 
windows_indep_sub$Year <- windows_indep_sub$Year+minYear-0.8
#Plot comparison
ggplot()+
  geom_pointrange(data = windows_indep, aes(x = Year, y = hD_50, ymin = hD_025, ymax = hD_975, color = "Aug 1 - Dec 1"))+
  geom_pointrange(data = windows_indep_sub, aes(x = Year, y = D_50, ymin = D_025, ymax = D_975, color = "Aug 20 - Nov 20"))+
  ggtitle("Indep. Normal, Poisson")+
  ylab("Date")+
  scale_y_continuous(breaks=c(213,244, 274, 305, 335),
                     labels=c("Aug 1", "Sept 1", "Oct 1", "Nov 1", "Dec 1"))+
  theme_bw()+
  scale_colour_manual("", values = c("Aug 1 - Dec 1"="navyblue", "Aug 20 - Nov 20"="darkorange2"))

ggsave(file = "Outputs/3_Plots/Indep_Poisson_Aug20_Nov20.png", width = 8, height = 5, units = "in")



#NB
JagsFit_Hier_NB_HQ <- jags.parallel(DataList_sub, model.file = Hier_NB_HQ, 
                                    n.chains =3, n.iter=20000, n.burnin = 5000, n.thin = 5, 
                                    parameters.to.save = c("rt_m", "rt_sd", "r", "q", "pred_abund", "rt_m_m", "rt_sd_m"))

All_Ests_Hier_NB_HQ <-   data.frame(JagsFit_Hier_NB_HQ$BUGSoutput$summary)
All_Ests_Hier_NB_HQ$Param <- row.names(All_Ests_Hier_NB_HQ) 
All_Ests_Hier_NB_HQ %>% filter(Rhat>1.01)

Summary_Hier_NB_HQ_sub <- Summarize_Data(JagsFit_Hier_NB_HQ, Catch_sub, Timing_Params = "rt_m[Year], rt_sd[Year], rt_m_m, rt_sd_m", Hier=T, ANorm=F)
saveRDS(Summary_Hier_NB_HQ_sub, "Outputs/1_ModelRuns/Summary_Hier_NB_HQ_sub.RDS")
rm(JagsFit_Hier_NB_HQ)

windows_nb_sub <- Summary_Hier_NB_HQ_sub$Yearly_Windows
#Recalc the year 
windows_nb_sub$Year <- windows_nb_sub$Year+minYear-0.8
#Plot comparison
ggplot()+
  geom_pointrange(data = windows_nb, aes(x = Year, y = hD_50, ymin = hD_025, ymax = hD_975, color = "Aug 1 - Dec 1"))+
  geom_pointrange(data = windows_nb_sub, aes(x = Year, y = D_50, ymin = D_025, ymax = D_975, color = "Aug 20 - Nov 20"))+
  ggtitle("Hier. Normal, Neg. Binomial")+
  ylab("Date")+
  scale_y_continuous(breaks=c(213,244, 274, 305, 335),
                     labels=c("Aug 1", "Sept 1", "Oct 1", "Nov 1", "Dec 1"))+
  theme_bw()+
  scale_colour_manual("", values = c("Aug 1 - Dec 1"="navyblue", "Aug 20 - Nov 20"="darkorange2"))

ggsave(file = "Outputs/3_Plots/Hier_NB_HQ_Aug20_Nov20.png", width = 8, height = 5, units = "in")


#Asym 
JagsFit_Hier_NB_ANorm <- jags.parallel(DataList_sub, model.file = Hier_NB_ANorm, 
                                       n.chains =3, n.iter=30000, n.burnin = 5000, n.thin = 5, 
                                       parameters.to.save = c("rt_m", "rt_sd", "rt_m_m", "rt_sd_m", "r", "q", "pred_abund"))
All_Ests_Hier_NB_ANorm <-   data.frame(JagsFit_Hier_NB_ANorm$BUGSoutput$summary)
All_Ests_Hier_NB_ANorm$Param <- row.names(All_Ests_Hier_NB_ANorm) 
All_Ests_Hier_NB_ANorm %>% filter(Rhat>1.01) %>% filter(!str_detect(Param, "^pred_abund"))


Summary_Hier_NB_ANorm_sub <- Summarize_Data(JagsFit_Hier_NB_ANorm, Catch_sub, Hier=T, ANorm=T,
                                            Timing_Params = "rt_m[Year], rt_sd[Side,Year], rt_m_m, rt_sd_m[Side]")
saveRDS(Summary_Hier_NB_ANorm_sub, "Outputs/1_ModelRuns/Summary_Hier_NB_ANorm_sub.RDS")
rm(JagsFit_Hier_NB_ANorm)

windows_asym_sub <- Summary_Hier_NB_ANorm_sub$Yearly_Windows
#Recalc the year 
windows_asym_sub$Year <- windows_asym_sub$Year+minYear-0.8 #trick to get the offset in the plot
#Plot comparison
ggplot()+
  geom_pointrange(data = windows_asym, aes(x = Year, y = hD_50, ymin = hD_025, ymax = hD_975, color = "Aug 1 - Dec 1"))+
  geom_pointrange(data = windows_asym_sub, aes(x = Year, y = D_50, ymin = D_025, ymax = D_975, color = "Aug 20 - Nov 20"))+
  ggtitle("Hier. Asym. Normal, Neg. Binomial")+
  ylab("Date")+
  scale_y_continuous(breaks=c(213,244, 274, 305, 335),
                     labels=c("Aug 1", "Sept 1", "Oct 1", "Nov 1", "Dec 1"))+
  theme_bw()+
  scale_colour_manual("", values = c("Aug 1 - Dec 1"="navyblue", "Aug 20 - Nov 20"="darkorange2"))

ggsave(file = "Outputs/3_Plots/Hier_NB_ANorm_Aug20_Nov20.png", width = 8, height = 5, units = "in")


#Get global windows 
Global1_sub <- Summary_Indep_Poisson_sub$Global_Windows %>% mutate(Model = "Indep_Poisson") %>% select(-starts_with("."))
Global2_sub <- Summary_Hier_NB_HQ_sub$Global_Windows %>% mutate(Model = "Hier_NB_HQ") %>% select(-Year) %>% select(-starts_with("."))
Global3_sub <- Summary_Hier_NB_ANorm_sub$Global_Windows %>% mutate(Model = "Hier_NB_ANorm") %>% select(-Year) %>% select(-starts_with("."))
Global_Winds_sub <- rbind(Global3_sub, Global2_sub, Global1_sub) %>% mutate_if(is.numeric, round)
#Convert day to date for excel formatting
Global_Winds_Dates_sub <- Global_Winds_sub %>% mutate_if(is.numeric, toDate)
write.csv(Global_Winds_Dates_sub, file = "Outputs/2_DataOut/Global_Windows_Dates_Aug20_Nov20.csv", row.names = F)


#Plot comparison of windows 
colnames(Global_Winds_sub) <- paste0(colnames(Global_Winds_sub),"sub")
Global_Winds_sub_joined <- cbind(Global_Winds,Global_Winds_sub)


ggplot(Global_Winds_sub_joined)+
  geom_point(aes(x=hD_50, y = hD_50sub, shape = Model, color = "Median"), alpha = 0.7, size = 5)+
  geom_point(aes(x=hD_025, y = hD_025sub, shape = Model, color = "2.5%"), alpha = 0.7, size = 5)+
  geom_point(aes(x=hD_975, y = hD_975sub, shape = Model, color = "97.5%"), alpha = 0.7, size = 5)+
  geom_abline(intercept = 0, slope = 1, color = "darkgrey")+
  theme_bw()+
  scale_y_continuous(breaks=c(213,244, 274, 305, 335),
                     labels=c("Aug 1", "Sept 1", "Oct 1", "Nov 1", "Dec 1"))+
  scale_x_continuous(breaks=c(213,244, 274, 305, 335),
                     labels=c("Aug 1", "Sept 1", "Oct 1", "Nov 1", "Dec 1"))+
  ylab("Date from model with truncated dates")+
  xlab("Date from model with original dates")+
  ggtitle("Comparison of Global Windows for each model")+
  scale_colour_manual("Quantiles", values = c("2.5%" = "darkred", "Median"="darkorange", "97.5%"="darkblue"))
ggsave(file = "Outputs/3_Plots/global_windows_comparison_Aug20_Nov20.png", width = 6, height = 5, units = "in")




# #Plot all the years on a one-to-one line 
# ggplot()+
#   geom_point(aes(x=windows_asym$hD_50, y = windows_asym_sub$D_50), color = "darkorange")+
#   geom_point(aes(x=windows_asym$hD_025, y = windows_asym_sub$D_025), color = "darkred")+
#   geom_point(aes(x=windows_asym$hD_975, y = windows_asym_sub$D_975), color = "darkblue")+
#   geom_point(aes(x=Global_Winds$hD_50[Global_Winds$Model == "Hier_NB_ANorm"], y = Global_Winds_sub$hD_50[Global_Winds$Model == "Hier_NB_ANorm"]), color = "darkorange", shape = 13, size = 5)+
#   geom_point(aes(x=Global_Winds$hD_025[Global_Winds$Model == "Hier_NB_ANorm"], y = Global_Winds_sub$hD_025[Global_Winds$Model == "Hier_NB_ANorm"]), color = "darkred", shape = 13, size =5)+
#   geom_point(aes(x=Global_Winds$hD_975[Global_Winds$Model == "Hier_NB_ANorm"], y = Global_Winds_sub$hD_975[Global_Winds$Model == "Hier_NB_ANorm"]), color = "darkblue", shape = 13, size = 5)+
#   
#   geom_abline(intercept = 0,slope = 1, color = "darkgrey")+
#   theme_bw()+
#   scale_y_continuous(breaks=c(213,244, 274, 305, 335),
#                      labels=c("Aug 1", "Sept 1", "Oct 1", "Nov 1", "Dec 1"))+
#   scale_x_continuous(breaks=c(213,244, 274, 305, 335),
#                      labels=c("Aug 1", "Sept 1", "Oct 1", "Nov 1", "Dec 1"))+
#   ylab("Date")+
#   xlab("Date")+
#   ggtitle("Hier. Asym. Normal, Neg. Binomial")
# ggsave(file = "Outputs/3_Plots/Comparison_Aug20_Nov_20_Asym.png", width = 6, height = 4, units = "in")
# 
# #Plot the comparison of windows 
# ggplot()+
#   geom_point(aes(x=Global_Winds$hD_50, y = Global_Winds_sub$hD_50, color = "Median"),alpha = 0.7,size = 5)+
#   geom_point(aes(x=Global_Winds$hD_025, y = Global_Winds_sub$hD_025, color = "2.5%"),alpha = 0.7, size =5)+
#   geom_point(aes(x=Global_Winds$hD_975, y = Global_Winds_sub$hD_975, color = "97.5%"),alpha = 0.7, size = 5)+
#   geom_abline(intercept = 0,slope = 1, color = "darkgrey")+
#   theme_bw()+
#   scale_y_continuous(breaks=c(213,244, 274, 305, 335),
#                      labels=c("Aug 1", "Sept 1", "Oct 1", "Nov 1", "Dec 1"))+
#   scale_x_continuous(breaks=c(213,244, 274, 305, 335),
#                      labels=c("Aug 1", "Sept 1", "Oct 1", "Nov 1", "Dec 1"))+
#   ylab("Date from model with truncated dates")+
#   xlab("Date from model with original dates")+
#   ggtitle("Comparison of Global Windows for each model")+
#   scale_colour_manual("Quantiles", values = c("2.5%" = "darkred", "Median"="darkorange", "97.5%"="darkblue"))
# ggsave(file = "Outputs/3_Plots/global_windows_comparison_Aug20_Nov20.png", width = 6, height = 5, units = "in")






#Run the model with different return data
annual_return_alt <- read.csv("Data/annual_SH_return_to_Albion_CatchChanges.csv")

#Data inputs - tried to make more flexible
DataList_alt <- list ("n_years"=length(unique(Catch$Year)), 
                      "n_fisheries" = length(unique(Catch$Fishery)),
                      "n_obs" = dim(Catch)[[1]],
                      "years" = (Catch$Year-minYear+1),
                      "days" = (Catch$Day-minDay+1),
                      "fishery" = ifelse(Catch$Fishery == "CN_test", 1, 2), 
                      "obs_catch"=Catch$SH_Catch,
                      "annual_return"=annual_return_alt$N,
                      # also include prior specification
                      "rt_m_mu" = 280-minDay, "rt_m_sig" = 40, # mean date of Oct 7
                      "halft_s" = 15  , "halft_df" = 2, 
                      "logit_q_mean" = -3, "logit_q_tau" = 0.8)



#Indep 
JagsFit_Indep_Poisson <- jags.parallel(DataList_alt, model.file =Non_Hier_Poisson_Q, 
                                       n.chains =3, n.iter = 30000, n.burnin = 5000, n.thin = 5,
                                       parameters.to.save = c("rt_m", "rt_sd", "q", "pred_abund"))
# Check convergence
All_Ests <-   data.frame(JagsFit_Indep_Poisson$BUGSoutput$summary)
All_Ests$Param <- row.names(All_Ests) 
All_Ests %>% filter(Rhat>1.01)

Summary_Indep_Poisson_alt <- Summarize_Data(JagsFit_Indep_Poisson, Catch, Hier=F, ANorm=F,
                                            Timing_Params = "rt_m[Year], rt_sd[Year]")
saveRDS(Summary_Indep_Poisson_alt, "Outputs/1_ModelRuns/Summary_Indep_Poisson_alt.RDS")
rm(JagsFit_Indep_Poisson)

#Get windows
windows_indep_alt <- Summary_Indep_Poisson_alt$Yearly_Windows
#Recalc the year 
windows_indep_alt$Year <- windows_indep_alt$Year+minYear-0.8
#Plot comparison
ggplot()+
  geom_pointrange(data = windows_indep, aes(x = Year, y = hD_50, ymin = hD_025, ymax = hD_975, color = "Original Catch"))+
  geom_pointrange(data = windows_indep_alt, aes(x = Year, y = D_50, ymin = D_025, ymax = D_975, color = "Adjusted Catch"))+
  ggtitle("Indep. Normal, Poisson w/ adjusted catch")+
  ylab("Date")+
  scale_y_continuous(breaks=c(213,244, 274, 305, 335),
                     labels=c("Aug 1", "Sept 1", "Oct 1", "Nov 1", "Dec 1"))+
  theme_bw()+
  scale_colour_manual("", values = c("Original Catch"="navyblue", "Adjusted Catch"="darkorange2"))

ggsave(file = "Outputs/3_Plots/Indep_Poisson_adjusted_catch.png", width = 8, height = 5, units = "in")


#NB
JagsFit_Hier_NB_HQ <- jags.parallel(DataList_alt, model.file = Hier_NB_HQ, 
                                    n.chains =3, n.iter=30000, n.burnin = 5000, n.thin = 5, 
                                    parameters.to.save = c("rt_m", "rt_sd", "r", "q", "pred_abund", "rt_m_m", "rt_sd_m"))

All_Ests_Hier_NB_HQ <-   data.frame(JagsFit_Hier_NB_HQ$BUGSoutput$summary)
All_Ests_Hier_NB_HQ$Param <- row.names(All_Ests_Hier_NB_HQ) 
All_Ests_Hier_NB_HQ %>% filter(Rhat>1.01)

Summary_Hier_NB_HQ_alt <- Summarize_Data(JagsFit_Hier_NB_HQ, Catch, Timing_Params = "rt_m[Year], rt_sd[Year], rt_m_m, rt_sd_m", Hier=T, ANorm=F)
saveRDS(Summary_Hier_NB_HQ_alt, "Outputs/1_ModelRuns/Summary_Hier_NB_HQ_alt.RDS")
rm(JagsFit_Hier_NB_HQ)

windows_nb_alt <- Summary_Hier_NB_HQ_alt$Yearly_Windows
#Recalc the year 
windows_nb_alt$Year <- windows_nb_alt$Year+minYear-0.8

#Plot comparison
ggplot()+
  geom_pointrange(data = windows_nb, aes(x = Year, y = hD_50, ymin = hD_025, ymax = hD_975, color = "Original Catch"))+
  geom_pointrange(data = windows_nb_alt, aes(x = Year, y = D_50, ymin = D_025, ymax = D_975, color = "Adjusted Catch"))+
  ggtitle("Hier. Normal, Neg. Binomial w/ adjusted catch")+
  ylab("Date")+
  scale_y_continuous(breaks=c(213,244, 274, 305, 335),
                     labels=c("Aug 1", "Sept 1", "Oct 1", "Nov 1", "Dec 1"))+
  theme_bw()+
  scale_colour_manual("", values = c("Original Catch"="navyblue", "Adjusted Catch"="darkorange2"))
ggsave(file = "Outputs/3_Plots/Hier_NB_HQ_adjusted_catch.png", width = 8, height = 5, units = "in")

#Asym 
JagsFit_Hier_NB_ANorm <- jags.parallel(DataList_alt, model.file = Hier_NB_ANorm, 
                                       n.chains =3, n.iter=30000, n.burnin = 5000, n.thin = 5, 
                                       parameters.to.save = c("rt_m", "rt_sd", "rt_m_m", "rt_sd_m", "r", "q", "pred_abund"))
All_Ests_Hier_NB_ANorm <-   data.frame(JagsFit_Hier_NB_ANorm$BUGSoutput$summary)
All_Ests_Hier_NB_ANorm$Param <- row.names(All_Ests_Hier_NB_ANorm) 
All_Ests_Hier_NB_ANorm %>% filter(Rhat>1.01) %>% filter(!str_detect(Param, "^pred_abund"))


Summary_Hier_NB_ANorm_alt <- Summarize_Data(JagsFit_Hier_NB_ANorm, Catch, Hier=T, ANorm=T,
                                            Timing_Params = "rt_m[Year], rt_sd[Side,Year], rt_m_m, rt_sd_m[Side]")
saveRDS(Summary_Hier_NB_ANorm_alt, "Outputs/1_ModelRuns/Summary_Hier_NB_ANorm_alt.RDS")
rm(JagsFit_Hier_NB_ANorm)

windows_asym_alt <- Summary_Hier_NB_ANorm_alt$Yearly_Windows
#Recalc the year 
windows_asym_alt$Year <- windows_asym_alt$Year+minYear-0.8 #trick to get the offset in the plot
#Plot comparison
ggplot()+
  geom_pointrange(data = windows_asym, aes(x = Year, y = hD_50, ymin = hD_025, ymax = hD_975,color = "Original Catch"))+
  geom_pointrange(data = windows_asym_alt, aes(x = Year, y = D_50, ymin = D_025, ymax = D_975, color = "Adjusted Catch"))+
  ggtitle("Hier. Asym. Normal, Neg. Binomial w/ adjusted catch")+
  ylab("Date")+
  scale_y_continuous(breaks=c(213,244, 274, 305, 335),
                     labels=c("Aug 1", "Sept 1", "Oct 1", "Nov 1", "Dec 1"))+
  theme_bw()+
  scale_colour_manual("", values = c("Original Catch"="navyblue", "Adjusted Catch"="darkorange2"))
ggsave(file = "Outputs/3_Plots/Hier_NB_ANorm_adusted_catch.png", width = 8, height = 5, units = "in")


#Get global windows 
Global1_alt <- Summary_Indep_Poisson_alt$Global_Windows %>% mutate(Model = "Indep_Poisson") %>% select(-starts_with("."))
Global2_alt <- Summary_Hier_NB_HQ_alt$Global_Windows %>% mutate(Model = "Hier_NB_HQ") %>% select(-Year) %>% select(-starts_with("."))
Global3_alt <- Summary_Hier_NB_ANorm_alt$Global_Windows %>% mutate(Model = "Hier_NB_ANorm") %>% select(-Year) %>% select(-starts_with("."))
Global_Winds_alt <- rbind(Global3_alt, Global2_alt, Global1_alt) %>% mutate_if(is.numeric, round)
#Convert day to date for excel formatting
Global_Winds_Dates_alt <- Global_Winds_alt %>% mutate_if(is.numeric, toDate)
write.csv(Global_Winds_Dates_alt, file = "Outputs/2_DataOut/Global_Windows_Dates_adjusted_catch.csv", row.names = F)


#Plot comparison of windows 
colnames(Global_Winds_alt) <- paste0(colnames(Global_Winds_alt),"sub_alt")
Global_Winds_alt_joined <- cbind(Global_Winds,Global_Winds_alt)


ggplot(Global_Winds_alt_joined)+
  geom_point(aes(x=hD_50, y = hD_50sub_alt, shape = Model, color = "Median"), alpha = 0.7, size = 5)+
  geom_point(aes(x=hD_025, y = hD_025sub_alt, shape = Model, color = "2.5%"), alpha = 0.7, size = 5)+
  geom_point(aes(x=hD_975, y = hD_975sub_alt, shape = Model, color = "97.5%"), alpha = 0.7, size = 5)+
  geom_abline(intercept = 0, slope = 1, color = "darkgrey")+
  theme_bw()+
  scale_y_continuous(breaks=c(213,244, 274, 305, 335),
                     labels=c("Aug 1", "Sept 1", "Oct 1", "Nov 1", "Dec 1"))+
  scale_x_continuous(breaks=c(213,244, 274, 305, 335),
                     labels=c("Aug 1", "Sept 1", "Oct 1", "Nov 1", "Dec 1"))+
  ylab("Date from model with adjusted catch")+
  xlab("Date from model with original catch")+
  ggtitle("Comparison of Global Windows for each model")+
  scale_shape_manual(values = c(4, 16, 13)) +
  scale_colour_manual("Quantiles", values = c("2.5%" = "darkred", "Median"="darkorange", "97.5%"="darkblue"))
ggsave(file = "Outputs/3_Plots/global_windows_comparison_adjusted_catch.png", width = 6, height = 5, units = "in")

