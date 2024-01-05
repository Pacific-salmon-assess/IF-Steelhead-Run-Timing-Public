######################################################
### Functions for IF Steelhead Run-timing Models   ###
###  Brooke Davis, Anna Potapova                   ###
###  Started March 2023                            ###
######################################################


#-----------------------------------------------------------------------------
# Load catch data and put in usable form, including extended dates and other species catch
Load_Catch_Data <- function() {
  albion_annual_extended <- read.csv("Data/Albion SH Catch 1983-2022_long_Aug 1 to Nov 30.csv")
  #Note: for now the data for non-SH in 1983 is missing so it will be NA
  albion_annual_other <- read.csv("Data/Albion non-SH Catch 1984-2022_long.csv") %>%
    filter(Simple.Net.Config != "Albion Chum Non-Assessment Net", Simple.Net.Config!= "Albion Chin VMN")
  
  #Aggregate data for by Day and reshape to match format from SH_functions script
  # Year | Day | Fishery | Catch
  #Add column for chinook vs coho
  albion_annual_extended <- albion_annual_extended %>%
    filter(Simple.Net.Config != "Albion Chum Non-Assessment Net", Simple.Net.Config!= "Albion Chin VMN")
  albion_annual_extended$Fishery <- NA
  albion_annual_extended$Fishery[grep("Chin", albion_annual_extended$Simple.Net.Config)] <- "CN_test" #Chinook test fishery
  albion_annual_extended$Fishery[grep("Chum", albion_annual_extended$Simple.Net.Config)] <- "CM_test" #Chum test fishery
  
  #Aggregate sum of catch per day, per year, per fishery
  Catch <- aggregate(CATCH_QTY ~ Year + DOY + Fishery, data = albion_annual_extended, sum)
  
  #Add column for chinook vs coho to the other species data 
  albion_annual_other$Fishery <- NA
  albion_annual_other$Fishery[grep("Chin", albion_annual_other$Simple.Net.Config)] <- "CN_test" #Chinook test fishery
  albion_annual_other$Fishery[grep("Chum", albion_annual_other$Simple.Net.Config)] <- "CM_test" #Chum test fishery
  #Aggregate the other species
  albion_annual_other <- aggregate(CATCH_QTY ~ Year + DOY + Fishery, data = albion_annual_other, sum)
  
  #Merge the two 
  Catch <- merge(Catch, albion_annual_other, by = c("Year", "DOY", "Fishery"), all = T)
  Catch <- rename(Catch, Day = DOY, SH_Catch = CATCH_QTY.x, OTHER_Catch = CATCH_QTY.y)
  Catch <- Catch %>% filter(Day <= 335) 
  Catch
}

#-----------------------------------------------------------------------------
#  JAGS Models
#-----------------------------------------------------------------------------
#JAGS model for non-hierarchical Poisson -- should be similar to Bison 2021 model
Non_Hier_Poisson_Q <- function(){
  
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
    rt_m[y] ~ dnorm(rt_m_mu,1/rt_m_sig^2)
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

#---------------------------------------------------------------------------

# Hier timing, NB error model

Hier_NB_HQ <- function(){
  
  #Likelihood
  # Global dists
  # Mean and sd around mean date global dist
  rt_m_m ~ dnorm(rt_m_mu, 1/rt_m_sig^2)
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

#----------------------------------------------------------
#Asymmetric normal, Hier timing, NB Obs model
Hier_NB_ANorm<- function(){
  
  # Global dists
  # Mean and sd around mean date global dist
  rt_m_m ~ dnorm(rt_m_mu, 1/rt_m_sig^2)
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
#----------------------------------------------------------
#Asymmetric normal, Hier timing, NB Obs model - test for changes in mean run-timing
Hier_NB_ANorm_TestMu <- function(){
  
  # Global dists
  # Intercept for linear rt_m
  rt_m_Int ~ dunif(-100, 100)
  rt_m_Slope ~ dnorm(0, 0.01)
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
    # mean now linear mod
    u[y] ~ dnorm(0, rt_m_tau)  # yearly random effect
    rt_m[y] <- rt_m_Int + rt_m_Slope*y + u[y]
    # sd from "global" shared dist
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


# data to make jags outputs usable

Summarize_Data <- function(JagsFit, Catch, Timing_Params, Dharma_draws = 1000, ANorm, Hier){
  
  # turn into mcmc list object
  JagsFit_MCMC<- as.mcmc(JagsFit)
  
  #Get the MinYear and MinDay from the Catch data 
  MinYear <- min(Catch$Year)
  MinDay <- min(Catch$Day)
  
  
  # Gather draws for timing params
    Date_Dist_Draws <- eval(parse(text = paste("tidybayes::spread_draws(JagsFit_MCMC,", Timing_Params, ")")))
    if("Year" %in% names(Date_Dist_Draws)){
    Date_Dist_Ints <- point_interval(Date_Dist_Draws) %>%
      # Put years and dates back to normal
      mutate(Year = Year + MinYear - 1) 
    } else {
      Date_Dist_Ints <- point_interval(Date_Dist_Draws)
    }
    
  # Get global windows, and uncertainty around them---------------------------------------------------
  # need to cut down number of years, since global params will be repeated for every year
    percs <- c(0.025, 0.05, 0.1, 0.5, 0.9, 0.95, 0.975)
    perc_names_h <- c("hD_025", "hD_05", "hD_10", "hD_50", "hD_90", "hD_95", "hD_975")
    perc_names <- c("D_025", "D_05", "D_10", "D_50", "D_90", "D_95", "D_975")
  # If Anorm spread so that have sd+ and sd- all in same row
    if(ANorm == T & Hier==T){
      # Get global windows
      Hier_Param_Draws <- Date_Dist_Draws %>%
        filter(Year == 1) %>%
        select(!c(rt_sd, rt_m)) %>%
        pivot_wider(values_from = rt_sd_m, names_from = Side, names_prefix = "rt_sd_m") 
      
      N_Draws <- dim(Hier_Param_Draws)[1]
      HWin <- matrix(nrow=N_Draws, ncol = length(percs))
      for(i in 1:N_Draws){
          HWin[i, ] <- qdnorm(p = percs, mean = Hier_Param_Draws$rt_m_m[i], sigma= as.numeric(Hier_Param_Draws[i,c("rt_sd_m1", "rt_sd_m2")])) + MinDay-1
      }
      HWin_DF <- as.data.frame(HWin)
      names(HWin_DF) <- perc_names_h
      Window_Draws <- cbind(Hier_Param_Draws, HWin_DF) %>% select(!c("rt_m_m", "rt_sd_m1", "rt_sd_m2", "Year"))
      Global_Windows <- point_interval(Window_Draws)
      
      # Now get yearly windows
      Yearly_Param_Draws <- Date_Dist_Draws %>%
        select(!c(rt_sd_m, rt_m_m)) %>%
        pivot_wider(values_from = rt_sd, names_from = Side, names_prefix = "rt_sd")
      N_Draws_Yearly <- dim(Yearly_Param_Draws)[1]
      Win <- matrix(nrow=N_Draws_Yearly, ncol = length(percs))
        for(i in 1:N_Draws_Yearly){
            Win[i, ] <- qdnorm(p = percs, mean = Yearly_Param_Draws$rt_m[i], sigma= as.numeric(Yearly_Param_Draws[i,c("rt_sd1", "rt_sd2")])) + MinDay-1
        }
      Win_DF <- as.data.frame(Win)
      names(Win_DF) <- perc_names
      Window_Draws <- cbind(Yearly_Param_Draws, Win_DF) %>% select(!c("rt_m", "rt_sd1", "rt_sd2"))
      Yearly_Windows <- point_interval(Window_Draws)
      
    } else if(ANorm==F) {
      # if normal, getting yearly windows the same
      N_Draws <- dim(Date_Dist_Draws)[1]
      Win <- matrix(nrow=N_Draws, ncol = length(percs))
      for(i in 1:N_Draws){
        for(j in 1:length(percs)){
          Win[i, ] <- qnorm(p = percs, mean = Date_Dist_Draws$rt_m[i],sd= Date_Dist_Draws$rt_sd[i]) + MinDay-1
        }
      }
      Win_DF <- as.data.frame(Win)
      names(Win_DF) <- perc_names
      Window_Draws <- cbind(Date_Dist_Draws, Win_DF) %>% select(-starts_with("rt")) 
      Yearly_Windows <- point_interval(Window_Draws) 
      
      # Get global Windows
        if(Hier==T){
          # do this with draws of hier params for Hier model
          Hier_Param_Draws <- Date_Dist_Draws %>%
            filter(Year == 1) %>%
            select(!c(rt_sd, rt_m)) 
          N_Draws <- dim(Hier_Param_Draws)[1]
          HWin <- matrix(nrow=N_Draws, ncol = length(percs))
          for(i in 1:N_Draws){
            HWin[i, ] <- qnorm(p = percs, mean = Hier_Param_Draws$rt_m_m[i], sd= Hier_Param_Draws$rt_sd_m[i]) + MinDay-1
          }
          HWin_DF <- as.data.frame(HWin)
          names(HWin_DF) <- perc_names_h
          Window_Draws <- cbind(Hier_Param_Draws, HWin_DF) %>% select(-starts_with("rt"))
          Global_Windows <- point_interval(Window_Draws)
        } else if(Hier==F){
          # for each draw take median across years
          Draw_Medians <- Window_Draws %>% rename_at(vars(perc_names), ~perc_names_h) %>%  group_by(.draw) %>% summarise_at(perc_names_h, median)
          # Could also do means across years? -- gives very wide windows
          #Draw_Means <- Window_Draws %>% rename_at(vars(perc_names), ~perc_names_h) %>%  group_by(.draw) %>% summarise_at(perc_names_h, mean)
          # now take median of those
          Global_Windows <- point_interval(Draw_Medians)
        }
      # for each draw get mean across years
      # do point_interval on those (should be 9000 long)
    }
     
  # Gather everything else we need - q, pred_Abund
    if(grepl("Poisson", deparse(substitute(JagsFit)))){
      print("Poisson")
      Pred_Draws <- tidybayes::spread_draws(JagsFit_MCMC, q[Fishery,Year], pred_abund[Day, Year], ndraws = Dharma_draws) %>%
        # Now get to get Pred_Catch 
        mutate(Pred_Catch = pred_abund*q) %>%
        # Nwo simulate to Post_Pred_Catch
        mutate(Year = Year + MinYear - 1) %>%
        mutate(Day = Day + MinDay -1)
      
      Pred_Draws$Post_Pred_Catch <- rpois(dim(Pred_Draws)[1], Pred_Draws$Pred_Catch)
      Pred_Draws$Fishery = ifelse(Pred_Draws$Fishery==1, "CN_test", "CM_test")
      
      Pred_Ints <- point_interval( Pred_Draws) %>%
        left_join(Catch[, c("Year", "Day", "Fishery", "SH_Catch")])
      
    } else if(grepl("NB", deparse(substitute(JagsFit)))){
  print("NB")
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
     
    } 
    
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
  
  q_draws <- ungroup(Pred_Draws) %>% select(q, .draw, Fishery, Year) %>% distinct(.draw, q, Fishery, Year)
    
  out <- list()
  out$Pred_Catch_Ints <- Pred_Ints
  out$Date_Dists_Ints <- Date_Dist_Ints
  out$q_draws <- q_draws
  out$sims <- sims
  out$Global_Windows <- Global_Windows
  out$Yearly_Windows <- Yearly_Windows
  out
}

#----------------------------------------------------------------------------
# use this to plot priors on catchability

dilogit <- function(p, mu, tau){
  
  x <- log(p/(1-p))
  
  jac <- (1/(1-p) + p/(1-p)^2)/(p/(1-p))
  
  ans <- dnorm(x, mean=mu, sd=1/sqrt(tau), 0) * jac
  
  return(ans)
  
}

#------------------------------------------------------------------------------
# Create matching plots for each model form
Create_Plots <- function(Summary_Data){
  
  
  #Grab the model name from the object - coerce the input to string, split at the underscore 
  #N=2 means it will only split into two pieces, use indexing to grab the second piece with the model name 
  Mod_Name <- str_split(deparse(substitute(Summary_Data)), "_", n = 2)[[1]][2]
  
  Model_Preds <- Summary_Data$Pred_Catch_Ints
  Obs_Abund <- Model_Preds %>% filter(is.na(SH_Catch)==F)
  nrows <- 3; ncols <- 3
  pdf(paste("Outputs/3_Plots/Pred_Abund_", Mod_Name, ".pdf", sep=""),
      width = 8, height = 11)
  for (i in 1:ceiling(length(unique(Model_Preds$Year))/(nrows*ncols))) {
    print(
      ggplot() +
        geom_pointrange(data = Obs_Abund, aes(x=Day, y=SH_Catch/q, ymax = SH_Catch/q.lower, ymin=SH_Catch/q.upper, color=Fishery), alpha=0.5) +
        geom_line(data = Model_Preds, aes(x=Day, y=pred_abund)) +
        geom_ribbon(data = Model_Preds, aes(x=Day, ymin = pred_abund.lower, ymax = pred_abund.upper), alpha = .15) +
        theme_light() +
        facet_wrap_paginate(~ Year, nrow = nrows, ncol=ncols, page = i) 
    )}
  dev.off()
  
  # Now plot posterior predictives for each fishery, on top of obs, to see if captures data well

  nrows <- 3; ncols <- 3
  pdf(paste("Outputs/3_Plots/Pred_Catch_", Mod_Name, ".pdf", sep=""),
      width = 8, height = 11)
  for (i in 1:ceiling(length(unique(Model_Preds))/(nrows*ncols))) {
    print(
      ggplot() +
        geom_point(data = Obs_Abund, aes(x=Day, y=SH_Catch, color=Fishery), alpha=0.5) +
        geom_line(data = Model_Preds, aes(x=Day, y=Post_Pred_Catch, color=Fishery)) +
        geom_ribbon(data = Model_Preds, aes(x=Day, ymin = Post_Pred_Catch.lower, ymax = Post_Pred_Catch.upper), alpha = .15) +
        theme_light() +
        facet_wrap_paginate(~ Year, nrow = nrows, ncol=ncols, page = i) 
    )}
  dev.off()
  
  #Plot mean dates 
    
    pdf(paste("Outputs/3_Plots/Mean_Dates_", Mod_Name, ".pdf", sep=""),  
        width = 8, height = 11)
    print(
      ggplot() +
        geom_pointrange(data = Summary_Data$Date_Dists_Ints, aes(x=Year, y=rt_m, ymin = rt_m.lower, ymax=rt_m.upper)) +
        theme_light()
    )
    dev.off()
  
  
  #Plot the DHARMa residuals 
  pdf(paste("Outputs/3_Plots/Dharma_", Mod_Name,".pdf", sep =""), width = 11, height = 8)
    plot(Summary_Data$sims)
  dev.off()
  
  # Add catchability plot back in, use Paul's fuction for logit-transformed normal
  
  # Plot q priors and draws
       q_ests <- Model_Preds %>% group_by(Fishery, Year) %>% summarise(q=mean(q), q.lower=mean(q.lower), q.upper=mean(q.upper))
       x = seq(0.0001, 0.1, by = 0.0001)
       prior_df <- data.frame(x = x,
                      prior_dens = dilogit(x, mu = -3, tau = 0.8) )
       #Plot the catchability draws over the prior
       pdf(paste("Outputs/3_Plots/Catchability_Prior_", Mod_Name,".pdf", sep =""), width = 8, height = 11)
       for (i in 1:ceiling(length(unique(Model_Preds$Year))/(nrows*ncols))) {
         print(
         ggplot(Summary_Data$q_draws)+
           geom_histogram(aes(x = q, y = after_stat(density), group = Fishery, fill = Fishery),colour = 1, alpha = 0.5, bins = 50)+
           geom_line(data = prior_df, aes(x, prior_dens)) +
           theme_bw()+
           geom_point(data = q_ests, aes(x=q, y = 0), size = 2.5)+
           geom_segment(data = q_ests, aes(x= q.lower, xend = q.upper, y = 0, yend = 0), linewidth = 1.25)+
           xlab("Catchability (q)") +
           facet_wrap_paginate(~ Year, nrow = nrows, ncol=ncols, page = i) 
         )}
       dev.off()
       
       
} # end plotting function

# Function to get asymmetric normal quantiles
qdnorm <- function(p,  mean, sigma ){
    z <- p
    r <- sigma[1]/(sigma[1] + sigma[2])
    z[p <= r] <- mean + sigma[1] * qnorm(0.5 * p[p <= r] * sum(sigma)/sigma[1], mean = 0, sd = 1)
    z[p > r] <- mean + sigma[2] * qnorm(0.5 * (sum(sigma) * (1 + p[p > r]) - 2 * sigma[1])/sigma[2], mean = 0, sd = 1)
    z
}
# and density
ddnorm <- function(x, mean, sigma){
    logC <- log(2) - log(sum(sigma)) -  0.5*log(2*pi)
    z <- (x - mean)
    z[z <=0 ] <- z[z <=0 ]/sigma[1]
    z[z > 0 ] <- z[z > 0 ]/sigma[2]
    ans <- logC - 0.5*z^2
    exp(ans)
}

# Inverse logit
InvLogit <- function(x){
   exp(x)/(1+exp(x))
}


