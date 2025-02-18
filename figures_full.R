###########################################
# Figures for urban rabies simulation     # 
# Full model: seroprevalence, immigration # 
# Model written in Julia                  #
# E.M. Beasley                            #
# Fall 2023                               #
###########################################

# Load packages ---------------------
library(pbmcapply) #Progress bar w/ETA
library(tidyverse)
library(viridis)
library(agricolae)
library(R0)
library(patchwork)

options(dplyr.summarise.inform = FALSE)
dir <- "./outputs" 

# Function: Calculate proportion of outbreaks eliminated --------
prop_elim <- function(){
  # Get names of files
  filenames <- list.files(path = dir, pattern = "*.csv")
  
  # Progress bar
  pb = progressBar(min = 1, max = length(filenames), initial = 1,
                   style = "ETA") 
  
  for(i in 1:length(filenames)){
    # Read 'em in
    testfile <- read.csv(paste(getwd(), 
                               str_replace(dir, ".", ""), "/",
                               filenames[i],sep = ""))
    
    # Calculate time in weeks & get time to first elimination
    prop_elim <- testfile %>%
      filter(elim == "True") %>%
      group_by(sero, type, disease, rate) %>%
      # rate=disease rate of immigrants
      filter(year >= 2) %>%
      select(rep,sero,type,disease,rate) %>%
      distinct() %>%
      summarise(prop = n()/50)
    
    if(i==1){
      prop_elim_frame <- prop_elim
    } else{
      prop_elim_frame <- rbind(prop_elim_frame, prop_elim)
    }
    
    setTxtProgressBar(pb,i)
  }
  return(prop_elim_frame)
}

prop_elim <- prop_elim()

# Figs: Proportion of outbreaks eliminated ---------
# All sims
ggplot(data=prop_elim, aes(x=factor(sero), y = prop,
                           color = factor(rate)))+
  geom_jitter()+
  #geom_boxplot(fill="lightgray")+
  geom_hline(yintercept=0.95, linetype="dashed") +
  labs(x = "Adult Vaccination Rate", 
       y = "Proportion Reaching Elimination")+
  theme_bw(base_size=14)+
  theme(panel.grid=element_blank())

# ggsave("./full_Figs/prop_elim_box.jpeg", width = 7, height = 5,
#        units = "in")

# Function for calculating time to first elimination --------
first_elim <- function(){
  # Get names of files
  filenames <- list.files(path = dir, pattern = "*.csv")

  # Progress bar
  pb = progressBar(min = 1, max = length(filenames), initial = 1,
                   style = "ETA") 
  
  for(i in 1:length(filenames)){
    # Read 'em in
    testfile <- read.csv(paste(getwd(), 
                               str_replace(dir, ".", ""), "/",
                               filenames[i],sep = ""))

    # Calculate time in weeks & get time to first elimination
    time_to_elim <- testfile %>%
      filter(elim == "True", year >= 2) %>%
      group_by(rep, sero, rate, type, disease) %>%
      mutate(nweek = ((year-1)*52) + week) %>%
      filter(nweek == min(nweek))
    
    if(i==1){
      first_elim_frame <- time_to_elim
    } else{
      first_elim_frame <- rbind(first_elim_frame, time_to_elim)
    }
    
    setTxtProgressBar(pb,i)
  }
  return(first_elim_frame)
}

first_elim_full <- first_elim()

# Time to elimination: figs ------------------
elim_sansbar <- first_elim_full %>%
  mutate(sero = factor(sero)) %>%
  filter(rate %in% c(1,5))
  
ggplot(data = elim_sansbar, aes(x=factor(sero),y=nweek,
                                fill = factor(rate)))+
  geom_boxplot()+
  # geom_jitter()+
  scale_fill_manual(values = c('lightgray', 'limegreen'))+
  labs(x = "Adult Vaccination Rate", 
       y = "Time to Elimination (Weeks)")+
  theme_bw()+
  theme(panel.grid=element_blank())

# ggsave(filename = "./full_Figs/time_elim_full.jpeg", width=6,
#        height=4, dpi=600, units="in")

imm_rate_elim <- first_elim_full %>%
  group_by(sero, rate) %>%
  summarise(center = median(nweek))

# dev.new(width = 80, height = 60, unit = "mm", res=600,
#         noRStudioGD=TRUE)

ggplot(data = imm_rate_elim, aes(x = factor(sero), 
                                 y = factor(rate),fill=center))+
  geom_tile()+
  scale_fill_viridis(name = "Weeks to Elimination",
                       option="B")+
  labs(x = "Vaccination Rate", y = "Immigration Rate")+
  theme_bw(base_size=16)+
  theme(panel.grid = element_blank())
# Immigration matters most when immunity is low

# ggsave(filename = "./full_Figs/time_elim_rate.jpeg", width=6,
#        height=4, dpi=600, units="in")

# interaction: proportion diseased immigrants
rate_interac <- first_elim_full %>%
  mutate(disease = factor(disease, levels = c("0","0.015", 
                                              "0.03", "0.045",
                                              "0.06")),
                          sero=factor(sero)) %>%
  group_by(sero, disease) %>%
  summarise(center = median(nweek))

ggplot(data = rate_interac, aes(x = sero, y = disease, 
                                fill = center))+
  geom_tile()+
  scale_fill_viridis(name = "Weeks to Elimination", option="B")+
  labs(x ="Vaccination Rate", y="Immigrant Infection Rate")+
  theme_bw(base_size=12)+
  theme(panel.grid=element_blank())
# Again: immigration matters when seroprevalence is low

# interaction: immigration type
type_interac <- first_elim_full %>%
  group_by(sero, type) %>%
  summarise(center = median(nweek))

ggplot(data = type_interac, aes(x = factor(sero), y = type, 
                                fill=center))+
  geom_tile()+
  scale_fill_viridis(name="Weeks to Elimination", option="B")+
  labs(x = "Vaccination Rate", y="Immigration Type")+
  theme_bw(base_size=12)+
  theme(panel.grid=element_blank())
# Immigration type doesn't seem to matter

# any interaction between immigrant prevalence & type?
im_interac <- first_elim_full %>%
  mutate(disease=factor(disease)) %>%
  group_by(sero, disease, type) %>%
  summarise(center = median(nweek))

ggplot(data=im_interac[im_interac$sero==0.2,], 
       aes(x = disease, y = type, fill=center))+
  geom_tile()+
  scale_fill_viridis_c()
# not really

# Disease/immigration interactions
im_interac <- first_elim_full %>%
  mutate(disease=factor(disease), rate = factor(rate)) %>%
  group_by(sero, disease, rate) %>%
  summarise(center = median(nweek))

interac_list <- list()

for(i in 1:length(c(0, 0.2, 0.4, 0.6))){
  samples <- c(0, 0.2, 0.4, 0.6)
  
  interac_list[[i]] <- ggplot(data=im_interac[im_interac$sero==samples[i],], 
        aes(x = disease, y = rate, fill=center))+
    geom_tile()+
    scale_fill_viridis_c(name = "Weeks to Elimination",
                         limits = c(100, 215))+
    labs(x = "Immigrant Disease Rate", y = "Immigration Rate",
         title = paste("Vaccination Rate = ", samples[i],
                       sep = ""))
}

((interac_list[[1]] + interac_list[[2]])/
  (interac_list[[3]] + interac_list[[4]])) +
  plot_layout(guides = 'collect')

# Prob of elimination at time t --------------------
first_elim_prob <- function(){
  # Get names of files
  filenames <- list.files(path = dir, pattern = "*.csv")
  
  # Progress bar
  pb = progressBar(min = 1, max = length(filenames), initial = 1,
                   style = "ETA") 
  
  for(i in 1:length(filenames)){
    # Read 'em in
    testfile <- read.csv(paste(getwd(), 
                               str_replace(dir, ".", ""), "/",
                               filenames[i],sep = "")) %>%
      filter(year >= 2) %>%
      mutate(nweek = ((year-1)*52) + week)
    
    recols <- testfile %>%
      group_by(rep, sero, rate, type, disease) %>%
      filter(elim == "False" & lag(elim) == "True")
    
    elims <- testfile %>%
      group_by(rep, sero, rate, type, disease) %>%
      filter(elim == "True" & lag(elim) == "False")
    
    recol.time <- rbind(elims, recols) %>%
      arrange(rep, nweek) %>%
      mutate(time = case_when(elim == "True" ~ nweek-lag(nweek),
                              TRUE ~ NA)) %>%
      filter(is.na(time) == F) %>%
      filter(time < 9)
    
    testfile <- suppressMessages(full_join(testfile, recol.time))
    
    testfile <- testfile %>%
      mutate(outbreak = case_when(elim=="True"~0,
                                  TRUE~1)) 
    
    indices <- which(is.na(testfile$time)==F)
    
    if(length(indices) > 0){
      for(j in 1:length(indices)){
        testfile$outbreak[(indices[j]-testfile$time[indices[j]]):(indices[j]-1)] <- 0
      }
    
    prop_nocases <- testfile %>%
      group_by(sero, rate, disease, type, nweek) %>%
      summarise(prop = (n()-sum(outbreak))/n())
  
    } else{
      prop_nocases <- testfile %>%
        group_by(sero, rate, disease, type, nweek) %>%
        summarise(prop = (n()-sum(outbreak))/n())
    }
    
    if(i==1){
      first_elim_frame <- prop_nocases
    } else{
      first_elim_frame <- rbind(first_elim_frame, prop_nocases)
    }
    
    setTxtProgressBar(pb,i)
  }
  return(first_elim_frame)
}

elim_prob_t <- first_elim_prob()

# Prob elim w/in time figs ------------
elim_prob_vax <- elim_prob_t %>%
  filter(rate %in% c(1,5), sero %in% c(0, 0.8)) %>%
  group_by(sero, nweek, rate) %>%
  summarise(mean_prop = mean(prop))

# dev.new(width = 80, height = 60, unit = "mm", res = 600,
#         noRStudioGD=TRUE)

ggplot(data=elim_prob_vax, aes(x = nweek, y = mean_prop,
                                color = factor(sero),
                               linetype = factor(rate)))+
  geom_line(linewidth = 1)+
  # geom_smooth()+
  geom_hline(yintercept = 0.95, linetype="dashed")+
  scale_color_viridis_d(end = 0.9, name = "Vaccination Rate")+
  scale_linetype(name = "Immigration Rate")+
  labs(x = "Week", y = "Probability of 0 Cases")+
  theme_bw(base_size=16)+
  theme(panel.grid = element_blank())

# Note: this figure accounts for recolonization!

# ggsave(filename = "./full_Figs/full_elim_meant.jpeg", width=8,
#        dpi=600, units="in", height = 6)

# Look @ other vars
imm_figs <- list()

for(i in 1:length(unique(elim_prob_t$sero))){
  subset <- elim_prob_t %>%
    filter(rate %in% c(1,5),
           sero == unique(elim_prob_t$sero)[i]) %>%
    group_by(rate, nweek) %>%
    summarise(mean_prop = mean(prop))
  
  imm_figs[[i]] <- ggplot(data=subset, aes(x = nweek, 
                                           y = mean_prop,
                                           color = factor(rate)))+
    geom_line(linewidth = 1)+
    # geom_smooth()+
    geom_hline(yintercept = 0.95, linetype="dashed")+
    scale_color_viridis_d(end = 0.9, name = "Immigration Rate")+
    labs(x = "Week", y = "Probability of 0 Cases")+
    theme_bw(base_size=16)+
    theme(panel.grid = element_blank())
}

# ggsave(plot = imm_figs[[1]],
#        filename = "./full_Figs/full_elim_vax0.jpeg", width=8,
#        dpi=600, units="in", height = 6)

# R0 calculation ------------------------
get_r0 <- function(){
  # Get names of files
  filenames <- list.files(path = dir, pattern = "*.csv")
  
  # Progress bar
  pb = progressBar(min = 1, max = length(filenames), initial = 1,
                   style = "ETA") 
  
  for(i in 1:length(filenames)){
    # Read 'em in
    testfile <- read.csv(paste(getwd(), 
                               str_replace(dir, ".", ""), "/",
                               filenames[i],sep = ""))
    
    # Calculate time in weeks & get time to first elimination
    time_to_elim <- testfile %>%
      filter(elim == "True", year >= 2) %>%
      group_by(rep, sero, rate, type, disease) %>%
      mutate(nweek = ((year-1)*52) + week) %>%
      filter(nweek == min(nweek))
    
    r0 <- data.frame(rep = integer(), sero = numeric(),
                     rate = numeric(), disease = numeric(),
                     type = character(), r.0 = numeric())
    
    for(j in 1:length(unique(time_to_elim$rep))){
      onerep <- testfile %>%
        filter(rep == unique(time_to_elim$rep)[j])
      
      first_elim <- filter(time_to_elim, 
                           rep == unique(time_to_elim$rep)[j])
      
      endpoint <- ifelse(nrow(first_elim)<1, 52*11, 
                         first_elim$nweek)
    
      test <- try(
        suppressMessages({
          estimate.R(epid = onerep$n_symptomatic[53:endpoint], 
                     GT=generation.time("gamma", c(4.5, 1)),
                     method = 'TD', begin = 1, 
                     end = endpoint-53, nsim = 100)}),
        silent= T
      )

      # Error somewhere in here:
      if(class(test) %in% 'try-error') {next} else {
        r.est <- suppressMessages({
          estimate.R(epid = onerep$n_symptomatic[53:endpoint], 
                     GT=generation.time("gamma", c(4.5, 1)),
                     method = 'TD', nsim = 1000)
        })
        
      }
      r.est <- r.est$estimates$TD$R 
      
      med.rest <- median(r.est[r.est > quantile(r.est, 0.025, 
                                                na.rm = T) &
                                 r.est < quantile(r.est, 
                                                  0.975, 
                                                  na.rm = T)])
      
      if(is.na(med.rest)==T) {next} else{
        
        row <- data.frame(rep=unique(onerep$rep), 
                          sero=unique(onerep$sero),
                          rate=unique(onerep$rate),
                          disease=unique(onerep$disease), 
                          type=unique(onerep$type), 
                          r.0=med.rest)
          
        r0 <- bind_rows(r0, row)
      }
    }
    
    if(i==1){
      r0.frame <- r0
    } else{
      r0.frame <- rbind(r0.frame, r0)
    }
    setTxtProgressBar(pb,i)
  }
  return(r0.frame)
}

r0 <- get_r0()

saveRDS(r0, "r0.rds")

# R0 figures ---------------
ggplot(data = r0, aes(x = factor(sero), y = r.0))+
  geom_boxplot(fill = 'lightgray')+
  labs(x = "Adult Vaccination Rate", 
       y = bquote("Median"~R[e]))+
  theme_bw(base_size=14)+
  theme(panel.grid=element_blank())

# ggsave(filename = "./full_Figs/full_re.jpeg", width=8,
#        dpi=600, units="in", height = 6)

# R0 by week -----------------

# Function to calculate # weeks rabies-free -------------
rabies_free <- function(){
  # Get names of files
  filenames <- list.files(path = dir, pattern = "*.csv")
  
  # Progress bar
  pb = progressBar(min = 1, max = length(filenames), initial = 1,
                   style = "ETA") 
  
  # turn off dplyr warning because it messes up the progress bar
  options(dplyr.summarise.inform = FALSE)
  
  for(i in 1:length(filenames)){
    # Read 'em in
    testfile <- read.csv(paste(getwd(), 
                               str_replace(dir, ".", ""), "/",
                               filenames[i],sep = ""))

    # Calculate time in weeks & get time to first elimination
    time_rabies_free <- testfile %>%
      filter(year >= 1) %>%
      group_by(rep, sero, rate, type, disease) %>%
      # rate=disease rate of immigrants
      filter(elim == "True") %>%
      mutate(nweek = ((year-1)*52) + week) %>%
      summarise(n_rabies_free = n())
    
    if(i==1){
      elim_frame <- time_rabies_free
    } else{
      elim_frame <- rbind(elim_frame, time_rabies_free)
    }
    
    setTxtProgressBar(pb,i)
  }
  return(elim_frame)
}

time_rabies_free <- rabies_free()

# Figures: Time rabies-free --------------
time_rabies_free <- time_rabies_free %>%
  filter(n_rabies_free > 0)

# Seroprevalence boxplot
ggplot(data = time_rabies_free, aes(x = factor(sero),
                                    y=n_rabies_free))+
  geom_jitter(fill="lightgray")+
  # geom_text(aes(label=groups, y=525))+
  labs(x="Adult Vaccination Rate", y="Weeks with 0 Cases")+
  theme_bw(base_size=12)+
  theme(panel.grid = element_blank())

# ggsave(filename = "./full_Figs/nweekfree_box.jpeg",
#        width = 6, height = 4, dpi= 600, units = "in")

# Did not save these figures since they're the inverse 
# of time to elimination

# Probabilty and duration of reinvasion ----------------------
reinfection <- function(){
  # Get names of files
  filenames <- list.files(path = dir, pattern = "*.csv")
  
  # Progress bar
  pb = progressBar(min = 1, max = length(filenames), initial = 1,
                   style = "ETA") 
  
  for(i in 1:length(filenames)){
    testfile <- read.csv(paste(getwd(), 
                               str_replace(dir, ".", ""), "/",
                               filenames[i],sep = "")) %>%
      filter(year >= 2) %>%
      mutate(nweek = ((year-1)*52) + week)
    
    recols <- testfile %>%
      group_by(rep, sero, rate, type, disease) %>%
      filter(elim == "False" & lag(elim) == "True")
    
    elims <- testfile %>%
      group_by(rep, sero, rate, type, disease) %>%
      filter(elim == "True" & lag(elim) == "False")
    
    recol.time <- rbind(elims, recols) %>%
      arrange(rep, nweek) %>%
      mutate(time = case_when(elim == "True" ~ nweek-lag(nweek),
                              TRUE ~ NA)) %>%
      filter(is.na(time) == F) %>%
      filter(time > 9)
    
    reinf_frame <- suppressMessages(full_join(testfile, 
                                           recol.time)) %>%
      dplyr::filter(is.na(time)==F) %>%
      dplyr::select(rep, sero, disease, rate, 
                    type, nweek, time) %>%
      mutate(prop = length(unique(.$rep))/
               length(unique(elims$rep)))

    if(i==1){
      reinf_frame_full <- reinf_frame
    } else{
      reinf_frame_full <- rbind(reinf_frame_full, reinf_frame)
    }
    setTxtProgressBar(pb,i)
  }
  return(reinf_frame_full)
}

reinf_outs <- reinfection()

# Figures: reinfection prob & length -------------------
reinf_outs <- reinf_outs %>%
  filter(disease != 0) %>%
  mutate(sero=factor(sero)) %>%
  # mutate(TimePeriod = case_when(between(week,19,41) ~
  #                                 "Juveniles w/Mom",
  #        TRUE ~ "Juveniles Independent"))
  mutate(InfStart = (nweek-time) %% 52) %>%
  mutate(TimePeriod = case_when(InfStart < 18 | 
                                  InfStart >= 43 ~ 
                                  "Post-dispersal",
                                between(InfStart,18,28) ~
                                  "Juveniles w/Mom",
                                between(InfStart,29,42) ~ 
                                  "Independent Juveniles"))

probs_condensed <- reinf_outs %>%
  filter(rate %in% c(1,5,10)) %>%
  dplyr::select(sero, type, disease, rate, prop) %>%
  distinct()

summary(lm(data=reinf_outs, prop~as.numeric(sero)+rate))

ggplot(data=probs_condensed, aes(x=factor(sero),y=prop,
                                 fill = factor(rate)))+
  geom_boxplot()+
  # scale_fill_manual(values = c("lightgray", "limegreen"))+
  scale_fill_viridis_d(end=0.9, name = "Immigration Rate")+
  labs(x = "Adult Vaccination Rate", 
       y = "Recolonization Probability")+
  theme_bw(base_size=12)+
  theme(panel.grid = element_blank())

# ggsave(filename = "./full_Figs/rinfprob_sero.jpeg",
#        width = 6, height = 4, dpi= 600, units = "in")

# Immigrant infection rate
rate_inter_rinf <- reinf_outs %>%
  mutate(disease=factor(disease), sero=factor(sero))%>%
  group_by(sero,disease)%>%
  summarise(med = median(prop))

ggplot(data=rate_inter_rinf, aes(x=sero,y=disease,
                                fill=med))+
  geom_tile()+
  scale_fill_viridis(name="Recolonization Probability",
                     option = "B")+
  labs(x="Seroprevalence", y="Immigrant Disease Rate")+
  theme_bw(base_size=12)+
  theme(panel.grid = element_blank())

# Immigration rate
rate_inter_rinf <- reinf_outs %>%
  mutate(rate=factor(rate), sero=factor(sero))%>%
  group_by(sero,rate)%>%
  summarise(med = median(prop))

ggplot(data=rate_inter_rinf, aes(x=sero,y=rate,
                                 fill=med))+
  geom_tile()+
  scale_fill_viridis(name="Recolonization Probability",
                     option = "B")+
  labs(x="Seroprevalence", y="Immigration Rate")+
  theme_bw(base_size=12)+
  theme(panel.grid = element_blank())

# Immigration type
type_inter_rinf <- reinf_outs %>%
  mutate(sero=factor(sero))%>%
  group_by(sero,type)%>%
  summarise(med = median(reinf_prob))

ggplot(data=type_inter_rinf, aes(x=sero,y=type,
                                 fill=med))+
  geom_tile()+
  scale_fill_viridis(name="Recolonization Probability",
                     option = "B")+
  labs(x="Adult Vaccination Rate", y="Immigration Type")+
  theme_bw(base_size=12)+
  theme(panel.grid = element_blank())

# Reinfection timing --------------------
reinf_timing <- reinf_outs %>%
  mutate(WeekPerPeriod = 
           case_when(TimePeriod == "Post-dispersal" ~ (54-42)+17,
                     TimePeriod == "Independent Juveniles" ~ 42-29,
                     TRUE ~ 20)) %>%
  group_by(TimePeriod, WeekPerPeriod) %>%
  summarise(ProbPerPeriod = (n()/nrow(.))) %>%
  distinct() %>%
  mutate(PropPerPeriod = ProbPerPeriod/WeekPerPeriod)

ggplot(data = reinf_timing, aes(x = TimePeriod, 
                                y = PropPerPeriod))+
  geom_col(fill = 'lightgray', color = 'black') +
  labs(x = "Time Period", y = "Adj. Probability of Reinvasion")+
  theme_bw(base_size = 14)+
  theme(panel.grid = element_blank())

# ggsave(filename = "./full_figs/reinf_by_time.jpeg", height = 5,
#        width = 8, units = "in", dpi = 600)

# Reinfection length ------------------
reinf_outs <- filter(reinf_outs, prop > 0)

# Test for associations between elimination time 
# and reinfection length
summary(lm(data = reinf_outs, time~(nweek-time)))

ggplot(data = reinf_outs, aes(x = (nweek-time), y = time))+
  geom_point()+
  geom_abline(intercept = 52*11, slope = -1, linetype = 'dashed')+
  labs(x = "Reinfection Start", y = "Reinfection Length") +
  theme_bw(base_size = 14)+
  theme(panel.grid = element_blank())

# ggsave(filename = "./full_Figs/rinflength_starttime.jpeg",
#        width = 6, height = 4, dpi= 600, units = "in")

ggplot(data=reinf_outs, aes(x=TimePeriod,y=time))+
  geom_boxplot(fill = 'lightgray', outlier.shape = NA)+
  # geom_jitter()+
  scale_fill_viridis_d(name = "Time Period", end = 0.9)+
  scale_y_continuous(limits = c(0,75))+
  labs(x = "Adult Vaccination Rate", 
       y = "Reinfection Length (Weeks)")+
  theme_bw(base_size=12)+
  theme(panel.grid = element_blank())

# ggsave(filename = "./full_Figs/rinflength_sero_juv.jpeg",
#        width = 6, height = 4, dpi= 600, units = "in")

# Disease rate
dis_inter_rlen <- reinf_outs %>%
  mutate(disease=factor(disease), sero=factor(sero))%>%
  group_by(sero,disease)%>%
  summarise(med = median(time))

ggplot(data=dis_inter_rlen, aes(x=sero,y=disease,
                                fill=med))+
  geom_tile()+
  scale_fill_viridis(name="Disease Rate",
                     option = "B")+
  labs(x="Seroprevalence", y="Immigrant Disease Rate")+
  theme_bw(base_size=12)+
  theme(panel.grid = element_blank())

dis_inter_rlen <- reinf_outs %>%
  mutate(disease=factor(disease), sero=factor(sero))%>%
  group_by(sero,disease, TimePeriod)%>%
  summarise(med = median(time))

ggplot(data=dis_inter_rlen, aes(x=TimePeriod,y=disease,
                                fill=med))+
  geom_tile()+
  scale_fill_viridis(name="Reinfection Length",
                     option = "B")+
  labs(x="Seroprevalence", y="Immigrant Disease Rate")+
  theme_bw(base_size=12)+
  theme(panel.grid = element_blank())

# Immigration rate
rate_inter_rinf <- reinf_outs %>%
  mutate(rate=factor(rate), sero=factor(sero))%>%
  group_by(sero,rate)%>%
  summarise(med = median(time))

ggplot(data=rate_inter_rinf, aes(x=sero,y=rate,
                                fill=med))+
  geom_tile()+
  scale_fill_viridis(name="Reinfection Length (Weeks)",
                     option = "B")+
  labs(x="Adult Vaccination Rate", y="Immigration Rate")+
  theme_bw(base_size=12)+
  theme(panel.grid = element_blank())

# ggsave(filename = "./full_Figs/rinflength_rate_vax.jpeg",
#               width = 6, height = 4, dpi= 600, units = "in")

# Immigration type
type_inter_rinf <- reinf_outs %>%
  mutate(sero=factor(sero))%>%
  group_by(sero,type)%>%
  summarise(med = median(reinf_length))

ggplot(data=type_inter_rinf, aes(x=sero,y=type,
                                 fill=med))+
  geom_tile()+
  scale_fill_viridis(name="Reinfection Length (Weeks)",
                     option = "B")+
  labs(x="Seroprevalence", y="Immigration Type")+
  theme_bw(base_size=12)+
  theme(panel.grid = element_blank())

# Total cases -----------------------
cases <- function(){
  # Get names of files
  filenames <- list.files(path = dir, pattern = "*.csv")
  
  # Progress bar
  pb = progressBar(min = 1, max = length(filenames), initial = 1,
                   style = "ETA") 
  
  for(i in 1:length(filenames)){
    # Read 'em in
    testfile <- read.csv(paste(getwd(), 
                               str_replace(dir, ".", ""), "/",
                               filenames[i],sep = ""))

    cases_frame <- testfile %>%
      mutate(nweek = ((year-1)*52) + week) %>%
      group_by(rep, sero, rate, type, disease) %>%
      summarise(ncases = sum(n_symptomatic))

    if(i==1){
      cases_frame_full <- cases_frame
    } else{
      cases_frame_full <- rbind(cases_frame_full, cases_frame)
    }
    
    setTxtProgressBar(pb,i)
  }
  return(cases_frame_full)
}

total_cases <- cases()

# Total cases: figs ---------------------
# sero
# dev.new(width = 80, height = 60, unit = "mm", res=600,
#         noRStudioGD=TRUE)

ggplot(data=total_cases, aes(x=factor(sero),y=ncases))+
  geom_boxplot(fill='lightgray')+
  # geom_text(aes(y=13500, label = groups))+
  labs(x = "Seroprevalence", y = "Total Cases")+
  theme_bw(base_size=16)+
  theme(panel.grid = element_blank())

# ggsave(filename = "./full_Figs/totalcases_sero.jpeg",
#        width = 6, height = 4, dpi= 600, units = "in")

# Disease rate
dis_inter_cases <- total_cases %>%
  mutate(disease=factor(disease), sero=factor(sero))%>%
  group_by(sero,disease) %>%
  summarise(med = median(ncases))

ggplot(data=dis_inter_cases, aes(x=sero,y=disease,
                                fill=med))+
  geom_tile()+
  scale_fill_viridis(name="Total Cases",
                     option = "B")+
  labs(x="Seroprevalence", y="Disease Rate")+
  theme_bw(base_size=12)+
  theme(panel.grid = element_blank())

# Immigration rate
rate_inter_cases <- total_cases %>%
  mutate(rate=factor(rate), sero=factor(sero))%>%
  group_by(sero,rate) %>%
  summarise(med = median(ncases))

ggplot(data=rate_inter_cases, aes(x=sero,y=rate,
                                 fill=med))+
  geom_tile()+
  scale_fill_viridis(name="Total Cases",
                     option = "B")+
  labs(x="Seroprevalence", y="Immigration Rate")+
  theme_bw(base_size=12)+
  theme(panel.grid = element_blank())

# Immigration type
type_inter_cases <- total_cases %>%
  mutate(sero=factor(sero))%>%
  group_by(sero,type) %>%
  summarise(med = median(ncases))

ggplot(data=type_inter_cases, aes(x=sero,y=type,
                                  fill=med))+
  geom_tile()+
  scale_fill_viridis(name="Total Cases",
                     option = "B")+
  labs(x="Seroprevalence", y="Immigration Type")+
  theme_bw(base_size=12)+
  theme(panel.grid = element_blank())

# Total cases after first elimination ------------------
reinf_cases <- function(){
  # Get names of files
  filenames <- list.files(path = dir, pattern = "*.csv")
  
  # Progress bar
  pb = progressBar(min = 1, max = length(filenames), initial = 1,
                   style = "ETA") 
  
  for(i in 1:length(filenames)){
    # Read 'em in
    testfile <- read.csv(paste(getwd(), 
                               str_replace(dir, ".", ""), "/",
                               filenames[i],sep = ""))
    
    # Calculate time in weeks & get time to first elimination
    time_to_elim <- testfile %>%
      group_by(rep, sero, rate, disease, type) %>%
      filter(year >= 2) %>%
      # rate=disease rate of immigrants
      filter(elim == "True") %>%
      mutate(nweek = ((year-1)*52) + week) %>%
      filter(nweek == min(nweek))
    
    f_elim <- data.frame(first_elim=time_to_elim$nweek,
                             rep=time_to_elim$rep)
    
    reinf_frame <- testfile %>%
      left_join(f_elim, by = "rep") %>%
      # rate=disease rate of immigrants
      mutate(nweek = ((year-1)*52) + week) %>%
      group_by(rep,sero, rate, type, disease) %>%
      filter(nweek > first_elim) %>%
      filter(elim != "True") %>%
      summarise(ncases = sum(n_symptomatic))
    
    if(i==1){
      reinf_frame_full <- reinf_frame
    } else{
      reinf_frame_full <- rbind(reinf_frame_full, reinf_frame)
    }
    
    setTxtProgressBar(pb,i)
  }
  return(reinf_frame_full)
}

reinf_case_frame <- reinf_cases()

# Cases post-reinfection: Figures -------------------
# sero
reinf_case_frame <- reinf_case_frame %>%
  mutate(sero=factor(sero))

ggplot(data=reinf_case_frame, aes(x=sero,y=ncases))+
  geom_boxplot(fill='lightgray')+
  # geom_text(aes(y=1300, label = groups))+
  labs(x = "Seroprevalence", y = "Cases Post-Reinvasion")+
  theme_bw(base_size=12)+
  theme(panel.grid = element_blank())

# ggsave(filename = "./full_Figs/reinfcases_sero.jpeg",
#        width = 6, height = 4, dpi= 600, units = "in")

# Disease Rate
dis_inter_rcases <- reinf_case_frame %>%
  mutate(disease=factor(disease), sero=factor(sero))%>%
  group_by(sero,disease) %>%
  summarise(med = median(ncases))

ggplot(data=dis_inter_rcases, aes(x=sero,y=disease,
                                 fill=med))+
  geom_tile()+
  scale_fill_viridis(name="Cases Post-Reinvasion",
                     option = "B")+
  labs(x="Seroprevalence", y="Disease Rate")+
  theme_bw(base_size=12)+
  theme(panel.grid = element_blank())

# Immigration rate
rate_inter_rcases <- reinf_case_frame %>%
  mutate(rate=factor(rate), sero=factor(sero))%>%
  group_by(sero,rate) %>%
  summarise(med = median(ncases))

ggplot(data=rate_inter_rcases, aes(x=sero,y=rate,
                                  fill=med))+
  geom_tile()+
  scale_fill_viridis(name="Cases Post-Reinvasion",
                     option = "B")+
  labs(x="Seroprevalence", y="Immigration Rate")+
  theme_bw(base_size=12)+
  theme(panel.grid = element_blank())

# Cases per week ----------------------
cases_per_week <- function(metric){
  # Get names of files
  filenames <- list.files(path = dir, pattern = "*.csv")
  
  # Progress bar
  pb = progressBar(min = 1, max = length(filenames), initial = 1,
                   style = "ETA") 
  
  for(i in 1:length(filenames)){
    # Read 'em in
    testfile <- read.csv(paste(getwd(), 
                               str_replace(dir, ".", ""), "/",
                               filenames[i],sep = ""))
    
    cases_frame <- testfile %>%
      mutate(nweek = ((year-1)*52) + week) %>%
      group_by(sero, rate, type, disease, nweek)
    
    if(metric=="mean"){
      cases_frame <- cases_frame %>%
        summarise(mean_cases = mean(n_symptomatic))
    }else if(metric == "max"){
      cases_frame <- cases_frame %>%
        summarise(max_cases = max(n_symptomatic))
    }else if(metric == 'median'){
      cases_frame <- cases_frame %>%
        summarise(median_cases = median(n_symptomatic))
    }
    
    if(i==1){
      cases_frame_full <- cases_frame
    } else{
      cases_frame_full <- rbind(cases_frame_full, cases_frame)
    }
    
    setTxtProgressBar(pb,i)
  }
  return(cases_frame_full)
}

# Mean cases per week ----------------
meancase <- cases_per_week(metric="median")

meancase_seros <- meancase %>%
  filter(nweek > 52) %>%
  filter(disease == 0)

# dev.new(width = 80, height = 60, unit = "mm", res=600,
#         noRStudioGD=TRUE)

ggplot(data=meancase_seros, aes(x=nweek, y=median_cases, 
                                color = factor(sero)))+
  geom_line()+
  # geom_vline(xintercept=c(52*3+20, 52*4+20, 52*5+20, 52*6+20),
  #            linetype="dashed", size = 1)+
  scale_color_viridis_d(end = 0.9, name="Vaccination Rate")+
  labs(x = "Week", y = "Median Cases")+
  theme_bw(base_size=16)+
  theme(panel.grid = element_blank())

# ggsave(filename = "./full_Figs/meanwkcases_sero.jpeg",
#        width = 6, height = 4, dpi= 600, units = "in")

# Avg Prevalence -----------------------
prev <- function(metric){
  # Get names of files
  filenames <- list.files(path = "./Outputs", pattern = "*.csv")
  
  # Progress bar
  pb = progressBar(min = 1, max = length(filenames), initial = 1,
                   style = "ETA") 
  
  for(i in 1:length(filenames)){
    # Read 'em in
    testfile <- read.csv(paste(getwd(),"/Outputs/",filenames[i],
                               sep = ""))

    cases_frame <- testfile %>%
      mutate(nweek = ((year-1)*52) + week) %>%
      group_by(rep, sero, rate, type, disease, nweek) %>%
      summarise(prev = n_infected/total_pop) 
    
    if(metric=="mean"){
      cases_frame <- cases_frame %>%
        ungroup() %>%
        group_by(sero, rate, type, disease, nweek) %>%
        summarise(mean_prev = mean(prev))
    }else if(metric == "max"){
      cases_frame <- cases_frame %>%
        ungroup() %>%
        group_by(sero, rate, type, disease, nweek) %>%
        summarise(max_prev = max(prev))
    }else if(metric == "median"){
      cases_frame <- cases_frame %>%
        ungroup() %>%
        group_by(sero, rate, type, disease, nweek) %>%
        summarise(median_prev = median(prev))
    }

    if(i==1){
      cases_frame_full <- cases_frame
    } else{
      cases_frame_full <- rbind(cases_frame_full, cases_frame)
    }
    
    setTxtProgressBar(pb,i)
  }
  return(cases_frame_full)
}

prev.mean <- prev(metric = "median")

# Prevalence figs -------------------
prev.sero <- prev.mean %>%
  filter(nweek >= 52) %>%
  filter(disease == 0.03) %>%
  group_by(sero, nweek) %>%
  summarise(med = median(median_prev))

# Sero only
ggplot(data=prev.sero, aes(x=nweek, y=med, 
                               color = factor(sero)))+
  geom_line()+
  scale_color_viridis_d(end=0.9, name = "Vaccination Rate")+
  labs(x = "Weeks", y = "Rabies Prevalence")+
  theme_bw(base_size = 12)+
  theme(panel.grid = element_blank())

# ggsave(filename = "./full_Figs/meanwkprev_sero.jpeg",
#        width = 6, height = 4, dpi= 600, units = "in")

# Weekly population -----------------------------
get_weekly_pop <- function(){
  # Get names of files
  filenames <- list.files(path = "./Outputs", pattern = "*.csv")
  
  # Progress bar
  pb = progressBar(min = 1, max = length(filenames), initial = 1,
                   style = "ETA") 
  
  for(i in 1:length(filenames)){
    # Read 'em in
    testfile <- read.csv(paste(getwd(),"/Outputs/",filenames[i],
                               sep = ""))
    
    pop_frame <- testfile %>%
      mutate(nweek = ((year-1)*52) + week) %>%
      group_by(sero, rate, type, disease, nweek) %>%
      summarise(mean_pop = mean(total_pop))
    
    if(i==1){
      pop_frame_full <- pop_frame
    } else{
      pop_frame_full <- rbind(pop_frame_full, pop_frame)
    }
    
    setTxtProgressBar(pb,i)
  }
  return(pop_frame_full)
}

weekly_pop <- get_weekly_pop()

# Population figs -------------------------
pop.sero <- weekly_pop %>%
  # filter(disease == 0, rate == 5)
  filter(disease == 0.015) %>%
  group_by(sero, nweek) %>%
  summarise(med_pop = median(mean_pop))

# Sero only, full sim
ggplot(data=pop.sero, aes(x=nweek, y=med_pop, 
                           color = factor(sero)))+
  geom_line()+
  geom_vline(xintercept=53, linetype="dashed")+
  scale_color_viridis_d(end=0.9, name = "Vaccination Rate")+
  labs(x = "Weeks", y = "Population Size")+
  theme_bw(base_size = 14)+
  theme(panel.grid = element_blank())

# ggsave(filename = "./full_Figs/medtotalpop.jpeg",
#        width = 6, height = 4, dpi= 600, units = "in")
