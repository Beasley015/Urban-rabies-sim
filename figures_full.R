###########################################
# Figures for urban rabies simulation     # 
# Full model: seroprevalence, immigration,#
# and landscape barriers                  #
# Model written in Julia                  #
# E.M. Beasley                            #
# Fall 2023                               #
###########################################

# Load packages ---------------------
library(pbmcapply) #Progress bar w/ETA
library(tidyverse)
library(viridis)
library(quantreg)

options(dplyr.summarise.inform = FALSE)

# Function for calculating time to first elimination --------
first_elim <- function(){
  # Get names of files
  filenames <- list.files(path = "./outs", pattern = "*.csv")

  # Progress bar
  pb = progressBar(min = 1, max = length(filenames), initial = 1,
                   style = "ETA") 
  
  for(i in 1:length(filenames)){
    # Read 'em in
    testfile <- read.csv(paste(getwd(),"/outs/",filenames[i],sep = ""))

    # add seroprevalence of immigrants because I forgot
    imm_sero <- str_extract(filenames[i], "(im_ser)(\\d*\\.\\d*)", 
                        group = 2)

    # Calculate time in weeks & get time to first elimination
    time_to_elim <- testfile %>%
      mutate(im_sero = imm_sero) %>%
      group_by(rep, sero, type, rate, im_sero, barrier_val) %>%
      # rate=disease rate of immigrants
      filter(n_infected == 0 & n_symptomatic == 0) %>%
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

# Time to elimination: prelim figs ------------------
# VERY prelim check
summary(aov(data = first_elim_full, formula = nweek~factor(sero)+factor(rate)+im_sero+factor(barrier_val)+type))
# seroprevalence, barriers, and immigration type may be important

# any interactions look promising?
summary(lm(data = first_elim_full, formula = nweek~sero+barrier_val+type+sero*barrier_val+sero*type+barrier_val*type))
# seroprevalence may be interacting with immigration type
# barrier value appers to be independent

# is the effect of seroprevalence similar to prelim model?
# filter out 1 barrier value & 1 immigration type, then plot
elim_sansbar_type <- first_elim_full %>%
  filter(barrier_val == 0 & type == "propagule")

ggplot(data = elim_sansbar_type, aes(x=factor(sero),y=nweek))+
  geom_boxplot(fill='lightgray')+
  labs(x = "Seroprevalence", y = "Time to Elimination (Weeks)")+
  theme_bw()+
  theme(panel.grid=element_blank())

# ggsave(filename = "./full_Figs/sero_elim_box.jpeg", width=6,
#        height=4, dpi=600, units="in")
# looks similar, try quantile regression again

summary(rq(data = elim_sansbar_type, 
           formula = nweek~sero,
           tau = c(0.05, 0.25, 0.5, 0.75, 0.95)))

rq.frame <- tibble(slopes = c(-22.5, -25, -38.75, -33, -177.5),
                   inters = c(75, 85, 106, 131, 265.5),
                   quants = c("5%", "25%", "50%", "75%", "95%"))%>%
  mutate(quants = factor(quants, levels = quants))

ggplot(data = elim_sansbar_type, aes(x = sero, y = nweek))+
  geom_point()+
  geom_abline(data = rq.frame, aes(intercept=inters, slope = slopes,
                                   color = quants))+
  labs(x = "Seroprevalence", y = "Outbreak Duration (Weeks)")+
  scale_color_viridis_d(end = 0.9, name = "Quantile")+
  theme_bw()+
  theme(panel.grid = element_blank())

# ggsave(filename = "./full_Figs/sero_elim_quantreg.jpeg", width=6,
#        height=4, dpi=600, units="in")

# Look at barrier strength
elim_sanssero_type <- first_elim_full %>%
  filter(sero == 0.8 & type == "propagule")

ggplot(data=elim_sanssero_type, aes(x = factor(barrier_val),
                                    y = nweek))+
  geom_boxplot()
# tried this figure at every sero value, 
# barrier strength doesn't seem to matter

# look at immigration type
# set barrier
bar_set <- first_elim_full %>%
  filter(barrier_val == 2)

ggplot(data = bar_set, aes(x=factor(sero), y=nweek, fill=type))+
  geom_boxplot()

#OVERALL: MAIN DRIVER OF TIME TO FIRST ELIMINATION IS SEROPREVALENCE

# Function to calculate # weeks rabies-free -------------
rabies_free <- function(){
  # Get names of files
  filenames <- list.files(path = "./outs", pattern = "*.csv")
  
  # Progress bar
  pb = progressBar(min = 1, max = length(filenames), initial = 1,
                   style = "ETA") 
  
  # turn off dplyr warning because it messes up the progress bar
  options(dplyr.summarise.inform = FALSE)
  
  for(i in 1:length(filenames)){
    # Read 'em in
    testfile <- read.csv(paste(getwd(),"/outs/",filenames[i],
                               sep = ""))
    
    # add seroprevalence of immigrants because I forgot
    imm_sero <- str_extract(filenames[i], "(im_ser)(\\d*\\.\\d*)", 
                            group = 2)
    
    # Calculate time in weeks & get time to first elimination
    time_rabies_free <- testfile %>%
      mutate(im_sero = imm_sero) %>%
      group_by(rep, sero, type, rate, im_sero, barrier_val) %>%
      # rate=disease rate of immigrants
      filter(n_infected == 0 & n_symptomatic == 0) %>%
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

# Prelim figures --------------
# VERY prelim check again
summary(lm(data = time_rabies_free, formula = n_rabies_free~sero+rate+as.numeric(im_sero)+barrier_val+type))
# about the same as time to elimination

# interactions?
summary(lm(data = time_rabies_free, 
           formula = n_rabies_free~sero+barrier_val+sero*barrier_val))
# nope

# Boxplot
ggplot(data = time_rabies_free, 
       aes(x = factor(sero),y=n_rabies_free))+
  geom_boxplot(fill="lightgray")+
  labs(x="Seroprevalence", y="Weeks with 0 Cases")+
  theme_bw()+
  theme(panel.grid=element_blank())

ggsave(filename = "./full_Figs/rabies_free_box.jpeg", width=6,
       height=4, dpi=600, units="in")

ggplot(data = time_rabies_free, 
       aes(x = factor(sero),y=n_rabies_free,
           fill=factor(barrier_val)))+
  geom_boxplot()
# same deal, seroprevalence matters the most and affects extreme outcomes

# Probabilty and duration of reinvasion ----------------------
reinfection <- function(){
  # Get names of files
  filenames <- list.files(path = "./outs", pattern = "*.csv")
  
  # Progress bar
  pb = progressBar(min = 1, max = length(filenames), initial = 1,
                   style = "ETA") 
  
  for(i in 1:length(filenames)){
    # Read 'em in
    testfile <- read.csv(paste(getwd(),"/outs/",filenames[i],sep = ""))
    
    # add seroprevalence of immigrants because I forgot
    imm_sero <- str_extract(filenames[i], "(im_ser)(\\d*\\.\\d*)", 
                            group = 2)
    
    # Calculate time in weeks & get time to first elimination
    time_to_elim <- testfile %>%
      mutate(im_sero = imm_sero) %>%
      group_by(rep, sero, type, rate, im_sero, barrier_val) %>%
      # rate=disease rate of immigrants
      filter(n_infected == 0 & n_symptomatic == 0) %>%
      mutate(nweek = ((year-1)*52) + week) %>%
      filter(nweek == min(nweek))
    
    first_elim <- data.frame(elim=time_to_elim$nweek,
                             rep=time_to_elim$rep)
    
    reinf_frame <- testfile %>%
      left_join(first_elim, by = "rep") %>%
      mutate(im_sero = imm_sero) %>%
      # rate=disease rate of immigrants
      mutate(nweek = ((year-1)*52) + week) %>%
      group_by(rep,sero, type, rate, im_sero, barrier_val) %>%
      filter(nweek > elim) %>%
      filter(n_infected > 0 | n_symptomatic > 0) %>%
      summarise(reinf_length = n()) %>%
      filter(reinf_length > 5) %>%
      ungroup() %>%
      mutate(reinf_prob = length(unique(rep))/50)
    
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
# quick look
summary(aov(reinf_prob~factor(sero)+factor(rate)+im_sero+factor(barrier_val)+type,data=reinf_outs))

ggplot(data=reinf_outs, aes(x=factor(sero),y=reinf_prob))+
  geom_boxplot(fill="lightgray")+
  labs(x="Seroprevalence", y="Probability of Reinvasion")+
  theme_bw()+
  theme(panel.grid=element_blank())

ggsave(filename = "./full_Figs/reinf_box.jpeg", width=6,
       height=4, dpi=600, units="in")

reinf_mean <- reinf_outs %>%
  group_by(sero, barrier_val, rate, im_sero) %>%
  summarise(mean_prob = mean(reinf_prob))

ggplot(data=reinf_outs, aes(x=factor(sero),y=reinf_length))+
  geom_boxplot(fill="lightgray")+
  labs(x="Seroprevalence", y="Length of Reinvasion (Weeks)")+
  theme_bw()+
  theme(panel.grid=element_blank())

# ggsave(filename = "./full_Figs/reinf_length.jpeg", width=6,
#        height=4, dpi=600, units="in")

# Total cases -----------------------
cases <- function(){
  # Get names of files
  filenames <- list.files(path = "./outs", pattern = "*.csv")
  
  # Progress bar
  pb = progressBar(min = 1, max = length(filenames), initial = 1,
                   style = "ETA") 
  
  for(i in 1:length(filenames)){
    # Read 'em in
    testfile <- read.csv(paste(getwd(),"/outs/",filenames[i],
                               sep = ""))
    
    # add seroprevalence of immigrants because I forgot
    imm_sero <- str_extract(filenames[i], "(im_ser)(\\d*\\.\\d*)", 
                            group = 2)
    
    cases_frame <- testfile %>%
      mutate(im_sero = imm_sero) %>%
      # rate=disease rate of immigrants
      mutate(nweek = ((year-1)*52) + week) %>%
      group_by(rep,sero, type, rate, im_sero, barrier_val) %>%
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
summary(aov(data=total_cases, ncases~factor(sero)+factor(rate)+
              factor(barrier_val)+im_sero+type))

cases_smol <- total_cases %>%
  filter(type == "propagule" & im_sero == "0.4" & rate == 0.05)

ggplot(data=cases_smol, aes(x=factor(sero), y = ncases,
                            fill=factor(barrier_val)))+
  geom_boxplot()+
  scale_fill_viridis_d(end=0.9, name = "Barrier Strength")+
  labs(x = "Seroprevalence", y = "Total Cases")+
  theme_bw()+
  theme(panel.grid = element_blank())

# ggsave(filename = "./full_Figs/total_cases.jpeg", width=6,
#        height=4, dpi=600, units="in")

# Total cases after first elimination ------------------
reinf_cases <- function(){
  # Get names of files
  filenames <- list.files(path = "./outs", pattern = "*.csv")
  
  # Progress bar
  pb = progressBar(min = 1, max = length(filenames), initial = 1,
                   style = "ETA") 
  
  for(i in 1:length(filenames)){
    # Read 'em in
    testfile <- read.csv(paste(getwd(),"/outs/",filenames[i],sep = ""))
    
    # add seroprevalence of immigrants because I forgot
    imm_sero <- str_extract(filenames[i], "(im_ser)(\\d*\\.\\d*)", 
                            group = 2)
    
    # Calculate time in weeks & get time to first elimination
    time_to_elim <- testfile %>%
      mutate(im_sero = imm_sero) %>%
      group_by(rep, sero, type, rate, im_sero, barrier_val) %>%
      # rate=disease rate of immigrants
      filter(n_infected == 0 & n_symptomatic == 0) %>%
      mutate(nweek = ((year-1)*52) + week) %>%
      filter(nweek == min(nweek))
    
    first_elim <- data.frame(elim=time_to_elim$nweek,
                             rep=time_to_elim$rep)
    
    reinf_frame <- testfile %>%
      left_join(first_elim, by = "rep") %>%
      mutate(im_sero = imm_sero) %>%
      # rate=disease rate of immigrants
      mutate(nweek = ((year-1)*52) + week) %>%
      group_by(rep,sero, type, rate, im_sero, barrier_val) %>%
      filter(nweek > elim) %>%
      filter(n_infected > 0 | n_symptomatic > 0) %>%
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
summary(aov(data=reinf_case_frame, ncases~factor(sero)+factor(rate)+
              factor(barrier_val)+im_sero+type))

reinf_cases_smol <- reinf_case_frame %>%
  filter(barrier_val == 2 & type == "propagule")

ggplot(data = reinf_case_frame, aes(x=factor(sero), y=ncases))+
  geom_boxplot(fill="lightgray")+
  labs(x="Seroprevalence",y="Cases Post-Elimination")+
  theme_bw()+
  theme(panel.grid = element_blank())

# ggsave(filename = "./full_Figs/cases_post_elim.jpeg", width=6,
#        height=4, dpi=600, units="in")

# Cases per week ----------------------
cases_per_week <- function(metric){
  # Get names of files
  filenames <- list.files(path = "./outs", pattern = "*.csv")
  
  # Progress bar
  pb = progressBar(min = 1, max = length(filenames), initial = 1,
                   style = "ETA") 
  
  for(i in 1:length(filenames)){
    # Read 'em in
    testfile <- read.csv(paste(getwd(),"/outs/",filenames[i],
                               sep = ""))
    
    # add seroprevalence of immigrants because I forgot
    imm_sero <- str_extract(filenames[i], "(im_ser)(\\d*\\.\\d*)", 
                            group = 2)
    
    cases_frame <- testfile %>%
      mutate(im_sero = imm_sero) %>%
      # rate=disease rate of immigrants
      mutate(nweek = ((year-1)*52) + week) %>%
      group_by(sero, type, rate, im_sero, barrier_val,nweek)
    
    if(metric=="mean"){
      cases_frame <- cases_frame %>%
        summarise(mean_cases = mean(n_symptomatic))
    }else if(metric == "max"){
      cases_frame <- cases_frame %>%
        summarise(max_cases = max(n_symptomatic))
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
meancase <- cases_per_week(metric="mean")

summary(aov(data=meancase, formula=mean_cases~factor(sero)+
              factor(barrier_val)+factor(rate)+im_sero+type+nweek))

# filter out 1 barrier value
meancase_seros <- meancase %>%
  filter(barrier_val == 2)

ggplot(data=meancase_seros, aes(x=nweek, y=mean_cases, 
                                color = factor(sero)))+
  geom_line()+
  geom_vline(xintercept=20, linetype="dashed")

# filter out 1 seroprevalence value
meancase_barrier <- meancase %>%
  filter(sero == 0.8)

ggplot(data=meancase_barrier, aes(x=nweek, y=mean_cases,
                                color = factor(barrier_val)))+
  geom_line()+
  geom_vline(xintercept=20, linetype = "dashed")

# Look after 100 weeks
meancase_seros_late <- meancase_seros %>%
  filter(nweek > 100)

ggplot(data=meancase_seros_late, aes(x=nweek, y=mean_cases, 
                                color = factor(sero)))+
  geom_line()

meancase_barrier_late <- meancase_barrier %>%
  filter(nweek > 100)

ggplot(data=meancase_barrier_late, aes(x=nweek, y=mean_cases, 
                                     color = factor(barrier_val)))+
  geom_line()

# Maximum cases per week ------------------------
maxcase <- cases_per_week(metric="max")

summary(aov(data=maxcase, formula=max_cases~factor(sero)+
              factor(barrier_val)+factor(rate)+im_sero+type+nweek))

# filter out 1 barrier value
maxcase_seros <- maxcase %>%
  filter(barrier_val == 2)

ggplot(data=maxcase_seros, aes(x=nweek, y=max_cases, 
                                color = factor(sero)))+
  geom_line()+
  geom_vline(xintercept=20, linetype="dashed")

# filter out 1 seroprevalence value
maxcase_barrier <- maxcase %>%
  filter(sero == 0.8)

ggplot(data=maxcase_barrier, aes(x=nweek, y=max_cases,
                                  color = factor(barrier_val)))+
  geom_line()+
  geom_vline(xintercept=20, linetype = "dashed")

# Look after 100 weeks
maxcase_seros_late <- maxcase_seros %>%
  filter(nweek > 100)

ggplot(data=maxcase_seros_late, aes(x=nweek, y=max_cases, 
                                     color = factor(sero)))+
  geom_line()

maxcase_barrier_late <- maxcase_barrier %>%
  filter(nweek > 100)

ggplot(data=maxcase_barrier_late, aes(x=nweek, y=max_cases, 
                                       color = factor(barrier_val)))+
  geom_line()
