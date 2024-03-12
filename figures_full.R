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
library(agricolae)

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
      group_by(rep, sero, type, rate, im_sero, barrier) %>%
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
# gross figs to start
ggplot(data=first_elim_full, aes(x=actual_sero, y=nweek,
                                 color=factor(rate)))+
  geom_point()

ggplot(data = first_elim_full, aes(x=factor(sero), y=nweek,
                                   fill=type))+
  geom_boxplot()

ggplot(data = first_elim_full, aes(x=factor(sero), y=nweek,
                                   fill=factor(barrier)))+
  geom_boxplot()

# is the effect of seroprevalence similar to prelim model?
# filter out 1 barrier value & 1 immigration type, then plot
elim_sansbar_type <- first_elim_full %>%
  filter(barrier == 3 & type == "propagule") %>%
  mutate(sero = factor(sero))

sero_nweek <- aov(elim_sansbar_type, formula=nweek~sero)
sero_nweek_HSD <- HSD.test(sero_nweek, 
                           trt = "sero")$groups %>%
  mutate(sero = rownames(.)) %>%
  select(-nweek)

elim_sansbar_type <- elim_sansbar_type %>%
  left_join(sero_nweek_HSD, by = "sero")
  
ggplot(data = elim_sansbar_type, aes(x=factor(sero),y=nweek))+
  geom_boxplot(fill='lightgray')+
  geom_text(aes(y = 500, label = groups))+
  labs(x = "Seroprevalence", y = "Time to Elimination (Weeks)")+
  theme_bw()+
  theme(panel.grid=element_blank())

# ggsave(filename = "./full_Figs/sero_elim.jpeg", width=6,
#        height=4, dpi=600, units="in")

bar_interac <- first_elim_full %>%
  group_by(sero, barrier) %>%
  summarise(center = median(nweek))

ggplot(data = bar_interac, aes(x = factor(sero), y = barrier,
                               fill=center))+
  geom_tile()+
  scale_fill_viridis(name = "Weeks to Elimination",
                       option="B")+
  labs(x = "Seroprevalence", y = "Barrier Strength")+
  theme_bw(base_size=12)+
  theme(panel.grid = element_blank())

# ggsave(filename = "./full_Figs/sero_barrier_interac.jpeg",
#        width = 6, height = 4, dpi= 600, units = "in")

# interaction: type
type_interac <- first_elim_full %>%
  group_by(sero, type) %>%
  summarise(center = median(nweek))

ggplot(data = type_interac, aes(x = sero, y = center, 
                                color = type))+
  geom_line()

# interaction: proportion diseased immigrants
rate_interac <- first_elim_full %>%
  mutate(rate = factor(rate), sero=factor(sero)) %>%
  group_by(sero, rate) %>%
  summarise(center = median(nweek))

ggplot(data = rate_interac, aes(x = sero, y = rate, 
                                fill = center))+
  geom_tile()+
  scale_fill_viridis(name = "Weeks to Elimination", option="B")+
  labs(x ="Seroprevalence", y="Immigrant Infection Rate")+
  theme_bw(base_size=12)+
  theme(panel.grid=element_blank())
# 
# ggsave(filename = "./full_Figs/sero_iminfec_interac.jpeg",
#        width = 6, height = 4, dpi= 600, units = "in")

# interaction: proportion vaccinated immigrants
vax_interac <- first_elim_full %>%
  group_by(sero, im_sero) %>%
  summarise(center = median(nweek))

ggplot(data = vax_interac, aes(x = factor(sero), y = im_sero, 
                                fill=center))+
  geom_tile()+
  scale_fill_viridis(name="Weeks to Elimination", option="B")+
  labs(x = "Seroprevalence", y="Immigrant Seroprevalence")+
  theme_bw(base_size=12)+
  theme(panel.grid=element_blank())

# ggsave(filename = "./full_Figs/sero_imsero_interac.jpeg",
#        width = 6, height = 4, dpi= 600, units = "in")

#OVERALL: MAIN DRIVER OF TIME TO FIRST ELIMINATION IS SEROPREVALENCE
#Longest duration at 0.2 seroprevalence
#0.6 is actually a good sero threshold for this metric
#Weird interactions with other variables at 0.2 seroprevalence

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
      group_by(rep, sero, type, rate, im_sero, barrier) %>%
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

# Figures: Time rabies-free --------------
# Start with interactions this time
# Barrier:
bar_inter_rfree <- time_rabies_free %>%
  mutate(barrier=factor(barrier), sero=factor(sero))%>%
  group_by(sero,barrier)%>%
  summarise(med = median(n_rabies_free))

ggplot(data=bar_inter_rfree, aes(x=sero,y=barrier,
                                 fill=med))+
  geom_tile()+
  scale_fill_viridis(end=0.9, name="Weeks Rabies-Free")+
  labs(x="Seroprevalence", y="Barrier Strength")+
  theme_bw(base_size=12)+
  theme(panel.grid = element_blank())

# ggsave(filename = "./full_Figs/nweekfree_barrier_interac.jpeg",
#        width = 6, height = 4, dpi= 600, units = "in")

# Immigration type
type_inter_rfree <- time_rabies_free %>%
  group_by(sero,type)%>%
  summarise(med = median(n_rabies_free))

ggplot(data=type_inter_rfree, aes(x=factor(sero),y=med,
                                 color=type))+
  geom_line(size=1.5)+
  scale_color_viridis_d(end=0.9, name="Immigration Type")+
  labs(x="Seroprevalence", y="Median Rabies-Free Weeks")+
  theme_bw(base_size=12)+
  theme(panel.grid = element_blank())

# Proportion diseased immigrants
iminfec_inter_rfree <- time_rabies_free %>%
  mutate(rate=factor(rate))%>%
  group_by(sero,rate)%>%
  summarise(med = median(n_rabies_free))

ggplot(data=iminfec_inter_rfree, aes(x=factor(sero),y=rate,
                                 fill=med))+
  geom_tile()+
  scale_fill_viridis(name="Weeks Rabies-Free", option="D")+
  labs(x="Seroprevalence", y="% Infected Immigrants")+
  theme_bw(base_size=12)+
  theme(panel.grid = element_blank())

# ggsave(filename = "./full_Figs/nweekfree_iminfec_interac.jpeg",
#        width = 6, height = 4, dpi= 600, units = "in")

imsero_inter_rfree <- time_rabies_free %>%
  group_by(sero,im_sero)%>%
  summarise(med = median(n_rabies_free))

ggplot(data=imsero_inter_rfree, aes(x=factor(sero),y=im_sero,
                                 fill=med))+
  geom_tile()+
  scale_fill_viridis(name="Weeks Rabies-Free", option="D")+
  labs(x="Seroprevalence", y="Immigrant Seroprevalence")+
  theme_bw(base_size=12)+
  theme(panel.grid = element_blank())

# ggsave(filename = "./full_Figs/nweekfree_imsero_interac.jpeg",
#        width = 6, height = 4, dpi= 600, units = "in")

# Seroprevalence boxplot
rfree_sero <- time_rabies_free %>%
  filter(barrier==3 & type=="propagule") 

rfree_aov <- aov(data=rfree_sero, n_rabies_free~factor(sero))
rfreeHSD <- HSD.test(rfree_aov, "factor(sero)")$group %>%
  mutate(sero=rownames(.)) %>%
  select(-n_rabies_free)

rfree_sero <- rfree_sero %>%
  mutate(sero=factor(sero)) %>%
  left_join(rfreeHSD, by = "sero")

ggplot(data = rfree_sero, aes(x = factor(sero),y=n_rabies_free))+
  geom_boxplot(fill="lightgray")+
  geom_text(aes(label=groups, y=525))+
  labs(x="Seroprevalence", y="Weeks with 0 Cases")+
  theme_bw(base_size=12)+
  theme(panel.grid = element_blank())

# ggsave(filename = "./full_Figs/nweekfree_box.jpeg",
#        width = 6, height = 4, dpi= 600, units = "in")

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
      group_by(rep, sero, type, rate, im_sero, barrier) %>%
      # rate=disease rate of immigrants
      filter(n_infected == 0 & n_symptomatic == 0)
      mutate(nweek = ((year-1)*52) + week) %>%
      filter(nweek == min(nweek))
    
    first_elim <- data.frame(elim=time_to_elim$nweek,
                             rep=time_to_elim$rep)
    
    reinf_frame <- testfile %>%
      left_join(first_elim, by = "rep") %>%
      mutate(im_sero = imm_sero) %>%
      # rate=disease rate of immigrants
      mutate(nweek = ((year-1)*52) + week) %>%
      group_by(rep,sero, type, rate, im_sero, barrier) %>%
      filter(nweek > elim) %>%
      filter(n_infected > 0 | n_symptomatic > 0) %>%
      summarise(reinf_length = n()) %>%
      filter(reinf_length > 0) %>%
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

ggplot(data=reinf_outs, aes(x=factor(sero),y=reinf_prob))+
  geom_boxplot()

# ggsave(filename = "./full_Figs/nweekfree_barrier_interac.jpeg",
#        width = 6, height = 4, dpi= 600, units = "in")





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
