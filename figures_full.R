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
dir <- "./outputs" 

# Function: Calculate proportion of outbreaks eliminated ----------
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
ggplot(data=prop_elim, aes(x=factor(sero), y = prop))+
  geom_boxplot(fill="lightgray")+
  geom_hline(yintercept=0.95, linetype="dashed") +
  labs(x = "Seroprevalence", 
       y = "Proportion Reaching Elimination")+
  theme_bw(base_size=14)+
  theme(panel.grid=element_blank())

# Interesting...

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
      mutate(nweek = ((year-2)*52) + week) %>%
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
  mutate(sero = factor(sero))

sero_nweek <- aov(elim_sansbar, formula=nweek~sero)
sero_nweek_HSD <- HSD.test(sero_nweek, 
                           trt = "sero")$groups %>%
  mutate(sero = rownames(.)) %>%
  select(-nweek)

elim_sansbar <- elim_sansbar %>%
  left_join(sero_nweek_HSD, by = "sero")
  
ggplot(data = elim_sansbar, aes(x=factor(sero),y=nweek))+
  geom_violin(fill='lightgray')+
  # geom_text(aes(y = 350, label = groups))+
  labs(x = "Seroprevalence", y = "Time to Elimination (Weeks)")+
  theme_bw()+
  theme(panel.grid=element_blank())

# ggsave(filename = "./full_Figs/sero_elim_full.jpeg", width=6,
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

# interaction: proportion diseased immigrants
rate_interac <- first_elim_full %>%
  mutate(disease = factor(disease), sero=factor(sero)) %>%
  group_by(sero, disease) %>%
  summarise(center = median(nweek))

ggplot(data = rate_interac, aes(x = sero, y = disease, 
                                fill = center))+
  geom_tile()+
  scale_fill_viridis(name = "Weeks to Elimination", option="B")+
  labs(x ="Vaccination Rate", y="Immigrant Infection Rate")+
  theme_bw(base_size=12)+
  theme(panel.grid=element_blank())

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

# any interaction between immigrant prevalence & type?
im_interac <- first_elim_full %>%
  mutate(disease=factor(disease)) %>%
  group_by(sero, disease, type) %>%
  summarise(center = median(nweek))

ggplot(data=im_interac[im_interac$sero==0.2,], 
       aes(x = disease, y = type, fill=center))+
  geom_tile()+
  scale_fill_viridis_c()

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
                               filenames[i],sep = ""))
    
    colnames=logical(length=10)
    for(yr in 1:11){colnames[yr]=paste("yr",yr,sep="")}

    # Calculate time in weeks & get time to first elimination
    time_to_elim <- testfile %>%
      group_by(rep, sero, type, rate, disease) %>%
      # rate=disease rate of immigrants
      filter(year >= 2) %>%
      filter(elim == "True") %>%
      mutate(nweek = ((year-1)*52) + week) %>%
      filter(nweek == min(nweek)) %>%
      ungroup()
    
    yrvec <- do.call(cbind,imap(colnames, 
                    .f=function(x,i){df = data.frame(x=rep(i*52,
                              length(unique(time_to_elim$rep))))}))
    
    time_to_elim <- cbind(time_to_elim, yrvec)
    colnames(time_to_elim)[14:24] <- colnames
    
    time_to_elim <- time_to_elim %>%
      mutate(across(yr1:yr11, 
                    .fns = function(x) length(which(nweek<x)),
                    .names = "{.col}_elim")) %>%
      select(-(yr1:yr11)) %>%
      select(-(rep:week)) %>%
      select(-(total_pop:actual_sero),-nweek) %>%
      distinct() %>%
      mutate(across(yr1_elim:yr11_elim, function(x) x/50))
    
    if(i==1){
      first_elim_frame <- time_to_elim
    } else{
      first_elim_frame <- rbind(first_elim_frame, time_to_elim)
    }
    
    setTxtProgressBar(pb,i)
  }
  return(first_elim_frame)
}

elim_prob_t <- first_elim_prob()

# Prob elim w/in time figs ------------
elim_prob_smol <- elim_prob_t %>%
  pivot_longer(cols = yr2_elim:yr11_elim, names_to = "year", 
               values_to = "prob") %>%
  mutate(year = as.numeric(str_extract(year, "(\\d)+"))) %>%
  group_by(year, sero) %>%
  mutate(year = year-1) %>%
  summarise(mean_prob = mean(prob))

# dev.new(width = 80, height = 60, unit = "mm", res = 600,
#         noRStudioGD=TRUE)

ggplot(data=elim_prob_smol, aes(x = year, y = mean_prob,
                                color = factor(sero)))+
  geom_line(size = 1)+
  geom_hline(yintercept = 0.95, linetype="dashed")+
  scale_color_viridis_d(end = 0.9, name = "Seroprevalence")+
  lims(x = c(1, 10))+
  labs(x = "Year", y = "Mean Elimination Probability")+
  theme_bw(base_size=16)+
  theme(panel.grid = element_blank())

# ggsave(filename = "./full_Figs/alt_sero_elim_meant.jpeg", width=8,
#        dpi=600, units="in", height = 6)

# Look at interactions: barrier
elim_prob_bar <- elim_prob_t %>%
  pivot_longer(cols = yr4_elim:yr14_elim, names_to = "year", 
               values_to = "prob") %>%
  mutate(year = as.numeric(str_extract(year, "(\\d)+"))) %>%
  mutate(year=year-4) %>%
  group_by(year, sero, barrier) %>%
  summarise(mean_prob = mean(prob))

ggplot(data = elim_prob_bar[elim_prob_bar$year==2,], 
       aes(x=factor(barrier), y=factor(sero), fill=mean_prob))+
  geom_tile()+
  scale_fill_viridis_c()

# Look at interactions: immigrant infection
elim_prob_rate <- elim_prob_t %>%
  pivot_longer(cols = yr4_elim:yr14_elim, names_to = "year", 
               values_to = "prob") %>%
  mutate(year = as.numeric(str_extract(year, "(\\d)+"))) %>%
  mutate(year=year-4) %>%
  group_by(year, sero, rate) %>%
  summarise(mean_prob = mean(prob))

ggplot(data = elim_prob_rate[elim_prob_rate$year==9,], 
       aes(x=factor(rate), y=factor(sero), fill=mean_prob))+
  geom_tile()+
  scale_fill_viridis_c()

# Try immigration type
elim_prob_sero <- elim_prob_t %>%
  pivot_longer(cols = yr4_elim:yr14_elim, names_to = "year", 
               values_to = "prob") %>%
  mutate(year = as.numeric(str_extract(year, "(\\d)+"))) %>%
  mutate(year=year-4) %>%
  group_by(year, sero, type) %>%
  summarise(mean_prob = mean(prob))

ggplot(data = elim_prob_sero[elim_prob_sero$year==5,], 
       aes(x=type, y=factor(sero), fill=mean_prob))+
  geom_tile()+
  scale_fill_viridis_c()

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
    
    # pull immigration type from file name
    type <- str_match(filenames[i], "type(.*?)\\.")[2]

    # Calculate time in weeks & get time to first elimination
    time_rabies_free <- testfile %>%
      mutate(type=type) %>%
      filter(year >= 4) %>%
      group_by(rep, sero, rate, type, barrier) %>%
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
  mutate(barrier=factor(barrier, levels = 0:4), sero=factor(sero))%>%
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

type_inter_rfree <- time_rabies_free %>%
  group_by(sero,type)%>%
  summarise(med = median(n_rabies_free))

ggplot(data=type_inter_rfree, aes(x=factor(sero),y=type,
                                 fill=med))+
  geom_tile()+
  scale_fill_viridis(name="Weeks Rabies-Free", option="D")+
  labs(x="Seroprevalence", y="Immigration Type")+
  theme_bw(base_size=12)+
  theme(panel.grid = element_blank())

# ggsave(filename = "./full_Figs/nweekfree_imsero_interac.jpeg",
#        width = 6, height = 4, dpi= 600, units = "in")

# Seroprevalence boxplot
rfree_sero <- time_rabies_free %>%
  filter(barrier==0) 

rfree_aov <- aov(data=rfree_sero, n_rabies_free~factor(sero))
rfreeHSD <- HSD.test(rfree_aov, "factor(sero)")$group %>%
  mutate(sero=rownames(.)) %>%
  select(-n_rabies_free)

rfree_sero <- rfree_sero %>%
  mutate(sero=factor(sero)) %>%
  left_join(rfreeHSD, by = "sero")

ggplot(data = rfree_sero, aes(x = factor(sero),y=n_rabies_free))+
  geom_boxplot(fill="lightgray")+
  # geom_text(aes(label=groups, y=525))+
  labs(x="Seroprevalence", y="Weeks with 0 Cases")+
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
    # Read 'em in
    testfile <- read.csv(paste(getwd(), 
                               str_replace(dir, ".", ""), "/",
                               filenames[i],sep = ""))
    
    # pull immigration type from file name
    type <- str_match(filenames[i], "type(.*?)\\.")[2]
    
    # Calculate time in weeks & get time to first elimination
    time_to_elim <- testfile %>%
      mutate(type=type) %>%
      filter(year >= 4) %>%
      group_by(rep, sero, rate, type, barrier) %>%
      # rate=disease rate of immigrants
      filter(n_infected == 0 & n_symptomatic == 0) %>%
      mutate(nweek = ((year-1)*52) + week) %>%
      filter(nweek == min(nweek))
    
    first_elim <- data.frame(elim=time_to_elim$nweek,
                             rep=time_to_elim$rep)
    
    reinf_frame <- testfile %>%
      mutate(type=type) %>%
      left_join(first_elim, by = "rep") %>%
      filter(year >= 4) %>%
      # rate=disease rate of immigrants
      mutate(nweek = ((year-1)*52) + week) %>%
      group_by(rep, barrier, sero, rate, type) %>%
      filter(nweek > elim) %>%
      filter(n_infected > 1) %>%
      mutate(reinf_length = n()) %>%
      filter(nweek==max(nweek)) %>%
      ungroup() %>%
      group_by(rep, sero, barrier, rate, type) %>%
      select(-c(n_infected, n_symptomatic,total_pop,actual_sero,
                year)) %>%
      filter(reinf_length == max(reinf_length)) %>%
      filter(reinf_length > 1) %>%
      ungroup() %>%
      mutate(reinf_prob = length(unique(rep))/nrow(first_elim),
             time_to_reinf = nweek-elim)
    
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
  mutate(elim=elim-(3*52+1)) %>%
  mutate(sero=factor(sero)) %>%
  mutate(TimePeriod = case_when(between(week,19,41) ~
                                  "Juveniles w/Mom",
         TRUE ~ "Juveniles Independent"))
  # mutate(TimePeriod = case_when(week < 20 ~ "PreBirth",
  #                               between(week,19,41) ~
  #                                 "PreIndependence",
  #                               week > 40 ~ "Juveniles Dispersed"))

ggplot(data=reinf_outs, aes(x=factor(sero),y=reinf_prob,
                            fill=TimePeriod))+
  geom_boxplot()+
  # scale_fill_manual(values = c("lightgray", "limegreen"))+
  scale_fill_viridis_d(end=0.9)+
  labs(x = "Seroprevalence", y = "Recolonization Probability")+
  theme_bw(base_size=12)+
  theme(panel.grid = element_blank())

# ggsave(filename = "./full_Figs/rinfprob_sero.jpeg",
#        width = 6, height = 4, dpi= 600, units = "in")

# Barrier:
bar_inter_rinf <- reinf_outs %>%
  mutate(barrier=factor(barrier), sero=factor(sero))%>%
  group_by(sero,barrier)%>%
  summarise(med = median(reinf_prob))

ggplot(data=bar_inter_rinf, aes(x=sero,y=barrier,
                                 fill=med))+
  geom_tile()+
  scale_fill_viridis(name="Recolonization Probability",
                     option = "B")+
  labs(x="Seroprevalence", y="Barrier Strength")+
  theme_bw(base_size=12)+
  theme(panel.grid = element_blank())

# Immigrant infection rate
rate_inter_rinf <- reinf_outs %>%
  mutate(rate=factor(rate), sero=factor(sero))%>%
  group_by(sero,rate)%>%
  summarise(med = median(reinf_prob))

ggplot(data=rate_inter_rinf, aes(x=sero,y=rate,
                                fill=med))+
  geom_tile()+
  scale_fill_viridis(name="Recolonization Probability",
                     option = "B")+
  labs(x="Seroprevalence", y="Immigrant Disease Rate")+
  theme_bw(base_size=12)+
  theme(panel.grid = element_blank())

type_inter_rinf <- reinf_outs %>%
  mutate(sero=factor(sero))%>%
  group_by(sero,type)%>%
  summarise(med = median(reinf_prob))

ggplot(data=type_inter_rinf, aes(x=sero,y=type,
                                 fill=med))+
  geom_tile()+
  scale_fill_viridis(name="Recolonization Probability",
                     option = "B")+
  labs(x="Seroprevalence", y="Immigration Type")+
  theme_bw(base_size=12)+
  theme(panel.grid = element_blank())

# Reinfection length ------------------
ggplot(data=reinf_outs, aes(x=factor(sero),y=reinf_length,
                            fill = TimePeriod))+
  geom_boxplot()+
  scale_fill_manual(values = c("lightgray", "limegreen"))+
  labs(x = "Seroprevalence", y = "Reinfection Length (Weeks)")+
  theme_bw(base_size=12)+
  theme(panel.grid = element_blank())

# ggsave(filename = "./full_Figs/rinflength_sero.jpeg",
#        width = 6, height = 4, dpi= 600, units = "in")

summary(aov(data=reinf_outs, formula = reinf_length ~ TimePeriod + sero))

# Barrier
bar_inter_rlen <- reinf_outs %>%
  mutate(barrier=factor(barrier), sero=factor(sero))%>%
  group_by(sero,barrier)%>%
  summarise(med = median(reinf_length))

ggplot(data=bar_inter_rlen, aes(x=sero,y=barrier,
                                fill=med))+
  geom_tile()+
  scale_fill_viridis(name="Reinfection Length (Weeks)",
                     option = "B")+
  labs(x="Seroprevalence", y="Barrier Strength")+
  theme_bw(base_size=12)+
  theme(panel.grid = element_blank())

# Immigrant disease rate
rate_inter_rinf <- reinf_outs %>%
  mutate(rate=factor(rate), sero=factor(sero))%>%
  group_by(sero,rate)%>%
  summarise(med = median(reinf_length))

ggplot(data=rate_inter_rinf, aes(x=sero,y=rate,
                                fill=med))+
  geom_tile()+
  scale_fill_viridis(name="Reinfection Length (Weeks)",
                     option = "B")+
  labs(x="Seroprevalence", y="Immigrant Disease Rate")+
  theme_bw(base_size=12)+
  theme(panel.grid = element_blank())

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

# Weeks until reinfection (con't here)----------------------
ggplot(data=reinf_outs, aes(x=factor(sero),y=time_to_reinf))+
  geom_boxplot(fill='lightgray')+
  labs(x = "Seroprevalence", y = "Time to Reinfection (Weeks)")+
  theme_bw(base_size=12)+
  theme(panel.grid = element_blank())

# ggsave(filename = "./full_Figs/time_to_reinf.jpeg",
#        width = 6, height = 4, dpi= 600, units = "in")

ggplot(data = reinf_outs, aes(x=elim, y = time_to_reinf))+
  geom_point(aes(color=factor(sero)))+
  scale_color_viridis_d()+
  labs(x="Time of Elimination (Weeks)", 
       y="Time to Reinfection (Weeks)")+
  theme_bw()+
  theme(panel.grid=element_blank())

# ggsave(filename = "./full_Figs/time_by_elim.jpeg",
#        width = 6, height = 4, dpi= 600, units = "in")

# Reinfection by time -----------------

  
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
    
    # pull immigration type from file name
    type <- str_match(filenames[i], "type(.*?)\\.")[2]

    cases_frame <- testfile %>%
      mutate(type=type) %>%
      # rate=disease rate of immigrants
      mutate(nweek = ((year-1)*52) + week) %>%
      group_by(rep,sero, rate, type, barrier) %>%
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
cases_aov <- aov(data=total_cases, ncases~factor(sero))
casesHSD <- HSD.test(cases_aov, "factor(sero)")$group %>%
  mutate(sero=rownames(.)) %>%
  select(-ncases)

total_cases <- total_cases %>%
  mutate(sero=factor(sero)) %>%
  left_join(casesHSD, by = "sero")

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

# Barrier
bar_inter_cases <- total_cases %>%
  mutate(barrier=factor(barrier), sero=factor(sero))%>%
  group_by(sero,barrier) %>%
  summarise(med = median(ncases))

ggplot(data=bar_inter_cases, aes(x=sero,y=barrier,
                                fill=med))+
  geom_tile()+
  scale_fill_viridis(name="Total Cases",
                     option = "B")+
  labs(x="Seroprevalence", y="Barrier Strength")+
  theme_bw(base_size=12)+
  theme(panel.grid = element_blank())

# Immigrant infections
rate_inter_cases <- total_cases %>%
  mutate(rate=factor(rate), sero=factor(sero))%>%
  group_by(sero,rate) %>%
  summarise(med = median(ncases))

ggplot(data=rate_inter_cases, aes(x=sero,y=rate,
                                 fill=med))+
  geom_tile()+
  scale_fill_viridis(name="Total Cases",
                     option = "B")+
  labs(x="Seroprevalence", y="Immigrant Infection Rate")+
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
  
  # pull immigration type from file name
  type <- str_match(filenames[i], "type(.*?)\\.")[2]
  
  for(i in 1:length(filenames)){
    # Read 'em in
    testfile <- read.csv(paste(getwd(), 
                               str_replace(dir, ".", ""), "/",
                               filenames[i],sep = ""))
    
    # Calculate time in weeks & get time to first elimination
    time_to_elim <- testfile %>%
      group_by(rep, sero, rate, barrier) %>%
      filter(year >= 4) %>%
      # rate=disease rate of immigrants
      filter(n_infected == 0 & n_symptomatic == 0) %>%
      mutate(nweek = ((year-1)*52) + week) %>%
      filter(nweek == min(nweek))
    
    first_elim <- data.frame(elim=time_to_elim$nweek,
                             rep=time_to_elim$rep)
    
    reinf_frame <- testfile %>%
      mutate(type=type) %>%
      left_join(first_elim, by = "rep") %>%
      # rate=disease rate of immigrants
      mutate(nweek = ((year-1)*52) + week) %>%
      group_by(rep,sero, rate, type, barrier) %>%
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
# sero
rcases_aov <- aov(data=reinf_case_frame, ncases~factor(sero))
rcasesHSD <- HSD.test(rcases_aov, "factor(sero)")$group %>%
  mutate(sero=rownames(.)) %>%
  select(-ncases)

reinf_case_frame <- reinf_case_frame %>%
  mutate(sero=factor(sero)) %>%
  left_join(rcasesHSD, by = "sero")

ggplot(data=reinf_case_frame, aes(x=factor(sero),y=ncases))+
  geom_boxplot(fill='lightgray')+
  # geom_text(aes(y=1300, label = groups))+
  labs(x = "Seroprevalence", y = "Cases Post-Reinvasion")+
  theme_bw(base_size=12)+
  theme(panel.grid = element_blank())

# ggsave(filename = "./full_Figs/reinfcases_sero.jpeg",
#        width = 6, height = 4, dpi= 600, units = "in")

# Barrier
bar_inter_rcases <- reinf_case_frame %>%
  mutate(barrier=factor(barrier), sero=factor(sero))%>%
  group_by(sero,barrier) %>%
  summarise(med = median(ncases))

ggplot(data=bar_inter_rcases, aes(x=sero,y=barrier,
                                 fill=med))+
  geom_tile()+
  scale_fill_viridis(name="Cases Post-Reinvasion",
                     option = "B")+
  labs(x="Seroprevalence", y="Barrier Strength")+
  theme_bw(base_size=12)+
  theme(panel.grid = element_blank())

# Immigrant infections
rate_inter_rcases <- reinf_case_frame %>%
  mutate(rate=factor(rate), sero=factor(sero))%>%
  group_by(sero,rate) %>%
  summarise(med = median(ncases))

ggplot(data=rate_inter_rcases, aes(x=sero,y=rate,
                                  fill=med))+
  geom_tile()+
  scale_fill_viridis(name="Cases Post-Reinvasion",
                     option = "B")+
  labs(x="Seroprevalence", y="Immigrant Infection Rate")+
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
    
    # pull immigration type from file name
    type <- str_match(filenames[i], "type(.*?)\\.")[2]
    
    cases_frame <- testfile %>%
      mutate(type=type) %>%
      mutate(nweek = ((year-1)*52) + week) %>%
      group_by(sero, rate, type, barrier,nweek)
    
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

# filter out 1 barrier value
meancase_seros <- meancase %>%
  filter(nweek >= (52*3)) %>%
  filter(barrier == 2)

# dev.new(width = 80, height = 60, unit = "mm", res=600,
#         noRStudioGD=TRUE)

ggplot(data=meancase_seros, aes(x=nweek, y=mean_cases, 
                                color = factor(sero)))+
  geom_line()+
  geom_vline(xintercept=c(52*3+20, 52*4+20, 52*5+20, 52*6+20),
             linetype="dashed", size = 1)+
  scale_color_viridis_d(end = 0.9, name="Seroprevalence")+
  labs(x = "Week", y = "Mean Cases")+
  theme_bw(base_size=16)+
  theme(panel.grid = element_blank())

# ggsave(filename = "./full_Figs/meanwkcases_sero.jpeg",
#        width = 6, height = 4, dpi= 600, units = "in")

# Maximum cases per week ------------------------
maxcase <- cases_per_week(metric="max")

# filter out 1 barrier value
maxcase_seros <- maxcase %>%
  filter(barrier == 2)

ggplot(data=maxcase_seros, aes(x=nweek, y=max_cases, 
                                color = factor(sero)))+
  geom_line()+
  geom_vline(xintercept=20, linetype="dashed")+
  scale_color_viridis_d(end=0.9, name = "Seroprevalence")+
  labs(x = "Weeks", y = "Maximum Cases")+
  theme_bw(base_size = 12)+
  theme(panel.grid = element_blank())

# ggsave(filename = "./full_Figs/maxwkcases_sero.jpeg",
#        width = 6, height = 4, dpi= 600, units = "in")

ggplot(data=maxcase_seros, aes(x=nweek, y=max_cases, 
                               color = factor(sero)))+
  geom_line()+
  lims(x = c(52, 500), y= c(0,75))+
  geom_vline(xintercept=72, linetype="dashed")+
  geom_vline(xintercept=124, linetype="dashed")+
  scale_color_viridis_d(end=0.9, name = "Seroprevalence")+
  labs(x = "Weeks", y = "Maximum Cases")+
  theme_bw(base_size = 12)+
  theme(panel.grid = element_blank())

# ggsave(filename = "./full_Figs/maxwkcases_sero1yr.jpeg",
#        width = 6, height = 4, dpi= 600, units = "in")

# Avg Prevalence -----------------------
prev <- function(metric){
  # Get names of files
  filenames <- list.files(path = "./outs", pattern = "*.csv")
  
  # Progress bar
  pb = progressBar(min = 1, max = length(filenames), initial = 1,
                   style = "ETA") 
  
  for(i in 1:length(filenames)){
    # Read 'em in
    testfile <- read.csv(paste(getwd(),"/outs/",filenames[i],
                               sep = ""))
    
    # pull immigration type from file name
    type <- str_match(filenames[i], "type(.*?)\\.")[2]

    cases_frame <- testfile %>%
      mutate(type=type) %>%
      mutate(nweek = ((year-1)*52) + week) %>%
      group_by(sero, rate, type, barrier, nweek) %>%
      summarise(prev = n_infected/total_pop) 
    
    if(metric=="mean"){
      cases_frame <- cases_frame %>%
        summarise(mean_prev = mean(prev))
    }else if(metric == "max"){
      cases_frame <- cases_frame %>%
        summarise(max_prev = max(prev))
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

prev.mean <- prev(metric = "mean")
prev.max <- prev(metric = "max")

# Prevalence figs -------------------
prev.sero <- prev.mean %>%
  filter(barrier == 2) %>%
  filter(nweek >= 3*52)

# Sero only
ggplot(data=prev.sero, aes(x=nweek, y=mean_prev, 
                               color = factor(sero)))+
  geom_line()+
  scale_color_viridis_d(end=0.9, name = "Seroprevalence")+
  labs(x = "Weeks", y = "Rabies Prevalence")+
  theme_bw(base_size = 12)+
  theme(panel.grid = element_blank())

# ggsave(filename = "./full_Figs/meanwkprev_sero.jpeg",
#        width = 6, height = 4, dpi= 600, units = "in")

# Weekly population -----------------------------
get_weekly_pop <- function(){
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
    
    # pull immigration type from file name
    type <- str_match(filenames[i], "type(.*?)\\.")[2]
    
    pop_frame <- testfile %>%
      mutate(type=type) %>%
      mutate(nweek = ((year-1)*52) + week) %>%
      group_by(sero, rate, type, barrier,nweek) %>%
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
  filter(barrier == 2)

# Sero only, full sim
ggplot(data=pop.sero, aes(x=nweek, y=mean_pop, 
                           color = factor(sero)))+
  geom_line()+
  geom_vline(xintercept=52*3+1, linetype="dashed")+
  scale_color_viridis_d(end=0.9, name = "Seroprevalence")+
  labs(x = "Weeks", y = "Population Size")+
  theme_bw(base_size = 12)+
  theme(panel.grid = element_blank())

# ggsave(filename = "./full_Figs/meantotalpop.jpeg",
#        width = 6, height = 4, dpi= 600, units = "in")
