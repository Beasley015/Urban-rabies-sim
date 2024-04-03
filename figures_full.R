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
# dir <- "./outs" # juveniles reproduce
dir <- "./outs2" # juveniles do not reproduce

# Function for calculating time to first elimination --------
first_elim <- function(){
  # Get names of files
  filenames <- list.files(path = dir, pattern = "*.csv")

  # Progress bar
  pb = progressBar(min = 1, max = length(filenames), initial = 1,
                   style = "ETA") 
  
  for(i in 1:length(filenames)){
    # Read 'em in
    testfile <- read.csv(paste(getwd(),"/outs/",filenames[i],sep = ""))

    # Calculate time in weeks & get time to first elimination
    time_to_elim <- testfile %>%
      group_by(rep, sero, rate, im_sero, barrier) %>%
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

# Time to elimination: figs ------------------
# is the effect of seroprevalence similar to prelim model?
# filter out 1 barrier value then plot
elim_sansbar <- first_elim_full %>%
  # filter(barrier == 0) %>%
  mutate(sero = factor(sero))

sero_nweek <- aov(elim_sansbar, formula=nweek~sero)
sero_nweek_HSD <- HSD.test(sero_nweek, 
                           trt = "sero")$groups %>%
  mutate(sero = rownames(.)) %>%
  select(-nweek)

elim_sansbar <- elim_sansbar %>%
  left_join(sero_nweek_HSD, by = "sero")
  
ggplot(data = elim_sansbar, aes(x=factor(sero),y=nweek))+
  geom_boxplot(fill='lightgray')+
  geom_text(aes(y = 550, label = groups))+
  labs(x = "Seroprevalence", y = "Time to Elimination (Weeks)")+
  theme_bw()+
  theme(panel.grid=element_blank())

# ggsave(filename = "./full_Figs/sero_elim_bar0.jpeg", width=6,
#        height=4, dpi=600, units="in")

bar_interac <- first_elim_full %>%
  group_by(sero, barrier) %>%
  summarise(center = median(nweek))

# dev.new(width = 80, height = 60, unit = "mm", res=600,
#         noRStudioGD=TRUE)

ggplot(data = bar_interac, aes(x = factor(sero), y = barrier,
                               fill=center))+
  geom_tile()+
  scale_fill_viridis(name = "Weeks to Elimination",
                       option="B")+
  labs(x = "Seroprevalence", y = "Barrier Strength")+
  theme_bw(base_size=16)+
  theme(panel.grid = element_blank())

# ggsave(filename = "./full_Figs/sero_barrier_interac.jpeg",
#        width = 6, height = 4, dpi= 600, units = "in")

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

# Prob of elimination at time t --------------------
first_elim_prob <- function(){
  # Get names of files
  filenames <- list.files(path = "./outs", pattern = "*.csv")
  
  # Progress bar
  pb = progressBar(min = 1, max = length(filenames), initial = 1,
                   style = "ETA") 
  
  for(i in 1:length(filenames)){
    # Read 'em in
    testfile <- read.csv(paste(getwd(),"/outs/",filenames[i],
                               sep = ""))
    
    colnames=logical(length=10)
    for(yr in 1:10){colnames[yr]=paste("yr",yr,sep="")}

    # Calculate time in weeks & get time to first elimination
    time_to_elim <- testfile %>%
      group_by(rep, sero, rate, im_sero, barrier) %>%
      # rate=disease rate of immigrants
      filter(n_infected == 0 & n_symptomatic == 0) %>%
      mutate(nweek = ((year-1)*52) + week) %>%
      filter(nweek == min(nweek)) %>%
      ungroup()
    
    yrvec <- do.call(cbind,imap(colnames, 
                    .f=function(x,i){df = data.frame(x=rep(i*52,
                              length(unique(time_to_elim$rep))))}))
    
    time_to_elim <- cbind(time_to_elim, yrvec)
    colnames(time_to_elim)[13:22] <- colnames
    
    time_to_elim <- time_to_elim %>%
      mutate(across(yr1:yr10, 
                    .fns = function(x) length(which(nweek<x)),
                    .names = "{.col}_elim")) %>%
      select(-(yr1:yr10)) %>%
      select(-(rep:week)) %>%
      select(-(total_pop:actual_sero),-nweek) %>%
      distinct() %>%
      mutate(across(yr1_elim:yr10_elim, function(x) x/50))
    
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
  pivot_longer(cols = yr1_elim:yr10_elim, names_to = "year", 
               values_to = "prob") %>%
  mutate(year = as.numeric(str_extract(year, "(\\d)+"))) %>%
  group_by(year, sero) %>%
  summarise(mean_prob = mean(prob))

# dev.new(width = 80, height = 60, unit = "mm", res = 600,
#         noRStudioGD=TRUE)

ggplot(data=elim_prob_smol, aes(x = year, y = mean_prob,
                                color = factor(sero)))+
  geom_line(size = 1)+
  geom_hline(yintercept = 0.95, linetype="dashed")+
  scale_color_viridis_d(end = 0.9, name = "Seroprevalence")+
  lims(x = c(1,8))+
  labs(x = "Year", y = "Mean Elimination Probability")+
  theme_bw(base_size=16)+
  theme(panel.grid = element_blank())

# ggsave(filename = "./full_Figs/sero_elim_meant.jpeg", width=8,
#        dpi=600, units="in", height = 6)

ggplot(data=elim_prob_smol, aes(x = year, y = mean_prob,
                                color = factor(sero)))+
  geom_line(size = 1)+
  geom_hline(yintercept = 0.95, linetype="dashed")+
  scale_color_viridis_d(end = 0.9, name = "Seroprevalence")+
  lims(x = c(3,8), y = c(0.8, 1))+
  labs(x = "Year", y = "Mean Elimination Probability")+
  theme_bw(base_size=12)+
  theme(panel.grid = element_blank())

# ggsave(filename = "./full_Figs/sero_elim_meant_zoom.jpeg", width=6,
#        height=4, dpi=600, units="in")

# Look at interactions: barrier
elim_prob_bar <- elim_prob_t %>%
  pivot_longer(cols = yr1_elim:yr10_elim, names_to = "year", 
               values_to = "prob") %>%
  mutate(year = as.numeric(str_extract(year, "(\\d)+"))) %>%
  group_by(year, sero, barrier) %>%
  summarise(mean_prob = mean(prob))

ggplot(data = elim_prob_bar[elim_prob_bar$year==6,], 
       aes(x=factor(barrier), y=factor(sero), fill=mean_prob))+
  geom_tile()
# nothing interesting, really

# Look at interactions: immigrant infection
elim_prob_rate <- elim_prob_t %>%
  pivot_longer(cols = yr1_elim:yr10_elim, names_to = "year", 
               values_to = "prob") %>%
  mutate(year = as.numeric(str_extract(year, "(\\d)+"))) %>%
  group_by(year, sero, rate) %>%
  summarise(mean_prob = mean(prob))

ggplot(data = elim_prob_rate[elim_prob_rate$year==6,], 
       aes(x=factor(rate), y=factor(sero), fill=mean_prob))+
  geom_tile()
# Nothing here either

# Try immigrant seroprevalence rates
elim_prob_sero <- elim_prob_t %>%
  pivot_longer(cols = yr1_elim:yr10_elim, names_to = "year", 
               values_to = "prob") %>%
  mutate(year = as.numeric(str_extract(year, "(\\d)+"))) %>%
  group_by(year, sero, im_sero) %>%
  summarise(mean_prob = mean(prob))

ggplot(data = elim_prob_sero[elim_prob_sero$year==6,], 
       aes(x=im_sero, y=factor(sero), fill=mean_prob))+
  geom_tile()

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

    # Calculate time in weeks & get time to first elimination
    time_rabies_free <- testfile %>%
      group_by(rep, sero, rate, im_sero, barrier) %>%
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
  geom_text(aes(label=groups, y=525))+
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
  filenames <- list.files(path = "./outs", pattern = "*.csv")
  
  # Progress bar
  pb = progressBar(min = 1, max = length(filenames), initial = 1,
                   style = "ETA") 
  
  for(i in 1:length(filenames)){
    # Read 'em in
    testfile <- read.csv(paste(getwd(),"/outs/",filenames[i],
                               sep = ""))
    
    # Calculate time in weeks & get time to first elimination
    time_to_elim <- testfile %>%
      group_by(rep, sero, rate, im_sero, barrier) %>%
      # rate=disease rate of immigrants
      filter(n_infected == 0 & n_symptomatic == 0) %>%
      mutate(nweek = ((year-1)*52) + week) %>%
      filter(nweek == min(nweek))
    
    first_elim <- data.frame(elim=time_to_elim$nweek,
                             rep=time_to_elim$rep)
    
    reinf_frame <- testfile %>%
      left_join(first_elim, by = "rep") %>%
      # rate=disease rate of immigrants
      mutate(nweek = ((year-1)*52) + week) %>%
      group_by(rep,sero, rate, im_sero, barrier) %>%
      filter(nweek > elim) %>%
      filter(n_infected > 0 | n_symptomatic > 0) %>%
      summarise(reinf_length = n()) %>%
      filter(reinf_length > 1) %>%
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
rinf_aov <- aov(data=reinf_outs, reinf_prob~factor(sero))
rinfHSD <- HSD.test(rinf_aov, "factor(sero)")$group %>%
  mutate(sero=rownames(.)) %>%
  select(-reinf_prob)

reinf_outs <- reinf_outs %>%
  mutate(sero=factor(sero)) %>%
  left_join(rinfHSD, by = "sero")

ggplot(data=reinf_outs, aes(x=factor(sero),y=reinf_prob))+
  geom_boxplot(fill='lightgray')+
  geom_text(aes(y=0.2, label = groups))+
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

# ggsave(filename = "./full_Figs/rinfprob_bar.jpeg",
#        width = 6, height = 4, dpi= 600, units = "in")

# Immigrant vax rate
imsero_inter_rinf <- reinf_outs %>%
  mutate(sero=factor(sero))%>%
  group_by(im_sero,sero)%>%
  summarise(med = median(reinf_prob))

ggplot(data=imsero_inter_rinf, aes(x=sero,y=im_sero,
                                fill=med))+
  geom_tile()+
  scale_fill_viridis(name="Recolonization Probability",
                     option = "B")+
  labs(x="Seroprevalence", y="Immigrant Vaccination Rate")+
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

# ggsave(filename = "./full_Figs/rinfprob_imdis.jpeg",
#        width = 6, height = 4, dpi= 600, units = "in")

# Reinfection length ------------------
# Sero only
rinf_aov <- aov(data=reinf_outs, reinf_length~factor(sero))
rinfHSD <- HSD.test(rinf_aov, "factor(sero)")$group %>%
  mutate(sero=rownames(.)) %>%
  select(-reinf_length)

reinf_outs <- reinf_outs %>%
  mutate(sero=factor(sero)) %>%
  left_join(rinfHSD, by = "sero")

ggplot(data=reinf_outs, aes(x=factor(sero),y=reinf_length))+
  geom_boxplot(fill='lightgray')+
  # geom_text(aes(y=200, label = groups.y))+
  labs(x = "Seroprevalence", y = "Reinfection Length (Weeks)")+
  theme_bw(base_size=12)+
  theme(panel.grid = element_blank())

# ggsave(filename = "./full_Figs/rinflength_sero.jpeg",
#        width = 6, height = 4, dpi= 600, units = "in")

# Add quantile regression?

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

# Immigrant vax rate
vax_inter_rinf <- reinf_outs %>%
  mutate(sero=factor(sero))%>%
  group_by(sero,im_sero)%>%
  summarise(med = median(reinf_length))

ggplot(data=vax_inter_rinf, aes(x=sero,y=im_sero,
                                fill=med))+
  geom_tile()+
  scale_fill_viridis(name="Reinfection Length (Weeks)",
                     option = "B")+
  labs(x="Seroprevalence", y="Immigrant Vaccination Rate")+
  theme_bw(base_size=12)+
  theme(panel.grid = element_blank())

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

    cases_frame <- testfile %>%
      # rate=disease rate of immigrants
      mutate(nweek = ((year-1)*52) + week) %>%
      group_by(rep,sero, rate, im_sero, barrier) %>%
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

# Immigrant vax rate
vax_inter_cases <- total_cases %>%
  mutate(sero=factor(sero))%>%
  group_by(sero,im_sero) %>%
  summarise(med = median(ncases))

ggplot(data=vax_inter_cases, aes(x=sero,y=im_sero,
                                  fill=med))+
  geom_tile()+
  scale_fill_viridis(name="Total Cases",
                     option = "B")+
  labs(x="Seroprevalence", y="Immigrant Vaccination Rate")+
  theme_bw(base_size=12)+
  theme(panel.grid = element_blank())

# Total cases after first elimination ------------------
reinf_cases <- function(){
  # Get names of files
  filenames <- list.files(path = "./outs", pattern = "*.csv")
  
  # Progress bar
  pb = progressBar(min = 1, max = length(filenames), initial = 1,
                   style = "ETA") 
  
  for(i in 1:length(filenames)){
    # Read 'em in
    testfile <- read.csv(paste(getwd(),"/outs/",filenames[i],
                               sep = ""))
    
    # Calculate time in weeks & get time to first elimination
    time_to_elim <- testfile %>%
      group_by(rep, sero, rate, im_sero, barrier) %>%
      # rate=disease rate of immigrants
      filter(n_infected == 0 & n_symptomatic == 0) %>%
      mutate(nweek = ((year-1)*52) + week) %>%
      filter(nweek == min(nweek))
    
    first_elim <- data.frame(elim=time_to_elim$nweek,
                             rep=time_to_elim$rep)
    
    reinf_frame <- testfile %>%
      left_join(first_elim, by = "rep") %>%
      # rate=disease rate of immigrants
      mutate(nweek = ((year-1)*52) + week) %>%
      group_by(rep,sero, rate, im_sero, barrier) %>%
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

# Immigrant vaccination
vax_inter_rcases <- reinf_case_frame %>%
  mutate(sero=factor(sero))%>%
  group_by(sero,im_sero) %>%
  summarise(med = median(ncases))

ggplot(data=vax_inter_rcases, aes(x=sero,y=im_sero,
                                  fill=med))+
  geom_tile()+
  scale_fill_viridis(name="Cases Post-Reinvasion",
                     option = "B")+
  labs(x="Seroprevalence", y="Immigrant Vaccination Rate")+
  theme_bw(base_size=12)+
  theme(panel.grid = element_blank())

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
    
    cases_frame <- testfile %>%
      # rate=disease rate of immigrants
      mutate(nweek = ((year-1)*52) + week) %>%
      group_by(sero, rate, im_sero, barrier,nweek)
    
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
  filter(barrier == 2)

dev.new(width = 80, height = 60, unit = "mm", res=600,
        noRStudioGD=TRUE)

ggplot(data=meancase_seros, aes(x=nweek, y=mean_cases, 
                                color = factor(sero)))+
  geom_line()+
  geom_vline(xintercept=c(20, 52+20, (52*2+20), 52*3+20), 
             linetype="dashed", size = 1)+
  scale_color_viridis_d(end = 0.9, name="Seroprevalence")+
  labs(x = "Week", y = "Mean Cases")+
  theme_bw(base_size=16)+
  theme(panel.grid = element_blank())

# ggsave(filename = "./full_Figs/meanwkcases_sero.jpeg",
#        width = 6, height = 4, dpi= 600, units = "in")

# At what time step does it hit the maximum?
meancase_max_t <- meancase_seros %>%
  group_by(sero, rate, im_sero) %>%
  select(-barrier) %>%
  filter(mean_cases == max(mean_cases))

ggplot(data=meancase_max_t, aes(x=factor(sero), y=nweek))+
  geom_boxplot(fill='lightgray')+
  labs(x="Seroprevalence", y="Week with Max Cases")+
  theme_bw(base_size=12)+
  theme(panel.grid = element_blank())

# ggsave(filename = "./full_Figs/TimeMaxCases.jpeg", width=6,
#        height=4, dpi=600, units="in")

# What about after the initial peak?
meancase_max_1yr <- meancase_seros %>%
  group_by(sero, rate, im_sero) %>%
  select(-barrier) %>%
  filter(nweek >= 72) %>%
  filter(mean_cases == max(mean_cases))

ggplot(data=meancase_max_1yr, aes(x=factor(sero), y=nweek))+
  geom_boxplot(fill='lightgray')+
  labs(x="Seroprevalence", y="Week with Max Cases")+
  theme_bw(base_size=12)+
  theme(panel.grid = element_blank())

# ggsave(filename = "./full_Figs/TimeMaxCases_1yr.jpeg", width=6,
#        height=4, dpi=600, units="in")

# Zoom in
ggplot(data=meancase_seros, aes(x=nweek, y=mean_cases, 
                                color = factor(sero)))+
  geom_line()+
  geom_vline(xintercept=72, linetype="dashed")+
  geom_vline(xintercept = 124, linetype = "dashed")+
  lims(x=c(52,500), y=c(0,100))+
  scale_color_viridis_d(end = 0.9, name="Seroprevalence")+
  labs(x = "Week", y = "Mean Cases")+
  theme_bw(base_size=12)+
  theme(panel.grid = element_blank())

# ggsave(filename = "./full_Figs/meanwkcases_sero1yr.jpeg",
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

    cases_frame <- testfile %>%
      # rate=disease rate of immigrants
      mutate(nweek = ((year-1)*52) + week) %>%
      group_by(sero, rate, im_sero, barrier, nweek) %>%
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
  filter(barrier == 2)

# Sero only
ggplot(data=prev.sero, aes(x=nweek, y=mean_prev, 
                               color = factor(sero)))+
  geom_line()+
  geom_vline(xintercept=20, linetype="dashed")+
  scale_color_viridis_d(end=0.9, name = "Seroprevalence")+
  labs(x = "Weeks", y = "Rabies Prevalence")+
  theme_bw(base_size = 12)+
  theme(panel.grid = element_blank())

# ggsave(filename = "./full_Figs/meanwkprev_sero.jpeg",
#        width = 6, height = 4, dpi= 600, units = "in")

ggplot(data=prev.sero, aes(x=nweek, y=mean_prev, 
                           color = factor(sero)))+
  geom_line()+
  geom_vline(xintercept=72, linetype="dashed")+
  geom_vline(xintercept=124, linetype="dashed")+
  lims(x=c(52,500), y=c(0,0.05))+
  scale_color_viridis_d(end=0.9, name = "Seroprevalence")+
  labs(x = "Weeks", y = "Rabies Prevalence")+
  theme_bw(base_size = 12)+
  theme(panel.grid = element_blank())

# ggsave(filename = "./full_Figs/meanwkprev_sero1yr.jpeg",
#        width = 6, height = 4, dpi= 600, units = "in")
