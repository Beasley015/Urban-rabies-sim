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
dir <- "./Outputs" 

# Function: Check immunity rates ------------------
immune <- function(){
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
    
    # Calculate mean immunity per week
    immun_rate <- testfile %>%
      mutate(nweek = (year-1)*52 + week) %>%
      select(rep,sero,type,disease,rate,nweek,actual_sero) %>%
      group_by(sero, type, disease, rate, nweek) %>%
      summarise(mean_vax_rate = mean(actual_sero))
    
    if(i==1){
      vax_frame <- immun_rate
    } else{
      vax_frame <- rbind(vax_frame, immun_rate)
    }
    
    setTxtProgressBar(pb,i)
  }
  return(vax_frame)
}

vax_frame <- immune()

ggplot(data = vax_frame, aes(x = nweek, y = mean_vax_rate,
                             color = factor(sero)))+
  geom_line()+
  labs(x = "Week", y = "Immunity Rate")+
  scale_color_viridis_d(end = 0.9, name = "Adult Vax Rate")+
  theme_bw(base_size=14)+
  theme(panel.grid=element_blank())

# ggsave("./full_Figs/vax_rates.jpeg", width = 7, height = 5,
#        units = "in")

# Function: Calculate proportion eliminated --------
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

prop_elim <- prop_elim() %>%
  mutate(type=case_when(type=="wave"~"Seasonal",
                        TRUE ~ "Continuous"))

# Figs: Proportion of outbreaks eliminated ---------
# Quick glance
summary(lm(data=prop_elim, prop~sero+rate+disease+type))

# dev.new(width = 80, height = 60, unit = "mm", res=600,
#         noRStudioGD=TRUE)

# All sims: vax & immigration rate
ggplot(data=prop_elim, aes(x=factor(sero), 
                                           y = prop,
                               color = factor(rate),
                               group = factor(rate)))+
  geom_point()+
  geom_smooth(method = 'lm', se = F)+
  #geom_boxplot(fill="lightgray")+
  scale_color_viridis_d(end = 0.9, name = "Weekly Immigrants")+
  geom_hline(yintercept=0.95, linetype="dashed") +
  labs(x = "Adult Vaccination Rate", 
       y = "Proportion Reaching Elimination")+
  theme_bw(base_size=16)+
  theme(panel.grid=element_blank())

# dev.off()

prop_means <- prop_elim %>%
  group_by(sero,rate) %>%
  summarise(mean.prop = mean(prop))

ggplot(data=prop_elim, aes(x=factor(sero), y = prop,
                           color = factor(disease),
                           group = factor(disease)))+
  geom_point()+
  geom_smooth(method = 'lm', se = F)+
  #geom_boxplot(fill="lightgray")+
  geom_hline(yintercept=0.95, linetype="dashed") +
  scale_color_viridis_d(end = 0.9, name = "Disease Rate")+
  labs(x = "Adult Vaccination Rate", 
       y = "Proportion Reaching Elimination")+
  theme_bw(base_size=14)+
  theme(panel.grid=element_blank())

# ggsave("./full_Figs/prop_elim_box.jpeg", width = 7, 
#        height = 5, units = "in")

ggplot(data = prop_elim, aes(x = factor(sero), y = prop,
                             fill = type))+
  geom_boxplot()+
  scale_fill_manual(values=c('limegreen', 'darkgray'),
                    name="Immigration Type")+
  labs(x = "Adult Vaccination Rate", 
       y = "Proportion Eliminated")+
  theme_bw(base_size=12)+
  theme(panel.grid = element_blank())

# ggsave("./full_Figs/prop_elim_box_type.jpeg", width = 7,
#        height = 5, units = "in")

# Immigration vars interaction
im_vars_elim <- prop_elim %>%
  group_by(rate, disease) %>%
  summarise(mean.prop = mean(prop))

ggplot(data = prop_elim, aes(x = factor(rate), 
                                y = prop,
                                fill = factor(disease)))+
  geom_boxplot()+
  scale_fill_viridis_d(name="Immigrant Prevalence")+
  labs(x = "Expected Weekly Immigrants", 
       y = "Proportion Eliminated")+
  theme_bw()+
  theme(panel.grid=element_blank())

# ggsave("./full_Figs/prop_elim_imvars.jpeg",
#        width = 7, height = 5, units = "in")

# Look at immigration type
elim.typesero <- prop_elim %>%
  group_by(type, sero) %>%
  summarise(mean.prop = mean(prop))
  
ggplot(data=prop_elim, aes(x=factor(sero), y=prop, 
                               fill=type))+
  geom_boxplot()+
  scale_fill_viridis_d(name = "Immigration Type", end = 0.9)+
  labs(x = "Adult Vaccination Rate", 
       y = "Proportion Eliminated")+
  theme_bw(base_size = 12)+
  theme(panel.grid=element_blank())
  
# ggsave("./full_Figs/prop_elim_heatmap_imtype.jpeg",
#        width = 7.5, height = 3, units = "in")

# Rate & type 
elim.typerate <- prop_elim %>%
  group_by(type, rate) %>%
  summarise(mean.prop = mean(prop))

ggplot(data=elim.typerate, aes(x=factor(rate), y=type, 
                               fill=mean.prop))+
  geom_tile()+
  scale_fill_viridis(name = "Proportion Eliminated")+
  labs(x = "Immigration Rate", y = "Immigration Type")+
  theme_bw(base_size = 12)+
  theme(panel.grid=element_blank())

# ggsave("./full_Figs/prop_elim_heatmap_imratetype.jpeg",
#        width = 7, height = 5, units = "in")

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
  mutate(nweek = nweek-52, years = nweek/52)

# Quick look:
summary(lm(data=elim_sansbar,
           years~as.numeric(sero)+rate+disease+type))

elim_sansbar %>%
  group_by(sero) %>%
  summarise(med = median(years))

quicklook <- aov(data=elim_sansbar, years~sero)
HSD.test(quicklook, trt = 'sero', console=T)

# dev.new(width = 80, height = 60, unit = "mm", res=600,
#         noRStudioGD=TRUE)
  
ggplot(data = elim_sansbar, aes(x=factor(sero),y=years,
                                fill = factor(rate)))+
  geom_boxplot()+#fill = 'lightgray')+
  # geom_jitter()+
  scale_fill_viridis_d(name='Weekly Immigrants')+
  labs(x = "Adult Vaccination Rate", 
       y = "Time to Elimination (Years)")+
  theme_bw(base_size=16)+
  theme(panel.grid=element_blank())

imm_rate_elim <- first_elim_full %>%
  group_by(sero, rate) %>%
  summarise(center = median(nweek))

ggplot(data = imm_rate_elim, aes(x = factor(sero), 
                                 y = factor(rate),fill=center))+
  geom_tile()+
  scale_fill_viridis(name = "Weeks to Elimination",
                       option="B")+
  labs(x = "Adult Vaccination Rate", y = "Expected Weekly Immigrants")+
  theme_bw(base_size=12)+
  theme(panel.grid = element_blank())

# ggsave(filename = "./full_Figs/time_elim_rate.jpeg",
#        width=7.5, height=4, dpi=600, units="in")

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
# not really

# Disease/immigration interactions
im_interac <- first_elim_full %>%
  mutate(disease=factor(disease, 
                        levels = c(0,0.015,0.03,0.045,0.6)),
         rate = factor(rate)) %>%
  group_by(sero, disease, rate) %>%
  summarise(center = median(nweek, na.rm=T))

interac_list <- list()

for(i in 1:length(c(0, 0.2, 0.4, 0.6))){
  samples <- c(0, 0.2, 0.4, 0.6)
  
  interac_list[[i]] <- ggplot(data=im_interac[im_interac$sero==samples[i],], 
        aes(x = disease, y = rate, fill=center))+
    geom_tile()+
    scale_fill_viridis_c(name = "Weeks to Elimination",
                         limits = c(100,300))+
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
      filter(time <= 10)
    
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
  # filter(rate %in% c(1,2), sero %in% c(0, 0.8)) %>%
  group_by(sero, nweek) %>%
  summarise(mean_prop = mean(prop))

# dev.new(width = 80, height = 60, unit = "mm", res = 600,
#         noRStudioGD=TRUE)

ggplot(data=elim_prob_vax, aes(x = nweek, y = mean_prop,
                                color = factor(sero)))+#,
                               # linetype = factor(rate)))+
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

# R0 calculation ------------------------
get_r0 <- function(estimate = "R0"){
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
    
    for(j in 1:length(unique(testfile$rep))){
      if(estimate == "R0"){
        onerep <- testfile %>%
          filter(rep == unique(testfile$rep)[j])
      
        first_elim <- filter(time_to_elim, 
                           rep == unique(testfile$rep)[j])
      
        start <- min(which(onerep$n_symptomatic>0))
      
        endpoint <- ifelse(nrow(first_elim)<1, 52*11, 
                          first_elim$nweek)
      
        tests <- try(est.R0.SB(epid = onerep$n_symptomatic[start:endpoint],
                             GT=generation.time("gamma", c(4.5, 1))),
                   silent=T)
      
        if(class(tests) %in% 'try-error' == T) {next} else{
          r0.test <- est.R0.SB(epid = onerep$n_symptomatic[start:endpoint],
                            GT=generation.time("gamma", c(4.5, 1)))
          
          re.est <- median(r0.test$R)
        }
        
      } else if(estimate == "Re"){
        onerep <- testfile %>%
          filter(rep == unique(testfile$rep)[j])
        
        first_elim <- filter(time_to_elim, 
                             rep == unique(testfile$rep)[j])
        
        start <- min(which(onerep$n_symptomatic>0))
        
        endpoint <- ifelse(nrow(first_elim)<1, 52*11, 
                           first_elim$nweek)
        
        tests <- try(est.R0.SB(epid = onerep$n_symptomatic[start:endpoint],
                               GT=generation.time("gamma", c(4.5, 1))),
                     silent=T)
        
        if(class(tests) %in% 'try-error' == T) {next} else{
          r0.test <- est.R0.SB(epid = onerep$n_symptomatic[start:endpoint],
                               GT=generation.time("gamma", c(4.5, 1)))
        }
          
          re.est <- median(r0.test$R)*
            (1-onerep$actual_sero[start:endpoint])
          
          re.est <- mean(re.est)

        
      } else(print("Error: estimate must be R0 or Re"))
      
      row <- data.frame(rep=unique(onerep$rep), 
                          sero=unique(onerep$sero),
                          rate=unique(onerep$rate),
                          disease=unique(onerep$disease), 
                          type=unique(onerep$type), 
                          r.0=re.est)
          
        r0 <- bind_rows(r0, row)
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
re <- get_r0(estimate = "Re")

# R0 figures ---------------
r0.grp <- r0 %>%
  group_by(sero) %>%
  summarise(med = median(r.0), mean = mean(r.0))

ggplot(data = r0, aes(x = factor(sero), y = r.0))+
  geom_boxplot(fill = 'lightgray')+
  labs(x = "Adult Vaccination Rate", 
       y = bquote(R[0]))+
  theme_bw(base_size=14)+
  theme(panel.grid=element_blank())

mean(r0$mean) #ESTIMATED R0: 1.38

# ggsave(filename = "./full_Figs/full_re.jpeg", width=8,
#        dpi=600, units="in", height = 6)

# Re figures ------------
re.grp <- re %>%
  group_by(sero) %>%
  summarise(med = median(r.0), mean = mean(r.0))

# when vax = 0%: mean re = 1.14
# when vax = 80%: mean re = 0.51

ggplot(data = re, aes(x = factor(sero), y = r.0))+
  geom_boxplot(fill = 'lightgray')+
  labs(x = "Adult Vaccination Rate", 
       y = bquote(R[e]))+
  theme_bw(base_size=14)+
  theme(panel.grid=element_blank())

# ggsave(filename = "./full_figs/re_estimated.jpeg", height = 4,
#        width = 6, units = "in")

# Time between elimination & recolonization ---------
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

    # Get time to first elimination
    time_to_elim <- testfile %>%
      filter(elim == "True", year >= 2) %>%
      group_by(rep, sero, rate, type, disease) %>%
      mutate(nweek = ((year-1)*52) + week) %>%
      filter(nweek == min(nweek))
    
    # Get time step of first reinfection
    recols <- testfile %>%
      mutate(nweek = ((year-1)*52) + week) %>%
      filter(nweek > 53) %>%
      group_by(rep, sero, rate, type, disease) %>%
      filter(elim == "False" & lag(elim) == "True")
    
    elims <- testfile %>%
      mutate(nweek = ((year-1)*52) + week) %>%
      group_by(rep, sero, rate, type, disease) %>%
      filter(elim == "True" & lag(elim) == "False")
    
    recol.time <- rbind(elims, recols) %>%
      arrange(rep, nweek) %>%
      mutate(time = case_when(elim == "True" ~ nweek-lag(nweek),
                              TRUE ~ NA)) %>%
      filter(time >= 10 | is.na(time)==T)
    
    # Skip if no reinfections
    if(nrow(recol.time)==0){next}
    
    # Elim time - reinfection time
    time_rabies_free <- recol.time %>%
      mutate(time.no.rabies = 
               case_when(elim=="False" ~ nweek-lag(nweek),
                         TRUE ~ NA)) %>%
      filter(elim == "False") %>%
      group_by(rep) %>%
      filter(nweek == min(nweek))
    
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
# Quick look
summary(lm(data=time_rabies_free, 
           time.no.rabies~sero+rate+disease+type))

# Seroprevalence boxplot
ggplot(data = time_rabies_free, aes(x = factor(sero),
                                    y=time.no.rabies))+
  geom_boxplot(fill="lightgray")+
  # geom_text(aes(label=groups, y=525))+
  labs(x="Adult Vaccination Rate", y="Weeks Until Reinfection")+
  theme_bw(base_size=12)+
  theme(panel.grid = element_blank())

# ggsave(filename = "./full_Figs/nweekfree_box.jpeg",
#        width = 6, height = 4, dpi= 600, units = "in")

im.vars <- time_rabies_free %>%
  group_by(rate, disease) %>%
  summarise(median = median(time.no.rabies), 
            mean = mean(time.no.rabies))

ggplot(data=im.vars, aes(x = factor(rate), y = factor(disease),
                         fill = median))+
  geom_tile()+
  scale_fill_viridis(direction = -1, option = 'B',
                     name = "Weeks Until Recolonization")+
  labs(x = "Weekly Immigrants", y = "Immigrant Disease Rate")+
  theme_bw(base_size = 12)

# ggsave(filename = "./full_Figs/nweekfree_heat_imvars.jpeg",
#        width = 6, height = 4, dpi= 600, units = "in")

imtype <- time_rabies_free %>%
  mutate(type=case_when(type=="wave"~"Seasonal",
                          TRUE~"Continuous")) %>%
  group_by(rate, type) %>%
  summarise(median = median(time.no.rabies), 
            mean = mean(time.no.rabies))

ggplot(data=imtype, aes(x = factor(rate), y = type, 
                        fill = median))+
  geom_tile()+
  scale_fill_viridis(name = "Weeks Until Reinvasion")+
  labs(x = "Expected Weekly Immigrants", y = "Immigration Type")+
  theme_bw(base_size=12)

# ggsave(filename = "./full_Figs/nweekfree_heat_type.jpeg",
#        width = 6, height = 4, dpi= 600, units = "in")

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
    
    if(unique(testfile$disease == 0)){next}
    
    recols <- testfile %>%
      group_by(rep, sero, rate, type, disease) %>%
      filter(elim == "False" & lag(elim) == "True")
    
    elims <- testfile %>%
      group_by(rep, sero, rate, type, disease) %>%
      filter(elim == "True" & lag(elim) == "False")
    
    first_elim <- elims %>%
      filter(nweek == min(nweek)) %>%
      mutate(first = nweek) %>%
      select(rep,sero,disease,rate,type,first,
             actual_sero,total_pop) 
    
    recol.time <- rbind(elims, recols) %>%
      arrange(rep, nweek) %>%
      mutate(time = case_when(elim == "True" ~ nweek-lag(nweek),
                              TRUE ~ NA)) %>%
      filter(is.na(time) == F) %>%
      filter(time >= 52)
    
    reinf_frame <- suppressMessages(full_join(testfile, 
                                           recol.time)) %>%
      dplyr::filter(is.na(time)==F) %>%
      left_join(y = first_elim, 
                by = c('rep', 'sero', 'disease', 
                       'rate', 'type')) %>%
      dplyr::select(rep, sero, disease, rate, 
                    type, nweek, time, first, actual_sero.x,
                    total_pop.x) %>%
      mutate(prop = length(unique(.$rep))/
               length(unique(elims$rep))) %>%
      group_by(rep) %>%
      mutate(weeks_elim = nweek-lag(nweek)) %>%
      mutate(weeks_elim = 
               case_when(is.na(weeks_elim)==T ~ nweek-first,
                         TRUE ~ weeks_elim)) %>%
      mutate(weekly_prop = n()/sum(weeks_elim))

    if(exists("reinf_frame_full")==F){
      reinf_frame_full <- reinf_frame
    } else{
      reinf_frame_full <- rbind(reinf_frame_full, reinf_frame)
    }
    setTxtProgressBar(pb,i)
  }
  return(reinf_frame_full)
}

reinf_outs <- reinfection()

# Figures: reinfection prob -----------
reinf_outs <- reinf_outs %>%
  mutate(sero=factor(sero)) %>%
  mutate(InfStart = (nweek-time) %% 52) %>%
  mutate(TimePeriod = case_when(InfStart < 18 | 
                                  InfStart >= 43 ~ 
                                  "Post-dispersal",
                                between(InfStart,18,28) ~
                                  "Juveniles w/Mom",
                                between(InfStart,29,42) ~ 
                                  "Independent Juveniles"))

probs_condensed <- reinf_outs %>%
  ungroup() %>%
  mutate(type=case_when(type=="wave"~"Seasonal",
                        TRUE ~ "Continuous")) %>%
  dplyr::select(sero, type, disease, rate, prop) %>%
  distinct()

weekly_probs_condensed <- reinf_outs %>%
  dplyr::select(rep, sero, type, disease, rate, prop,
                weekly_prop) %>%
  distinct() %>%
  mutate(type=case_when(type=="wave"~"Seasonal",
                   TRUE ~ "Continuous"))

summary(lm(data=probs_condensed,
           prop~as.numeric(sero)+rate+disease+type))
summary(lm(data=weekly_probs_condensed,
           weekly_prop~as.numeric(sero)+rate+disease+type+
             as.numeric(sero)*rate))

reinf_imms <- ggplot(data=weekly_probs_condensed, 
                      aes(x=factor(rate), y=weekly_prop,
                          fill = factor(disease)))+
  geom_boxplot(outlier.shape=NA)+
  scale_fill_viridis_d(end=0.9, 
                       name = "Immigrant Disease Rate")+
  scale_y_continuous(limits = c(0, 0.016))+
  labs(x = "Expected Weekly Immigrants", 
       y = "Weekly Recolonization Probability")+
  theme_bw(base_size=12)+
  theme(panel.grid = element_blank())

# ggsave(filename = "./full_Figs/rinfprob_imms.jpeg",
#        width = 6, height = 4, dpi= 600, units = "in")

reinf_type <- ggplot(data=weekly_probs_condensed, 
                     aes(x=factor(rate), y = weekly_prop,
                         fill=type))+
  geom_boxplot(outlier.shape=NA)+
  scale_fill_viridis_d(end=0.9,
                    name = "Immigration Type")+
  scale_y_continuous(limits = c(0, 0.016))+
  labs(x= "Expected Weekly Immigrants", 
       y = "Weekly Recolonization Probability")+
  theme_bw(base_size=12)+
  theme(panel.grid=element_blank())

recol_vax <- ggplot(data=weekly_probs_condensed, 
       aes(x=factor(sero), y = weekly_prop,
           fill=factor(rate)))+
  geom_boxplot(outlier.shape=NA)+
  scale_fill_viridis_d(end=0.9,
                       name = "Weekly Immigrants")+
  scale_y_continuous(limits = c(0, 0.016))+
  labs(x= "Adult Vaccination Rate", 
       y = "Weekly Recolonization Probability")+
  theme_bw(base_size=12)+
  theme(panel.grid=element_blank())

dev.new(width = 160, height = 120, unit = "mm", res=600)

((reinf_imms | reinf_type)/recol_vax)+
  plot_annotation(tag_levels = "a")

# ggsave(filename = "./full_Figs/rinfprob_imms_outliers.jpeg",
#        width=12.5, height = 4.5, unit='in')

ggplot(data=probs_condensed, aes(x=factor(sero), y = prop,
                                 fill=factor(rate)))+
  geom_boxplot()+
  # scale_fill_manual(values = c('limegreen', 'lightgray'),
  #                   name = "Immigration Type")+
  scale_fill_viridis_d()+
  labs(x= "Expected Weekly Immigrants", 
       y = "Reinvasion Probability")+
  theme_bw(base_size=12)+
  theme(panel.grid=element_blank())
  
weekly_probs_med <- weekly_probs_condensed %>%
  group_by(sero, rate) %>%
  summarise(med = median(weekly_prop), mean = mean(weekly_prop))

ggplot(data=weekly_probs_med, aes(x=sero,
                                        y=factor(rate),
                                 fill = mean))+
  geom_tile()+
  scale_fill_viridis(option = "B", direction = -1,
                       name = "Weekly Recol. Probability")+
  labs(x = "Adult Vaccination Rate", 
       y = "Expected Weekly Immigrants")+
  theme_bw(base_size=12)+
  theme(panel.grid = element_blank())

# ggsave(filename = "./full_Figs/rinfprob_heat_ratevax.jpeg",
#        width = 7.5, height = 3.5, dpi= 600, units = "in")

# Immigration type
type_inter_rinf <- reinf_outs %>%
  mutate(sero=factor(sero))%>%
  group_by(sero,type)%>%
  summarise(med = median(weekly_prop))

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
  filter(type=="propagule") %>%
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

# ggsave(filename = "./full_figs/reinf_by_time.jpeg", 
#        height = 5, width = 8, units = "in", dpi = 600)

# Reinfection length ------------------
reinf_outs <- filter(reinf_outs, prop > 0)

# Test for associations between elimination time 
# and reinfection length
summary(lm(data = reinf_outs, time~(nweek-time)))

ggplot(data=reinf_outs, aes(x = factor(sero), y = time,
                            fill=factor(rate)))+
  geom_boxplot(outlier.shape = NA)

dev.new(width = 80, height = 60, unit = "mm", res=600,
        noRStudioGD=TRUE)

ggplot(data=reinf_outs, aes(x=TimePeriod,y=time))+
  geom_boxplot(fill = 'lightgray', outlier.shape = NA)+
  # geom_boxplot(fill = 'lightgray')+
  # geom_boxplot(aes(fill = factor(disease)), outlier.shape=NA)+
  # scale_fill_viridis_d(name = "Time Period", end = 0.9)+
  scale_y_continuous(limits = c(0,75))+
  labs(x = "Adult Vaccination Rate", 
       y = "Recolonization Length (Weeks)")+
  theme_bw(base_size=12)+
  theme(panel.grid = element_blank())

# ggsave(filename = "./full_Figs/rinflength.jpeg",
#        width = 6, height = 4, dpi= 600, units = "in")

# Quick look at vax rates
reinf_outs %>%
  select(rep, sero, disease, rate, type, time) %>%
  group_by(sero) %>%
  group_map(~mutate(., n_reinf = n()), .keep = TRUE) %>%
  bind_rows() %>%
  filter(time > 52) %>%
  group_by(sero) %>%
  group_map(~mutate(., n_long = n()), .keep = TRUE) %>%
  bind_rows() %>%
  group_by(sero) %>%
  mutate(perc_long = n_long/n_reinf) %>%
  mutate(max = max(time)) %>%
  select(sero, perc_long, max) %>%
  distinct()

rq025 <- rq(time ~ as.numeric(sero), data = reinf_outs, 
            tau = 0.025)
rq25 <- rq(time ~ as.numeric(sero), data = reinf_outs, tau = 0.25)
rq5 <- rq(time ~ as.numeric(sero), data = reinf_outs, tau = 0.5)
rq75 <- rq(time ~ as.numeric(sero), data = reinf_outs, tau = 0.75)
rq975 <- rq(time ~ as.numeric(sero), data = reinf_outs, 
            tau = 0.975)
rq1 <- rq(time ~ as.numeric(sero), data = reinf_outs, 
            tau = 0.999)

reinf_outs_quant <- reinf_outs %>%
  mutate(quantile = 
           case_when(time < quantile(time, 0.025) ~ "0",
                     time >= quantile(time,0.025) &
                       time < quantile(time,0.25) ~ "0.025",
                     time >= quantile(time,0.25) &
                       time < quantile(time, 0.55) ~ "0.25",
                     time >= quantile(time,0.5) &
                       time < quantile(time, 0.75) ~ "0.5",
                     time >= quantile(time, 0.75) &
                       time < quantile(time, 0.975) ~ "0.75",
                     TRUE ~ "1")) 

colors <- viridis(6, end = 0.9, direction = -1)

ggplot(data=reinf_outs_quant, aes(x=sero, y=time, 
                                  color = quantile))+
  geom_point(alpha=0.1)+
  scale_color_viridis_d(end = 0.9, direction = -1,
                        name = "Quantile")+
  guides(color=guide_legend(override.aes=list(alpha=1)))+
  # scale_color_manual(values = c('limegreen', 'black'))+
  geom_abline(intercept=coef(rq025)[1], slope=coef(rq025)[2],
              colour = colors[1], linewidth=1)+
  geom_abline(intercept=coef(rq25)[1], slope=coef(rq25)[2],
              colour = colors[2],linewidth=1)+
  geom_abline(intercept=coef(rq5)[1], slope=coef(rq5)[2],
              colour = colors[3],linewidth=1)+
  geom_abline(intercept=coef(rq75)[1], slope=coef(rq75)[2],
              colour = colors[4], linewidth=1)+
  geom_abline(intercept=coef(rq975)[1], slope=coef(rq975)[2],
              colour = colors[5],linewidth=1)+
  geom_abline(intercept=coef(rq1)[1], slope=coef(rq1)[2],
              colour = colors[6],linewidth=1)+
  labs(x = "Adult Vaccination Rate", y = "Reinfection Length")+
  theme_bw(base_size=12)+
  theme(panel.grid = element_blank())

# ggsave(filename = "./full_Figs/rinflength_vax.jpeg",
#        width = 6, height = 4, dpi= 600, units = "in")
  
# Disease rate
reinf_outs %>%
  group_by(disease) %>%
  summarise(mean = mean(time), median=median(time))

# Immigration rate
imm.vars <- reinf_outs %>%
  group_by(rate,sero) %>%
  summarise(mean = mean(time), median=median(time))

ggplot(data=imm.vars, aes(x=factor(sero),y=factor(rate),
                            fill=mean))+
  geom_tile()+
  scale_fill_viridis()

# ggsave(filename = "./full_Figs/rinflength_rate_vax.jpeg",
#               width = 6, height = 4, dpi= 600, units = "in")

# Immigration type
type_inter_rinf <- reinf_outs %>%
  mutate(sero=factor(sero))%>%
  group_by(sero,type)%>%
  summarise(med = median(time))

ggplot(data=type_inter_rinf, aes(x=sero,y=type,
                                 fill=med))+
  geom_tile()+
  scale_fill_viridis(name="Reinfection Length (Weeks)",
                     option = "B")+
  labs(x="Seroprevalence", y="Immigration Type")+
  theme_bw(base_size=12)+
  theme(panel.grid = element_blank())

# Reinf length ~ assorted vars -------------
# Length ~ time rabies free
ggplot(data=reinf_outs, aes(x = weeks_elim, y = time,
                            color = factor(rate)))+
  geom_point(alpha = 0.5)+
  geom_smooth(se = F)

# Length ~ pop size at outbreak start
ggplot(data=reinf_outs, aes(x = total_pop.x, y = time,
                            color = factor(rate)))+
  geom_point(alpha = 0.5)+
  geom_smooth(se = F)

# Length ~ vax rate at outbreak start
ggplot(data=reinf_outs, aes(x = actual_sero.x, y = time,
                            color = factor(type)))+
  geom_point(alpha = 0.5)+
  geom_smooth(se = F)

# Time to next birth pulse
pulses <- seq(20, 572, by = 52)

time_to_bp <- reinf_outs %>%
  ungroup() %>%
  filter(type == "propagule" & nweek < 540) %>%
  rowwise() %>%
  mutate(nxt_bp = min(pulses[which(pulses > nweek)])-nweek)

ggplot(data=time_to_bp, aes(x = nxt_bp, y = time,
                            color = factor(rate)))+
  geom_point()+
  geom_smooth(se = F)
  
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
reinf_case_frame.vax <- reinf_case_frame %>%
  mutate(sero=factor(sero), rate=factor(rate)) %>%
  group_by(sero, rate) %>%
  summarise(mean=mean(ncases), median=median(ncases))

ggplot(data=reinf_case_frame, aes(x=factor(sero),y=ncases,
                                  fill=factor(rate)))+
  geom_boxplot(outlier.shape=NA)+
  scale_fill_viridis_d(name="Weekly Immigrants")+
  labs(x = "Adult Vaccination Rate", 
       y = "Cases Post-Recolonization")+
  scale_y_continuous(limits = c(0,2000))+
  theme_bw(base_size=12)+
  theme(panel.grid = element_blank())

# ggsave(filename = "./full_Figs/reinfcases_box.jpeg",
#        width = 6, height = 4, dpi= 600, units = "in")

# Closer look at low imm rate
reinf_lowimm <- reinf_case_frame.vax %>%
  filter(rate %in% c(1,2,3))

ggplot(reinf_lowimm, aes(x = factor(sero), y = factor(rate),
                         fill=median))+
  geom_tile()+
  scale_fill_viridis(name = "Median Cases")+
  labs(x = "Adult Vaccination Rate", 
       y = "Expected Weekly Immigrants")+
  theme_bw(base_size=14)

# ggsave(filename = "./full_Figs/reinfcases_smolheat.jpeg",
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
  # filter(disease == 0.03) %>%
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
