###########################################
# Figures for urban rabies simulation     # 
# E.M. Beasley                            #
# Fall 2023                               #
###########################################

# Load packages ---------------------
library(tidyverse)
library(viridis)

# Load data --------------------
setwd("c:/users/beasl/documents/urban-rabies-sim")
sero_outs <- read.csv("outputs.csv")

# Clean data -----------------------------
# Add seroprevalence column
nrep = 50
nyear = 10
nweek = 52

seros = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8)

sero_outs$sero = rep(seros, each = nrep*nyear*nweek)

sero_time <- sero_outs %>%
  mutate(sero = factor(sero)) %>%
  group_by(rep, sero) %>%
  filter(n_infected == 0 & n_symptomatic == 0) %>%
  mutate(nweek = ((year-1)*52) + week) %>%
  filter(nweek == min(nweek))

ggplot(data = sero_time, aes(x = sero, y = nweek))+
  geom_boxplot(fill = "lightgray")+
  labs(x = "Seroprevalence", y = "Outbreak Duration (Weeks)")+
  theme_bw()+
  theme(panel.grid = element_blank())

ggsave(filename = "outbreakduration.jpeg", width = 5, height = 4,
       units = "in", dpi = 600)

test <- aov(data=sero_time, nweek~sero)
TukeyHSD(test)

# % of reps with outbreaks > 100 weeks
long_outbreaks <- sero_time %>%
  filter(nweek >= 100) %>%
  group_by(sero) %>%
  summarise(count = n())

ggplot(data = long_outbreaks, aes(x = sero, y = count))+
  geom_col(fill = "lightgray", color = "black")+
  labs(x = "Seroprevalence", y = "Outbreaks > 100 Weeks")+
  theme_bw()+
  theme(panel.grid = element_blank())

ggsave(filename = "longoutbreaks.jpeg", width = 5, height = 4,
       units = 'in', dpi = 600)

# max weekly cases
infecs <- sero_outs %>%
  mutate(sero = factor(sero)) %>%
  group_by(rep, sero) %>%
  filter(n_infected != 0 | n_symptomatic != 0) %>%
  group_by(sero, rep) %>%
  summarise(max = max(n_infected))

ggplot(data = infecs, aes(x = sero, y = max))+
  geom_boxplot(fill = "lightgray")+
  labs(x = "Seroprevalence", y = "Maximum Weekly Infections")+
  theme_bw()+
  theme(panel.grid = element_blank())

ggsave(filename = "maxcases.jpeg", width = 5, height = 4, 
       units = 'in', dpi = 600)

# mean cases per week
infecs_time <- sero_outs %>%
  mutate(sero = factor(sero)) %>%
  group_by(rep, sero) %>%
  filter(n_infected != 0 | n_symptomatic != 0) %>%
  mutate(nweek = ((year-1)*52) + week) %>%
  group_by(sero, nweek) %>%
  summarise(mean = mean(n_infected))

ggplot(data = infecs_time, aes(x = nweek, y = mean, color = sero))+
  geom_line()+
  scale_color_viridis_d(name = "Seroprevalence")+
  labs(x = "Week", y = "Mean Weekly Infections")+
  theme_bw()+
  theme(panel.grid = element_blank())

ggsave(filename = "casesperweek.jpeg", width = 5, height = 4, 
       units = 'in', dpi = 600)

# max mortality
sero_smol <- sero_outs %>%
  mutate(sero = factor(sero)) %>%
  group_by(rep, sero) %>%
  filter(n_infected != 0 | n_symptomatic != 0)

sero_smol$deaths <- NA

get_mortality <- function(dat){
  for(i in 1:nrow(dat)){
    if(dat$year[i] == 1 & dat$week[i] == 1){
      dat$deaths[i] <- NA
    }
    else{
      dat$deaths[i] <- dat$total_pop[i-1]-dat$total_pop[i]
    }
  }
  return(dat)
}

sero_smol <- get_mortality(sero_smol)

sero_death <- sero_smol %>%
  group_by(rep, sero) %>%
  summarise(maxdeath = max(deaths, na.rm = T), 
            meddeath = median(deaths, na.rm = T))

ggplot(sero_death, aes(x = sero, y = maxdeath))+
  geom_boxplot(fill = "lightgray")+
  labs(x = "Seroprevalence", y = "Maximum weekly deaths")+
  theme_bw()+
  theme(panel.grid = element_blank())

ggsave(filename = "maxdeaths.jpeg", width = 5, height = 4,
       units = 'in', dpi = 600)

sero_death_time <- sero_smol %>%
  filter(n_infected != 0 | n_symptomatic != 0) %>%
  filter(deaths > 0) %>%
  mutate(nweek = ((year-1)*52) + week) %>%
  group_by(sero, nweek) %>%
  summarise(mean = mean(deaths, na.rm = T))

ggplot(sero_death_time, aes(x = nweek, y = mean, color=sero))+
  geom_line()+
  scale_color_viridis_d(name = "Seroprevalence")+
  labs(x = "Week", y = "Mean Deaths")+
  theme_bw()+
  theme(panel.grid = element_blank())

ggsave(filename = "deathsperweek.jpeg", width = 5, height = 4,
       units = 'in', dpi = 600)
