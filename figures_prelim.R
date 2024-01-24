###########################################
# Figures for urban rabies simulation     # 
# E.M. Beasley                            #
# Fall 2023                               #
###########################################

# Load packages ---------------------
library(tidyverse)
library(viridis)
library(quantreg)

# Load data --------------------
setwd("c:/users/beasl/documents/urban-rabies-sim")
sero_outs <- read.csv("outputs_prelim.csv")
sero_extended <- read.csv("outputs_extended.csv")

# Clean data -----------------------------
# Add seroprevalence column
nrep <- 50
nyear <- 10
nweek <- 52

seros <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8)
seros_extended <- c(0.05, 0.1, 0.15, 0.55, 0.65, 0.75)

sero_outs$sero <- rep(seros, each = nrep*nyear*nweek)
sero_extended$sero <- rep(seros_extended, each = nrep*nyear*nweek)

sero_outs <- rbind(sero_outs, sero_extended)

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

# ggsave(filename = "outbreakduration.jpeg", width = 5, height = 4,
#        units = "in", dpi = 600)

summary(lm(data = sero_time, 
           formula = nweek~as.numeric(as.character(sero))))

summary(rq(data = sero_time, 
           formula = nweek~as.numeric(as.character(sero)),
           tau = c(0.05, 0.25, 0.5, 0.75, 0.95)))

rq.frame <- tibble(slopes = c(-20, -24.44, -38, -29.33, -140),
                      inters = c(75, 84.33, 105.8, 127, 248),
                      quants = c("5%", "25%", "50%", "75%", "95%"))%>%
  mutate(quants = factor(quants, levels = quants))

ggplot(data = sero_time, aes(x = as.numeric(as.character(sero)), 
                             y = nweek))+
  geom_point()+
  geom_abline(data = rq.frame, aes(intercept=inters, slope = slopes,
                                   color = quants))+
  labs(x = "Seroprevalence", y = "Outbreak Duration (Weeks)")+
  scale_color_viridis_d(end = 0.9, name = "Quantile")+
  theme_bw()+
  theme(panel.grid = element_blank())

# ggsave(filename = "quantileoubreakduration.jpeg", width = 5, 
#        height = 4, units = "in", dpi = 600)
  
# max weekly cases
infecs <- sero_outs %>%
  # mutate(sero = factor(sero)) %>%
  group_by(rep, sero) %>%
  filter(n_infected != 0 | n_symptomatic != 0) %>%
  group_by(sero, rep) %>%
  summarise(max = max(n_infected))

summary(lm(max~sero, data = infecs))

ggplot(data = infecs, aes(x = factor(sero), y = max))+
  geom_boxplot(fill = "lightgray")+
  labs(x = "Seroprevalence", y = "Maximum Weekly Infections")+
  theme_bw()+
  theme(panel.grid = element_blank())

# ggsave(filename = "maxcases.jpeg", width = 5, height = 4, 
#        units = 'in', dpi = 600)

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

birth_pulse <- c(20, (20+52), 20+(2*52), 20+(3*52), 20+(4*52),
                 20+(5*52), 20+(6*52), 20+(7*52), 20+(8*52))

juvie_dispersal <- c(50+(2*52), 50+(3*52))

ggplot(sero_death_time, aes(x = nweek, y = mean, color=sero))+
  geom_line()+
  geom_vline(xintercept = birth_pulse, linetype = "dashed",
             alpha = 0.5)+
  geom_vline(xintercept = juvie_dispersal, color = "red",
             linetype = "dashed", alpha = 0.5)+
  scale_color_viridis_d(name = "Seroprevalence")+
  labs(x = "Week", y = "Mean Deaths")+
  theme_bw()+
  theme(panel.grid = element_blank())

# ggsave(filename = "deathsperweek_juvie_dispersal.jpeg", width = 5, 
#        height = 4, units = 'in', dpi = 600)

