library(tidyverse)
library(viridis)
library(DescTools)

setwd("./ParamSensitivity")

# Population size ---------------------
pop <- read.csv("kmax10.csv") %>%
  select(rep, year, week, total_pop, a_mort, j_mort) %>%
  mutate(nweek = ((year-1)*52)+week)

ggplot(data=pop, aes(x = nweek, y = total_pop, 
                     color = factor(rep)))+
  geom_line()+
  scale_color_viridis_d(name = "Rep", end = 0.9)+
  facet_grid(rows = vars(a_mort), cols=vars(j_mort))+
  labs(x = "Week", y = "Population Size", title = "Max K = 40")+
  theme_bw(base_size = 12)+
  theme(panel.grid = element_blank())

# ggsave("maxcc_10.jpeg", width = 7, height = 4, units = "in")

# Compare July and October population sizes ------------
seasonal <- pop %>%
  filter(week %in% c(28:31, 41:44)) %>%
  mutate(season = case_when(week %in% c(28:31) == T ~ "Summer",
                            TRUE ~ "Fall"))

seasonal.summary <- seasonal %>%
  group_by(a_mort, j_mort, season) %>%
  summarise(mean = mean(total_pop), min = min(total_pop),
            max = max(total_pop))

ggplot(data = seasonal, aes(x = season, y = total_pop))+
  geom_boxplot(fill = 'lightgray')+
  facet_grid(rows = vars(a_mort), cols = vars(j_mort))+
  labs(x = "Season", y = "Carrying Capacity", title = "Max K = 40")+
  theme_bw()+
  theme(panel.grid = element_blank())

# ggsave("seasonal_pop_10.jpeg", width = 5, height = 4, units = "in")

# Transmission Rates -----------------------
# Note: pop immunity set at 0, immigration rate low, no imm disease
dis <- read.csv("disease_test.csv") %>%
  select(rep, year, week, total_pop, n_infected, n_symptomatic, elim,
         l1, l2) %>%
  mutate(nweek = ((year-1)*52)+week)

# Compare proportion of outbreaks eliminated
prop_eliminated <- dis %>%
  filter(year >= 5, elim == "True") %>%
  select(rep, l1, l2) %>%
  group_by(l1, l2) %>%
  distinct() %>%
  summarise(prop = n()/5)

unique(prop_eliminated$prop) #ok not bad

# Some combos had 0% elimination, add to data
combos <- expand_grid(prop_eliminated$l1, prop_eliminated$l2)
colnames(combos) <- c("l1", "l2")

prop_eliminated <- prop_eliminated %>%
  right_join(combos, by = c("l1", "l2")) %>%
  distinct() %>%
  mutate(prop = case_when(is.na(prop) == T ~ 0,
                          TRUE ~ prop))

# Proportion eliminated: figs
ggplot(prop_eliminated, aes(x = factor(l1), 
                            y = factor(l2),
                            fill = factor(prop)))+
  geom_tile()+
  scale_fill_viridis_d(name = "Proportion eliminated")+
  labs(x = "Within-cell transmission", y = "Home range transmission")+
  theme_bw()

# ggsave(filename = "propelimheat.jpeg", width = 5, height = 4,
#        units = "in")

# Weeks to elimination
time_to_elim <- dis %>%
  group_by(l1, l2) %>%
  filter(year >= 5, elim == "True") %>%
  filter(nweek == min(nweek)) %>%
  right_join(combos, by = c("l1", "l2")) %>%
  distinct() %>%
  mutate(l1 = factor(l1), 
         l2 = factor(l2))

ggplot(data = time_to_elim, aes(x = l1, y = nweek))+
  geom_boxplot(fill = "lightgray")+
  labs(x = "Within-cell transmission", y = "Week")+
  theme_bw()+
  theme(panel.grid = element_blank())

# ggsave(filename = "celltransmissionbox.jpeg", width = 5, height = 4,
#        units = "in")

ggplot(data = time_to_elim, aes(x = l2, y = nweek))+
  geom_boxplot(fill = "lightgray")+
  labs(x = "Home Range Transmission", y = "Week")+
  theme_bw()+
  theme(panel.grid = element_blank())

# ggsave(filename = "hrtransmissionbox.jpeg", width = 5, height = 4,
#        units = "in")

ggplot(data = time_to_elim, aes(x = l1, y = l2,
                                fill = nweek))+
  geom_tile()+
  scale_fill_viridis_c(name = "Week")+
  labs(x = "Within-cell transmission", y = "Home range transmission")+
  theme_bw()

# ggsave(filename = "transmissionheatmap.jpeg", width = 5, height = 4,
#        units = "in")

# Lower home range transmission = outbreaks last longer
# Do sims with high transmission just burn out?

# Check with cases per week
ggplot(data = dis, aes(x = nweek, y = n_symptomatic, 
                       color = factor(rep)))+
  geom_line()+
  facet_grid(rows = vars(l1), cols = vars(l2))+
  scale_color_viridis_d(end = 0.9, name = "Rep")+
  xlim(c((52*4)+1, 800))+
  theme_bw()+
  theme(panel.grid=element_blank())

# ggsave(filename = "facet_disease.jpeg", width = 10, height = 9,
#        units = "in")

# looks like home range transmission influences initial peak
# Check mean & max number of weekly cases

mean_cases <- dis %>%
  filter(nweek >=250) %>%
  group_by(rep, homerange, cell) %>%
  summarise(mean.cases = mean(n_symptomatic)) %>%
  mutate(cell = factor(cell), 
         homerange = factor(homerange))

ggplot(data=mean_cases, aes(x = cell, y = homerange,
                           fill = mean.cases))+
  geom_tile()+
  scale_fill_viridis(name = "Mean Weekly Cases")+
  labs(x = "Within-cell transmission", y = "Home range transmission")+
  theme_bw()

# ggsave(filename = "meancases_heatmap.jpeg", width = 5, height = 4,
#        units = "in")

# Zoom in on best parameter combos
smol03 <- dis %>%
  filter(cell == 0.03 & homerange == 0.001)

smol025 <- dis %>%
  filter(cell == 0.025 & homerange == 0.001)

ggplot(data=smol025, aes(x = nweek, y = n_symptomatic, 
                        color = factor(rep)))+
  geom_line()+
  scale_color_viridis_d(end = 0.9, name = "Rep")+
  xlim(c((52*4)+1, 800))+
  labs(x = "Week", y = "Symptomatic Cases")+
  theme_bw(base_size = 14)+
  theme(panel.grid = element_blank())

# ggsave(filename = "cell025_weekly.jpeg", width = 6, height = 4, 
#        units = "in")

# Demonstrate that a smaller peak can lead to more realistic cases
ggplot(data = filter(smol025, rep == 5), 
       aes(x = nweek, y = n_symptomatic))+
  geom_line()+
  xlim(c((52*4)+1, 800))+
  labs(x = "Week", y = "Symptomatic Cases")+
  theme_bw(base_size = 14)+
  theme(panel.grid = element_blank())

# ggsave(filename = "small_peak.jpeg", width = 6, height = 4,
#        units = "in")

ggplot(data=mean_cases, aes(x = cell, y = mean.cases))+
  geom_boxplot(fill = "lightgray")+
  labs(x = "Within-cell transmission", y = "Mean Weekly Cases")+
  theme_bw()+
  theme(panel.grid = element_blank())

# ggsave(filename = "meancases_cell_box.jpeg", width = 5, height = 4,
#        units = "in")

ggplot(data=mean_cases, aes(x = homerange, y = mean.cases))+
  geom_boxplot(fill = "lightgray")+
  labs(x = "Home range transmission", y = "Mean Weekly Cases")+
  theme_bw()+
  theme(panel.grid = element_blank())

# ggsave(filename = "meancases_hr_box.jpeg", width = 5, height = 4,
#        units = "in")

max_cases <- dis %>%
  group_by(rep, homerange, cell) %>%
  summarise(max.cases = max(n_symptomatic)) %>%
  mutate(cell = factor(cell), 
         homerange = factor(homerange))

ggplot(data=max_cases, aes(x = cell, y = max.cases))+
  geom_boxplot(fill = "lightgray")+
  labs(x = "Within-cell transmission", y = "Max Weekly Cases")+
  theme_bw()+
  theme(panel.grid = element_blank())

# ggsave(filename = "maxcases_cell_box.jpeg", width = 5, height = 4,
#        units = "in")

ggplot(data=max_cases, aes(x = homerange, y = max.cases))+
  geom_boxplot(fill = "lightgray")+
  labs(x = "Home range transmission", y = "Max Weekly Cases")+
  theme_bw()+
  theme(panel.grid = element_blank())

# ggsave(filename = "maxcases_hr_box.jpeg", width = 5, height = 4,
#        units = "in")

ggplot(data=max_cases, aes(x = cell, y = homerange,
                           fill = max.cases))+
  geom_tile()+
  scale_fill_viridis(name = "Max Weekly Cases")+
  labs(x = "Within-cell transmission", y = "Home range transmission")+
  theme_bw()

# ggsave(filename = "maxcases_heatmap.jpeg", width = 5, height = 4,
#        units = "in")

dis_pop <- dis %>%
  group_by(cell, homerange, nweek) %>%
  summarise(mean_pop = mean(total_pop))

ggplot(data = dis_pop, aes (x = nweek, y = mean_pop, 
                        color = factor(cell)))+
  geom_line()+
  scale_color_viridis_d(name = "Within-cell transmission",
                        end = 0.9)+
  facet_grid(rows = vars(homerange))+
  theme_bw()+
  theme(panel.grid.minor = element_blank())

# ggsave(filename = "pop_direct.jpeg", width = 8, height = 6,
#        units = "in")

ggplot(data = dis_pop, aes (x = nweek, y = mean_pop, 
                            color = factor(homerange)))+
  geom_line()+
  scale_color_viridis_d(name = "Home range transmission",
                        end = 0.9)+
  facet_grid(rows = vars(cell))+
  theme_bw()+
  theme(panel.grid.minor = element_blank())

# ggsave(filename = "pop_indirect.jpeg", width = 8, height = 6,
#        units = "in")

# Messing with # of starting cases --------------
start <- read.csv("disease_low_initialization.csv") %>%
  select(rep, year, week, total_pop, n_infected, n_symptomatic, elim,
         cell, homerange) %>%
  mutate(nweek = ((year-1)*52)+week)

# Proportion of outbreaks eliminated
prop_elim <- start %>%
  filter(year >= 5, elim == "True") %>%
  select(rep, cell, homerange) %>%
  group_by(cell, homerange) %>%
  distinct() %>%
  summarise(prop = n()/5)

unique(prop_elim$prop) #still good

# Weekly cases
ggplot(data = start, aes(x = nweek, y = n_symptomatic,
                             color = factor(rep)))+
  geom_line()+
  facet_grid(cols = vars(cell))+
  scale_color_viridis_d(end = 0.9, name = "Rep")+
  xlim(c((52*4)+1), 800)+
  theme_bw()+
  theme(panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank())

# ggsave("weekly_cases_start15.jpeg", width = 6, height = 4,
#        units = "in")

# Low within-cell prob ------------
dis <- read.csv("disease02.csv") %>%
  select(rep, year, week, total_pop, n_infected, n_symptomatic, elim) %>%
  mutate(nweek = ((year-1)*52)+week)

# Compare proportion of outbreaks eliminated
dis %>%
  filter(year >= 5, elim == "True") %>%
  select(rep) %>%
  distinct() %>%
  summarise(prop = n()/5)
#ok none were eliminated

# Check weekly # of cases
ggplot(data = dis, aes(x = nweek, y = n_symptomatic, 
                       color = factor(rep)))+
  geom_line()+
  lims(x = c((52*4)+1, 800))+
  scale_color_viridis_d(end = 0.9, name = "Rep")+
  theme_bw()+
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())

# ggsave(filename = "cellprob02.jpeg", width = 6, height = 4,
#        units = "in")
