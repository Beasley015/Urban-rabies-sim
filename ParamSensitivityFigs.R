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

# ggsave("seasonal_pop_10.jpeg", width = 5, height = 4, 
#        units = "in")

# Transmission Rates -----------------------
# Note: pop immunity set at 0, immigration rate low, no imm disease
dis <- read.csv("disease_test.csv") %>%
  select(rep, year, week, total_pop, n_infected, n_symptomatic, elim,
         l1, l2) %>%
  mutate(nweek = ((year-1)*52)+week)

# Compare proportion of outbreaks eliminated
prop_eliminated <- dis %>%
  filter(year >= 2, elim == "True") %>%
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
  labs(x = "Within-cell transmission", y = "Between-cell transmission")+
  theme_bw()

# ggsave(filename = "propelimheat.jpeg", width = 5, height = 4,
#        units = "in")

# Weeks to elimination
time_to_elim <- dis %>%
  group_by(l1, l2) %>%
  filter(year >= 2, elim == "True") %>%
  filter(nweek == min(nweek)) %>%
  # right_join(combos, by = c("l1", "l2")) %>%
  distinct() %>%
  mutate(l1 = factor(l1), 
         l2 = factor(l2))

ggplot(data = time_to_elim, aes(x = l1, y = nweek))+
  geom_boxplot(fill = "lightgray")+
  labs(x = "Within-cell transmission", y = "Week")+
  theme_bw()+
  theme(panel.grid = element_blank())

# ggsave(filename = "celltransmissionbox.jpeg", width = 5,
#        height = 4, units = "in")

ggplot(data = time_to_elim, aes(x = l2, y = nweek))+
  geom_boxplot(fill = "lightgray")+
  labs(x = "Home Range Transmission", y = "Week")+
  theme_bw()+
  theme(panel.grid = element_blank())

# ggsave(filename = "hrtransmissionbox.jpeg", width = 5, 
#        height = 4, units = "in")

ggplot(data = time_to_elim, aes(x = l1, y = l2,
                                fill = nweek))+
  geom_tile()+
  scale_fill_viridis_c(name = "Week")+
  labs(x = "Within-cell transmission", y = "Home range transmission")+
  theme_bw()

# ggsave(filename = "transmissionheatmap.jpeg", width = 5,
#        height = 4, units = "in")

# Check with cases per week
ggplot(data = dis, aes(x = nweek, y = n_symptomatic, 
                       color = factor(rep)))+
  geom_line()+
  geom_hline(yintercept = 50, linetype = "dashed")+
  facet_grid(rows = vars(l1), cols = vars(l2))+
  scale_color_viridis_d(end = 0.9, name = "Rep")+
  xlim(c(52+1, 600))+
  theme_bw()+
  theme(panel.grid=element_blank())

# ggsave(filename = "facet_disease.jpeg", width = 10, height = 9,
#        units = "in")

mean_cases <- dis %>%
  filter(nweek > 52, elim == "False") %>%
  group_by(rep, l2, l1) %>%
  summarise(mean.cases = mean(n_symptomatic)) %>%
  mutate(l1 = factor(l1), 
         l2 = factor(l2))

ggplot(data=mean_cases, aes(x = l1, y = l2,
                           fill = mean.cases))+
  geom_tile()+
  scale_fill_viridis(name = "Mean Weekly Cases")+
  labs(x = "Within-cell transmission", y = "Home range transmission")+
  theme_bw()

# ggsave(filename = "meancases_heatmap.jpeg", width = 5,
#        height = 4, units = "in")

# Population sizes with disease ------------------
dis_pop <- dis %>%
  group_by(l1, l2, nweek) %>%
  summarise(mean_pop = mean(total_pop))

ggplot(data = dis_pop, aes (x = nweek, y = mean_pop, 
                        color = factor(l1)))+
  geom_line()+
  geom_vline(xintercept = 53, linetype = 'dashed')+
  scale_color_viridis_d(name = "Within-cell transmission",
                        end = 0.9)+
  facet_grid(rows = vars(l2))+
  theme_bw()+
  theme(panel.grid.minor = element_blank())

# ggsave(filename = "pop_direct.jpeg", width = 8, height = 6,
#        units = "in")

ggplot(data = dis_pop, aes (x = nweek, y = mean_pop, 
                            color = factor(l2)))+
  geom_line()+
  geom_vline(xintercept = 53, linetype = 'dashed')+
  scale_color_viridis_d(name = "Between-cell transmission",
                        end = 0.9)+
  labs(x = "Week", y = "Mean Population Size")+
  facet_grid(rows = vars(l1))+
  theme_bw()+
  theme(panel.grid.minor = element_blank())

# ggsave(filename = "pop_indirect.jpeg", width = 8, height = 6,
#        units = "in")

# Disease: fewer params, more reps per param -----------
dis <- read.csv("disease_test_smol.csv") %>%
  select(rep, year, week, total_pop, n_infected, n_symptomatic,
         elim, l1, l2) %>%
  mutate(nweek = ((year-1)*52)+week)

# Compare proportion of outbreaks eliminated
prop_eliminated <- dis %>%
  filter(year >= 2, elim == "True") %>%
  select(rep, l1, l2) %>%
  group_by(l1, l2) %>%
  distinct() %>%
  summarise(prop = n()/20)

unique(prop_eliminated$prop) #ok not bad

prop_eliminated <- prop_eliminated %>%
  # right_join(combos, by = c("l1", "l2")) %>%
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

# ggsave(filename = "propelimheat_smol.jpeg", width = 5,
#        height = 4, units = "in")

# Weeks to elimination
time_to_elim <- dis %>%
  group_by(l1, l2) %>%
  filter(year >= 2, elim == "True") %>%
  filter(nweek == min(nweek)) %>%
  # right_join(combos, by = c("l1", "l2")) %>%
  distinct() %>%
  mutate(l1 = factor(l1), 
         l2 = factor(l2))

ggplot(data = time_to_elim, aes(x = l1, y = nweek))+
  geom_boxplot(fill = "lightgray")+
  labs(x = "Within-cell transmission", y = "Week")+
  theme_bw()+
  theme(panel.grid = element_blank())

# ggsave(filename = "celltransmissionbox_smol.jpeg", width = 5,
#        height = 4, units = "in")

ggplot(data = time_to_elim, aes(x = l2, y = nweek))+
  geom_boxplot(fill = "lightgray")+
  labs(x = "Home Range Transmission", y = "Week")+
  theme_bw()+
  theme(panel.grid = element_blank())

# ggsave(filename = "hrtransmissionbox_smol.jpeg", width = 5,
#        height = 4,units = "in")

ggplot(data = time_to_elim, aes(x = l1, y = l2,
                                fill = nweek))+
  geom_tile()+
  scale_fill_viridis_c(name = "Week")+
  labs(x = "Within-cell transmission", y = "Home range transmission")+
  theme_bw()

# ggsave(filename = "transmissionheatmap_smol.jpeg", width = 5,
#        height = 4, units = "in")

# Check with cases per week
ggplot(data = dis, aes(x = nweek, y = n_symptomatic, 
                       color = factor(rep)))+
  geom_line()+
  geom_hline(yintercept = 50, linetype = "dashed")+
  facet_grid(rows = vars(l1), cols = vars(l2))+
  scale_color_viridis_d(end = 0.9, name = "Rep")+
  xlim(c(53, 600))+
  theme_bw()+
  theme(panel.grid=element_blank())

# ggsave(filename = "facet_disease.jpeg", width = 10, height = 9,
#        units = "in")

mean_cases <- dis %>%
  filter(nweek > 52, elim == "False") %>%
  group_by(rep, l2, l1) %>%
  summarise(mean.cases = mean(n_symptomatic)) %>%
  mutate(l1 = factor(l1), 
         l2 = factor(l2))

ggplot(data=mean_cases, aes(x = l1, y = l2,
                            fill = mean.cases))+
  geom_tile()+
  scale_fill_viridis(name = "Mean Weekly Cases")+
  labs(x = "Within-cell transmission", y = "Home range transmission")+
  theme_bw()

# ggsave(filename = "meancases_heatmap_smol.jpeg", width = 5,
#        height = 4,units = "in")

# Birth pulse ---------
births <- pop %>%
  filter(week %in% c(17,18))

end.year <- pop %>%
  filter(week == 52) %>%
  group_by(rep, a_mort, j_mort) %>%
  mutate(diff = total_pop - lag(total_pop, 
                                default = first(total_pop),
                                order_by=year)) %>%
  summarize(mean_growth = mean(diff/total_pop)) %>%
  filter(a_mort == 0.005 & j_mort == 0.02)

summary(end.year$mean_growth)

# R-naught -------------------
dis <- read.csv("disease_test.csv") %>%
  select(rep, year, week, total_pop, n_infected, n_symptomatic,
         elim, l1, l2) %>%
  mutate(nweek = ((year-1)*52)+week) %>%
  filter(l1 == 0.03 & l2 == 0.002) %>%
  filter(year > 1)

r0.list <- list()
for(i in 1:length(unique(dis$rep))){
  test <- filter(dis, rep==i)  

  r0.list[[i]] <- estimate.R(epid = test$n_symptomatic, 
             GT=generation.time("gamma", c(4.5, 1)),
             pop.size = test$total_pop,
             methods = c('ML', 'EG'))
}
print(r0.list)
