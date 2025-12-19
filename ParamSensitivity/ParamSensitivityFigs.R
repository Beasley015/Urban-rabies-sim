library(tidyverse)
library(patchwork)
library(viridis)
library(DescTools)
library(R0)

setwd("./ParamSensitivity")

# Carrying capacity: max K, adult mortality, juvie mortality ---------
kmort <- read.csv("K_mortality_sensitivity.csv") %>%
  select(rep, year, week, total_pop, a_mort, j_mort, maxK) %>%
  mutate(nweek = ((year-1)*52)+week)

seasonal <- kmort %>%
  filter(week %in% c(28:31, 41:44)) %>%
  mutate(season = case_when(week %in% c(28:31) == T ~ "Summer",
                            TRUE ~ "Fall"))

seasonal.summary <- seasonal %>%
  group_by(maxK, a_mort, j_mort, season) %>%
  summarise(mean_dens = mean(total_pop)/625, min = min(total_pop)/625,
            max = max(total_pop)/625)

# summer empirical range 7.75–27.25, mean 15
# fall empirical range 6.75–15, mean 11

likely.combos <- seasonal.summary %>%
  mutate(likely = case_when(season=="Summer" & mean_dens>12 & mean_dens<18 ~ "YES",
                            season=="Fall" & mean_dens>8 & mean_dens<14 ~ "YES",
                            TRUE ~ "NO")) %>%
  filter(likely == "YES") 
  # Drop combos that are not realistic in both seasons:
  
best.combos <- likely.combos %>%
  group_by(maxK, a_mort, j_mort) %>%
  summarise(nlikely = n()) %>%
  filter(nlikely == 2) %>%
  select(-nlikely)
  
# Filter best combos from full output
kmort.best <- best.combos %>%
  left_join(kmort, by = c('maxK', 'a_mort', 'j_mort'))

# Make figures: seasonal comparisons
season.allreps <- seasonal %>%
  group_by(maxK, a_mort, j_mort, season, rep) %>%
  summarise(mean_dens = mean(total_pop)/625, min = min(total_pop)/625,
            max = max(total_pop)/625)

best.season <- best.combos %>%
  left_join(season.allreps, by = c('maxK', 'a_mort', 'j_mort'))

dens_fig_list <- list()
for(i in 1:nrow(best.combos)){
  row <- best.combos[i,]
  
  test <- filter(best.season, maxK==row[1], a_mort==row[2], j_mort==row[3])
  
  dens_fig_list[[i]] <- ggplot(data=test, aes(x = season, y = mean_dens))+
    geom_boxplot(fill = 'lightgray')+
    geom_segment(aes(x = 0.65, xend = 1.35, y = 11, yend = 11),
                 linewidth = 1.5, linetype = "dashed")+
    geom_segment(aes(x = 1.65, xend = 2.35, y = 15, yend = 15),
                 linewidth = 1.5, linetype = "dashed")+
    labs(x = "Season", y = "Mean Density",
         title = paste("Max K = ", row[1], ", Adult Mortality = ", row[2],
                       ", Juvenile Mortality = ", row[3], sep = ""))+
    theme_bw(base_size = 12)+
    theme(panel.grid = element_blank()) 
}

layout <- "
AABB
CCDD
EEFF
"

dens_fig_list[[1]]+dens_fig_list[[2]]+dens_fig_list[[3]]+dens_fig_list[[4]]+
  dens_fig_list[[5]]+dens_fig_list[[6]]+
  plot_layout(design=layout)

ggsave(filename="density_sens.jpeg", width = 14, height = 10, units = "in")

# Make figures: total population
pop_fig_list <- list()
for(i in 1:nrow(best.combos)){
  row <- best.combos[i,]
  
  test <- filter(kmort.best, maxK==row[1], a_mort==row[2], j_mort==row[3])
  
  pop_fig_list[[i]] <- ggplot(data=test, aes(x = nweek, y = total_pop))+
    stat_summary(geom = "ribbon", fun.data = mean_cl_normal, fill = 'lightgray')+
    stat_summary(geom = "line", fun = mean)+
    labs(x = "Week", y = "Population Size",
         title = paste("Max K = ", row[1], ", Adult Mortality = ", row[2],
                       ", Juvenile Mortality = ", row[3], sep = ""))+
    theme_bw(base_size = 12)+
    theme(panel.grid=element_blank()) 
}

pop_fig_list[[3]]/pop_fig_list[[5]]/pop_fig_list[[6]]

# ggsave(filename="K_sens.jpeg", width = 14, height = 7, units = "in")

# Transmission Rates: lambda1 Wide Sweep -----------
# Combine files into 1, if needed
# filenames <- list.files(pattern = "broadsweep.csv")
# 
# file.list <- list()
# for(i in 1:length(filenames)){
#   file.list[[i]] <- read.csv(filenames[i])
# }
# 
# lambda1.frame <- do.call(rbind, file.list)
# 
# write.csv(lambda1.frame, "l1_sens_wide.csv")
# file.remove(filenames)

# Read in data
dis.wide <- read.csv("l1_sens_wide.csv") %>%
  select(rep, year, week, total_pop, n_infected, 
         n_symptomatic, elim, lambda1) %>%
  mutate(nweek = ((year-1)*52)+week)

# proportion eliminated
prop_eliminated <- dis.wide %>%
  filter(year >= 2, elim == "True") %>%
  select(rep, lambda1) %>%
  group_by(lambda1) %>%
  distinct() %>%
  summarise(prop = n()/10)
# 0.0224 has lowest elimination rate at 60%

unique(prop_eliminated$prop)

# Time to elimination
time_to_elim <- dis.wide %>%
  group_by(lambda1, rep) %>%
  filter(year >= 2, elim == "True") %>%
  filter(nweek == min(nweek)) %>%
  # Filter trials with no persistence:
  filter(lambda1 > 0.005 & lambda1 < 0.13) %>% 
  distinct() %>%
  mutate(lambda1 = factor(lambda1))

all_combos <- expand_grid(unique(time_to_elim$lambda1), 
                          unique(time_to_elim$rep))
colnames(all_combos) <- c("lambda1", "rep")

time_to_elim <- time_to_elim %>%
  right_join(all_combos, by=c("rep", "lambda1")) %>%
  mutate(nweek = case_when(is.na(nweek)==T ~ 52*11,
                           TRUE ~ nweek))

time_to_elim %>%
  ungroup() %>%
  group_by(lambda1) %>%
  summarise(meantime = median(nweek))

ggplot(data = time_to_elim, aes(x = lambda1, y = nweek))+
  geom_boxplot(fill = "lightgray")+
  labs(x = "Transmission Rate", y = "Week of Elimination")+
  theme_bw(base_size = 14)+
  theme(panel.grid = element_blank())

ggsave("weekelim_wide.jpeg", width = 8, height = 6,
       units = "in")

# Weekly cases
mean_cases <- dis.wide %>%
  filter(nweek > 70, elim == "False") %>%
  filter(lambda1 > 0.005 & lambda1 < 0.13) %>% 
  mutate(lambda1 = factor(lambda1, levels = sort(unique(lambda1)))) %>%
  group_by(rep, lambda1) %>%
  summarise(mean.cases = median(n_symptomatic))

mean_cases %>%
  ungroup() %>%
  group_by(lambda1) %>%
  summarise(meancases = median(mean.cases))

ggplot(data=mean_cases, aes(x=lambda1, y = mean.cases))+
  geom_boxplot(fill = 'lightgray')+
  labs(x = "Transmission Rate", y = "Median Weekly Cases")+
  theme_bw(base_size = 14)+
  theme(panel.grid = element_blank())

# ggsave("medcase_wide.jpeg", width = 8, height = 6,
#        units = "in")

# R_e calculation
dis.wide <- dis.wide %>%
  filter(lambda1 > 0.005 & lambda1 < 0.13)

first_elim <- dis.wide %>%
  filter(year >= 2) %>%
  filter(elim == "True") %>%
  group_by(rep, lambda1) %>%
  summarise(first = min(nweek))

r.list.l1 <- list()
for(i in 1:length(unique(dis.wide$rep))){
  for(j in 1:length(unique(dis.wide$lambda1))){
    test <- filter(dis.wide, rep==i & 
                     lambda1==unique(dis.wide$lambda1)[j]) %>%
      filter(year >= 2)
    
    elim_test <- filter(first_elim, 
                        rep==i & lambda1==unique(dis.wide$lambda1)[j]) 
    
    if(nrow(elim_test) != 1){next}
    
    start <- as.numeric(min(which(test$n_symptomatic>0)))
    end <- elim_test$first-53
    
    re.test <- try(re <- estimate.R(epid = 
                                      test$n_symptomatic[start:end], 
                                    GT=generation.time("gamma", c(4.5, 1)),
                                    method = 'TD', nsim = 1000))
    
    if(class(re.test) %in% 'try-error') {next} else {
      re <- estimate.R(epid = 
                         test$n_symptomatic[start:end], 
                       GT=generation.time("gamma", c(4.5, 1)),
                       method = 'TD', nsim = 1000)
    }
    
    vec <- c(unique(test$rep), unique(test$lambda1),
             median(re$estimates$TD$R))
    
    len <- length(r.list.l1)
    r.list.l1[[len+1]] <- vec
  }
}

re.df <- as.data.frame(do.call(rbind, r.list.l1))
colnames(re.df) <- c("rep", "lambda1", "Re")

ggplot(data = re.df, aes(x = factor(lambda1), y = Re))+
  geom_boxplot(fill = 'lightgray')+
  geom_hline(yintercept = 1.13, linetype = "dashed",
             linewidth = 1)+
  labs(x = expression(lambda[1]), y = expression(R[e]))+
  theme_bw(base_size = 14)+
  theme(panel.grid = element_blank())

# ggsave("re_l1_wide.jpeg", width = 8, height = 6, units = "in")

# Transmission Rates: l2 wide sweep -------------
# Combine files 
# filenames <- list.files(pattern = "lambda2")
# 
# file.list <- list()
# for(i in 1:length(filenames)){
#   file.list[[i]] <- read.csv(filenames[i])
# }
# 
# lambda2.frame <- do.call(rbind, file.list)
# 
# write.csv(lambda2.frame, "l2_sens_wide.csv")
# file.remove(filenames)

# Read in data
dis.wide <- read.csv("l2_sens_wide.csv") %>%
  select(rep, year, week, total_pop, n_infected, 
         n_symptomatic, elim, lambda2) %>%
  mutate(nweek = ((year-1)*52)+week)

# proportion eliminated
prop_eliminated <- dis.wide %>%
  filter(year >= 2, elim == "True") %>%
  select(rep, lambda2) %>%
  group_by(lambda2) %>%
  distinct() %>%
  summarise(prop = n()/10)
# Elimination ranges from 20% to 70%

unique(prop_eliminated$prop)

# Time to elimination
time_to_elim <- dis.wide %>%
  group_by(lambda2, rep) %>%
  filter(year >= 2, elim == "True") %>%
  filter(nweek == min(nweek)) %>%
  distinct() %>%
  mutate(lambda2 = factor(lambda2))

all_combos <- expand_grid(unique(time_to_elim$lambda2), 
                          unique(time_to_elim$rep))
colnames(all_combos) <- c("lambda2", "rep")

time_to_elim <- time_to_elim %>%
  right_join(all_combos, by=c("rep", "lambda2")) %>%
  mutate(nweek = case_when(is.na(nweek)==T ~ 52*11,
                           TRUE ~ nweek))

time_to_elim %>%
  ungroup() %>%
  group_by(lambda2) %>%
  summarise(meantime = median(nweek))

ggplot(data = time_to_elim, aes(x = lambda2, y = nweek))+
  geom_boxplot(fill = "lightgray")+
  labs(x = "Transmission Rate", y = "Week of Elimination")+
  theme_bw(base_size = 14)+
  theme(panel.grid = element_blank())

ggsave("weekelim_wide_l2.jpeg", width = 8, height = 6,
       units = "in")
# Lots of variability

# Weekly cases
mean_cases <- dis.wide %>%
  filter(nweek > 70, elim == "False") %>%
  mutate(lambda2 = factor(lambda2, levels = sort(unique(lambda2)))) %>%
  group_by(rep, lambda2) %>%
  summarise(mean.cases = median(n_symptomatic))

mean_cases %>%
  ungroup() %>%
  group_by(lambda2) %>%
  summarise(meancases = median(mean.cases))
# < 0.002 or 0.0152 are most promising

ggplot(data=mean_cases, aes(x=lambda2, y = mean.cases))+
  geom_boxplot(fill = 'lightgray')+
  labs(x = "Transmission Rate", y = "Median Weekly Cases")+
  theme_bw(base_size = 14)+
  theme(panel.grid = element_blank())

# ggsave("medcase_wide_l2.jpeg", width = 8, height = 6,
#        units = "in")

# R_e calculation
first_elim <- dis.wide %>%
  filter(year >= 2) %>%
  filter(elim == "True") %>%
  group_by(rep, lambda2) %>%
  summarise(first = min(nweek)) %>%
  mutate(lambda2 = factor(lambda2))

r.list.l2 <- list()
for(i in 1:length(unique(dis.wide$rep))){
  for(j in 1:length(unique(dis.wide$lambda2))){
    test <- filter(dis.wide, rep==i & 
                     lambda2==unique(dis.wide$lambda2)[j]) %>%
      filter(year >= 2)
    
    elim_test <- filter(first_elim, 
                        rep==i & lambda2==unique(dis.wide$lambda2)[j]) 
    
    if(nrow(elim_test) != 1){next}
    
    start <- as.numeric(min(which(test$n_symptomatic>0)))
    end <- elim_test$first-53
    
    re.test <- try(re <- estimate.R(epid = 
                                      test$n_symptomatic[start:end], 
                                    GT=generation.time("gamma", c(4.5, 1)),
                                    method = 'TD', nsim = 1000))
    
    if(class(re.test) %in% 'try-error') {next} else {
      re <- estimate.R(epid = 
                         test$n_symptomatic[start:end], 
                       GT=generation.time("gamma", c(4.5, 1)),
                       method = 'TD', nsim = 1000)
    }
    
    vec <- c(unique(test$rep), unique(test$lambda2),
             median(re$estimates$TD$R))
    
    len <- length(r.list.l2)
    r.list.l2[[len+1]] <- vec
  }
}

re.df <- as.data.frame(do.call(rbind, r.list.l2))
colnames(re.df) <- c("rep", "lambda2", "Re")

ggplot(data = re.df, aes(x = factor(lambda2), y = Re))+
  geom_boxplot(fill = 'lightgray')+
  geom_hline(yintercept = 1.13, linetype = "dashed",
             linewidth = 1)+
  labs(x = expression(lambda[2]), y = expression(R[e]))+
  theme_bw(base_size = 14)+
  theme(panel.grid = element_blank())

# ggsave("re_l2_wide.jpeg", width = 8, height = 6, units = "in")

# Glance at median Re
re.df %>%
  group_by(lambda2) %>%
  summarise(re = median(Re))
# All are *slightly* too high... but the 0.004-0.008 range is closest

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

unique(prop_eliminated$prop)

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
  group_by(l1, rep, l2) %>%
  filter(year >= 2, elim == "True") %>%
  filter(nweek == min(nweek)) %>%
  # right_join(combos, by = c("l1", "l2")) %>%
  distinct() %>%
  mutate(l1 = factor(l1), 
         l2 = factor(l2, levels = c(0.007, 0.008, 0.009, 0.01,
                                    0.011)))

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
  filter(nweek > 100, elim == "False") %>%
  group_by(rep, l1,l2) %>%
  summarise(mean.cases = median(n_symptomatic)) %>%
  mutate(l1 = factor(l1),l2 = factor(l2))

ggplot(data=mean_cases, aes(x=l1, y = mean.cases))+
  geom_boxplot()

ggplot(data=mean_cases, aes(x = l1, y = l2,
                           fill = mean.cases))+
  geom_tile()+
  scale_fill_viridis(name = "Median Weekly Cases")+
  labs(x = "Within-cell transmission", y = "Home range transmission")+
  theme_bw()

# ggsave(filename = "meancases_heatmap.jpeg", width = 5,
#        height = 4, units = "in")

# Population sizes with disease ------------------
dis_pop <- dis %>%
  group_by(nweek,l1) %>%
  summarise(mean_pop = mean(total_pop))

ggplot(data = dis_pop, aes (x = nweek, y = mean_pop))+
  geom_line()+
  geom_vline(xintercept = 53, linetype = 'dashed')+
  # scale_color_viridis_d(name = "Within-cell transmission",
  #                       end = 0.9)+
  # facet_grid(rows = vars(l2))+
  labs(x="Week", y = "Mean Population Size")+
  theme_bw(base_size = 14)+
  theme(panel.grid = element_blank())

# ggsave(filename = "pop_direct.jpeg", width = 8, height = 6,
#        units = "in")

ggplot(data = dis_pop, aes (x = nweek, y = mean_pop, 
                            color = factor(l2)))+
  geom_line()+
  geom_vline(xintercept = 53, linetype = 'dashed')+
  scale_color_viridis_d(name = "Between-cell transmission",
                        end = 0.9)+
  labs(x = "Week", y = "Mean Population Size")+
  # facet_grid(rows = vars(l1))+
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
  group_by(l1, l2, rep) %>%
  filter(year >= 2, elim == "True") %>%
  filter(nweek == min(nweek)) %>%
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

time_to_elim %>%
  group_by(l2) %>%
  summarise(median = median(nweek))

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
  filter(nweek > 100, elim == "False") %>%
  group_by(rep, l2, l1) %>%
  summarise(mean.cases = mean(n_symptomatic)) %>%
  mutate(l1 = factor(l1), 
         l2 = factor(l2))

ggplot(data=mean_cases, aes(x = l2, y = mean.cases))+
  geom_boxplot(fill='lightgray')+
  scale_fill_viridis(name = "Mean Weekly Cases")+
  labs(x = "Within-cell transmission", y = "Home range transmission")+
  theme_bw()+
  theme(panel.grid=element_blank())

# ggsave(filename = "meancases_heatmap_smol.jpeg", width = 5,
#        height = 4,units = "in")

# One transmission value, just to look -----------------
dis <- read.csv("disease_test_035.csv") %>%
  select(rep, year, week, total_pop, n_infected, n_symptomatic,
         elim) %>%
  mutate(nweek = ((year-1)*52)+week)

# Proportion eliminated
dis %>%
  filter(year >= 2, elim == "True") %>%
  select(rep) %>%
  distinct() %>%
  summarise(prop = n()/20) #90%: not great, no worse than others

time_to_elim <- dis %>%
  group_by(rep) %>%
  filter(year >= 2, elim == "True") %>%
  filter(nweek == min(nweek)) 

summary(time_to_elim$nweek) # Mean 215, again no worse than others
# median is slightly lower than with an l2

# Birth pulse ---------
births <- pop %>%
  filter(week %in% c(17,18))

end.year <- pop %>%
  filter(week == 52) %>%
  group_by(rep, amort, jmort) %>%
  mutate(diff = total_pop - lag(total_pop, 
                                default = first(total_pop),
                                order_by=year)) %>%
  summarize(mean_growth = mean(diff/total_pop)) %>%
  filter(amort == 0.005 & jmort == 0.02)

summary(end.year$mean_growth)

# R-naught -------------------
dis <- read.csv("disease_test_l2wide.csv") %>%
  select(rep, year, week, total_pop, n_infected, 
         n_symptomatic, elim, l2) %>%
  mutate(nweek = ((year-1)*52)+week) %>%
  filter(year > 1)

dis <- read.csv("disease_test_smol.csv") %>%
  select(rep, year, week, total_pop, n_infected, n_symptomatic,
         elim, l2) %>%
  mutate(nweek = ((year-1)*52)+week) %>%
  # filter(l1 == 0.09 & l2 == 0.01) %>%
  filter(year > 1)

first_elim <- dis %>%
  filter(elim == "True") %>%
  group_by(rep, l2) %>%
  summarise(first = min(nweek))

r0.list <- list()
for(i in 1:length(unique(dis$rep))){
  # for(j in 1:length(unique(dis$l1))){
    for(k in 1:length(unique(dis$l2))){
      test <- filter(dis, rep==i & l2==unique(dis$l2)[k]) 
      elim_test <- filter(first_elim, 
                          rep==i & l2==unique(dis$l2)[k]) 
    
      if(nrow(elim_test) != 1){next}
  
      start <- as.numeric(min(which(test$n_symptomatic>0)))
      end <- elim_test$first-53

      r0.test <- try(r0 <- estimate.R(epid = 
                               test$n_symptomatic[start:end], 
                  GT=generation.time("gamma", c(4.5, 1)),
                  method = 'TD', nsim = 1000))
    
      if(class(r0.test) %in% 'try-error') {next} else {
            r0 <- estimate.R(epid = 
                               test$n_symptomatic[start:end], 
                  GT=generation.time("gamma", c(4.5, 1)),
                  method = 'SB', nsim = 1000)
      }
    
      vec <- c(unique(test$rep), unique(test$l2),
              median(r0$estimates$TD$R))
             
      len <- length(r0.list)
      r0.list[[len+1]] <- vec
     }
   # }
}

r0.df <- as.data.frame(do.call(rbind, r0.list))
colnames(r0.df) <- c("rep", "l1", "R0")

r0.df <- r0.df %>%
  group_by(rep, l1) %>%
  filter(R0 < quantile(.$R0, 0.975) && 
           Re > quantile(.$R0, 0.025))

r0.df %>%
  ungroup() %>%
  group_by(l1) %>%
  summarise(mean.r0 = median(R0))

ggplot(data = r0.df, aes(x = factor(l1), y = R0))+
  geom_boxplot(fill='lightgray')+
  labs(x = "Core Transmission", y = bquote(R[0]))+
  theme_bw(base_size = 14)+
  theme(panel.grid=element_blank())

# ggsave(filename="R_0_sensitivity.jpeg", width = 5,
# height = 3.5, units = 'in')

# Epicurves, just for fun ----------------------
curve.df <- dis %>%
  filter(l1 == 0.035, l2 ==0.007) %>%
  mutate(nweek = (year-1)*52+week) %>%
  select(rep, nweek, n_symptomatic)


ggplot(data = curve.df, aes(x=nweek, y = n_symptomatic,
                            group = rep, color = factor(rep)))+
  geom_smooth(se = F, span = 0.2)+
  scale_y_continuous(limits = c(0, NA))+
  scale_x_continuous(limits = c(NA, 350))+
  scale_color_viridis_d(end = 0.9)+
  labs(x = "Week", y = "Cases")+
  theme_bw(base_size = 15)+
  theme(panel.grid=element_blank(), legend.position = "None")

ggsave(filename="example_epicurve.jpeg", width = 5, height = 4, 
       units = "in")
