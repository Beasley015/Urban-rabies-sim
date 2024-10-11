library(tidyverse)
library(viridis)

setwd("./ParamSensitivity")

# Population size ---------------------
pop <- read.csv("popsize10.csv") %>%
  select(rep, year, week, total_pop) %>%
  mutate(nweek = ((year-1)*52)+week)

ggplot(data=pop, aes(x = nweek, y = total_pop, 
                     color = factor(rep)))+
  geom_line()+
  scale_color_viridis_d(name = "Rep", end = 0.9)+
  labs(x = "Week", y = "Population Size")+
  theme_bw(base_size = 12)+
  theme(panel.grid = element_blank())

# ggsave("maxcc_10.jpeg", width = 5, height = 4, units = "in")

# Land cover test ---------------
# Do raccoons aggregate in favorable habitat?

for(i in 1:5){
  # Read in land cover data
  land <- read.csv(paste("land", i, ".csv", sep = ""), header = F)
  land <- land[-c(1:5, 56:60), -c(1:5, 56:60)]

  # Tabulate land cover types 
  landtab <- table(factor(unlist(land), 
                          levels = unique(unlist(land))))

  # Combine last 5 habitats because they are rare
  landtab[4] <- sum(landtab[4:length(landtab)])
  landtab <- landtab[1:4]
  names(landtab) <- c("forest", "developed", "pasture",
                      "other")

  # Convert to proportions
  landtab <- landtab/sum(landtab)
  
  # Load raccoon data and label some columns
  landprops <- as.data.frame(t(read.csv("raccoon_land_props.csv")))
  colnames(landprops)[1:6] <- c("rep", "year", "week","forest",
                              "developed", "pasture")

  # Combine remaining habitats into "other" category
  landprops <- landprops %>%
    filter(rep == i) %>%
    mutate("Other" = V7+V8+V9+V10+V11) %>%
    dplyr::select(rep:pasture, Other) %>%
    summarise_at(.vars = vars(forest:Other), .funs = mean)

  # Compare with bar plots
  land.compare <- data.frame(Habitat = colnames(landprops),
                             Actual = as.vector(landtab),
                             Used = as.numeric(landprops[1,])) %>%
    pivot_longer(cols = -Habitat, names_to = "type", 
                 values_to = "prop")

  ggplot(data = land.compare, aes(x = Habitat, y = prop, 
                                  fill = type))+
    geom_col(position = "dodge", color = "black") +
    scale_fill_manual(values = c("lightgray", "limegreen"), 
                      name = "")+
    labs(y = "Proportion")+
    theme_bw(base_size=12)+
    theme(legend.title = , panel.grid = element_blank())

  ggsave(filename = paste("landprop", i, ".jpeg"), width = 6,
         height = 4, units = "in")
}

# Mortality rates from carrying capacity -------------------
# Read in simulation results
ccmort <- read.csv("carrying_capacity_mortality12.csv") %>%
  select(rep, year, week, total_pop) %>%
  mutate(nweek = ((year-1)*52)+week)

# Add mortality probs because I forgot to save them
a = c(0.005, 0.01, 0.015)
j = c(0.015, 0.02, 0.025)

mort.probs <- expand.grid(j,a)
colnames(mort.probs) <- c("j_mort", "a_mort")

mort.full <- mort.probs[rep(row.names(mort.probs), 
                            each = max(ccmort$nweek)),]

ccmort <- cbind(ccmort, mort.full)

# Full, faceted plot
ggplot(data=ccmort, aes(x = nweek, y = total_pop, 
                        color = factor(rep)))+
  geom_line()+
  facet_grid(rows = vars(a_mort), cols=vars(j_mort))+
  scale_color_viridis_d(end = 0.9, name = "Rep")+
  # geom_vline(xintercept = ((0:9)*52)+20, linetype = "dashed")+
  # geom_vline(xintercept = ((0:9)*52)+43, linetype = "dotted")+
  labs(x = "Week", y = "Population Size")+
  theme_bw(base_size = 24)+
  theme(panel.grid = element_blank())

# ggsave("full_cc_max12.jpeg", width=12, height=10, units = "in")

# Look at panel(s) that are most realistic

ccmort.005.02 <- ccmort %>%
  filter(a_mort == 0.005 & j_mort == 0.02)

ccmort.005.025 <- ccmort %>%
  filter(a_mort == 0.005 & j_mort == 0.025)

ccmort.01.02 <- ccmort %>%
  filter(a_mort == 0.01 & j_mort == 0.02)

ggplot(data=ccmort.01.02, aes(x=nweek, y = total_pop, 
                          color = factor(rep)))+
  geom_line()+
  scale_color_viridis_d(end = 0.9, name = "Rep")+
  geom_hline(yintercept = 6000, linetype = "dashed")+
  labs(x = "Week", y = "Population Size")+
  theme_bw()+
  theme(panel.grid = element_blank())

# ggsave("smol_panel_max12_01_02.jpeg", width=5, height=4, 
#        units = "in")

# Transmission Rates -----------------------
# Note: pop immunity set at 0, immigration rate low, no imm disease
dis <- read.csv("disease_rates.csv") %>%
  select(rep, year, week, total_pop, n_infected, n_symptomatic, elim,
         cell, homerange) %>%
  mutate(nweek = ((year-1)*52)+week)

# Compare proportion of outbreaks eliminated
prop_eliminated <- dis %>%
  filter(year >= 5, elim == "True") %>%
  select(rep, cell, homerange) %>%
  group_by(cell, homerange) %>%
  distinct() %>%
  summarise(prop = n()/5)

unique(prop_eliminated$prop) #ok not bad

# Proportion eliminated: figs
ggplot(prop_eliminated, aes(x = factor(cell), 
                            y = factor(homerange),
                            fill = factor(prop)))+
  geom_tile()+
  scale_fill_viridis_d(name = "Proportion eliminated")+
  labs(x = "Within-cell transmission", y = "Home range transmission")+
  theme_bw()

# ggsave(filename = "propelimheat.jpeg", width = 5, height = 4,
#        units = "in")

# Weeks to elimination
time_to_elim <- dis %>%
  group_by(cell, homerange) %>%
  filter(year >= 5, elim == "True") %>%
  filter(nweek == min(nweek)) %>%
  mutate(cell = factor(cell), 
         homerange = factor(homerange))

ggplot(data = time_to_elim, aes(x = cell, y = nweek))+
  geom_boxplot(fill = "lightgray")+
  labs(x = "Within-cell transmission", y = "Week")+
  theme_bw()+
  theme(panel.grid = element_blank())

# ggsave(filename = "celltransmissionbox.jpeg", width = 5, height = 4,
#        units = "in")

ggplot(data = time_to_elim, aes(x = homerange, y = nweek))+
  geom_boxplot(fill = "lightgray")+
  labs(x = "Home Range Transmission", y = "Week")+
  theme_bw()+
  theme(panel.grid = element_blank())

# ggsave(filename = "hrtransmissionbox.jpeg", width = 5, height = 4,
#        units = "in")

ggplot(data = time_to_elim, aes(x = cell, y = homerange,
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
  facet_grid(rows = vars(cell), cols = vars(homerange))+
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

# ggsave("weekly_cases_start20.jpeg", width = 6, height = 4,
#        units = "in")
