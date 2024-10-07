library(tidyverse)
library(viridis)

setwd("./ParamSensitivity")

# Population size ---------------------
pop <- read.csv("popsize.csv") %>%
  select(rep, year, week, total_pop) %>%
  mutate(nweek = ((year-1)*52)+week)

ggplot(data=pop, aes(x = nweek, y = total_pop, 
                     color = factor(rep)))+
  geom_line()+
  scale_color_viridis_d(name = "Rep", end = 0.9)+
  labs(x = "Week", y = "Population Size")+
  theme_bw(base_size = 12)+
  theme(panel.grid = element_blank())

# ggsave("maxcc_12.jpeg", width = 5, height = 4, units = "in")

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
ccmort <- read.csv("carrying_capacity_mortality.csv") %>%
  select(rep, year, week, total_pop, ac_mort, jc_mort) %>%
  mutate(nweek = ((year-1)*52)+week)

# Full, faceted plot
ggplot(data=ccmort, aes(x = nweek, y = total_pop, color = factor(rep)))+
  geom_line()+
  facet_grid(rows = vars(ac_mort), cols=vars(jc_mort))+
  scale_color_viridis_d(end = 0.9, name = "Rep")+
  # geom_vline(xintercept = ((0:9)*52)+20, linetype = "dashed")+
  # geom_vline(xintercept = ((0:9)*52)+43, linetype = "dotted")+
  labs(x = "Week", y = "Population Size")+
  theme_bw()+
  theme(panel.grid = element_blank())

# ggsave("full_cc_birthrate95.jpeg", width=10, height=9, units = "in")

# Look at a particularly interesting panel, 
# where effects of juvie dispersal are visible

juv.disp <- ccmort %>%
  filter(ac_mort == 0.005 & jc_mort == 0.02)

ggplot(data=juv.disp, aes(x=nweek, y = total_pop, color = factor(rep)))+
  geom_line()+
  scale_color_viridis_d(end = 0.9, name = "Rep")+
  geom_vline(xintercept = ((0:9)*52)+43, linetype = "dashed", alpha = 0.5)+
  labs(x = "Week", y = "Population Size")+
  theme_bw()+
  theme(panel.grid = element_blank())

# ggsave("juvie_disperse.jpeg", width=5, height=4, units = "in")

jmort_only <- ccmort %>%
  filter(ac_mort == 0.005) %>%
  group_by(jc_mort, nweek) %>%
  summarise(mean_pop = mean(total_pop))

ggplot(data=jmort_only, aes(x = nweek, y = mean_pop, 
                            color = factor(jc_mort)))+
  geom_line()+
  scale_color_viridis_d(end = 0.9, name = "Juvenile Mortality")+
  labs(x = "Week", y = "Mean Population Size")+
  theme_bw()+
  theme(panel.grid = element_blank())

# ggsave("cc_adult005_birthrate95.jpeg", width=5, height=3, units = "in")

# Transmission Rates -----------------------
# Note: pop immunity set at 0, immigration rate low, no imm disease
dis <- read.csv("disease_test.csv") %>%
  select(rep, year, week, total_pop, n_infected, n_symptomatic,
         direct_prob, indirect_prob, elim) %>%
  mutate(nweek = ((year-1)*52)+week)

# Compare proportion of outbreaks eliminated
prop_eliminated <- dis %>%
  filter(year >= 5, elim == "True") %>%
  select(rep, direct_prob, indirect_prob) %>%
  group_by(direct_prob, indirect_prob) %>%
  distinct() %>%
  summarise(prop = n()/5)

unique(prop_eliminated$prop) #ok not bad

# Proportion eliminated: figs
ggplot(prop_eliminated, aes(x = factor(direct_prob), 
                            y = factor(indirect_prob),
                            fill = factor(prop)))+
  geom_tile()+
  scale_fill_viridis_d(name = "Proportion eliminated")+
  labs(x = "Within-cell transmission", y = "Home range transmission")+
  theme_bw()

# ggsave(filename = "propelimheat.jpeg", width = 5, height = 4,
#        units = "in")

# Weeks to elimination
time_to_elim <- dis %>%
  group_by(direct_prob, indirect_prob) %>%
  filter(year >= 5, elim == "True") %>%
  filter(nweek == min(nweek)) %>%
  mutate(direct_prob = factor(direct_prob), 
         indirect_prob = factor(indirect_prob))

ggplot(data = time_to_elim, aes(x = direct_prob, y = nweek))+
  geom_boxplot(fill = "lightgray")+
  labs(x = "Within-cell transmission", y = "Week")+
  theme_bw()+
  theme(panel.grid = element_blank())

# ggsave(filename = "celltransmissionbox.jpeg", width = 5, height = 4,
#        units = "in")

ggplot(data = time_to_elim, aes(x = indirect_prob, y = nweek))+
  geom_boxplot(fill = "lightgray")+
  labs(x = "Home Range Transmission", y = "Week")+
  theme_bw()+
  theme(panel.grid = element_blank())

# ggsave(filename = "hrtransmissionbox.jpeg", width = 5, height = 4,
#        units = "in")

ggplot(data = time_to_elim, aes(x = direct_prob, y = indirect_prob,
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
  facet_grid(rows = vars(direct_prob), cols = vars(indirect_prob))+
  scale_color_viridis_d(end = 0.9, name = "Rep")+
  xlim(c(52*4+1, 800))+
  theme_bw()+
  theme(panel.grid=element_blank())

# ggsave(filename = "facet_disease.jpeg", width = 10, height = 9,
#        units = "in")

# looks like home range transmission influences initial peak
# Check max number of weekly cases

max_cases <- dis %>%
  group_by(rep, indirect_prob, direct_prob) %>%
  summarise(max.cases = max(n_symptomatic)) %>%
  mutate(direct_prob = factor(direct_prob), 
         indirect_prob = factor(indirect_prob))

ggplot(data=max_cases, aes(x = direct_prob, y = max.cases))+
  geom_boxplot(fill = "lightgray")+
  labs(x = "Within-cell transmission", y = "Max Weekly Cases")+
  theme_bw()+
  theme(panel.grid = element_blank())

# ggsave(filename = "maxcases_cell_box.jpeg", width = 5, height = 4,
#        units = "in")

ggplot(data=max_cases, aes(x = indirect_prob, y = max.cases))+
  geom_boxplot(fill = "lightgray")+
  labs(x = "Home range transmission", y = "Max Weekly Cases")+
  theme_bw()+
  theme(panel.grid = element_blank())

# ggsave(filename = "maxcases_hr_box.jpeg", width = 5, height = 4,
#        units = "in")

ggplot(data=max_cases, aes(x = direct_prob, y = indirect_prob,
                           fill = max.cases))+
  geom_tile()+
  scale_fill_viridis(name = "Max Weekly Cases")+
  labs(x = "Within-cell transmission", y = "Home range transmission")+
  theme_bw()

# ggsave(filename = "maxcases_heatmap.jpeg", width = 5, height = 4,
#        units = "in")

# Look at lower values of home-range transmission
# To check weekly case numbers

smol <- dis %>%
  filter(indirect_prob < 0.01)

ggplot(data = smol, aes(x = nweek, y = n_symptomatic, 
                       color = factor(rep)))+
  geom_line()+
  facet_grid(rows = vars(direct_prob), cols = vars(indirect_prob))+
  scale_color_viridis_d(end = 0.9, name = "Rep")+
  xlim(c(52*4+1, 800))+
  theme_bw()+
  theme(panel.grid=element_blank())

# ggsave("weekly_cases_smol.jpeg", width = 6, height = 5, 
#        units = "in")

dis_pop <- dis %>%
  group_by(direct_prob, indirect_prob, nweek) %>%
  summarise(mean_pop = mean(total_pop))

ggplot(data = dis_pop, aes (x = nweek, y = mean_pop, 
                        color = factor(direct_prob)))+
  geom_line()+
  scale_color_viridis_d(name = "Within-cell transmission",
                        end = 0.9)+
  facet_grid(rows = vars(indirect_prob))+
  theme_bw()+
  theme(panel.grid.minor = element_blank())

# ggsave(filename = "pop_direct.jpeg", width = 8, height = 6,
#        units = "in")

ggplot(data = dis_pop, aes (x = nweek, y = mean_pop, 
                            color = factor(indirect_prob)))+
  geom_line()+
  scale_color_viridis_d(name = "Home range transmission",
                        end = 0.9)+
  facet_grid(rows = vars(direct_prob))+
  theme_bw()+
  theme(panel.grid.minor = element_blank())

# ggsave(filename = "pop_indirect.jpeg", width = 8, height = 6,
#        units = "in")
