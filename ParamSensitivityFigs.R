library(tidyverse)
library(viridis)
library(DescTools)
library(R0)

setwd("./ParamSensitivity")

# Population size ---------------------
pop <- read.csv("kmax10.csv") %>%
  select(rep, year, week, total_pop, amort, jmort) %>%
  mutate(nweek = ((year-1)*52)+week)

ggplot(data=pop, aes(x = nweek, y = total_pop, 
                     color = factor(rep)))+
  geom_line()+
  scale_color_viridis_d(name = "Rep", end = 0.9)+
  facet_grid(rows = vars(amort), cols=vars(jmort))+
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
  group_by(amort, jmort, season) %>%
  summarise(mean = mean(total_pop), min = min(total_pop),
            max = max(total_pop))

ggplot(data = seasonal, aes(x = season, y = total_pop))+
  geom_boxplot(fill = 'lightgray')+
  facet_grid(rows = vars(amort), cols = vars(jmort))+
  labs(x = "Season", y = "Carrying Capacity", title = "Max K = 40")+
  theme_bw()+
  theme(panel.grid = element_blank())

# ggsave("seasonal_pop_10.jpeg", width = 5, height = 4,
#        units = "in")

# Transmission Rates: Wide Sweep -----------
# dis.wide <- read.csv("disease_test_widenet.csv") %>%
dis.wide <- read.csv("disease_test_l1.csv") %>%
  select(rep, year, week, total_pop, n_infected, 
         n_symptomatic, elim, l1) %>%
  mutate(nweek = ((year-1)*52)+week)

# proportion eliminated
prop_eliminated <- dis.wide %>%
  filter(year >= 2, elim == "True") %>%
  select(rep, l1) %>%
  group_by(l1) %>%
  distinct() %>%
  summarise(prop = n()/10)

unique(prop_eliminated$prop)

# Time to elimination
time_to_elim <- dis.wide %>%
  group_by(l1, rep) %>%
  filter(year >= 2, elim == "True") %>%
  filter(nweek == min(nweek)) %>%
  distinct() %>%
  mutate(l1 = factor(l1))

time_to_elim %>%
  ungroup() %>%
  group_by(l1) %>%
  summarise(meantime = median(nweek))

elim_l1 <- ggplot(data = time_to_elim, aes(x = l1, y = nweek))+
  geom_boxplot(fill = "lightgray")+
  labs(x = "Transmission Rate", y = "Week of Elimination")+
  theme_bw()+
  theme(panel.grid = element_blank())

(elim_wide/elim_l1)+
  plot_annotation(tag_levels = "a")

# ggsave("weekelim_wide.jpeg", width = 10, height = 4,
#        units = "in")

mean_cases <- dis.wide %>%
  filter(nweek > 70, elim == "False") %>%
  mutate(l1 = factor(l1, levels = sort(unique(l1)))) %>%
  group_by(rep, l1) %>%
  summarise(mean.cases = median(n_symptomatic))

mean_cases %>%
  ungroup() %>%
  group_by(l1) %>%
  summarise(meancases = median(mean.cases))

ggplot(data=mean_cases, aes(x=l1, y = mean.cases))+
  geom_boxplot(fill = 'lightgray')+
  labs(x = "Transmission Rate", y = "Median Weekly Cases")+
  theme_bw(base_size = 14)+
  theme(panel.grid = element_blank())

# ggsave("medcase_wide_10rep.jpeg", width = 6, height = 4,
#        units = "in")

# Transmission Rates: l2 wide sweep -------------
dis.wide <- read.csv("disease_test_l2wide.csv") %>%
  select(rep, year, week, total_pop, n_infected, 
         n_symptomatic, elim, l2) %>%
  mutate(nweek = ((year-1)*52)+week)

prop_eliminated <- dis.wide %>%
  filter(year >= 2, elim == "True") %>%
  select(rep, l2) %>%
  group_by(l2) %>%
  distinct() %>%
  summarise(prop = n()/10)

unique(prop_eliminated$prop)

# Time to elimination
time_to_elim <- dis.wide %>%
  group_by(l2, rep) %>%
  filter(year >= 2, elim == "True") %>%
  filter(nweek == min(nweek)) %>%
  distinct() %>%
  mutate(l2 = factor(l2))

time_to_elim %>%
  ungroup() %>%
  group_by(l2) %>%
  summarise(mean=median(nweek))

ggplot(data = time_to_elim, aes(x = l2, y = nweek))+
  geom_boxplot(fill = "lightgray")+
  labs(x = "Transmission Rate", y = "Week of Elimination")+
  theme_bw()+
  theme(panel.grid = element_blank())

# ggsave("weekelim_wide_l2.jpeg", width = 6, height = 4,
#        units = "in")

mean_cases <- dis.wide %>%
  filter(nweek > 70, elim == "False") %>%
  mutate(l2 = factor(l2, levels = sort(unique(l2)))) %>%
  group_by(rep, l2) %>%
  summarise(mean.cases = median(n_symptomatic))

mean_cases %>%
  group_by(l2) %>%
  summarise(median = median(mean.cases))

ggplot(data=mean_cases, aes(x=l2, y = mean.cases))+
  geom_boxplot(fill = 'lightgray')+
  labs(x = "Transmission Rate", y = "Median Weekly Cases")+
  theme_bw(base_size = 14)+
  theme(panel.grid = element_blank())

# ggsave("medcase_wide_l2.jpeg", width = 6, height = 4,
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
                  method = 'TD', nsim = 1000)
      }
    
      vec <- c(unique(test$rep), unique(test$l2),
              median(r0$estimates$TD$R))
             
      len <- length(r0.list)
      r0.list[[len+1]] <- vec
     }
   # }
}

r0.df <- as.data.frame(do.call(rbind, r0.list))
colnames(r0.df) <- c("rep", "l1", "Re")

r0.df <- r0.df %>%
  group_by(rep, l1) %>%
  filter(Re < quantile(.$Re, 0.975) && Re > quantile(.$Re, 0.025))

r0.df %>%
  ungroup() %>%
  group_by(l1) %>%
  summarise(mean.re = median(Re))

ggplot(data = r0.df, aes(x = factor(l1), y = Re))+
  geom_boxplot(fill='lightgray')+
  labs(x = "Core Transmission", y = bquote(R[e]))+
  theme_bw(base_size = 14)+
  theme(panel.grid=element_blank())

# ggsave(filename="R_e_sensitivity.jpeg", width = 5, height = 3.5,
#        units = 'in')

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
