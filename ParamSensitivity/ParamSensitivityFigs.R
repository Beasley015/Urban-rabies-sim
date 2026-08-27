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
  row1 <- as.vector(row[1])[[1]]
  row2 <- as.vector(row[2])[[1]]
  row3 <- as.vector(row[3])[[1]]
  
  test <- filter(best.season, maxK==row[1], a_mort==row[2], j_mort==row[3])
  
  dens_fig_list[[i]] <- ggplot(data=test, aes(x = season, y = mean_dens))+
    geom_boxplot(fill = 'lightgray')+
    geom_segment(aes(x = 0.65, xend = 1.35, y = 11, yend = 11),
                 linewidth = 1.5, linetype = "dashed")+
    geom_segment(aes(x = 1.65, xend = 2.35, y = 15, yend = 15),
                 linewidth = 1.5, linetype = "dashed")+
    labs(x = "Season", y = "Mean Density",
         title = bquote(expr = ~ K[max] == .(row1) ~ ", Adult Mortality =" ~
                          .(row2) ~ ", Juvenile Mortality =" ~ .(row3),
                        where = globalenv()))+
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
  row1 <- as.vector(row[1])[[1]]
  row2 <- as.vector(row[2])[[1]]
  row3 <- as.vector(row[3])[[1]]
  
  test <- filter(kmort.best, maxK==row[1], a_mort==row[2], j_mort==row[3])
  
  pop_fig_list[[i]] <- ggplot(data=test, aes(x = nweek, y = total_pop))+
    stat_summary(geom = "ribbon", fun.data = mean_cl_normal, fill = 'lightgray')+
    stat_summary(geom = "line", fun = mean)+
    labs(x = "Week", y = "Population Size",
         title = bquote(expr = ~ K[max] == .(row1) ~ ", Adult Mortality =" ~
                          .(row2) ~ ", Juvenile Mortality =" ~ .(row3),
                        where = globalenv()))+
    theme_bw(base_size = 12)+
    theme(panel.grid=element_blank()) 
}

pop_fig_list[[3]]/pop_fig_list[[5]]/pop_fig_list[[6]]

# ggsave(filename="K_sens.jpeg", width = 14, height = 7, units = "in")

# Kmax: Variance tests --------------
rep5 <- vector("list", 20)
rep5 <- lapply(rep5, function(x) x <- sample(unique(seasonal$rep), 
                                     size = 5, replace = F))

rep5.var <- tibble()
for(i in 1:20){
  subst <- seasonal %>%
    filter(rep %in% rep5[[i]],
           year == 5) %>%
    group_by(a_mort, j_mort, season) %>%
    summarise(vnce = var(total_pop)) %>%
    mutate(num.rep = 5)
  
  rep5.var <- bind_rows(rep5.var, subst)
}

rep10 <- list()

rep15 <- list()

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
  summarise(prop = n()/20)
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

# all_combos <- expand_grid(unique(time_to_elim$lambda1), 
#                           unique(time_to_elim$rep))
# colnames(all_combos) <- c("lambda1", "rep")

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
  labs(x = bquote("Transmission Rate (" ~ lambda  [1] ~ ")"), 
       y = "Week of Elimination")+
  theme_bw(base_size = 14)+
  theme(panel.grid = element_blank())

ggsave("weekelim_wide.jpeg", width = 10, height = 6,
       units = "in", dpi = 300)

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
  labs(x = bquote("Transmission Rate (" ~ lambda[1] ~ ")"), 
       y = "Median Weekly Cases")+
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
  labs(x = bquote("Transmission Rate (" ~ lambda[1] ~ ")"),
       y = expression(R[e]))+
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
  summarise(prop = n()/20)
# Elimination ranges from 20% to 70%

unique(prop_eliminated$prop)

# Time to elimination
time_to_elim <- dis.wide %>%
  group_by(lambda2, rep) %>%
  filter(year >= 2, elim == "True") %>%
  filter(nweek == min(nweek)) %>%
  distinct() %>%
  mutate(lambda2 = factor(lambda2))

# all_combos <- expand_grid(unique(time_to_elim$lambda2), 
#                           unique(time_to_elim$rep))
# colnames(all_combos) <- c("lambda2", "rep")

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
  labs(x = bquote("Transmission Rate (" ~ lambda[2] ~ ")"), 
       y = "Week of Elimination")+
  theme_bw(base_size = 12)+
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
  labs(x = bquote("Transmission Rate (" ~ lambda[2] ~ ")"), y = "Median Weekly Cases")+
  theme_bw(base_size = 12)+
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
  labs(x = bquote("Transmission Rate (" ~ lambda[2] ~ ")"), y = expression(R[e]))+
  theme_bw(base_size = 12)+
  theme(panel.grid = element_blank())

# ggsave("re_l2_wide.jpeg", width = 8, height = 6, units = "in")

# Glance at median Re
re.df %>%
  group_by(lambda2) %>%
  summarise(re = median(Re))
# 0.004ish looks like the winner

# Transmission Rates: full sweep -----------------------
# Combine files 
# filenames <- list.files(pattern = "l1\\d+")
# 
# file.list <- list()
# for(i in 1:length(filenames)){
#   file.list[[i]] <- read.csv(filenames[i])
# }
# 
# lambda.frame <- do.call(rbind, file.list)
# 
# write.csv(lambda.frame, "lambda_full.csv")
# file.remove(filenames)

# Read in data
dis.wide <- read.csv("lambda_full.csv") %>%
  select(rep, year, week, total_pop, n_infected, 
         n_symptomatic, elim, lambda1, lambda2) %>%
  mutate(nweek = ((year-1)*52)+week)

# proportion eliminated
prop_eliminated <- dis.wide %>%
  filter(year >= 2, elim == "True") %>%
  select(rep,lambda1,lambda2) %>%
  group_by(lambda1,lambda2) %>%
  distinct() %>%
  summarise(prop = n()/20)
# Elimination ranges from 30% to 70%

unique(prop_eliminated$prop)

ggplot(data=prop_eliminated, aes(x=factor(lambda1), y=factor(lambda2),
                                 fill=prop))+
  geom_tile()+
  scale_fill_viridis(name = "Proportion Eliminated")+
  labs(x = "Lambda 1", y = "Lambda 2")+
  theme_bw(base_size=12)
# no clear pattern

# ggsave("full_lambda_elim.jpeg", height = 5, width = 7, units = "in")

# Time to elimination
time_to_elim <- dis.wide %>%
  group_by(lambda1, lambda2, rep) %>%
  filter(year >= 2, elim == "True") %>%
  filter(nweek == min(nweek)) %>%
  distinct() %>%
  mutate(lambda2 = factor(lambda2), lambda1=factor(lambda1))

# all_combos <- expand_grid(unique(time_to_elim$lambda1), 
#                           unique(time_to_elim$lambda2), 
#                           unique(time_to_elim$rep))
# colnames(all_combos) <- c("lambda1", "lambda2", "rep")

time_to_elim <- time_to_elim %>%
  # right_join(all_combos, by=c("rep", "lambda1", "lambda2")) %>%
  mutate(nweek = case_when(is.na(nweek)==T ~ 52*11,
                           TRUE ~ nweek)) %>%
  group_by(lambda1, lambda2) %>%
  summarise(meantime = median(nweek))

ggplot(data = time_to_elim, aes(x = lambda1, y = lambda2,
                                fill=meantime))+
  geom_tile()+
  scale_fill_viridis(name = "Week of Elimination")+
  labs(x = expression(lambda[1]), y = expression(lambda[2]))+
  theme_bw(base_size = 12)+
  theme(panel.grid = element_blank())
# Top left corner eliminates very quickly (low l1 and high l2)

# Start filtering values
poss.combos <- time_to_elim %>%
  filter(meantime > 52*7) # want rabies to last at least half of the sim

# ggsave("TimeToElim_full.jpeg", width = 8, height = 6,
#        units = "in")

# Weekly cases
mean_cases <- dis.wide %>%
  filter(nweek > 70, elim == "False") %>%
  mutate(lambda1 = factor(lambda1, levels = sort(unique(lambda1))), 
         lambda2 = factor(lambda2, levels = sort(unique(lambda2)))) %>%
  group_by(rep, lambda1, lambda2) %>%
  summarise(mean.cases = median(n_symptomatic)) %>%
  ungroup() %>%
  group_by(lambda1, lambda2) %>%
  summarise(meancases = median(mean.cases))

ggplot(data=mean_cases, aes(x=lambda1, y = lambda2,
                            fill = meancases))+
  geom_tile()+
  scale_fill_viridis(name = "Median Weekly Cases")+
  labs(x = expression(lambda[1]), y = expression(lambda[2]))+
  theme_bw(base_size = 12)+
  theme(panel.grid = element_blank())
# Highest l2 generally has too many

poss.combos <- poss.combos %>%
  left_join(mean_cases, by = c("lambda1", "lambda2")) %>%
  filter(meancases < 10) # remove cases that are too high

# ggsave("medcase_full.jpeg", width = 8, height = 6,
#        units = "in")

# Re calculation
dis <- dis.wide %>%
  filter(year > 1)

first_elim <- dis %>%
  filter(elim == "True") %>%
  group_by(rep, lambda1, lambda2) %>%
  summarise(first = min(nweek))

r0.list <- list()
for(i in 1:length(unique(dis$rep))){
  for(j in 1:length(unique(dis$lambda1))){
    for(k in 1:length(unique(dis$lambda2))){
      test <- filter(dis, rep==i & lambda1==unique(dis$lambda1)[j] & lambda2==unique(dis$lambda2)[k]) 
      elim_test <- filter(first_elim, 
                        rep==i & lambda1==unique(dis$lambda1)[j] & lambda2==unique(dis$lambda2)[k]) 
    
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
    
      vec <- c(unique(test$rep), unique(test$lambda1), unique(test$lambda2), median(r0$estimates$TD$R))
    
      len <- length(r0.list)
      r0.list[[len+1]] <- vec
    }
  }
}

r0.df <- as.data.frame(do.call(rbind, r0.list))
colnames(r0.df) <- c("rep", "lambda1", "lambda2", "Re")

r0.df <- r0.df %>%
  group_by(lambda1, lambda2) #%>%
  # filter(Re < quantile(Re, 0.975) & Re > quantile(Re, 0.025))

re.means <- r0.df %>%
  ungroup() %>%
  group_by(lambda1, lambda2) %>%
  summarise(mean.re = median(Re)) %>%
  mutate(lambda1=factor(lambda1), lambda2=factor(lambda2))

poss.combos <- poss.combos %>%
  left_join(re.means, by=c("lambda1", "lambda2")) %>%
  filter(mean.re < 1.3)

ggplot(data=re.means, aes(x=lambda1, y = lambda2, fill = mean.re))+
  geom_tile()+
  scale_fill_viridis(name=expression(R[e]))+
  labs(x = expression(lambda[1]), y = expression(lambda[2]))+
  theme_bw(base_size=12)

# ggsave("re_full.jpeg", width = 8, height = 6,
#        units = "in")

re.tab <-filter(re.means, mean.re > 1 & mean.re < 1.25)
# write.table(re.tab, file="re_table.csv", sep = ",")

# write.table(poss.combos, "poss.combos.csv")

# Population sizes with disease ------------------
dis_pop <- dis %>%
  group_by(nweek,lambda1, lambda2) %>%
  mutate(lambda1=factor(lambda1), lambda2=factor(lambda2)) %>%
  filter(lambda1 %in% c(0.015, 0.025) & lambda2 %in% c(0.002, 0.005)) %>%
  filter((lambda1 == 0.015 & lambda2 == 0.002) | (lambda1 == 0.025 & lambda2 == 0.005))

ggplot(data = dis_pop, aes (x = nweek, y = total_pop,
                            color = lambda1))+
  stat_summary(geom = "line", fun = mean)+
  geom_vline(xintercept = 53, linetype = 'dashed')+
  scale_color_viridis_d(name = "Within-cell transmission",
                        end = 0.9)+
  labs(x="Week", y = "Mean Population Size")+
  theme_bw(base_size = 14)+
  theme(panel.grid = element_blank())

# ggsave(filename = "pop_change.jpeg", width = 8, height = 6,
#        units = "in")

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

# Starting cases ------------------
# Read it in
stcase <- read.csv("starting_cases.csv") %>%
  filter(year > 1) %>%
  select(rep, year, week, starting_cases, total_pop, n_infected,
         n_symptomatic, elim) %>%
  mutate(nweek = ((year-1)*52)+week)

# Differences in elimination probability?
stcase.elim <- stcase %>%
  filter(elim == "True") %>%
  select(rep,starting_cases) %>%
  group_by(starting_cases) %>%
  distinct() %>%
  summarise(prop = n()/20)
# Not really- 40-55% probability with no clear pattern

# Time to elimination
stcase.time <- stcase %>%
  group_by(starting_cases, rep) %>%
  filter(elim == "True") %>%
  filter(nweek == min(nweek)) %>%
  distinct() %>%
  mutate(starting_cases=factor(starting_cases))

ggplot(stcase.time, aes(x = starting_cases, y = nweek))+
  geom_boxplot(fill = 'lightgray')+
  labs(x = "Starting Cases", y = "Week of Elimination")+
  theme_bw(base_size = 12)+
  theme(panel.grid = element_blank())

summary(aov(nweek~starting_cases, data=stcase.time))

# Weekly cases, max and median
stcase.cases <- stcase %>%
  select(starting_cases, rep, n_symptomatic) %>%
  mutate(starting_cases=factor(starting_cases)) %>%
  group_by(starting_cases, rep) %>%
  summarise(medcase = median(n_symptomatic), maxcase = max(n_symptomatic))
  
ggplot(stcase.cases, aes(x = starting_cases, y=medcase))+
  geom_boxplot(fill = 'lightgray')+
  labs(x = "Starting Cases", y = "Median Weekly Cases")+
  theme_bw(base_size = 12)+
  theme(panel.grid = element_blank())

ggplot(stcase.cases, aes(x = starting_cases, y=maxcase))+
  geom_boxplot(fill = 'lightgray')+
  labs(x = "Starting Cases", y = "Maximum Weekly Cases")+
  theme_bw(base_size = 12)+
  theme(panel.grid = element_blank())

maxcase.mod <- aov(maxcase~starting_cases, data=stcase.cases)
TukeyHSD(maxcase.mod)
# 10 & 5 are equivalent; 20 and 40 lead to increase

maxcase.week <- stcase %>%
  select(starting_cases, rep, n_symptomatic, nweek) %>%
  mutate(starting_cases=factor(starting_cases)) %>%
  group_by(starting_cases, rep) %>%
  filter(n_symptomatic==max(n_symptomatic))

ggplot(data = maxcase.week, aes(x = starting_cases, y = nweek))+
  geom_boxplot(fill='lightgray')+
  labs(x = "Starting Cases", y = "Week of Max Cases")+
  theme_bw(base_size=12)+
  theme(panel.grid = element_blank())

# Re
first_elim <- stcase %>%
  filter(year >= 2) %>%
  filter(elim == "True") %>%
  group_by(rep, starting_cases) %>%
  summarise(first = min(nweek))

r.list <- list()
for(i in 1:length(unique(stcase$rep))){
  for(j in 1:length(unique(stcase$starting_cases))){
    test <- filter(stcase, rep==i & 
                     starting_cases==unique(stcase$starting_cases)[j]) %>%
      filter(year >= 2)
    
    elim_test <- filter(first_elim, 
                        rep==i & 
                          starting_cases==unique(stcase$starting_cases)[j]) 
    
    start <- as.numeric(min(which(test$n_symptomatic>0)))
    end <- ifelse(nrow(elim_test) != 1,
                  52*11,
                  elim_test$first-53)
    
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
    
    vec <- c(unique(test$rep), unique(test$starting_cases),
             median(re$estimates$TD$R))
    
    len <- length(r.list)
    r.list[[len+1]] <- vec
  }
}

re.df <- as.data.frame(do.call(rbind, r.list))
colnames(re.df) <- c("rep", "starting_cases", "Re")

ggplot(data = re.df, aes(x = factor(starting_cases), y = Re))+
  geom_boxplot(fill = 'lightgray')+
  geom_hline(yintercept = 1.13, linetype = "dashed",
             linewidth = 1)+
  labs(x = "Starting Cases",
       y = expression(R[e]))+
  theme_bw(base_size = 12)+
  theme(panel.grid = element_blank())

summary(aov(Re~starting_cases, data=re.df))
re.mod <- aov(Re~factor(starting_cases), data=re.df) 
TukeyHSD(re.mod)
# it's a gradient, slowly increasing

# Landscape size --------------------
# Read in outputs
landsims <- read.csv("land_size.csv") %>%
  filter(year > 1) %>%
  select(rep, year, week, land_size, total_pop, n_infected,
         n_symptomatic, elim) %>%
  mutate(nweek = ((year-1)*52)+week)

# Differences in elimination probability
landsims.elim <- landsims %>%
  filter(elim == "True") %>%
  select(rep,land_size) %>%
  group_by(land_size) %>%
  distinct() %>%
  summarise(prop = n()/20)
# Strong decreases in elim probability until 80x80

# Prevalence
land.prev <- landsims %>%
  filter(elim == "False") %>%
  select(rep, land_size, nweek, n_symptomatic, total_pop) %>%
  mutate(prev = n_symptomatic/total_pop) %>%
  mutate(land_size = case_when(land_size == 40 ~ "40x40",
                               land_size == 60 ~ "60x60",
                               land_size == 80 ~ "80x80",
                               TRUE ~ "100x100")) %>%
  mutate(land_size = factor(land_size, 
                            levels = c("40x40", "60x60", "80x80", 
                                       "100x100"))) %>%
  group_by(rep, land_size) %>%
  summarise(med_prev = median(prev))

ggplot(land.prev, aes(x = factor(land_size), y = med_prev))+
  geom_boxplot(fill = 'lightgray') +
  labs(x = "Landscape size (# cells)", y = "Median Rabies Prevalence")+
  theme_bw(base_size = 12)+
  theme(panel.grid = element_blank())
  
