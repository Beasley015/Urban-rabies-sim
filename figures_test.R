library(tidyverse)

trans5 <- read.csv("test_5.csv") %>%
  mutate(transmission=5)

trans3 <- read.csv("test_3.csv") %>%
  mutate(transmission=3)

trans2.5 <- read.csv("test_25.csv") %>%
  mutate(transmission=2.5)

trans2.75 <- read.csv("test_275.csv") %>%
  mutate(transmission=2.75)

trans2 <- read.csv("test_2.csv") %>%
  mutate(transmission=2)

alltest <- rbind(trans5, trans3, trans2.75, trans2.5, trans2) %>%
  mutate(sero=factor(sero), transmission=factor(transmission)) %>%
  mutate(nweek = ((year-1)*52) + week)
# year 4 week 1 = week 157

# eliminated
elim <- alltest %>%
  filter(year >= 4) %>%
  filter(n_infected == 0 & n_symptomatic == 0) %>%
  group_by(sero, transmission) %>%
  select(sero, transmission, rep) %>%
  distinct() %>%
  summarise(count = n())

ggplot(data=elim, aes(x=sero, y=count, fill=transmission))+
  geom_col(position="dodge")

# Weeks to elimination
elim <- alltest %>%
  filter(year >= 4) %>%
  filter(n_infected == 0 & n_symptomatic == 0) %>%
  group_by(sero, transmission, rep) %>%
  filter(nweek == min(nweek))

ggplot(data=elim, aes(x=sero, y=nweek, fill=transmission))+
  geom_boxplot()

# Mean cases per week
perweek <- alltest %>%
  group_by(sero, transmission, nweek) %>%
  summarise(mean_cases = mean(n_symptomatic))

ggplot(data=perweek[perweek$sero==0,], aes(x=nweek, y=mean_cases,
                                           color=transmission))+
  geom_line()

perweek_3 <- alltest %>%
  filter(transmission == 3 & sero %in% c(0, 0.2, 0.4, 0.6, 0.8)) %>%
  group_by(sero, nweek) %>%
  summarise(mean_cases = mean(n_symptomatic))

ggplot(data=perweek_3, aes(x=nweek, y = mean_cases, color=sero))+
  geom_line()

perweek_3_4 <- alltest %>%
  filter(transmission == 3 & sero==0.6)


ggplot(data=perweek_3_4, aes(x=nweek, y = n_symptomatic, 
                             color = factor(rep)))+
  geom_line()


perweek_275 <- alltest %>%
  filter(transmission == 2.75 & sero %in% c(0, 0.1, 0.3, 0.5, 0.7)) %>%
  group_by(sero, nweek) %>%
  summarise(mean_cases = mean(n_symptomatic))

ggplot(data=perweek_275, aes(x=nweek, y = mean_cases, color=sero))+
  geom_line()

perweek_275_4 <- alltest %>%
  filter(transmission == 2.75 & sero==0.4)


ggplot(data=perweek_275_4, aes(x=nweek, y = n_symptomatic, 
                             color = factor(rep)))+
  geom_line()

#Population size
weekpop <- alltest %>%
  group_by(sero, transmission, nweek) %>%
  summarise(mean_pop = mean(total_pop))

ggplot(data=weekpop[weekpop$sero==0,], aes(x=nweek, y=mean_pop,
                                             color=transmission))+
  geom_line()+
  labs(x = "Week", y = "Population Size")+
  theme_bw(base_size=14)+
  theme(panel.grid=element_blank())

ggsave(filename = "testpop.jpeg", width = 8, height = 6,
       units = "in")

# Full sims, 3% transmission & reduced repro -----------
full3 <- read.csv("test_full3.csv") %>%
  mutate(nweek = ((year-1)*52) + week) %>%
  mutate(trans="3%")

full35 <- read.csv("test_full35.csv") %>%
  mutate(nweek = ((year-1)*52) + week) %>%
  mutate(trans="3.5%")

full4 <- read.csv("test_full4.csv") %>%
  mutate(nweek = ((year-1)*52) + week) %>%
  mutate(trans="4%")

all_full <-rbind(full3, full35, full4)

nelim <- all_full %>%
  filter(year >= 4) %>%
  filter(n_infected == 0 & n_symptomatic == 0) %>%
  select(rep, trans) %>%
  distinct() %>%
  group_by(trans) %>%
  summarise(perc.elim = n()/50)

nelim

# weeks to elimination
nweek0 <- all_full %>%
  filter(year >= 4) %>%
  filter(n_infected == 0 & n_symptomatic == 0) %>%
  group_by(rep, trans) %>%
  filter(nweek == min(nweek))

ggplot(data=nweek0, aes(x=factor(trans), y = nweek))+
  geom_boxplot()

# cases per week
cases0 <- all_full %>%
  group_by(nweek, trans) %>%
  summarise(mean_cases = mean(n_symptomatic))

ggplot(data=cases0, aes(x = nweek, y = mean_cases, 
                        color=factor(trans)))+
  geom_line()

# population size
pop0 <- all_full %>%
  group_by(nweek, trans) %>%
  summarise(mean_pop = mean(total_pop))

ggplot(data=pop0, aes(x = nweek, y = mean_pop, color=factor(trans)))+
  geom_line()

