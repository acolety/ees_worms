# Earthworm abundance and diversity
# EES 4th year field course 2023
# Jess, Caitlin, Simon, Cory, Heather, Imara, Alex
# 11-12 Sept 2023 

# libraries ----
library(skimr)
library(tidyverse)
library(devtools)
devtools::install_github("cmartin/ggConvexHull")
library(ggConvexHull)
library(lme4)
library(MuMIn)
library(MASS)
library(ggpubr)

# colors ----
"#8BB26C" # Anecic
"#225235" # Endogeic
"#45ADC7" # Epigeic
"#A9BFC4" # Juveniles
"#49d423" # mown
"#6224d8" #unmown

# loading data ----
worms <- read.csv("data/earthworm_data.csv")

# viewing and cleaning data ----
view(worms)
skim(worms)

worms_sub <- filter(worms, habitat %in% c("unmowed_grassland",
                                          "mowed_grassland"))  # removing woodland data

# raw data visualization ----
(rough <- ggplot(worms_sub, aes(x = habitat, y = earthworm_count)) +
   geom_boxplot() +
   theme_classic())  # any difference between worm counts in 2 habitats?

(rough2 <- ggplot(worms_sub, aes(x = soil_ph, y = earthworm_count, color = habitat)) +
    geom_point() +
    theme_classic())  # pH impact on earthworm count?

# Data wrangling ----

(worms_sub <- worms_sub %>% 
   mutate(epigeic_prop = epigeic / adults, 
          anecic_prop = anecic / adults,
          endogeic_prop = endogeic / adults,
          adult_prop = adults / earthworm_count, 
          juv_prop = juveniles / earthworm_count))  # adding columns for ecotype 
                                                    # and age proportions

(comm <- worms_sub %>% 
    mutate(nonepi_prop = anecic_prop + endogeic_prop, 
           nonepi = anecic + endogeic) %>%  # creating data subset with anecic 
                                            # and endogeic worms grouped
    dplyr::select(-c(date, soil_t1, soil_t2, soil_t3, soil_t4, soil_m1, soil_m2, soil_m3,
              soil_m4, anecic_prop, endogeic_prop, adult_prop, juv_prop)))

(comm_long <- comm %>% 
    pivot_longer(cols = c('epigeic', 'endogeic', 'anecic', 'juveniles'),
                 names_to = "ecotype", 
                 values_to = "count"))  # changing data to long format


## checking distributions of ecotypes
### checking proportion distributions
hist(worms_sub$epigeic_prop, breaks = 30)
hist(worms_sub$anecic_prop, breaks = 30)
hist(worms_sub$endogeic_prop, breaks = 30)  # looks normal enough 


# modelling ----

## modelling worm abundance as impacted by habitat

abund_mod <- glm(earthworm_count ~ habitat, data = worms_sub, family = poisson)
summary(abund_mod)
plot(abund_mod)  # significant difference between worm abundance in 2 habitats, due 
                 # to juvenile abundance differences
139.31/28 # 4.97, so very overdispersed

## using negative binomial instead
abund_nb <- glm.nb(earthworm_count ~ habitat, data = worms_sub)
summary(abund_nb)  # no difference in abundance between habitats
plot(abund_nb)
30.69/28 # 1.09, results with trepidation but not super overdispersed


## community structure questions

### juveniles 

juv_hab_mod <- glm(juveniles ~ habitat, data = worms_sub, family = poisson)
summary(juv_hab_mod)
plot(huv_hab_mod)  # checking just abundance for juveniles
120.97/28  # really high

juv_hab_nb <- glm.nb(juveniles ~ habitat, data = worms_sub)
summary(juv_hab_nb)
plot(juv_hab_nb)
31.571/28 # 1.13

### adult proportions

adult_prop_hab <- lm(adult_prop ~ habitat, data = worms_sub)
summary(adult_prop_hab)
plot(adult_prop_hab)

### epigeic vs nonepi proportions

epi_prop_hab <- lm(epigeic_prop ~ habitat, data = worms_sub)
summary(epi_prop_hab)
plot(epi_prop_hab)

### ecotypes

epi_hab_model <- glm(epigeic ~ habitat, data = worms_sub, family = poisson)
summary(epi_hab_model)  # overdispersed
plot(epi_hab_model)

epi_hab_nb <- glm.nb(epigeic ~ habitat, data = worms_sub)
summary(epi_hab_nb)

non_epi_hab <- glm(nonepi ~ habitat, data = comm, family = poisson)
summary(non_epi_hab)  # overdispersed
plot(non_epi_hab)

non_epi_nb <- glm.nb(nonepi ~ habitat, data = comm)
summary(non_epi_nb)


### looking at environmental variables
#### scaling  moisture and pH data

worms_sub <- worms_sub %>% 
  mutate(scaled_m = scale(soil_m_avg), 
         scaled_ph = scale(soil_ph))  # scaled data still has variance/resid issues

#### boxcox transformation
m_boxcox <- boxcox(lm(soil_m_avg ~ 1, data = worms_sub))
lambda <- m_boxcox$x[which.max(m_boxcox$y)]
lambda

ph_boxcox <- boxcox(lm(soil_ph ~ 1, data = worms_sub))
lambda_ph <- ph_boxcox$x[which.max(ph_boxcox$y)]
lambda_ph

worms_sub <- worms_sub %>% 
  mutate(boxcox_m = (1/sqrt(soil_m_avg)), 
         boxcox_ph = (1/((soil_ph)^2)))  # applying appropriate transformations


### do enviro var differ bt habitats
temp_hab <- aov(soil_t_avg ~ habitat, data = worms_sub)
summary(temp_hab) # no sig diff bt soil temps bt habitats
plot(temp_hab)  # all looks good

moisture_habitat <- aov(boxcox_m ~ habitat, data = worms_sub)
summary(moisture_habitat)  # no sig diff in moisture bt habitats
plot(moisture_habitat)
shapiro.test(resid(moisture_habitat))  # all good
bartlett.test(boxcox_m ~ habitat, data = worms_sub)  # variance not homogeneous, but good w boxcox

ph_habitat <- aov(boxcox_ph ~ habitat, data = worms_sub)
summary(ph_habitat)  # sig diff
plot(moisture_habitat)
shapiro.test(resid(ph_habitat))  # not good, doesn't pass even with transformation
bartlett.test(boxcox_ph ~ habitat, data = worms_sub)  # good 

### trying Mann-Whitney U test for pH instead
ph_habitat_MWU <- wilcox.test(soil_ph ~ habitat, data = worms_sub) 
      # significant, error abt ties is bc of identical values
(ph_habitat_MWU)
qnorm(ph_habitat_MWU$p.value)  # getting z value

#### determining which has higher pH
worms_sub %>% 
  group_by(habitat) %>% 
  summarise(mean(soil_ph)) %>% 
  ungroup()  # unmown has the higher pH
  

#### modelling environmental factors against worms
temp_abund_model <- glm(soil_t_avg ~ earthworm_count, data = worms_sub)
summary(temp_abund_model) # No sig diff bt worm abundances with
plot(temp_abund_model)
shapiro.test(resid(temp_abund_model)) # pass

temp_epigeic <- glm.nb(epigeic ~ soil_t_avg, data = worms_sub)
summary(temp_epigeic)  # significant
plot(temp_epigeic)
shapiro.test(resid(temp_epigeic)) # normal enough

ph_mod <- glm(log(soil_ph) ~ earthworm_count, data = worms_sub)
summary(ph_mod) # soil pH unrelated to general earthworm abund
plot(ph_mod)
shapiro.test(resid(ph_mod))  # initially not normally distrib, logged fixed

soil_moisture <- glm.nb(earthworm_count ~ soil_m_avg, data = worms_sub)
summary(soil_moisture)  # unrelated
plot(soil_moisture)
shapiro.test(resid(soil_moisture)) #good

# plotting community structure by UG vs MG ----

## plotting ecotype proportions by count per plot
(community <- ggplot(comm, aes(x = epigeic, y = nonepi, colour = habitat))     +
   geom_point()                                                               + 
   ggConvexHull::geom_convexhull(alpha = 0.3, aes(fill = habitat))            +
   theme_classic())                                                           +
  xlab("Mean epigeic worms per plot")                       +
  ylab("Mean non-epigeic worms per plot")                   +
  scale_fill_manual(name = "Habitat", values=c("#49d423", "#6224d8"), 
                    labels = c("Mown grassland", "Unmown grassland"))        +
  scale_colour_manual(name = "Habitat", values=c("#49d423", "#6224d8"), 
                      labels = c("Mown grassland", "Unmown grassland"))


# bar chart for community structure by ecotype and abundances ----

### creating smaller dataframe w standard errors
comm_long_errors <- comm_long %>% 
  group_by(habitat, ecotype) %>% 
  summarise(mean = mean(count), se = (sd(count)/sqrt(length(count)))) %>% 
  ungroup()

### bar chart with new data frame
(community_bar_err <- ggplot(comm_long_errors, aes(x = habitat, y = mean, 
                                                   fill = ecotype))            +
    geom_bar(stat = "identity", position = "dodge")                            +
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se, y = mean), 
                  width=0.25, position = position_dodge(0.9))                  +
    theme_classic())                                                           +
  xlab("Habitat")                                                              + 
  ylab("Average number of worms per plot")                                     +
  scale_fill_manual(name = "Ecotype", values = c("#8BB26C", "#225235", 
                                                 "#45ADC7", "#A9BFC4"), 
                    labels = c("Anecic", "Endogeic", "Epigeic", "Juveniles"))  +
  scale_x_discrete(labels = c("Mowed grassland", "Unmowed grassland"))         +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)))

comm_long_errors2 <- comm_long %>% 
  group_by(habitat) %>% 
  mutate(mean = mean(count), total_count = sum(count)) %>% 
  mutate(se = sd(count) / sqrt(length(count))) %>%  # can't get se right
  dplyr::select(c(habitat, plot, ecotype, mean, se, total_count)) %>% 
  ungroup() 

(community_bar_2 <- ggplot(comm_long_errors2, aes(x = habitat, y = mean, 
                                                   fill = ecotype))            +
    geom_bar(stat = "identity", aes(x = habitat, y = mean), width = 0.7)       +
    theme_classic())                                                           +
  xlab("Habitat")                                                              + 
  ylab("Number of worms")                                                      +
  scale_fill_manual(name = "Ecotype", values = c("#8BB26C", "#225235", 
                                                 "#45ADC7", "#A9BFC4"), 
                    labels = c("Anecic", "Endogeic", "Epigeic", "Juveniles"))  +
  scale_x_discrete(labels = c("Mown grassland", "Unmown grassland"))           +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)))


## boxplot to check
(community_box <- ggplot(comm_long, aes(x = habitat, y = count, 
                                        fill = ecotype))                       +
    geom_boxplot()                                                             +
    theme_classic())

# Plotting habitat worm abundance against soil pH ----

(habitat_ph <- ggplot(worms, aes(x = habitat, y = soil_ph, color = habitat))   +
   geom_boxplot()                                                              +
   scale_colour_manual(name = "Habitat", values=c("#225235", "#45ADC7", 
                                                  "#8BB26C"), 
                       labels = c("Mixed Woodland", "Mowed grassland",
                                  "Unmowed grassland"))                        +
   labs(x = "Habitat", y = "Soil pH")                                          +
   scale_x_discrete(labels = c("Mixed woodland", "Mowed grassland", 
                               "Unmowed grassland"))                           +
   theme_classic())

# plotting worm abundance against soil variables ----
(abund_x_pH <- ggplot(worms_sub, aes(x = soil_ph, y = earthworm_count))        +
    geom_point()                                                               +
    theme_classic()                                                            +
    labs(y = "Count per plot", x = "Soil pH", title = "(a) pH"))

(abund_x_moisture <- ggplot(worms_sub, aes(x = soil_m_avg, 
                                           y = earthworm_count))               +
    geom_point()                                                               +
    theme_classic()                                                            +
    labs(y = "Count per plot", x = "Soil moisture (%)", title = "(b) moisture"))

(abund_x_temp <- ggplot(worms_sub, aes(x = soil_t_avg, y = earthworm_count))   +
    geom_point()                                                               +
    theme_classic()                                                            +
    labs(y = "Count per plot", x = "Soil temperature (°C)", 
         title = "(c) overall abundance and temperature"))

(epi_x_temp <- ggplot(worms_sub, aes(x = soil_t_avg, y = epigeic))             +
    geom_point()                                                               +
    geom_smooth(method = "glm.nb", se = F, color = "#45ADC7")                  +
    theme_classic()                                                            +
    labs(y = "Epigeic count per plot", x = "Soil temperature (°C)", 
         title = "(d) epigeic abundance and temperature"))


(abund_x_soil <- ggarrange(abund_x_pH, abund_x_moisture, abund_x_temp, 
                           epi_x_temp, ncol = 2, nrow = 2))

# saving figures ----
ggsave("img/community_struc_bar.png", plot = community_bar_err, width = 7,
       height = 5)
ggsave("img/community_struc_scatter.png", plot = community, width = 7, height = 5)
ggsave("img/habitat_ph.png", plot = habitat_ph, width = 7, height = 5)
ggsave("img/community_bar_2.png", plot = community_bar_2, width = 7, height = 5)
ggsave("img/abund_and_soil.png", plot = abund_x_soil, width = 9, height = 7)
