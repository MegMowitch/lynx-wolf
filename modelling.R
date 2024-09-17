library(dplyr)
library(vegan)
library(readr)
library(tidyverse)
library(car)
library(moments)
library(e1071)

#HELLO! Welcome to the modelling script. :) The code is not supposed to be 
#used A to Z, just browse and scan, and choose what you need.
#For me the order was standarization -> Kruskal Wallis -> modelling

setwd("your-working-directory")
read_csv("your-csv-with-covariates")
species_count <- read_csv("your-csv-with-trapping-rates")

#displaying the data to make sure it's ok

head(species_count)

#in my project I wanted to check what impacts the lynx/wolf trapping rates
# (canlup_trap and lynx_trap)
#I chose landscape/human-presence-associated covariates like: dist-fence
#(distance to the border fence) or dist_building (distance to the nearest
#building based on the data acquired with QGIS)

# Data standarization

data_BF <- species_count %>%
  mutate(
    human_trap_std = scale(human_trap),
    dist_fence_std = scale(dist_fence),
    dist_road_std = scale(dist_road),
    dist_paved_std = scale(dist_paved),
    dist_building_std = scale(dist_building),
    dist_narewka_std = scale(dist_narewka),
    dist_bialowieza_std = scale(dist_bialowieza),
    dist_hajn_std = scale(dist_hajn),
    vulp_trap_std = scale(vulp_trap),
    canlup_trap_std=scale(canlup_tr),
    lynx_trap_std=scale(lynx_trap)
  )

head(species_count)

# Correlations between variables

## 1. graphic analysis of the canlup_tr(human_trap) 

ggplot(species_count, aes(x = human_trap)) +
  geom_point(aes(y = canlup_tr), color = "blue") +
  geom_point(aes(y = lynx_trap), color = "red") +
  labs(title = "Trapping rate of wolf and lynx vs. human trapping rate",
       x = "dist_fence", y = "trapping rate") +
  theme_minimal()

ggplot(species_count, aes(x = dist_fence)) +
  geom_point(aes(y = canlup_tr), color = "blue") +
  geom_point(aes(y = lynx_trap), color = "red") +
  labs(title = "Trapping rate of wolf and lynx vs. distance to fence",
       x = "dist_fence", y = "trapping rate") +
  theme_minimal()

species_count %>% 
  filter(human_trap <10) %>% 
  ggplot(aes(x= lynx_trap, y = canlup_tr)) +
  geom_smooth()+
  geom_smooth(method = "lm")+
  geom_jitter()

species_count %>% 
  ggplot(aes(x= dist_building, y = canlup_tr)) +
  geom_smooth()+
  geom_smooth(method = "lm")+
  geom_jitter()

species_count %>% 
  ggplot(aes(x= dist_hajn, y = canlup_tr)) +
  geom_smooth()+
  geom_smooth(method = "lm")+
  geom_jitter()

species_count %>% 
  filter(human_trap <10) %>% 
  ggplot(aes(x= human_trap, y = lynx_trap)) +
  geom_smooth()+
  geom_smooth(method = "lm")+
  geom_jitter()

species_count %>% 
  ggplot(aes(x= dist_building, y = lynx_trap)) +
  geom_smooth()+
  geom_smooth(method = "lm")+
  geom_jitter()

species_count %>% 
  ggplot(aes(x= dist_hajn, y = lynx_trap)) +
  geom_smooth()+
  geom_smooth(method = "lm")+
  geom_jitter()


## 2. correlation matrix

cor_matrix <- cor(data_BF %>% select(human_trap_std, dist_fence_std, 
              dist_road_std, dist_paved_std, dist_building_std,
              dist_narewka_std, dist_bialowieza_std, dist_hajn_std, 
              lynx_trap_std, canlup_trap_std
))

print(cor_matrix)

# multiple linear regression y=a1x1 + a2x2...+b  

model_wolf <- lm(canlup_tr ~ human_trap_std + dist_fence_std + dist_road_std
                 + dist_paved_std + dist_building_std + dist_narewka_std
                 + dist_bialowieza_std + dist_hajn_std + lynx_trap_std, data = data_BF)

summary(model_wolf) 

plot(model_wolf)

vif_values_model_wolf <- vif(model_wolf)
print(vif_values_model_wolf)
#there are some colinear variables = too bad :(
#ANYWAY WOLF SEEMS TO BE TOTALLY INDIFFERENT TO HUMANS
#WHICH IS NOT BREAKING NEWS.. TEMPORAL NICHE PARTITIONING
#IS GOOD ENOUGH, SPATIAL VARIABLES DON'T BOTHER THE WOLF


# model B for trying different variables - TEST THEM, it's a 
# trial-and-error version

model_wolf_B <- lm(canlup_tr ~ human_trap_std + dist_building_std
                   + dist_hajn_std + lynx_trap_std, data = data_BF)

summary(model_wolf_B)

plot(model_wolf_B)

vif_values_wolfB <- vif(model_wolf_B)
print(vif_values_wolfB)

# LYNXXX TIME same - multiple regression for lynx 

model_lynx <- lm(lynx_trap ~ human_trap_std + dist_fence_std + dist_road_std
                 + dist_paved_std + dist_building_std + dist_narewka_std
                 + dist_bialowieza_std + dist_hajn_std, data = data_BF)
summary(model_lynx)

plot(model_lynx)

vif_values_lynx <- vif(model_lynx)
print(vif_values_lynx)
#dist_fence and dist_narewka over 5
#I leave the code for all var here, although the results point out it's not 
#the best variables' collection

# model B for trying different variables - trial and error version

model_lynx_B <- lm(lynx_trap ~ human_trap_std + canlup_trap_std +
                     + dist_building_std + dist_hajn_std, data = data_BF)
summary(model_lynx_B)

plot(model_lynx_B)

vif_values_lynxB <- vif(model_lynx_B)
print(vif_values_lynxB)

# THIS IS A valid model :))) PERFECT MODEL :D 
#LYNX DOES NOT LIKE HUMANS, BUT LIKES BUILDINGS VICCINITY FOR BETTER AMBUSH
#OF COURSE, IT'S A CAT... 
model_lynx_C <- lm(lynx_trap ~  canlup_trap_std
                + dist_building_std + vulp_trap_std, data = data_BF)
summary(model_lynx_C)

plot(model_lynx_C)

#checking VIF

vif_values_lynxC <- vif(model_lynx_C)
print(vif_values_lynxC)

#VIF is excellent :)) 

# NOW checking out the linear mixed effects models
library(lme4)

length(unique(data_BF$site))
nrow(data_BF)

# categorizing the data

summary(data_BF$human_trap)
n_categories <- 5 #5 categories :) 

data_BF$human_trap_category <- cut(data_BF$human_trap, 
                                   breaks = quantile(data_BF$human_trap, 
              probs = seq(0, 1, length.out = n_categories + 1), na.rm = TRUE),
                                   include.lowest = TRUE, 
                                   labels = FALSE)

head(data_BF)

table(data_BF$human_trap_category)
write.csv(data_BF, "species_count_model.csv", row.names = FALSE)
#It would be a pity to lose all this... 

# THE MODEL finally

if(!require(nlme)) install.packages("nlme")

library(nlme)

mixed_lynx_model <- lme(fixed = lynx_trap ~ dist_building_std + canlup_trap_std
                        + dist_hajn_std  + human_trap_std, 
                        random = ~ 1 | human_trap_category, data = data_BF)
summary(mixed_lynx_model)

# Predictions for groups
predict(mixed_lynx_model, level = 0:1)

#for wolf
library(nlme)

mixed_wolf_model <- lme(fixed = canlup_tr ~ dist_building_std + lynx_trap_std
                        + dist_hajn_std  + human_trap_std, 
                        random = ~ 1 | human_trap_category, data = data_BF)
summary(mixed_wolf_model)

# Predictions for groups
predict(mixed_wolf_model, level = 0:1)

#KRUSKAL WALLIS TEST :) just to check whether those trapping rates are 
#in fact different

# skewness for lynx trapping rate

skewness(data_BF$lynx_trap)

#acceptable skewness

skewness(data_BF$canlup_tr)

#skewness high for the wolf data :((( but will check the test anyway

data_BF$lynx_trap <- as.numeric(data_BF$lynx_trap)
data_BF$canlup_tr <- as.numeric(data_BF$canlup_tr)

n <- nrow(data_BF)

# checking that length
length_lynx <- length(data_BF$lynx_trap)
length_wolf <- length(data_BF$canlup_tr)

print(paste("Length of lynx_trap:", length_lynx))
print(paste("Length of canlup_tr:", length_wolf))

if (length_lynx != length_wolf) {
  stop("Lenghts of lynx_trap i canlup_tr are different.")
}

# grouping and creating a vector
trapping_rate <- c(data_BF$lynx_trap, data_BF$canlup_tr)
species <- factor(rep(c("lynx", "wolf"), each = length_lynx))

# are we good length-wise?
print(paste("Trapping rate length:", length(trapping_rate)))
print(paste("Sepcies length:", length(species)))

# variance check of trapping between species
levene_test <- leveneTest(trapping_rate ~ species)
print(levene_test)

# Kruskal-Wallis

lynx_data <- unlist(data["lynx_trap"])
wolf_data <- unlist(data["canlup_tr"])

kruskal_test <- kruskal.test(list(lynx = data_BF$lynx_trap, 
                                  wolf = data_BF$canlup_tr))
print(kruskal_test)

# the trapping rate is DIFFERENT ! Also breaking news :))) 
