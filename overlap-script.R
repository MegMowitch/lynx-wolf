# Analysing carnivore community of Bialowieza forest - temporal activity patterns
# Contact: Meg Mowitch, megmowitch@gmail.com

library(readr)
library(lubridate)
library(dplyr)
library(overlap)
library(ggplot2)
library(aod)
library(lmtest)

setwd("your-working-directory")
dat <- read_csv("overlap_animal_site_ready.csv") ## YES! the final file from 
#the database-prep script

#adding 'month' column
dat$month_number <- month(mdy_hms(dat$time))
dat$month <- month.name[dat$month_number]
write.csv(dat, "overlap_animal_site_ready.csv", row.names = FALSE)

# correcting date format
dat$time <- strptime(dat$time, format = "%m/%d/%Y %H:%M:%S")
# transforming into radians
dat$radiany <- (as.numeric(format(dat$time, "%H")) + 
                  as.numeric(format(dat$time, "%M"))/60 + 
                  as.numeric(format(dat$time, "%S"))/3600) / 24 * 2 * pi
head(dat)

# analysis: lynx vs wolf
lynx <- dat[dat$species == "Lynx lynx",]$radiany
wolf <- dat[dat$species == "Canis lupus",]$radiany
apex <- dat[dat$species %in% c("Lynx lynx", "Canis lupus"),]$radiany

densityPlot(lynx, rug=TRUE, main="Lynx temporal activity")
densityPlot(wolf, rug=TRUE, main="Wolf temporal activity")
densityPlot(apex, rug=TRUE)

overlapPlot(lynx, wolf, rug=TRUE, main="Lynx and wolf temporal activity")

legend('topright', legend = c(expression(italic("Lynx lynx")), 
                              expression(italic("Canis lupus"))), 
                              lty = c(1, 2), col = c(1, 4), bty = 'n')
legend('topright', legend = c(expression(italic("Canis lupus"))), 
       lty = c(1, 2), col = c(1, 4), bty = 'n')

########### dla człowieka / for human

human_dat <- read_csv("human_overlap_site_ready.csv")

#adding 'month' column
human_dat$month_number <- month(mdy_hms(human_dat$time))
human_dat$month <- month.name[human_dat$month_number]
write.csv(dat, "human_overlap_site_ready.csv", row.names = FALSE)

#transforming date to dif format and radians
human_dat$time <- strptime(human_dat$time, format = "%m/%d/%Y %H:%M:%S")

human_dat$radiany <- (as.numeric(format(human_dat$time, "%H")) + 
                  as.numeric(format(human_dat$time, "%M"))/60 + 
                  as.numeric(format(human_dat$time, "%S"))/3600) / 24 * 2 * pi
head(human_dat)

# human activity pattern
human <- human_dat[human_dat$human == "Homo sapiens",]$radiany


densityPlot(human, rug=TRUE, main="Human temporal activity")


## overlap plot
overlapPlot(lynx, human, rug=TRUE, main="Lynx and human temporal activity")
overlapPlot(wolf, human, rug=TRUE, main="Wolf and human temporal activity")
overlapPlot(apex, human, rug=TRUE, main="Apex predators and human temporal activity")

#legends
legend('topright', legend = c(expression(italic("Lynx lynx")), 
                              expression(italic("Homo sapiens"))), 
       lty = c(1, 2), col = c(1, 4), bty = 'n')

legend('topright', legend = c(expression(italic("Canis lupus")), 
                              expression(italic("Homo sapiens"))), 
       lty = c(1, 2), col = c(1, 4), bty = 'n')

legend('topright', legend = c("Apex predators", expression(italic("Homo sapiens"))), 
       lty = c(1, 2), col = c(1, 4), bty = 'n')


############ Współczynnik nakładania się (coefficient of overlap)
lynx_wolf_est <- overlapEst(lynx, wolf, type="Dhat1")
lynx_wolf_est
lynx_human_est <- overlapEst(lynx, human, type="Dhat1")
lynx_human_est
wolf_human_est <- overlapEst(wolf, human, type="Dhat4")
wolf_human_est
apex_human_est <- overlapEst(apex, human, type="Dhat4")
apex_human_est

#BOOTSTRAP
lynx_boot <- resample(lynx, 1000)
wolf_boot <- resample(wolf, 1000)
human_boot <- resample(human, 1000)
apex_boot <- resample(apex, 1000)

lynx_wolf_boot <- bootEst(lynx_boot, wolf_boot, type="Dhat1")
(BSmean <- mean(lynx_wolf_boot))
#bootstrapped mean overlap
print(BSmean)
print(lynx_wolf_boot)
#use 'basic0' output as your CI
lynx_wolf_ci <- bootCI(lynx_wolf_est, lynx_wolf_boot, conf = 0.95)

lynx_human_boot <- bootEst(lynx_boot, human_boot, type="Dhat1")
(BSmean <- mean(lynx_human_boot))
#bootstrapped mean overlap
print(BSmean)
print(lynx_human_boot)
#use 'basic0' output as your CI
lynx_human_ci <- bootCI(lynx_human_est, lynx_human_boot, conf = 0.95)

wolf_human_boot <- bootEst(wolf_boot, human_boot, type="Dhat4")
(BSmean <- mean(wolf_human_boot))
#bootstrapped mean overlap
print(BSmean)
print(wolf_human_boot)
#use 'basic0' output as your CI
wolf_human_ci <- bootCI(wolf_human_est, wolf_human_boot, conf = 0.95)

apex_human_boot <- bootEst(apex_boot, human_boot, type="Dhat4")
(BSmean <- mean(apex_human_boot))
#bootstrapped mean overlap
print(BSmean)
print(apex_human_boot)
#use 'basic0' output as your CI
apex_human_ci <- bootCI(apex_human_est, apex_human_boot, conf = 0.95)

#Wald test
D_hat <- wolf_human_est_site12
SE_D_hat <- sd(wolf_human_boot_12) / sqrt(length(wolf_human_boot_12))
SE_D_hat
wald_result <- wald.test(b = D_hat, Sigma = SE_D_hat^2, Terms = 1)
summary(wald_result)
print(wald_result)

######### Whole section for the overlap per site/months!!!!!!!!!!!!!
#########
##1. ORDER in categories and checking out what we have - I am checking 
## categories 4 and 5 = max human presence

site_data <- read_csv("site_category.csv")

dat_animal <- read_csv("overlap_animal_site_ready.csv")

dat_animal$time <- strptime(dat_animal$time, format = "%m/%d/%Y %H:%M:%S")

dat_animal$radiany <- (as.numeric(format(dat_animal$time, "%H")) + 
                  as.numeric(format(dat_animal$time, "%M"))/60 + 
                  as.numeric(format(dat_animal$time, "%S"))/3600) / 24 * 2 * pi

head(dat_animal)

# lynx vs wolf
lynx <- dat_animal[dat_animal$species == "Lynx lynx",]$radiany
wolf <- dat_animal[dat_animal$species == "Canis lupus",]$radiany
apex <- dat_animal[dat_animal$species %in% c("Lynx lynx", "Canis lupus"),]$radiany

# filtering only 3 - 5 categories
wolf_site_3_5 <- dat_animal[dat_animal$species == "Canis lupus" & dat_animal$human_trap_category %in% c(3, 4, 5),]$radiany
lynx_site_3_5 <- dat_animal[dat_animal$species == "Lynx lynx" & dat_animal$human_trap_category %in% c(3, 4, 5),]$radiany
n_wolf <- length(wolf_site_3_5)
n_lynx <- length(lynx_site_3_5)
print(n_wolf)
print(n_lynx)

# filtering only 3 - 5 categories
wolf_site_1_2 <- dat_animal[dat_animal$species == "Canis lupus" & dat_animal$human_trap_category %in% c(1, 2),]$radiany
lynx_site_1_2 <- dat_animal[dat_animal$species == "Lynx lynx" & dat_animal$human_trap_category %in% c(1, 2),]$radiany
n_wolf <- length(wolf_site_1_2)
n_lynx <- length(lynx_site_1_2)
print(n_wolf)
print(n_lynx)

#checking how many species in categories
species_check <- dat %>%
  filter(species == "Canis lupus", month_number == 4)
spiecies_check_am <- nrow(species_check)
print(spiecies_check_am)


# densityPlot for cat 3 to 5
densityPlot(wolf_site_3_5, rug=TRUE, main="Wolf temporal activity (Categories 3 to 5)")
densityPlot(lynx_site_3_5, rug=TRUE, main="Lynx temporal activity (Categories 3 to 5)")


## overlap for cat 3 to 5
overlapPlot(lynx_site_3_5, wolf_site_3_5, rug=TRUE, main="Lynx and wolf temporal activity (cat. 3 to 5)")
## overlap for cat 1 to 2
overlapPlot(lynx_site_1_2, wolf_site_1_2, rug=TRUE, main="Lynx and wolf temporal activity (cat. 1 to 2)")

#legend
legend('topright', legend = c(expression(italic("Lynx lynx")), 
                              expression(italic("Canis lupus"))), 
       lty = c(1, 2), col = c(1, 4), bty = 'n')

## overlap for different months 
lynx_dec <- subset(dat_animal, species == "Lynx lynx" & month(time) == 12)$radiany
wolf_dec <- subset(dat_animal, species == "Canis lupus" & month(time) == 12)$radiany
lynx_jan <- subset(dat_animal, species == "Lynx lynx" & month(time) == 1)$radiany
wolf_jan <- subset(dat_animal, species == "Canis lupus" & month(time) == 1)$radiany
lynx_feb <- subset(dat_animal, species == "Lynx lynx" & month(time) == 2)$radiany
wolf_feb <- subset(dat_animal, species == "Canis lupus" & month(time) == 2)$radiany
ynx_mar <- subset(dat_animal, species == "Lynx lynx" & month(time) == 3)$radiany
wolf_mar <- subset(dat_animal, species == "Canis lupus" & month(time) == 3)$radiany
lynx_ap <- subset(dat_animal, species == "Lynx lynx" & month(time) == 4)$radiany
wolf_ap <- subset(dat_animal, species == "Canis lupus" & month(time) == 4)$radiany

# Tworzenie wykresów overlap dla stycznia i marca
overlapPlot(lynx_dec, wolf_dec, rug=TRUE, main="Lynx and wolf temporal activity - December")
overlapPlot(lynx_jan, wolf_jan, rug=TRUE, main="Lynx and wolf temporal activity - January")
overlapPlot(lynx_feb, wolf_feb, rug=TRUE, main="Lynx and wolf temporal activity - February")
overlapPlot(lynx_mar, wolf_mar, rug=TRUE, main="Lynx and wolf temporal activity - March")

#legend
legend('topright', legend = c(expression(italic("Lynx lynx")), 
                              expression(italic("Canis lupus"))), 
       lty = c(1, 2), col = c(1, 4), bty = 'n')


########### dla człowieka / for human

human_site_dat <- read_csv("human_overlap_site_ready.csv")
head(human_site_dat)

human_site_dat$time <- strptime(human_site_dat$time, format = "%m/%d/%Y %H:%M:%S")

human_site_dat$radiany <- (as.numeric(format(human_site_dat$time, "%H")) + 
                        as.numeric(format(human_site_dat$time, "%M"))/60 + 
                        as.numeric(format(human_site_dat$time, "%S"))/3600) / 24 * 2 * pi


head(human_site_dat)

# human activity pattern
human_site <- human_site_dat[human_site_dat$human == "Homo sapiens",]$radiany

densityPlot(human_site, rug=TRUE, main="Human temporal activity")

# filtering only 4 and 5 categories
human_site_1_2 <- human_site_dat[human_site_dat$human == "Homo sapiens" & human_site_dat$human_trap_category %in% c(1, 2),]$radiany
human_site_3_5 <- human_site_dat[human_site_dat$human == "Homo sapiens" & human_site_dat$human_trap_category %in% c(3, 4, 5),]$radiany

# densityPlot for cat 4 and 5
densityPlot(human_site_3_5, rug=TRUE, main="Human temporal activity (Categories 3 and 5)")

## overlap for cat 4 and 5
overlapPlot(lynx_site_3_5, human_site_3_5, rug=TRUE, main="Lynx and human temporal activity (cat. 3 to 5)")
overlapPlot(wolf_site_3_5, human_site_3_5, rug=TRUE, main="Wolf and human temporal activity (cat. 3 to 5)")
overlapPlot(wolf_site_1_2, human_site_1_2, rug=TRUE, main="Wolf and human temporal activity (cat. 1 and 2)")
overlapPlot(lynx_site_1_2, human_site_1_2, rug=TRUE, main="Lynx and human temporal activity (cat. 1 and 2)")

#overlap for particular months
human_dec <- subset(human_site_dat, human == "Homo sapiens" & month(time) == 12)$radiany
human_jan <- subset(human_site_dat, human == "Homo sapiens" & month(time) == 1)$radiany
human_feb <- subset(human_site_dat, human == "Homo sapiens" & month(time) == 2)$radiany
human_mar <- subset(human_site_dat, human == "Homo sapiens" & month(time) == 3)$radiany

# overlap for January and March
overlapPlot(lynx_dec, human_dec, rug=TRUE, main="Lynx and human temporal activity - December")
overlapPlot(lynx_jan, human_jan, rug=TRUE, main="Lynx and human temporal activity - January")
overlapPlot(lynx_feb, human_feb, rug=TRUE, main="Lynx and human temporal activity - February")
overlapPlot(lynx_mar, human_mar, rug=TRUE, main="Lynx and human temporal activity - March")
overlapPlot(wolf_dec, human_dec, rug=TRUE, main="Wolf and human temporal activity - December")
overlapPlot(wolf_jan, human_jan, rug=TRUE, main="Wolf and human temporal activity - January")
overlapPlot(wolf_feb, human_feb, rug=TRUE, main="Wolf and human temporal activity - February")
overlapPlot(wolf_mar, human_mar, rug=TRUE, main="Wolf and human temporal activity - March")


#Legend
legend('topright', legend = c(expression(italic("Lynx lynx")), 
                              expression(italic("Homo sapiens"))), 
       lty = c(1, 2), col = c(1, 4), bty = 'n')

legend('topright', legend = c(expression(italic("Canis lupus")), 
                              expression(italic("Homo sapiens"))), 
       lty = c(1, 2), col = c(1, 4), bty = 'n')

############ Współczynnik nakładania się (coefficient of overlap)

#site
lynx_wolf_est_site35 <- overlapEst(lynx_site_3_5, wolf_site_3_5, type="Dhat1")
lynx_wolf_est_site35
lynx_human_est_site35 <- overlapEst(lynx_site_3_5, human_site_3_5, type="Dhat1")
lynx_human_est_site35
wolf_human_est_site35 <- overlapEst(wolf_site_3_5, human_site_3_5, type="Dhat4")
wolf_human_est_site35
lynx_wolf_est_site12 <- overlapEst(lynx_site_1_2, wolf_site_1_2, type="Dhat1")
lynx_wolf_est_site12
lynx_human_est_site12 <- overlapEst(lynx_site_1_2, human_site_1_2, type="Dhat1")
lynx_human_est_site12
wolf_human_est_site12 <- overlapEst(wolf_site_1_2, human_site_1_2, type="Dhat4")
wolf_human_est_site12

#december
lynx_wolf_est_dec <- overlapEst(lynx_dec, wolf_dec, type="Dhat1")
lynx_wolf_est_dec
lynx_human_est_dec <- overlapEst(lynx_dec, human_dec, type="Dhat1")
lynx_human_est_dec
wolf_human_est_dec <- overlapEst(wolf_dec, human_dec, type="Dhat4")
wolf_human_est_dec

#january
lynx_wolf_est_jan <- overlapEst(lynx_jan, wolf_jan, type="Dhat1")
lynx_wolf_est_jan
lynx_human_est_jan <- overlapEst(lynx_jan, human_jan, type="Dhat1")
lynx_human_est_jan
wolf_human_est_jan <- overlapEst(wolf_jan, human_jan, type="Dhat4")
wolf_human_est_jan

#february
lynx_wolf_est_feb <- overlapEst(lynx_feb, wolf_feb, type="Dhat1")
lynx_wolf_est_feb
lynx_human_est_feb <- overlapEst(lynx_feb, human_feb, type="Dhat1")
lynx_human_est_feb
wolf_human_est_feb <- overlapEst(wolf_feb, human_feb, type="Dhat4")
wolf_human_est_feb

#march
lynx_wolf_est_mar <- overlapEst(lynx_mar, wolf_mar, type="Dhat1")
lynx_wolf_est_mar
lynx_human_est_mar <- overlapEst(lynx_mar, human_mar, type="Dhat1")
lynx_human_est_mar
wolf_human_est_mar <- overlapEst(wolf_mar, human_mar, type="Dhat4")
wolf_human_est_mar

#checking if results for each month differ from the general overlap
mean_overlap_months <- c(wolf_human_est_dec, wolf_human_est_jan, wolf_human_est_feb, wolf_human_est_mar)
mean_overlap_all <- wolf_human_est
t_result <- t.test(mean_overlap_months, mu = 0.707476, alternative = "two.sided", conf.level = 0.99)
print(t_result)
anova_result <- aov(mean_overlap_months ~ 1)
tukey_result <- TukeyHSD(anova_result)

#BOOTSTRAP
lynx_boot_35 <- resample(lynx_site_3_5, 1000)
wolf_boot_35 <- resample(wolf_site_3_5, 1000)
human_boot_35 <- resample(human_site_3_5, 1000)

lynx_boot_12 <- resample(lynx_site_1_2, 1000)
wolf_boot_12 <- resample(wolf_site_1_2, 1000)
human_boot_12 <- resample(human_site_1_2, 1000)

lynx_wolf_boot_12 <- bootEst(lynx_boot_12, wolf_boot_12, type="Dhat1")
(BSmean <- mean(lynx_wolf_boot_12))
print(BSmean)
lynx_human_boot_12 <- bootEst(lynx_boot_12, human_boot_12, type="Dhat1")
(BSmean <- mean(lynx_human_boot_12))
print(BSmean)
wolf_human_boot_12 <- bootEst(wolf_boot_12, human_boot_12, type="Dhat4")
(BSmean <- mean(wolf_human_boot_12))
print(BSmean)

lynx_wolf_boot_35 <- bootEst(lynx_boot_35, wolf_boot_35, type="Dhat1")
(BSmean <- mean(lynx_wolf_boot_35))
print(BSmean)
lynx_human_boot_35 <- bootEst(lynx_boot_35, human_boot_35, type="Dhat1")
(BSmean <- mean(lynx_human_boot_35))
print(BSmean)
wolf_human_boot_35 <- bootEst(wolf_boot_35, human_boot_35, type="Dhat4")
(BSmean <- mean(wolf_human_boot_35))
print(BSmean)

#use 'basic0' output as your CI
lynx_wolf_ci_12 <- bootCI(lynx_wolf_est_site12, lynx_wolf_boot_12, conf = 0.95)
lynx_wolf_ci_12
lynx_human_ci_12 <- bootCI(lynx_human_est_site12, lynx_human_boot_12, conf = 0.95)
lynx_human_ci_12
wolf_human_ci_12 <- bootCI(wolf_human_est_site12, wolf_human_boot_12, conf = 0.95)
wolf_human_ci_12

lynx_wolf_ci_35 <- bootCI(lynx_wolf_est_site35, lynx_wolf_boot_35, conf = 0.95)
lynx_wolf_ci_35
lynx_human_ci_35 <- bootCI(lynx_human_est_site35, lynx_human_boot_35, conf = 0.95)
lynx_human_ci_35
wolf_human_ci_35 <- bootCI(wolf_human_est_site35, wolf_human_boot_35, conf = 0.95)
wolf_human_ci_35



