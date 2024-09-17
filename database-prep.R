setwd("your-working-directory")

library(dplyr)
library(lubridate)
library(readr)
library(purrr)

###### DATABASE PREPARATION SPECIFICALLY FOR TRAPPER FILES
###### script for site overlap database preparation
###### entry file: TRAPPER output standard csv for camera trap monitoring
###### with columns: site, time, species, count
###### site = site no, time = date and exact time of record, species = 
###### name of the recorded species, count = number of records
###### !IMPORTANT! for human database I used TRAPPER's 'count', for animal
###### database that was 'count new' - those measures are not to be compared
###### cause their meaning is different
###### also - the rows need to be multiplied if n in count is >1, cause overlap
###### only accepts 1 record in a row
## to be used in fragments :)
## 1. HUMANS 
# read the csv file
dat <- read_csv("human_overlap_site.csv")
dat <- dat[dat$count != 0, ]
write_csv(dat, "human_overlap_site_acc.csv")
df <- read_csv("C:/Users/malha/Documents/thesis-lynx-2024/human_overlap_site_acc.csv")

# multiplying the rows - I need to have 1 human record per row so if n>1, they
# need to be multiplied
powiel_wiersze <- function(df) {
  df %>%
    pmap_dfr(function(...) {
      row <- tibble(...)
      count <- as.integer(row$count)  # projecting count to integer
      if (!is.na(count) && count > 1) {
        bind_rows(replicate(count, row, simplify = FALSE))
      } else {
        row
      }
    })
}

# multiplied verses for each record (count > 1)
df_multiplied <- powiel_wiersze(df)

#saving in csv
write_csv(df_multiplied, "your-file-directory-and-name")
df <- read_csv("your-file-directory-and-name")


##### now I am assigning categories to sites - I have both csvs ready, now I am
##### merging them
human_overlap <- read_csv("human_overlap_site_acc.csv")
site_category <- read_csv("site_category.csv")

merged_data <- left_join(human_overlap, site_category, by = "site")
head(merged_data)
write_csv(merged_data, "human_overlap_site_ready.csv")

filtered_data <- merged_data %>%
  filter(!is.na(human_trap_category))

head(filtered_data)

write_csv(filtered_data, "human_overlap_site_ready.csv")

#### same for lynx and wolf file

df <- read_csv("C:/nnnn/overlap_animal_site.csv")
dat <- read_csv("overlap_animal_site.csv")
dat <- dat[dat$count != 0, ]
write_csv(dat, "overlap_animal_site_acc.csv")
df <- read_csv("C:/nnnn/overlap_animal_site_acc.csv")

powiel_wiersze <- function(df) {
  df %>%
    pmap_dfr(function(...) {
      row <- tibble(...)
      count <- as.integer(row$count)  # projecting count to integer
      if (!is.na(count) && count > 1) {
        bind_rows(replicate(count, row, simplify = FALSE))
      } else {
        row
      }
    })
}

df_multiplied <- powiel_wiersze(df)

write_csv(df_multiplied, "C:/nnnn/thesis-lynx-2024/overlap_animal_site_acc.csv")

df <- read_csv("C:/nnnn/overlap_animal_site_acc.csv")


##### now I am assigning categories to sites
animal_overlap <- read_csv("overlap_animal_site_acc.csv")
site_category <- read_csv("site_category.csv")

# checking the assigned ranges
ranges <- site_data %>%
  group_by(human_trap_category) %>%
  summarise(min_human_trap = min(human_trap),
            max_human_trap = max(human_trap))
print("Ranges of human_trap for each human_trap_category:")
print(ranges)

merged_data <- left_join(animal_overlap, site_category, by = "site")
head(merged_data)
write_csv(merged_data, "overlap_animal_site_acc_with_category.csv")

filtered_data <- merged_data %>%
  filter(!is.na(human_trap_category))

head(filtered_data)

write_csv(filtered_data, "overlap_animal_site_ready.csv")
