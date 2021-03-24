## load libraries
library(tidyverse)
library(readxl)
library(socialmixr)

## load CoMix data
contact <- read_excel("20200327_comix_social_contacts.xlsx", sheet = "All_contacts_imputed") %>%
    select(-1)

## write out csv
write_csv(contact, "coMix_matrix.csv", col_names = FALSE)

## load POLYMOD data
contact <- contact_matrix(polymod, countries = "United Kingdom", 
    age.limits = c(0, 5, 18, 30, 40, 50, 60, 70), quiet = TRUE)

## write out csv
write_csv(as_tibble(contact$matrix), "POLYMOD_matrix.csv", col_names = FALSE)
