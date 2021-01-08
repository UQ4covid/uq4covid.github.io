## load libraries
library(tidyverse)
library(readxl)

## load data
contact <- read_excel("20200327_comix_social_contacts.xlsx", sheet = "All_contacts_imputed") %>%
    select(-1)

## write out csv
write_csv(contact, "contact_matrix.csv", col_names = FALSE)
