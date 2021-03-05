## load libraries
library(tidyverse)
library(readxl)

## load data
IFR <- read_csv("IFRImperial.csv", col_names = TRUE) %>%
    rename(ageCat = age) %>%
    mutate(ageCat = gsub("–", "-", ageCat)) %>%
    mutate(IFR = gsub("%", "", IFR)) %>%
    mutate(IFR = as.numeric(IFR) / 100) %>%
    select(ageCat, propD = IFR) %>%
    mutate(ageCat = ifelse(ageCat == "≥80", "80-89", ageCat)) %>%
    separate(ageCat, c("LB", "UB"), sep = "-", remove = FALSE) %>%
    mutate_at(vars(c("UB", "LB")), as.numeric) %>%
    mutate(ageMid = (LB + UB) / 2) %>%
    select(-LB, -UB) %>%
    mutate(data = "Imperial")

## load CDC data
hosp <- read_excel("CasesHospDeathInHospForUS_Update.xlsx") %>%
    rename(propH = `Estimated hospitalisation rate`, propD = `Overall death rate`) %>%
    select(ageCat = Age, propH, propD) %>%
    mutate(ageCat = gsub(" Years", "", ageCat)) %>%
    mutate(ageCat = ifelse(ageCat == "85+", "85-94", ageCat)) %>%
    separate(ageCat, c("LB", "UB"), sep = "-", remove = FALSE) %>%
    mutate_at(vars(c("UB", "LB")), ~gsub(" ", "", .)) %>%
    mutate_at(vars(c("UB", "LB")), as.numeric) %>%
    mutate(ageMid = (LB + UB) / 2) %>%
    select(-LB, -UB) %>%
    mutate(data = "CDC")

## extract IFR data
IFR <- select(hosp, -propH) %>%
    rbind(IFR)
hosp <- select(hosp, -propD)

## load Imperial data
hosp <- read_csv("hospitalisationImperial.csv", col_names = TRUE) %>%
    rename(ageCat = X1, propH = `Proportion of infected individuals hospitalised`) %>%
    select(ageCat, propH) %>%
    mutate(ageCat = gsub("–", "-", ageCat)) %>%
    mutate(ageCat = gsub(" years", "", ageCat)) %>%
    mutate(propH = gsub("·", ".", propH)) %>%
    separate(propH, c("propH", "CI"), sep = "% ") %>%
    mutate(propH = as.numeric(propH) / 100) %>%
    select(-CI) %>%
    mutate(ageCat = ifelse(ageCat == "≥80", "80-89", ageCat)) %>%
    separate(ageCat, c("LB", "UB"), sep = "-", remove = FALSE) %>%
    mutate_at(vars(c("UB", "LB")), as.numeric) %>%
    mutate(ageMid = (LB + UB) / 2) %>%
    select(-LB, -UB) %>%
    mutate(data = "Imperial") %>%
    rbind(hosp)

## save cleaned data
saveRDS(IFR, "IFR.rds")
saveRDS(hosp, "hosp.rds")

