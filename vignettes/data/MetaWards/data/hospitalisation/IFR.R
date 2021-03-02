## load libraries
library(tidyverse)
library(svglite)

## load data
IFR <- read_csv("IFRImperial.csv", col_names = TRUE) %>%
    rename(ageCat = age) %>%
    mutate(ageCat = gsub("–", "-", ageCat)) %>%
    mutate(IFR = gsub("%", "", IFR)) %>%
    mutate(IFR = as.numeric(IFR) / 100) %>%
    mutate(CI = gsub("\\(", "", CI)) %>%
    mutate(CI = gsub("\\)", "", CI)) %>%
    separate(CI, c("LCI", "UCI"), sep = "–") %>%
    mutate(LCI = as.numeric(LCI) / 100) %>%
    mutate(UCI = as.numeric(UCI) / 100)

## match data roughly to age-classes used in the model
IFR <- mutate(IFR, ageCat = ifelse(ageCat == "70-79", "70+", ageCat)) %>%
    mutate(ageCat = ifelse(ageCat == "≥80", "70+", ageCat)) %>%
    group_by(ageCat) %>%
    summarise_all(mean) %>%
    arrange(ageCat)

## extract mid-points of age-classes
IFR <- IFR %>%
    mutate(ageCat = ifelse(ageCat == "70+", "70-80", ageCat)) %>%
    separate(ageCat, c("LB", "UB"), sep = "-", remove = FALSE) %>%
    mutate_at(vars(c("UB", "LB")), as.numeric) %>%
    mutate(ageMid = (LB + UB) / 2) %>%
    dplyr::select(-LB, -UB)

## plot bounds
p <- ggplot(IFR, aes(x = ageMid)) +
    geom_point(aes(y = IFR)) +
    geom_ribbon(aes(ymin = LCI, ymax = UCI), alpha = 0.5) +
    xlab("Age") + ylab("IFR")
ggsave("IFR.pdf", p)
ggsave("IFR.svg", p)
ggsave("IFR.png", p)

## save cleaned data
saveRDS(IFR, "IFR.rds")

