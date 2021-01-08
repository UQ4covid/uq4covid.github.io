## load libraries
library(tidyverse)
library(fitdistrplus)

## load data
hosp <- read_csv("admissionToOutcomeLineList.csv", col_names = TRUE)

## summarise data
mutate_if(hosp, is.character, as.factor) %>%
    summary()

## match data roughly to age-classes used in the model
hosp <- mutate(hosp, ageCat = ifelse(ageCat == "70-79", "70+", ageCat)) %>%
    mutate(ageCat = ifelse(ageCat == "80-89", "70+", ageCat)) %>%
    mutate(ageCat = ifelse(ageCat == "90+", "70+", ageCat)) %>%
    arrange(ageCat)

## now plot time spent in hospital
hosp %>%
    filter(finaloutcome == "Death" | finaloutcome == "Discharged") %>%
    ggplot(aes(x = time, fill = finaloutcome)) +
        geom_density(alpha = 0.5, colour = NA) +
        facet_wrap(~ageCat, scales = "free_y")
hosp %>%
    filter(finaloutcome == "Death" | finaloutcome == "Discharged") %>%
    ggplot(aes(x = time)) +
    geom_density() +
    facet_wrap(~ageCat, scales = "free_y")

## fit some simple models to each age-class and ignore outcome
hosp <- hosp %>%
    filter(finaloutcome == "Death" | finaloutcome == "Discharged") %>%
    dplyr::select(time, ageCat) %>%
    group_by(ageCat) %>%
    nest() %>%
    mutate(exp = map(data, ~fitdist(.$time, "exp"))) %>%
    mutate(gamma = map(data, ~fitdist(.$time, "gamma"))) %>%
    mutate(expAIC = map_dbl(exp, "aic")) %>%
    mutate(gammaAIC = map_dbl(gamma, "aic")) %>%
    mutate(model = ifelse(expAIC < gammaAIC, exp, gamma)) %>%
    mutate(coefs = map(model, "estimate")) %>%
    mutate(preds = map(coefs, ~{
        t <- 0:100
        if(length(.) == 1) {
            out <- dexp(t, rate = .)
        } else {
            out <- dgamma(t, shape = .[1], rate = .[2])
        }
        tibble(time = t, out = out)
    }))

## fitted plots
dplyr::select(hosp, ageCat, data) %>%
    unnest(cols = data) %>%
    ggplot(aes(x = time)) +
        geom_density() +
        facet_wrap(~ageCat, scales = "free_y") +
        geom_line(aes(y = out), data = dplyr::select(hosp, ageCat, preds) %>%
              unnest(cols = preds), col = "red")
ggsave("bestFits.pdf")

## take a look at shape / scale
hosp$coefs

## all apart from last one are close to one scale
## so try to refit all using exponential
hosp <- hosp %>%
    mutate(coefs = map(exp, "estimate")) %>%
    mutate(preds = map(coefs, ~{
        t <- 0:100
        out <- dexp(t, rate = .)
        tibble(time = t, out = out)
    }))

## fitted plots
dplyr::select(hosp, ageCat, data) %>%
    unnest(cols = data) %>%
    ggplot(aes(x = time)) +
    geom_density() +
    facet_wrap(~ageCat, scales = "free_y") +
    geom_line(aes(y = out), data = dplyr::select(hosp, ageCat, preds) %>%
                  unnest(cols = preds), col = "red")
ggsave("expFits.pdf")

## take a look at rates
reduce(hosp$coefs, c)
