## load libraries
library(tidyverse)
library(Rcpp)
library(sf)
library(gganimate)
library(viridis)

## create output directory
dir.create("outputs")

## read in parameters, remove guff and reorder
pars <- read_delim("disease.dat", delim = " ") %>%
    rename(nu = `beta[1]`, nuA = `beta[6]`) %>%
    select(!c(starts_with("beta"), repeats, starts_with(".lock"), .p_home_weekend)) %>%
    select(nu, nuA, !output, output)

## read in contact matrix
contact <- read_csv("POLYMOD_matrix.csv", col_names = FALSE) %>%
    as.matrix()

## extract parameters for simulation   
pars <- select(slice(pars, 50), !output) %>%
    unlist()

## read in commuter data
EW19 <- read_delim("EW19.dat", delim = " ", col_names = FALSE)

## add seeding information
EW19 <- mutate(EW19, X4 = 0)
for(i in 1:10) {
    val <- 0
    while(val == 0) {
        seed <- sample(which(EW19[, 1] == 317), 1)
        if((EW19$X4[seed] + 1) <= EW19$X3[seed]) {
            EW19$X4[seed] <- EW19$X4[seed] + 1
            val <- 1
        }
    }
}
EW19 <- as.matrix(EW19)

## add age probabilities
age_probs <- read_csv("age_seeds.csv", col_names = FALSE)$X2

## set seed
set.seed(4578)

## try discrete-time model
sourceCpp("discreteStochModel.cpp")
disSims <- list()
# for(i in 1:50) {
    disSims[[1]] <- discreteStochModel(pars, 0, 100, age_probs, EW19, contact)
# }
stageNms <- map(c("S", "E", "A", "RA", "P", "Ione", "DI", "Itwo", "RI", "H", "RH", "DH"), ~paste0(., "_", 1:8)) %>%
    map(~map(., ~paste0(., "_", 1:max(EW19[, 1])))) %>%
    reduce(c) %>%
    reduce(c)
disSims <- map(disSims, ~as_tibble(.)) %>%
    bind_rows(.id = "rep") %>%
    set_names(c("rep", "t", stageNms)) %>%
    mutate(t = t + 1)

## spatial animation of simulation

## read in shapefile
lad19 <- st_read("LAD19_shapefile/LAD19_shapefile.shp")

## extract cases over time in each age-group
p <- select(disSims, t, starts_with("E_")) %>%
    filter(t <= 30) %>%
    pivot_longer(!t, values_to = "counts", names_to = "var") %>%
    separate(var, c("var", "age", "lad"), sep = "_") %>%
    select(!var) %>%
    mutate(age = as.numeric(age), lad = as.numeric(lad)) %>%
    group_by(lad) %>%
    mutate(tot = sum(counts)) %>%
    ungroup() %>%
    filter(tot > 0) %>%
    select(!tot) %>%
    group_by(lad, age) %>%
    mutate(tot = cumsum(counts)) %>%
    group_by(lad, t) %>%
    mutate(tot = sum(tot)) %>%
    ungroup() %>%
    mutate(counts = ifelse(tot == 0, NA, counts)) %>%
    select(!tot)
max_p <- max(p$counts, na.rm = TRUE)
p <- inner_join(lad19, p, by = c("objectid" = "lad")) %>%
    ggplot() +
        geom_sf(aes(fill = counts), colour = NA) +
        facet_wrap(~age) +
        scale_fill_viridis_c(limits = c(0, max_p)) +
        theme_bw()

## add transitions
p <- p + transition_time(t) + ggtitle("Day = {frame_time}")
spatial_gif <- animate(p, nframes = 30, fps = 1)
spatial_gif
anim_save("spatial.gif", spatial_gif)

## extract simulation closest to median
medRep <- pivot_longer(disSims, !c(rep, t), names_to = "var", values_to = "n") %>%
    group_by(t, var) %>%
    summarise(
        median = median(n),
        .groups = "drop"
    ) %>%
    inner_join(
        pivot_longer(disSims, !c(rep, t), names_to = "var", values_to = "n"),
        by = c("t", "var")
    ) %>%
    mutate(diff = (median - n)^2) %>%
    group_by(rep) %>%
    summarise(diff = sum(diff), .groups = "drop") %>%
    arrange(diff) %>%
    slice(2) %>%
    inner_join(disSims, by = "rep") %>%
    select(!c(diff, rep))

## apply underdispersion model to counts
obsScale <- 0.8
scaleFn <- function(count, obsScale) {
    ## create incidence over time
    inc <- diff(c(0, count))
    inc <- rep(1:length(count), times = inc)
    if(length(inc) > 1) {
        inc <- sample(inc, round(length(inc) * obsScale))
    }
    inc <- table(inc)
    count <- numeric(length(count))
    count[as.numeric(names(inc))] <- inc
    cumsum(count)
}
medRep <- mutate(medRep, across(starts_with("DI"), scaleFn, obsScale = obsScale, .names = "{.col}obs")) %>%
    mutate(across(starts_with("DH"), scaleFn, obsScale = obsScale, .names = "{.col}obs"))

## plot replicates
p <- pivot_longer(disSims, !c(rep, t), names_to = "var", values_to = "n") %>%
    group_by(t, var) %>%
    summarise(
        LCI = quantile(n, probs = 0.025),
        LQ = quantile(n, probs = 0.25),
        median = median(n),
        UQ = quantile(n, probs = 0.75),
        UCI = quantile(n, probs = 0.975),
        .groups = "drop"
    ) %>%
    mutate(age = gsub("[^0-9]", "", var)) %>%
    mutate(var = gsub("[^a-zA-Z]", "", var)) %>%
    mutate(var = gsub("one", "1", var)) %>%
    mutate(var = gsub("two", "2", var)) %>%
    ggplot(aes(x = t)) +
        geom_ribbon(aes(ymin = LCI, ymax = UCI), alpha = 0.5) +
        geom_ribbon(aes(ymin = LQ, ymax = UQ), alpha = 0.5) +
        geom_line(aes(y = median)) +
        geom_line(
            aes(y = n), 
            data = pivot_longer(select(medRep, !ends_with("obs")), !t, names_to = "var", values_to = "n") %>%
                mutate(age = gsub("[^0-9]", "", var)) %>%
                mutate(var = gsub("[^a-zA-Z]", "", var)) %>%
                mutate(var = gsub("one", "1", var)) %>%
                mutate(var = gsub("two", "2", var)),
            col = "red", linetype = "dashed"
        ) +
        geom_line(
            aes(y = n), 
            data = pivot_longer(select(medRep, t, ends_with("obs")), !t, names_to = "var", values_to = "n") %>%
                mutate(var = gsub("obs", "", var)) %>%
                mutate(age = gsub("[^0-9]", "", var)) %>%
                mutate(var = gsub("[^a-zA-Z]", "", var)) %>%
                mutate(var = gsub("one", "1", var)) %>%
                mutate(var = gsub("two", "2", var)),
            col = "blue", linetype = "dashed"
        ) +
        facet_grid(var ~ age, scales = "free") +
        xlab("Days") + 
        ylab("Counts")
ggsave("outputs/sims.pdf", p, width = 10, height = 10)

saveRDS(medRep, "outputs/disSims.rds")
saveRDS(pars, "outputs/pars.rds")
