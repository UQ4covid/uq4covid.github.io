## load libraries
library(tidyverse)
library(Rcpp)
library(sf)
library(gganimate)
library(viridis)
library(parallel)
library(abind)
library(patchwork)

## set seed
set.seed(4578)

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
pars <- select(slice(pars, 6), !output) %>%
    unlist()

## solution to round numbers preserving sum
## adapted from:
## https://stackoverflow.com/questions/32544646/round-vector-of-numerics-to-integer-while-preserving-their-sum
smart_round <- function(x) {
    y <- floor(x)
    indices <- tail(order(x - y), round(sum(x)) - sum(y))
    y[indices] <- y[indices] + 1
    y
}

## add age probabilities
ageProbs <- read_csv("age_seeds.csv", col_names = FALSE)$X2

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
## expand to deal with age-classes
u1 <- apply(EW19, 1, function(x, ageProbs) {
        u <- matrix(0, 12, length(ageProbs))
        u[2, ] <- smart_round(ageProbs * x[4])
        u[1, ] <- smart_round(ageProbs * x[3]) - u[2, ]
        list(u)
    }, ageProbs = ageProbs) %>%
    map(1) %>%
    abind(along = 3)
u1_moves <- as.matrix(EW19[, 1:2])
u1_day <- list()
u1_night <- list()
for(i in 1:max(u1_moves[, 1])) {
    temp <- u1[, , u1_moves[, 2] == i]
    u1_day[[i]] <- apply(temp, c(1, 2), sum)
    temp <- u1[, , u1_moves[, 1] == i]
    u1_night[[i]] <- apply(temp, c(1, 2), sum)
}
u1_day <- abind(u1_day, along = 3)
u1_night <- abind(u1_night, along = 3)
N_day <- apply(u1_day, c(2, 3), sum)
N_night <- apply(u1_night, c(2, 3), sum)

## write inputs out
saveRDS(u1, "outputs/u1.rds")
saveRDS(u1_moves, "outputs/u1_moves.rds")
saveRDS(u1_day, "outputs/u1_day.rds")
saveRDS(u1_night, "outputs/u1_night.rds")
saveRDS(N_day, "outputs/N_day.rds")
saveRDS(N_night, "outputs/N_night.rds")

## try discrete-time model
sourceCpp("discreteStochModel.cpp")
disSims <- mclapply(1:8, function(i) {
    discreteStochModel(pars, 0, 100, u1_moves, u1, u1_day, u1_night, N_day, N_night, contact)
}, mc.cores = 8)
stageNms <- map(c("S", "E", "A", "RA", "P", "Ione", "DI", "Itwo", "RI", "H", "RH", "DH"), ~paste0(., "_", 1:8)) %>%
    map(~map(., ~paste0(., "_", 1:max(EW19[, 1])))) %>%
    reduce(c) %>%
    reduce(c)
disSims <- map(disSims, ~as_tibble(.)) %>%
    bind_rows(.id = "rep") %>%
    set_names(c("rep", "t", stageNms)) %>%
    mutate(t = t + 1)

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

## plot replicates at national level
p <- pivot_longer(disSims, !c(rep, t), names_to = "var", values_to = "n") %>%
    mutate(var = gsub("_[0-9]*$", "", var)) %>%
    group_by(rep, t, var) %>%
    summarise(n = sum(n), .groups = "drop") %>%
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
                mutate(var = gsub("_[0-9]*$", "", var)) %>%
                group_by(t, var) %>%
                summarise(n = sum(n), .groups = "drop") %>%
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
                mutate(var = gsub("_[0-9]*$", "", var)) %>%
                group_by(t, var) %>%
                summarise(n = sum(n), .groups = "drop") %>%
                mutate(age = gsub("[^0-9]", "", var)) %>%
                mutate(var = gsub("[^a-zA-Z]", "", var)) %>%
                mutate(var = gsub("one", "1", var)) %>%
                mutate(var = gsub("two", "2", var)),
            col = "blue", linetype = "dashed"
        ) +
        facet_grid(var ~ age, scales = "free") +
        xlab("Days") + 
        ylab("Counts")
ggsave("outputs/simsNational.pdf", p, width = 10, height = 10)

## plot replicates at LAD level for LADs with largest epidemic load
p <- select(medRep, !ends_with("obs")) %>%
    pivot_longer(!t, names_to = "var", values_to = "n") %>%
    mutate(age = gsub('^(?:[^_]*_)(.*)', '\\1', var)) %>%
    mutate(LAD = gsub('^(?:[^_]*_)(.*)', '\\1', age)) %>%
    mutate(age = gsub('(.*)_[0-9]*', '\\1', age)) %>%
    mutate(var = gsub('^(.*)_[0-9]*_.*', '\\1', var)) %>%
    filter(var == "DI") %>%
    filter(t == max(t)) %>%
    group_by(LAD) %>%
    summarise(n = sum(n), .groups = "drop") %>%
    arrange(desc(n)) %>%
    slice(1:5) %>%
    select(-n) %>%
    inner_join(
        pivot_longer(disSims, !c(rep, t), names_to = "var", values_to = "n") %>%
        mutate(age = gsub('^(?:[^_]*_)(.*)', '\\1', var)) %>%
        mutate(LAD = gsub('^(?:[^_]*_)(.*)', '\\1', age)),
        by = "LAD"
    ) %>%
    mutate(age = gsub('(.*)_[0-9]*', '\\1', age)) %>%
    mutate(var = gsub('^(.*)_[0-9]*_.*', '\\1', var)) %>%
    group_by(t, var, age, LAD) %>%
    summarise(
        LCI = quantile(n, probs = 0.025),
        LQ = quantile(n, probs = 0.25),
        median = median(n),
        UQ = quantile(n, probs = 0.75),
        UCI = quantile(n, probs = 0.975),
        .groups = "drop"
    ) %>%
    mutate(var = gsub("one", "1", var)) %>%
    mutate(var = gsub("two", "2", var))
p1 <- list()
for(i in 1:5) {
    p1[[i]] <- filter(p, LAD == unique(p$LAD)[i]) %>%
        select(-LAD) %>%
        ggplot(aes(x = t)) +
        geom_ribbon(aes(ymin = LCI, ymax = UCI), alpha = 0.5) +
        geom_ribbon(aes(ymin = LQ, ymax = UQ), alpha = 0.5) +
        geom_line(aes(y = median)) +
        geom_line(
            aes(y = n), 
            data = pivot_longer(select(medRep, !ends_with("obs")), !t, names_to = "var", values_to = "n") %>%
                mutate(age = gsub('^(?:[^_]*_)(.*)', '\\1', var)) %>%
                mutate(LAD = gsub('^(?:[^_]*_)(.*)', '\\1', age)) %>%
                mutate(age = gsub('(.*)_[0-9]*', '\\1', age)) %>%
                mutate(var = gsub('^(.*)_[0-9]*_.*', '\\1', var)) %>%
                filter(LAD == unique(p$LAD)[i]) %>%
                mutate(var = gsub("one", "1", var)) %>%
                mutate(var = gsub("two", "2", var)), 
            col = "red", linetype = "dashed"
        ) +
        geom_line(
            aes(y = n), 
            data = pivot_longer(select(medRep, t, ends_with("obs")), !t, names_to = "var", values_to = "n") %>%
                mutate(var = gsub("obs", "", var)) %>%
                mutate(age = gsub('^(?:[^_]*_)(.*)', '\\1', var)) %>%
                mutate(LAD = gsub('^(?:[^_]*_)(.*)', '\\1', age)) %>%
                mutate(age = gsub('(.*)_[0-9]*', '\\1', age)) %>%
                mutate(var = gsub('^(.*)_[0-9]*_.*', '\\1', var)) %>%
                filter(LAD == unique(p$LAD)[i]) %>%
                mutate(var = gsub("one", "1", var)) %>%
                mutate(var = gsub("two", "2", var)),
            col = "blue", linetype = "dotted"
        ) +
        facet_grid(var ~ age, scales = "free") +
        xlab("Days") + 
        ylab("Counts") +
        ggtitle(paste0("LAD = ", unique(p$LAD)[i]))
}
p1 <- wrap_plots(p1, nrow = 3)
ggsave("outputs/simsTopLADs.pdf", p1, width = 20, height = 30)

## spatial animation of simulation

## read in shapefile
lad19 <- st_read("LAD19_shapefile/LAD19_shapefile.shp")

## extract cases over time in each age-group
p <- select(medRep, t, starts_with("DI_")) %>%
    select(!ends_with("obs")) %>%
    filter(t <= 50) %>%
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
spatial_gif <- animate(p, nframes = 50, fps = 1)
anim_save("outputs/simsSpatialDI.gif", spatial_gif)

## extract cases over time in each age-group
p <- select(medRep, t, starts_with("E")) %>%
    filter(t <= 50) %>%
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
spatial_gif <- animate(p, nframes = 50, fps = 1)
anim_save("outputs/simsSpatialE.gif", spatial_gif)

saveRDS(medRep, "outputs/disSims.rds")
saveRDS(pars, "outputs/pars.rds")
