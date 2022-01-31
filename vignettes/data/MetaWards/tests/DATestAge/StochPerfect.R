## load libraries
library(tidyverse)
library(Rcpp)

## create output directory
dir.create("outputs")

## source Skellam code
source("trSkellam.R")

## read in parameters, remove guff and reorder
pars <- read_delim("disease.dat", delim = " ") %>%
  rename(nu = `beta[1]`, nuA = `beta[6]`) %>%
  dplyr::select(!c(starts_with("beta"), repeats, starts_with(".lock"), .p_home_weekend)) %>%
  dplyr::select(nu, nuA, !output, output)

## read in contact matrix
contact <- read_csv("POLYMOD_matrix.csv", col_names = FALSE) %>%
  as.matrix()

## extract parameters for simulation   
pars <- dplyr::select(slice(pars, 6), !output) %>%
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

## set up number of initial individuals in each age-class
N <- 10000
N <- smart_round(read_csv("age_seeds.csv", col_names = FALSE)$X2 * N)
I0 <- smart_round(read_csv("age_seeds.csv", col_names = FALSE)$X2 * 1)
N <- N - I0

## set initial counts
u <- matrix(0, 12, 8)
u[1, ] <- N
u[2, ] <- I0

## set seed
set.seed(4578)

## try discrete-time model
sourceCpp("discreteStochModel.cpp")
#disSims <- list()
source("perfectModel.R")
disSims <- perfectModel(pars, C=contact, u=u, ndays = 150, npart=50, MD=TRUE, a_dis=0.05, b_dis = 0.05)
#for(i in 1:50) {
#  disSims[[i]] <- discreteStochModel(pars, 0, 150, u, contact)
#}

stageNms <- map(c("S", "E", "A", "RA", "P", "Ione", "DI", "Itwo", "RI", "H", "RH", "DH"), ~paste0(., 1:8)) %>%
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
  dplyr::select(!c(diff, rep))

###################################
#Old underdispersion model for counts
## apply underdispersion model to counts
#obsScale <- 0.8
#scaleFn <- function(count, obsScale) {
  ## create incidence over time
#  inc <- diff(c(0, count))
#  inc <- rep(1:length(count), times = inc)
#  if(length(inc) > 1) {
#    inc <- sample(inc, round(length(inc) * obsScale))
#  }
#  inc <- table(inc)
#  count <- numeric(length(count))
#  count[as.numeric(names(inc))] <- inc
#  cumsum(count)
#}
#medRep <- mutate(medRep, across(starts_with("DI"), scaleFn, obsScale = obsScale, .names = "{.col}obs")) %>%
#  mutate(across(starts_with("DH"), scaleFn, obsScale = obsScale, .names = "{.col}obs"))
############################################################

#Skellam noise model for observations
a1=0.01
a2=0.2
b=0.1
skelNoise <- function(count, a1, a2, b1, b2){
  ##create incidence over time
  inc <- diff(c(0, count))
  #Amend U to have first upper bound of 0 (e.g. H[t-1] to offset deaths come from hospitals at previous time)
  #U <- c(0, U)
  for(i in 1:length(inc)){
    inc[i] <- inc[i] + rtskellam(a1+b1*inc[i], a2+b2*inc[i], -inc[i])#, U[i]-inc[i])
  }
  cumsum(inc)
}
medRep <- mutate(medRep, across(starts_with("DI"), skelNoise, a1=a1, a2=a2, b1=b, b2=b, .names = "{.col}obs")) %>%
  mutate(across(starts_with("DH"), skelNoise, a1=a1, a2=a2, b1=b, b2=b, .names = "{.col}obs"))
###########################################################

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
    data = pivot_longer(dplyr::select(medRep, !ends_with("obs")), !t, names_to = "var", values_to = "n") %>%
      mutate(age = gsub("[^0-9]", "", var)) %>%
      mutate(var = gsub("[^a-zA-Z]", "", var)) %>%
      mutate(var = gsub("one", "1", var)) %>%
      mutate(var = gsub("two", "2", var)),
    col = "red", linetype = "dashed"
  ) +
  geom_line(
    aes(y = n), 
    data = pivot_longer(dplyr::select(medRep, t, ends_with("obs")), !t, names_to = "var", values_to = "n") %>%
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

