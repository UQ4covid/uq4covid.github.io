## load libraries
library(tidyverse)
library(Rcpp)
library(RcppArmadillo)
library(coda)
library(GGally)

## source full iFFBS MCMC sampler
sourceCpp("iFFBSFullMCMC.cpp")

## read in data and generate incidence curves
Npop <- 150
data <- readRDS("outputs/disSims.rds") %>%
    mutate(across(c(RA, RI, DI, RH, DH), ~. - lag(., default = 0))) %>%
    mutate(H = H - lag(H, default = 0) + DH + RH) %>%
    mutate(I2 = I2 - lag(I2, default = 0) + RI) %>%
    mutate(I1 = I1 - lag(I1, default = 0) + I2 + DI + H) %>%
    mutate(P = P - lag(P, default = 0) + I1) %>%
    select(t, RA, DI, RI, RH, DH, P)
DI_prime <- data$DI
DH_prime <- data$DH
RA_prime <- data$RA
RI_prime <- data$RI
RH_prime <- data$RH
P_prime <- data$P
pini <- 1 / (Npop * length(DH_prime))

## run iFFBS algorithm
set.seed(666)
post <- iFFBS(DI_prime, DH_prime, RA_prime, RI_prime, RH_prime, P_prime, 1:length(DI_prime), Npop, 1000, pini, pars = rep(NA, 11))

## extract current pars
currpars <- as.matrix(post)[, 1:11]
currpars <- currpars[nrow(currpars), ]

## restart algorithm at current parameter values
## to help adaptation
post <- iFFBS(DI_prime, DH_prime, RA_prime, RI_prime, RH_prime, P_prime, 1:length(DI_prime), Npop, 50000, pini, pars = currpars)

## extract current pars
currpars <- as.matrix(post)[, 1:11]
currpars <- currpars[nrow(currpars), ]

## restart algorithm at current parameter values
## to help adaptation
post <- iFFBS(DI_prime, DH_prime, RA_prime, RI_prime, RH_prime, P_prime, 1:length(DI_prime), Npop, 300000, pini, pars = currpars)

## save output
saveRDS(post, "outputs/iFFBSpost.rds")

## convert to correct format and add posterior for pA
colnames(post) <- c("pH", "pHD", "pI1", "pI1D", "pI1H", "pI2", "pP", "pE", "pEP", "beta", "betaA", apply(expand.grid(c("S", "E", "A", "RA", "P", "I1", "I2", "RI", "DI", "H", "RH", "DH"), 1:length(DI_prime)), 1, function(x) paste0(x[1], "_", as.numeric(x[2]))), "loglike")
post <- cbind(post[, 1:11], 1 - exp(-1.0 / (1.0 / (-log(1.0 - post[, "pP"])) + 1.0 / (-log(1.0 - post[, "pI1"])) + 1.0 / (-log(1.0 - post[, "pI2"])))), post[, ncol(post)], post[, -c(1:11, ncol(post))])
colnames(post)[12] <- "pA"
colnames(post) <- c("pH", "pHD", "pI1", "pI1D", "pI1H", "pI2", "pP", "pE", "pEP", "beta", "betaA", "pA", "loglike", apply(expand.grid(c("S", "E", "A", "RA", "P", "I1", "I2", "RI", "DI", "H", "RH", "DH"), 1:length(DI_prime)), 1, function(x) paste0(x[1], "_", as.numeric(x[2]))))
post <- as.mcmc(post)

## plot trace plots (and convert to raster to save memory)
## this uses Linux commands including ImageMagick and
## pdftk to do the conversions and bind back
## to a PDF (so that multiple pages can be combined)
png("iFFBSTraces%d.png")
plot(post[, 1:13])
dev.off()
system("for i in iFFBSTraces*.png; do convert \"$i\" \"$i\".pdf; done")
system("rm iFFBSTraces*.png")
system("pdftk *.png.pdf output outputs/iFFBSTraces.pdf")
system("rm iFFBSTraces*.png.pdf")

## save output
saveRDS(post, "outputs/iFFBSpost.rds")

## remove burnin
post <- window(post, start = 50000)

## plot posteriors
p <- as.matrix(post) %>%
    as_tibble() %>%
    select(1:12) %>%
    pivot_longer(everything(), names_to = "var", values_to = "post") %>%
    ggplot() +
        geom_density(aes(x = post)) +
        geom_point(aes(x = pars, y = 0), data = as.data.frame(t(readRDS("outputs/pars.rds"))) %>%
            pivot_longer(everything(), names_to = "var", values_to = "pars") %>%
            mutate(var = ifelse(var == "nuA", "betaA", var)) %>%
            mutate(var = ifelse(var == "nu", "beta", var))) +
        facet_wrap(~ var, scales = "free")
ggsave("outputs/iFFBSPost.pdf", p)

## joint posteriors  (convert to raster to save memory)     
p <- as.matrix(post) %>%
    as_tibble() %>%
    select(1:12) %>%
    ggpairs()
ggsave("iFFBSJointPost.png", p, width = 10, height = 10)
system("convert iFFBSJointPost.png outputs/iFFBSJointPost.pdf")
system("rm iFFBSJointPost.png")

## plot marginal densities over time
p1 <- as.matrix(window(post, thin = 50)) %>%
    as_tibble() %>%
    select(-c(1:13)) %>%
    summarise(across(everything(), list(
        median = ~median(.), 
        LCI = ~quantile(., probs = 0.025),
        LQ = ~quantile(., probs = 0.25),
        UQ = ~quantile(., probs = 0.75),
        UCI = ~quantile(., probs = 0.975)
    ))) %>%
    pivot_longer(everything(), names_to = "var", values_to = "value") %>%
    separate(var, c("var", "t", "estimate"), sep = "_") %>%
    pivot_wider(names_from = estimate, values_from = value) %>%
    mutate(t = as.numeric(t))
p2 <- readRDS("outputs/disSims.rds") %>%
    pivot_longer(!t, names_to = "var", values_to = "count")
p <- ggplot(p1) +
        geom_ribbon(aes(x = t, ymin = LCI, ymax = UCI), colour = NA, alpha = 0.5) +
        geom_ribbon(aes(x = t, ymin = LQ, ymax = UQ), colour = NA, alpha = 0.5) +
        geom_line(aes(x = t, y = count), data = p2, col = "red") +
        facet_wrap(~var, scales = "free") +
        xlab("Days") + 
        ylab("Counts")
ggsave("outputs/iFFBSStates.pdf", p)

