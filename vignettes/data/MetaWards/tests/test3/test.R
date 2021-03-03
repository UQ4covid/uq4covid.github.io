library(tidyverse)
library(patchwork)

## source reconstruct function
source("../../R_tools/dataTools.R")

## read in pars
pars <- read_delim("disease.dat", delim = " ") %>%
    select(ends_with("_1"), nu = `beta[1]`, nuA = `beta[6]`)
colnames(pars) <- gsub("_1", "", colnames(pars))
colnames(pars) <- gsub("\\.", "", colnames(pars))
pars <- pars %>%
    mutate_at(vars(pE, pP, pA, pI1, pI2, pH), ~ifelse(. == 1, 0.99, .)) %>%
    mutate(gammaE = -log(1 - pE)) %>%
    mutate(gammaP = -log(1 - pP)) %>%
    mutate(gammaA = -log(1 - pA)) %>%
    mutate(gammaI1 = -log(1 - pI1)) %>%
    mutate(gammaI2 = -log(1 - pI2)) %>%
    mutate(gammaH = -log(1 - pH)) %>%
    select(-pE, -pP, -pA, -pI1, -pI2, -pH) %>%
    unlist()

## set population size
N <- c(10011, 25694, 9985, 8654, 8322, 8654, 6823, 21857)

## set contact matrix
contact <- read_csv("contact_matrix.csv", col_names = FALSE) %>%
    as.matrix()

## check NGM
S0 <- c(N[1] - 10, N[-1])
R0 <- NGM(R0 = NA, nu = pars["nu"], C = contact, S0 = S0, N = N, 
    nuA = pars["nuA"], gammaE = pars["gammaE"], pEP = pars["pEP"], gammaA = pars["gammaA"], 
    gammaP = pars["gammaP"], gammaI1 = pars["gammaI1"], pI1H = pars["pI1H"], pI1D = pars["pI1D"], 
    gammaI2 = pars["gammaI2"])
K <- R0$K
R0 <- R0$R0
print(paste0("R0 = ", R0))

## reduce from large domain NGM to small domain NGM (Diekmann et al, 2010)
Kind <- !apply(K, 1, function(x) all(x == 0))
Ksmall <- matrix(0, nrow(K), sum(Kind))
Ksmall[1:sum(Kind), 1:sum(Kind)] <- diag(sum(Kind))
Ksmall <- t(Ksmall) %*% K %*% Ksmall

## solve balance equation for final size (from Andreasen 2011)
balance <- function(z, C) {
    sum(abs(z - (1 - exp(-C %*% z))))
}
A <- t(N * t(Ksmall / N))
finalsize <- optim(rep(0.5, length(N)), balance, C = A, control = list(maxit = 10000), method = "L-BFGS-B", lower = rep(0, length(N)), upper = rep(1, length(N)))
stopifnot(finalsize$value <= 0.0007 | finalsize$convergence == 0)
finalsize <- finalsize$par * N

## fit deterministic model
source("../detModel.R")

## fit discrete-time stochastic model
source("stochModel.R")

## loop over age classes
pdf("test.pdf", width = 10, height = 5)
rec <- list()
for(i in 1:8) {
    ## set up list for runs
    rec[[i]] <- list()
    ## loop over repeats
    for(j in 1:10) {
        ## establish connection
        system(paste0("bzip2 -dkf raw_outputs/testx", str_pad(j, 3, pad = "0"), "/age", i, ".db.bz2"))
        con <- DBI::dbConnect(RSQLite::SQLite(), paste0("raw_outputs/testx", str_pad(j, 3, pad = "0"), "/age", i, ".db"))

        ## extract data
        compact <- tbl(con, "compact") %>%
            collect()
        
        DBI::dbDisconnect(con)
        
        ## reconstruct counts from incidence
        rec[[i]][[j]] <- reconstruct(
                    compact$Einc, compact$Pinc, compact$I1inc, compact$I2inc, compact$RIinc, 
                    compact$DIinc, compact$Ainc, compact$RAinc,
                    compact$Hinc, compact$RHinc, compact$DHinc
                ) %>%
                magrittr::set_colnames(c(
                    "Einc", "E", "Pinc", "P", "I1inc", "I1", "I2inc", "I2", 
                    "RI", "DI",
                    "Ainc", "A", "RA",
                    "Hinc", "H", "RH", "DH"
                )) %>%
                as_tibble()
        rec[[i]][[j]]$day <- compact$day
    }
    
    ## collapse to mean and CIs
    rec[[i]] <- bind_rows(rec[[i]], .id = "rep") %>%
        select(-ends_with("inc")) %>%
        complete(rep, day = 1:max(.$day)) %>%
        group_by(rep) %>%
        mutate_at(-c(1, 2), ~ifelse(day == 1 & is.na(.), 0, .)) %>%
        fill(names(.)) %>%
        ungroup() %>%
        select(-rep) %>%
        gather(stage, count, -day) %>%
        group_by(day, stage) %>%
        nest() %>%
        mutate(data = map(data, ~{
            data.frame(
                mn = mean(.$count), 
                LCI = quantile(.$count, probs = 0.025), 
                UCI = quantile(.$count, probs = 0.975)
            )
        })) %>%
        unnest(cols = data) %>%
        gather(estimate, value, -day, -stage) %>%
        spread(stage, value) %>%
        select(day, estimate, E, P, I1, DI)
        
    ## extract discrete-time stochastic simulations
    tempStoch <- select(disSims, day = t, ends_with(as.character(i))) %>%
        set_names(gsub(paste0(i, "$"), "", colnames(.))) %>%
        gather(stage, count, -day) %>%
        group_by(day, stage) %>%
        nest() %>%
        mutate(data = map(data, ~{
            data.frame(
                mn = mean(.$count), 
                LCI = quantile(.$count, probs = 0.025), 
                UCI = quantile(.$count, probs = 0.975)
            )
        })) %>%
        unnest(cols = data) %>%
        gather(estimate, value, -day, -stage) %>%
        spread(stage, value) %>%
        select(day, estimate, E, P, I1, DI) %>%
        ungroup()
        
    ## extract deterministic run
    tempDet <- out %>%
        select(-starts_with("S")) %>%
        select(time, ends_with(as.character(i))) %>%
        mutate_all(as.numeric) %>%
        rename(day = time) %>%
        mutate(estimate = "mn") %>%
        set_names(gsub(paste0(as.character(i), "$"), "", colnames(.))) %>%
        select(day, estimate, E, P, I1, DI)
    
    ## bind models
    temp <- list(stoch =  tempStoch, det = tempDet, MW = rec[[i]]) %>%
        bind_rows(.id = "model") %>%
        gather(stage, count, -model, -estimate, -day) %>%
        spread(estimate, count)
        
    p1 <- ggplot(temp, aes(x = day)) +
        geom_ribbon(aes(ymin = LCI, ymax = UCI, fill = stage), data = filter(temp, model == "MW"), alpha = 0.3) + 
        geom_line(aes(y = mn, colour = stage, linetype = model)) +
        geom_hline(yintercept = finalsize[i], linetype = "dashed") +
        ggtitle(paste0("age", i))
    p2 <- ggplot(filter(temp, model != "det"), aes(x = day)) +
        geom_ribbon(aes(ymin = LCI, ymax = UCI, fill = stage), data = filter(temp, model == "MW"), alpha = 0.3) +
        geom_ribbon(aes(ymin = LCI, ymax = UCI, fill = stage), data = filter(temp, model == "stoch"), alpha = 0.3) +
        geom_line(aes(y = mn, colour = stage, linetype = model)) +
        ggtitle(paste0("age", i))
    print(p1 + p2)
}
dev.off()
