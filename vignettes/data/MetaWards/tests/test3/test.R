## R script to query database
library(tidyverse)

## source reconstruct function
source("../../R_tools/dataTools.R")

## fit deterministic model
source("detModel.R")

## loop over age classes
p <- list()
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
        select(day, estimate, E, P, I1, DI = RI)
        
    ## plot against deterministic model
    tempDet <- out %>%
        select(-starts_with("S")) %>%
        select(time, ends_with(as.character(i))) %>%
        mutate_all(as.numeric) %>%
        rename(day = time) %>%
        mutate(estimate = "mn")
    colnames(tempDet) <- gsub(paste0(as.character(i), "$"), "", colnames(tempDet))
    tempDet <- select(tempDet, day, estimate, E, P, I1, DI)
    
    ## bind deterministic and stochastic models
    temp <- list(det =  tempDet, MW = rec[[i]]) %>%
        bind_rows(.id = "model") %>%
        gather(stage, count, -model, -estimate, -day) %>%
        spread(estimate, count)
        
    p[[i]] <- ggplot(temp, aes(x = day)) +
        geom_line(aes(y = mn, colour = stage, linetype = model)) +
        geom_ribbon(aes(ymin = LCI, ymax = UCI, fill = stage), data = filter(temp, model == "MW"), alpha = 0.5)
    ggsave(paste0("test", i, ".pdf"), p[[i]])
}

## check NGM
print(paste0("NGM R0 = ", NGM(R0 = NA, nu = pars["nu"], C = contact, N = N, 
    nuA = pars["nuA"], gammaE = pars["gammaE"], pEA = pars["pEA"], gammaA = pars["gammaA"], 
    gammaP = pars["gammaP"], gammaI1 = pars["gammaI1"], pI1I2 = pars["pI1I2"], 
    gammaI2 = pars["gammaI2"])))
    
## check against SEIR without using age-specific NGM
F <- matrix(c(0, 0, 0, pars["nu"], 0, 0, pars["nu"], 0, 0), 3, 3)
V <- matrix(c(pars["gammaE"], -pars["gammaE"], 0, 0, pars["gammaP"], -pars["gammaP"], 0, 0, pars["gammaI1"]), 3, 3)
K <- F %*% solve(V)
print(paste0("Simple R0 = ", max(eigen(K)$values)))

## solve balance equation for final size
balance <- function(z, R0, epsilon) {
    abs(1 - z - (1 - epsilon) * exp(-R0 * z))
}
finalsize <- optimise(balance, interval = c(0, 1), R0 = 3, epsilon = 10 / (N[1] - 10))
finalsize <- finalsize$minimum * N[1]

## add final size to plot
p[[1]] <- p[[1]] + geom_hline(yintercept = finalsize, linetype = "dashed")  
ggsave(paste0("test1.pdf"), p[[1]]) 

## run SimBIID stochastic model
source("stochModel.R")
    
## wrangle outputs to correct format
tempDet <- out %>%
    select(-starts_with("S")) %>%
    select(time, ends_with("1")) %>%
    mutate_all(as.numeric) %>%
    rename(day = time) %>%
    mutate(estimate = "mn")
colnames(tempDet) <- gsub("1$", "", colnames(tempDet))
tempDet <- select(tempDet, day, estimate, E, P, I1, DI)
sims <- sims %>%
    pluck("runs") %>%
    select(-rep) %>%
    rename(day = t) %>%
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
disSims <- map(disSims, as_tibble) %>%
    bind_rows(.id = "rep") %>%
    rename(day = V1, S = V2, E = V3, P = V4, I1 = V5, DI = V6) %>%
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
    select(day, estimate, E, P, I1, DI) %>%
    ungroup()
##filter(disSims, estimate == "mn") %>% select(-estimate) %>% gather(stage, count, -day) %>% ggplot(aes(x = day, y = count, colour = stage)) + geom_line()

## bind model outputs and plot
temp <- list(det = tempDet, stoch = sims, disStoch = disSims, MW = rec[[1]]) %>%
    bind_rows(.id = "model") %>%
    gather(stage, count, -model, -estimate, -day) %>%
    spread(estimate, count)
p <- ggplot(temp, aes(x = day)) +
    geom_line(aes(y = mn, colour = stage, linetype = model)) +
    geom_ribbon(aes(ymin = LCI, ymax = UCI, fill = stage), data = filter(temp, model == "disStoch"), alpha = 0.5)
ggsave(paste0("test1nonMW.pdf"), p)

