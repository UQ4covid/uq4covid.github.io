## load libraries
library(tidyverse)
library(deSolve)
        
## write deterministic model
detModel <- function(time, state, parameters, contact, N) {
    
    with(as.list(c(state, parameters)), {
    
        states <- c(P1 + I11 + I21,
            P2 + I12 + I22,
            P3 + I13 + I23,
            P4 + I14 + I24,
            P5 + I15 + I25,
            P6 + I16 + I26,
            P7 + I17 + I27,
            P8 + I18 + I28)
        statesA <- c(A1, A2, A3, A4, A5, A6, A7, A8)
        
        N <- ifelse(N == 0, 1, N)
            
        beta <- nu * contact %*% ((states + nuA * statesA) / N)
        
        dS1 <- -beta[1] * S1
        dS2 <- -beta[2] * S2
        dS3 <- -beta[3] * S3
        dS4 <- -beta[4] * S4
        dS5 <- -beta[5] * S5
        dS6 <- -beta[6] * S6
        dS7 <- -beta[7] * S7
        dS8 <- -beta[8] * S8
        
        dE1 <- beta[1] * S1 - gammaE1 * E1
        dE2 <- beta[2] * S2 - gammaE2 * E2
        dE3 <- beta[3] * S3 - gammaE3 * E3
        dE4 <- beta[4] * S4 - gammaE4 * E4
        dE5 <- beta[5] * S5 - gammaE5 * E5
        dE6 <- beta[6] * S6 - gammaE6 * E6
        dE7 <- beta[7] * S7 - gammaE7 * E7
        dE8 <- beta[8] * S8 - gammaE8 * E8
        
        dA1 <- gammaE1 * (1 - pEP1) * E1 - gammaA1 * A1
        dA2 <- gammaE2 * (1 - pEP2) * E2 - gammaA2 * A2
        dA3 <- gammaE3 * (1 - pEP3) * E3 - gammaA3 * A3
        dA4 <- gammaE4 * (1 - pEP4) * E4 - gammaA4 * A4
        dA5 <- gammaE5 * (1 - pEP5) * E5 - gammaA5 * A5
        dA6 <- gammaE6 * (1 - pEP6) * E6 - gammaA6 * A6
        dA7 <- gammaE7 * (1 - pEP7) * E7 - gammaA7 * A7
        dA8 <- gammaE8 * (1 - pEP8) * E8 - gammaA8 * A8
        
        dRA1 <- gammaA1 * A1
        dRA2 <- gammaA2 * A2
        dRA3 <- gammaA3 * A3
        dRA4 <- gammaA4 * A4
        dRA5 <- gammaA5 * A5
        dRA6 <- gammaA6 * A6
        dRA7 <- gammaA7 * A7
        dRA8 <- gammaA8 * A8
        
        dP1 <- gammaE1 * pEP1 * E1 - gammaP1 * P1
        dP2 <- gammaE2 * pEP2 * E2 - gammaP2 * P2
        dP3 <- gammaE3 * pEP3 * E3 - gammaP3 * P3
        dP4 <- gammaE4 * pEP4 * E4 - gammaP4 * P4
        dP5 <- gammaE5 * pEP5 * E5 - gammaP5 * P5
        dP6 <- gammaE6 * pEP6 * E6 - gammaP6 * P6
        dP7 <- gammaE7 * pEP7 * E7 - gammaP7 * P7
        dP8 <- gammaE8 * pEP8 * E8 - gammaP8 * P8
        
        dI11 <- gammaP1 * P1 - gammaI11 * I11
        dI12 <- gammaP2 * P2 - gammaI12 * I12
        dI13 <- gammaP3 * P3 - gammaI13 * I13
        dI14 <- gammaP4 * P4 - gammaI14 * I14
        dI15 <- gammaP5 * P5 - gammaI15 * I15
        dI16 <- gammaP6 * P6 - gammaI16 * I16
        dI17 <- gammaP7 * P7 - gammaI17 * I17
        dI18 <- gammaP8 * P8 - gammaI18 * I18
        
        dI21 <- (1 - pI1H1 - pI1D1) * gammaI11 * I11 - gammaI21 * I21
        dI22 <- (1 - pI1H2 - pI1D2) * gammaI12 * I12 - gammaI22 * I22
        dI23 <- (1 - pI1H3 - pI1D3) * gammaI13 * I13 - gammaI23 * I23
        dI24 <- (1 - pI1H4 - pI1D4) * gammaI14 * I14 - gammaI24 * I24
        dI25 <- (1 - pI1H5 - pI1D5) * gammaI15 * I15 - gammaI25 * I25
        dI26 <- (1 - pI1H6 - pI1D6) * gammaI16 * I16 - gammaI26 * I26
        dI27 <- (1 - pI1H7 - pI1D7) * gammaI17 * I17 - gammaI27 * I27
        dI28 <- (1 - pI1H8 - pI1D8) * gammaI18 * I18 - gammaI28 * I28
        
        dH1 <- pI1H1 * gammaI11 * I11 - gammaH1 * H1
        dH2 <- pI1H2 * gammaI12 * I12 - gammaH2 * H2
        dH3 <- pI1H3 * gammaI13 * I13 - gammaH3 * H3
        dH4 <- pI1H4 * gammaI14 * I14 - gammaH4 * H4
        dH5 <- pI1H5 * gammaI15 * I15 - gammaH5 * H5
        dH6 <- pI1H6 * gammaI16 * I16 - gammaH6 * H6
        dH7 <- pI1H7 * gammaI17 * I17 - gammaH7 * H7
        dH8 <- pI1H8 * gammaI18 * I18 - gammaH8 * H8                 
        
        dDI1 <- pI1D1 * gammaI11 * I11
        dDI2 <- pI1D2 * gammaI12 * I12
        dDI3 <- pI1D3 * gammaI13 * I13
        dDI4 <- pI1D4 * gammaI14 * I14
        dDI5 <- pI1D5 * gammaI15 * I15
        dDI6 <- pI1D6 * gammaI16 * I16
        dDI7 <- pI1D7 * gammaI17 * I17
        dDI8 <- pI1D8 * gammaI18 * I18              
        
        dRI1 <- gammaI21 * I21
        dRI2 <- gammaI22 * I22
        dRI3 <- gammaI23 * I23
        dRI4 <- gammaI24 * I24
        dRI5 <- gammaI25 * I25
        dRI6 <- gammaI26 * I26
        dRI7 <- gammaI27 * I27
        dRI8 <- gammaI28 * I28 
        
        dRH1 <- (1 - pHD1) * gammaH1 * H1
        dRH2 <- (1 - pHD2) * gammaH2 * H2
        dRH3 <- (1 - pHD3) * gammaH3 * H3
        dRH4 <- (1 - pHD4) * gammaH4 * H4
        dRH5 <- (1 - pHD5) * gammaH5 * H5
        dRH6 <- (1 - pHD6) * gammaH6 * H6
        dRH7 <- (1 - pHD7) * gammaH7 * H7
        dRH8 <- (1 - pHD8) * gammaH8 * H8 
        
        dDH1 <- pHD1 * gammaH1 * H1
        dDH2 <- pHD2 * gammaH2 * H2
        dDH3 <- pHD3 * gammaH3 * H3
        dDH4 <- pHD4 * gammaH4 * H4
        dDH5 <- pHD5 * gammaH5 * H5
        dDH6 <- pHD6 * gammaH6 * H6
        dDH7 <- pHD7 * gammaH7 * H7
        dDH8 <- pHD8 * gammaH8 * H8
        
        return(list(c(
            dS1, dS2, dS3, dS4, dS5, dS6, dS7, dS8,
            dE1, dE2, dE3, dE4, dE5, dE6, dE7, dE8,
            dA1, dA2, dA3, dA4, dA5, dA6, dA7, dA8,
            dRA1, dRA2, dRA3, dRA4, dRA5, dRA6, dRA7, dRA8,
            dP1, dP2, dP3, dP4, dP5, dP6, dP7, dP8,
            dI11, dI12, dI13, dI14, dI15, dI16, dI17, dI18,
            dI21, dI22, dI23, dI24, dI25, dI26, dI27, dI28,
            dRI1, dRI2, dRI3, dRI4, dRI5, dRI6, dRI7, dRI8,
            dDI1, dDI2, dDI3, dDI4, dDI5, dDI6, dDI7, dDI8,
            dH1, dH2, dH3, dH4, dH5, dH6, dH7, dH8,
            dRH1, dRH2, dRH3, dRH4, dRH5, dRH6, dRH7, dRH8,
            dDH1, dDH2, dDH3, dDH4, dDH5, dDH6, dDH7, dDH8
        )))
    })
}

## set seeds
nseeds <- 10

## set initial conditions
inits <- c(
    S1 = N[1] - nseeds, S2 = N[2], S3 = N[3], S4 = N[4], S5 = N[5], S6 = N[6], S7 = N[7], S8 = N[8],
    E1 = nseeds, E2 = 0, E3 = 0, E4 = 0, E5 = 0, E6 = 0, E7 = 0, E8 = 0,
    A1 = 0, A2 = 0, A3 = 0, A4 = 0, A5 = 0, A6 = 0, A7 = 0, A8 = 0,
    RA1 = 0, RA2 = 0, RA3 = 0, RA4 = 0, RA5 = 0, RA6 = 0, RA7 = 0, RA8 = 0,
    P1 = 0, P2 = 0, P3 = 0, P4 = 0, P5 = 0, P6 = 0, P7 = 0, P8 = 0,
    I11 = 0, I12 = 0, I13 = 0, I14 = 0, I15 = 0, I16 = 0, I17 = 0, I18 = 0,
    I21 = 0, I22 = 0, I23 = 0, I24 = 0, I25 = 0, I26 = 0, I27 = 0, I28 = 0,
    RI1 = 0, RI2 = 0, RI3 = 0, RI4 = 0, RI5 = 0, RI6 = 0, RI7 = 0, RI8 = 0,
    DI1 = 0, DI2 = 0, DI3 = 0, DI4 = 0, DI5 = 0, DI6 = 0, DI7 = 0, DI8 = 0,
    H1 = 0, H2 = 0, H3 = 0, H4 = 0, H5 = 0, H6 = 0, H7 = 0, H8 = 0,
    RH1 = 0, RH2 = 0, RH3 = 0, RH4 = 0, RH5 = 0, RH6 = 0, RH7 = 0, RH8 = 0,
    DH1 = 0, DH2 = 0, DH3 = 0, DH4 = 0, DH5 = 0, DH6 = 0, DH7 = 0, DH8 = 0
)

## read in pars
pars <- read_delim("disease.dat", delim = " ") %>%
    select(-contains("beta"), -contains("lock_"), 
    nu = `beta[1]`, nuA = `beta[6]`, -repeats, -output)
colnames(pars) <- gsub("\\.", "", colnames(pars))
colnames(pars) <- gsub("_", "", colnames(pars))
pinds <- map(c("pE", "pP", "pA", "pI1", "pI2", "pH"), function(x, nms) {
        which(str_detect(nms, paste0(x, "[1-8]")))
    }, nms = colnames(pars)) %>%
    reduce(c)
pars <- mutate_at(pars, pinds, ~-log(1 - ifelse(. == 1, 0.99, .)))
colnames(pars)[pinds] <- gsub("p", "gamma", colnames(pars)[pinds])  
  
## set contact matrix
contact <- read_csv("contact_matrix.csv", col_names = FALSE) %>%
    as.matrix()
    
## set up time
times <- 1:200

## solve using ODE solver
out <- ode(y = inits, times = times, func = detModel, parms = pars, contact = contact, N = N) %>%
    as.matrix() %>%
    as_tibble()

### example plot    
#select(out, time, ends_with("1")) %>%
#    select(-starts_with("S")) %>%
#    select(which(apply(., 2, sum) != 0)) %>%
#    gather(stage, count, -time) %>%
#    ggplot(aes(x = time, y = count, colour = stage)) +
#        geom_line()
    
