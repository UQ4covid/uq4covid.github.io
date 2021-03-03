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
        
        dE1 <- beta[1] * S1 - gammaE * E1
        dE2 <- beta[2] * S2 - gammaE * E2
        dE3 <- beta[3] * S3 - gammaE * E3
        dE4 <- beta[4] * S4 - gammaE * E4
        dE5 <- beta[5] * S5 - gammaE * E5
        dE6 <- beta[6] * S6 - gammaE * E6
        dE7 <- beta[7] * S7 - gammaE * E7
        dE8 <- beta[8] * S8 - gammaE * E8
        
        dA1 <- gammaE * (1 - pEP) * E1 - gammaA * A1
        dA2 <- gammaE * (1 - pEP) * E2 - gammaA * A2
        dA3 <- gammaE * (1 - pEP) * E3 - gammaA * A3
        dA4 <- gammaE * (1 - pEP) * E4 - gammaA * A4
        dA5 <- gammaE * (1 - pEP) * E5 - gammaA * A5
        dA6 <- gammaE * (1 - pEP) * E6 - gammaA * A6
        dA7 <- gammaE * (1 - pEP) * E7 - gammaA * A7
        dA8 <- gammaE * (1 - pEP) * E8 - gammaA * A8
        
        dRA1 <- gammaA * A1
        dRA2 <- gammaA * A2
        dRA3 <- gammaA * A3
        dRA4 <- gammaA * A4
        dRA5 <- gammaA * A5
        dRA6 <- gammaA * A6
        dRA7 <- gammaA * A7
        dRA8 <- gammaA * A8
        
        dP1 <- gammaE * pEP * E1 - gammaP * P1
        dP2 <- gammaE * pEP * E2 - gammaP * P2
        dP3 <- gammaE * pEP * E3 - gammaP * P3
        dP4 <- gammaE * pEP * E4 - gammaP * P4
        dP5 <- gammaE * pEP * E5 - gammaP * P5
        dP6 <- gammaE * pEP * E6 - gammaP * P6
        dP7 <- gammaE * pEP * E7 - gammaP * P7
        dP8 <- gammaE * pEP * E8 - gammaP * P8
        
        dI11 <- gammaP * P1 - gammaI1 * I11
        dI12 <- gammaP * P2 - gammaI1 * I12
        dI13 <- gammaP * P3 - gammaI1 * I13
        dI14 <- gammaP * P4 - gammaI1 * I14
        dI15 <- gammaP * P5 - gammaI1 * I15
        dI16 <- gammaP * P6 - gammaI1 * I16
        dI17 <- gammaP * P7 - gammaI1 * I17
        dI18 <- gammaP * P8 - gammaI1 * I18
        
        dI21 <- (1 - pI1H - pI1D) * gammaI1 * I11 - gammaI2 * I21
        dI22 <- (1 - pI1H - pI1D) * gammaI1 * I12 - gammaI2 * I22
        dI23 <- (1 - pI1H - pI1D) * gammaI1 * I13 - gammaI2 * I23
        dI24 <- (1 - pI1H - pI1D) * gammaI1 * I14 - gammaI2 * I24
        dI25 <- (1 - pI1H - pI1D) * gammaI1 * I15 - gammaI2 * I25
        dI26 <- (1 - pI1H - pI1D) * gammaI1 * I16 - gammaI2 * I26
        dI27 <- (1 - pI1H - pI1D) * gammaI1 * I17 - gammaI2 * I27
        dI28 <- (1 - pI1H - pI1D) * gammaI1 * I18 - gammaI2 * I28
        
        dH1 <- pI1H * gammaI1 * I11 - gammaH * H1
        dH2 <- pI1H * gammaI1 * I12 - gammaH * H2
        dH3 <- pI1H * gammaI1 * I13 - gammaH * H3
        dH4 <- pI1H * gammaI1 * I14 - gammaH * H4
        dH5 <- pI1H * gammaI1 * I15 - gammaH * H5
        dH6 <- pI1H * gammaI1 * I16 - gammaH * H6
        dH7 <- pI1H * gammaI1 * I17 - gammaH * H7
        dH8 <- pI1H * gammaI1 * I18 - gammaH * H8                 
        
        dDI1 <- pI1D * gammaI1 * I11
        dDI2 <- pI1D * gammaI1 * I12
        dDI3 <- pI1D * gammaI1 * I13
        dDI4 <- pI1D * gammaI1 * I14
        dDI5 <- pI1D * gammaI1 * I15
        dDI6 <- pI1D * gammaI1 * I16
        dDI7 <- pI1D * gammaI1 * I17
        dDI8 <- pI1D * gammaI1 * I18              
        
        dRI1 <- gammaI2 * I21
        dRI2 <- gammaI2 * I22
        dRI3 <- gammaI2 * I23
        dRI4 <- gammaI2 * I24
        dRI5 <- gammaI2 * I25
        dRI6 <- gammaI2 * I26
        dRI7 <- gammaI2 * I27
        dRI8 <- gammaI2 * I28 
        
        dRH1 <- (1 - pHD) * gammaH * H1
        dRH2 <- (1 - pHD) * gammaH * H2
        dRH3 <- (1 - pHD) * gammaH * H3
        dRH4 <- (1 - pHD) * gammaH * H4
        dRH5 <- (1 - pHD) * gammaH * H5
        dRH6 <- (1 - pHD) * gammaH * H6
        dRH7 <- (1 - pHD) * gammaH * H7
        dRH8 <- (1 - pHD) * gammaH * H8 
        
        dDH1 <- pHD * gammaH * H1
        dDH2 <- pHD * gammaH * H2
        dDH3 <- pHD * gammaH * H3
        dDH4 <- pHD * gammaH * H4
        dDH5 <- pHD * gammaH * H5
        dDH6 <- pHD * gammaH * H6
        dDH7 <- pHD * gammaH * H7
        dDH8 <- pHD * gammaH * H8
        
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

## read in parameters
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
    
## set contact matrix
contact <- read_csv("contact_matrix.csv", col_names = FALSE) %>%
    as.matrix()
    
## set up time
times <- 1:150

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
    
