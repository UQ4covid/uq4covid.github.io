## log-sum-exp function
log_sum_exp <- function(x, mn = FALSE) {
    maxx <- max(x)
    y <- maxx + log(sum(exp(x - maxx)))
    if(mn) y <- y - log(length(x))
    y
}

## function to check counts (mainly useful for error checking)
## u: matrix of counts in order:
##       S, E, A, RA, P, I1, DI, I2, RI, H, RH, DH
## cu: matrix of cumulative counts
## N: population size
checkCounts <- function(u, cu, N) {
    
    ## check states match population size
    stopifnot(all(rowSums(u) == N))
    
    ## check states are positive
    stopifnot(all(u >= 0) & all(cu >= 0))
    
    ## check cumulative counts match up
    stopifnot(all((cu[, 1] + cu[, 2]) == N))
    stopifnot(all((cu[, 2] - rowSums(cu[, c(3, 5)])) >= 0))
    stopifnot(all((cu[, 3] - cu[, 4]) >= 0))
    stopifnot(all((cu[, 5] - cu[, 6]) >= 0))
    stopifnot(all((cu[, 6] - rowSums(cu[, c(7, 8, 10)])) >= 0))
    stopifnot(all((cu[, 8] - cu[, 9]) >= 0))
    stopifnot(all((cu[, 10] - rowSums(cu[, c(11, 12)])) >= 0))
    
    # ## check infective states and new incidence are valid
    ## KEEPING FOR POSTERITY BUT DON'T NEED DUE TO IMPORTS
    # Einc <- cu[, 2] - cuprev[, 2]
    # print(which(Einc > 0))
    # stopifnot(all(rowSums(uprev[Einc > 0, c(3, 5, 6, 8)]) > 0))
}

## run model for each set of design points in pars
## pars: matrix / data frame of parameters in order: 
##       nu, nuA, pE, pEP, pA, pP, pI1, pI1H, pI1D, pI2, pH, pHD
## data: data frame of data to fit to
## u:    vector of initial states in order:
##       S, E, A, RA, P, I1, DI, I2, RI, H, RH, DH
## ndays: the number of days to fit to
## npart: the number of particles
## MD:    turn model discrepancy process on (TRUE) or off (FALSE)
## obsScale: scaling parameter for Poisson observation process (see code)
## disScale: scaling parameter for model discrepancy variance (see code)
## whichSave: a design point from which to save the particles over time (useful just
##      to check some runs, but not returned if whichSave = NA)

PF <- function(pars, data, u, ndays, npart = 10, MD = TRUE, obsScale = 1, disScale = 1, whichSave = NA) {
    runs <- mclapply(1:nrow(pars), function(j, pars, u, npart, ndays, data, MD, obsScale, disScale, whichSave) {
        
        ## set initial log-likelihood
        ll <- 0
        
        ## set population size (used as simple check on processes)
        N <- sum(u)
        
        ## set up particles
        u <- matrix(rep(u, npart), npart, byrow = TRUE)
        
        ## capture cumulative counts (used in MD process)
        cu <- u
        
        ## set up proposed particles
        disSims <- matrix(NA, nrow(u), ncol(u))
        
        ## vector of particle weights
        weights <- numeric(npart)
        
        ## temporary list to store outputs if !is.na(whichSave)
        tempOut <- list()
        
        ## loop over time
        for(t in 1:ndays) {
            
            ## loop over particles
            for(i in 1:npart) {
                
                ## run model
                disSims[i, ] <- discreteStochModel(unlist(pars[j, ]), t - 1, t, u[i, ])[2, -1]
                
                ## cols: c("S", "E", "A", "RA", "P", "I1", "DI", "I2", "RI", "H", "RH", "DH")
                ##          1,   2,   3,   4,    5,   6,    7,    8,    9,    10,  11,   12
                # print(c(j, t, i))
                # if(j == 1 & t == 19 & i == 1) browser()
                
                ## adjust states according to model discrepancy
                if(MD) {
                    ## DH (MD on incidence)
                    DHinc <- disSims[i, 12] - u[i, 12]
                    DHinc <- DHinc + trSkellam(1 + disScale * DHinc, 1 + disScale * DHinc, -DHinc, u[i, 10] - DHinc)
                    disSims[i, 12] <- u[i, 12] + DHinc
                    cu[i, 12] <- cu[i, 12] + DHinc
                    
                    ## RH given DH (MD on icidence)
                    RHinc <- disSims[i, 11] - u[i, 11]
                    RHinc <- RHinc + trSkellam(1 + disScale * RHinc, 1 + disScale * RHinc, -RHinc, u[i, 10] - DHinc - RHinc)
                    disSims[i, 11] <- u[i, 11] + RHinc
                    cu[i, 11] <- cu[i, 11] + RHinc
                    
                    ## H
                    disSims[i, 10] <- disSims[i, 10] + trSkellam(1 + disScale * disSims[i, 10], 1 + disScale * disSims[i, 10], -disSims[i, 10], u[i, 10] + u[i, 6] - DHinc - RHinc - disSims[i, 10])  
                    Hinc <- disSims[i, 10] - u[i, 10] + DHinc + RHinc
                    cu[i, 10] <- cu[i, 10] + Hinc
                    
                    ## DI given H (MD on incidence)
                    DIinc <- disSims[i, 7] - u[i, 7]
                    DIinc <- DIinc + trSkellam(1 + disScale * DIinc, 1 + disScale * DIinc, -DIinc, u[i, 6] - Hinc - DIinc)
                    disSims[i, 7] <- u[i, 7] + DIinc
                    cu[i, 7] <- cu[i, 7] + DIinc            
                    
                    ## RI (MD on incidence)
                    RIinc <- disSims[i, 9] - u[i, 9]
                    RIinc <- RIinc + trSkellam(1 + disScale * RIinc, 1 + disScale * RIinc, -RIinc, u[i, 8] - RIinc)
                    disSims[i, 9] <- u[i, 9] + RIinc
                    cu[i, 9] <- cu[i, 9] + RIinc
                    
                    ## I2 given H, RI and DI
                    disSims[i, 8] <- disSims[i, 8] + trSkellam(1 + disScale * disSims[i, 8], 1 + disScale * disSims[i, 8], -disSims[i, 8], u[i, 8] + u[i, 6] - Hinc - DIinc - RIinc - disSims[i, 8])
                    I2inc <- disSims[i, 8] - u[i, 8] + RIinc
                    cu[i, 8] <- cu[i, 8] + I2inc
                    
                    ## I1 given later
                    disSims[i, 6] <- disSims[i, 6] + trSkellam(1 + disScale * disSims[i, 6], 1 + disScale * disSims[i, 6], -disSims[i, 6], u[i, 6] + u[i, 5] - I2inc - Hinc - DIinc - disSims[i, 6])
                    I1inc <- disSims[i, 6] - u[i, 6] + DIinc + Hinc + I2inc
                    cu[i, 6] <- cu[i, 6] + I1inc
                    
                    ## P given later
                    disSims[i, 5] <- disSims[i, 5] + trSkellam(1 + disScale * disSims[i, 5], 1 + disScale * disSims[i, 5], -disSims[i, 5], u[i, 5] + u[i, 2] - I1inc - disSims[i, 5])
                    Pinc <- disSims[i, 5] - u[i, 5] + I1inc
                    cu[i, 5] <- cu[i, 5] + Pinc
                    
                    ## RA (MD on incidence)
                    RAinc <- disSims[i, 4] - u[i, 4]
                    RAinc <- RAinc + trSkellam(1 + disScale * RAinc, 1 + disScale * RAinc, -RAinc, u[i, 3] - RAinc)
                    disSims[i, 4] <- u[i, 4] + RAinc
                    cu[i, 4] <- cu[i, 4] + RAinc
                    
                    ## A given later
                    disSims[i, 3] <- disSims[i, 3] + trSkellam(1 + disScale * disSims[i, 3], 1 + disScale * disSims[i, 3], -disSims[i, 3], u[i, 3] + u[i, 2] - Pinc - RAinc - disSims[i, 3])
                    Ainc <- disSims[i, 3] - u[i, 3] + RAinc
                    cu[i, 3] <- cu[i, 3] + Ainc
                    
                    ## E
                    disSims[i, 2] <- disSims[i, 2] + trSkellam(1 + disScale * disSims[i, 2], 1 + disScale * disSims[i, 2], -disSims[i, 2], u[i, 2] + u[i, 1] - Ainc - Pinc - disSims[i, 2])
                    Einc <- disSims[i, 2] - u[i, 2] + Ainc + Pinc
                    cu[i, 2] <- cu[i, 2] + Einc
                    
                    ## S
                    disSims[i, 1] <- u[i, 1] - Einc
                    cu[i, 1] <- cu[i, 1] - Einc
                } else {
                    ## DH
                    DHinc <- disSims[i, 12] - u[i, 12]
                    cu[i, 12] <- cu[i, 12] + DHinc
                    
                    ## RH 
                    RHinc <- disSims[i, 11] - u[i, 11]
                    cu[i, 11] <- cu[i, 11] + RHinc
                    
                    ## H
                    Hinc <- disSims[i, 10] - u[i, 10] + DHinc + RHinc
                    cu[i, 10] <- cu[i, 10] + Hinc
                    
                    ## DI 
                    DIinc <- disSims[i, 7] - u[i, 7]
                    cu[i, 7] <- cu[i, 7] + DIinc            
                    
                    ## RI
                    RIinc <- disSims[i, 9] - u[i, 9]
                    cu[i, 9] <- cu[i, 9] + RIinc
                    
                    ## I2 
                    I2inc <- disSims[i, 8] - u[i, 8] + RIinc
                    cu[i, 8] <- cu[i, 8] + I2inc
                    
                    ## I1 
                    I1inc <- disSims[i, 6] - u[i, 6] + DIinc + Hinc + I2inc
                    cu[i, 6] <- cu[i, 6] + I1inc
                    
                    ## P given later
                    Pinc <- disSims[i, 5] - u[i, 5] + I1inc
                    cu[i, 5] <- cu[i, 5] + Pinc
                    
                    ## RA
                    RAinc <- disSims[i, 4] - u[i, 4]
                    cu[i, 4] <- cu[i, 4] + RAinc
                    
                    ## A 
                    Ainc <- disSims[i, 3] - u[i, 3] + RAinc
                    cu[i, 3] <- cu[i, 3] + Ainc
                    
                    ## E
                    Einc <- disSims[i, 2] - u[i, 2] + Ainc + Pinc
                    cu[i, 2] <- cu[i, 2] + Einc
                    
                    ## S
                    disSims[i, 1] <- u[i, 1] - Einc
                    cu[i, 1] <- cu[i, 1] - Einc
                }
                
                ## checks
                if(!all(rowSums(disSims[1:i, , drop = FALSE]) == N)) {
                    stop("Error in discrepancies")
                }
                
                ## calculate log observation error weights
                obsDiffs <- c(
                    DIinc - (data$DIobs[t] - ifelse(t > 1, data$DIobs[t - 1], 0)),
                    DHinc - (data$DHobs[t] - ifelse(t > 1, data$DHobs[t - 1], 0))
                )  
                weights[i] <- sum(dpois(obsDiffs, 1 + obsScale * c(DIinc, DHinc), log = TRUE))
            }
            checkCounts(disSims, cu, N)
            tempOut[[t]] <- disSims
            
            ## calculate log-likelihood contribution
            ll <- ll + log_sum_exp(weights, mn = TRUE)
            
            ## if zero likelihood then return
            if(!is.finite(ll)) {
                return(ll)
            }
            
            ## normalise weights
            weights <- exp(weights - log_sum_exp(weights))
            
            ## resample
            inds <- apply(rmultinom(npart, 1, weights), 2, function(x) which(x == 1))
            u <- disSims[inds, ]
            cu <- cu[inds, ]
        }
        if(!is.na(whichSave)) {
            if(j == whichSave) {
                assign("sims", tempOut, pos = 1)
            }
        }
        ll
    }, pars = pars, u = u, npart = npart, ndays = ndays, data = data, MD = MD, obsScale = obsScale, disScale = disScale, whichSave = whichSave, mc.cores = detectCores())
    runs <- do.call("c", runs)
    runs
}
