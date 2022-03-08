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
  
  ## check states match population size in each age-class
  stopifnot(all((colSums(u) - N) == 0))
  
  ## check states are positive
  stopifnot(all(u >= 0) & all(cu >= 0))
  
  ## check cumulative counts match up
  cu <- t(cu)
  stopifnot(all((cu[, 1] + cu[, 2] - N) == 0))
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
##       each parameter except nu/nuA have entries for
##       each age-group e.g. nu, nuA, pE1, pE2, pE3 etc.
## C:    contact matrix for mixing between age-classes
## data: data frame of data to fit to, with columns
##       DI1, DI2, ..., DI8, DH1, ..., DH8
## u:    matrix of initial states (rows) in order:
##       S, E, A, RA, P, I1, DI, I2, RI, H, RH, DH
##       with age-groups along the columns
## ndays: the number of days to fit to
## npart: the number of particles
## MD:    turn model discrepancy process on (TRUE) or off (FALSE)
## obsScale: scaling parameter for Poisson observation process (see code)
## disScale: scaling parameter for model discrepancy variance (see code)
## whichSave: a design point from which to save the particles over time (useful just
##      to check some runs, but not returned if whichSave = NA)

#AMENDED PF FUNCTION TO RUN A SINGLE SIMUALTION, ADDING DISCREPANCY AT EACH TIME
#NEED TO REMOVE, OTHER DESIGN POINTS, PARTICLES (MAYBE), REWEIGHTING, ETC.
perfectModel <- function(pars, C, u, ndays, npart = 10, MD = TRUE, a_dis, b_dis) {
  
    ## set population size vector (used as simple check on processes)
    N <- colSums(u)
    
    ## set up particles
    nages <- ncol(u)
    u <- rep(list(u), npart)
    
    ## capture cumulative counts (used in MD process)
    cu <- u
    
    ## set up proposed particles
    disSims <- u
    
    ## vector of particle weights
    weights <- numeric(npart)
    
    ## temporary list to store outputs if !is.na(whichSave)
    tempOut <- list()
    
    ## loop over time
    for(t in 1:ndays) {
      
      ## loop over particles
      for(i in 1:npart) {
        
        ## run model
        temp <- discreteStochModel(unlist(pars), t - 1, t, u[[i]], C)[2, -1]
        disSims[[i]] <- matrix(temp, ncol = length(N), byrow = TRUE)
        
        ## cols: c("S", "E", "A", "RA", "P", "I1", "DI", "I2", "RI", "H", "RH", "DH")
        ##          1,   2,   3,   4,    5,   6,    7,    8,    9,    10,  11,   12
        # print(c(j, t, i))
        # if(j == 1 & t == 19 & i == 1) browser()
        
        ## set weights
        #weights[i] <- 0
        
        ## loop over age-classes
        for(j in 1:nages) {
          ## adjust states according to model discrepancy
          if(MD) {
            ## DH (MD on incidence)
            DHinc <- disSims[[i]][12, j] - u[[i]][12, j]
            DHinc <- DHinc + rtskellam(1, a_dis + b_dis * DHinc, a_dis + b_dis * DHinc, -DHinc, u[[i]][10, j] - DHinc)
            disSims[[i]][12, j] <- u[[i]][12, j] + DHinc
            cu[[i]][12, j] <- cu[[i]][12, j] + DHinc
            
            ## RH given DH (MD on icidence)
            RHinc <- disSims[[i]][11, j] - u[[i]][11, j]
            RHinc <- RHinc + rtskellam(1, a_dis + b_dis * RHinc, a_dis + b_dis * RHinc, -RHinc, u[[i]][10, j] - DHinc - RHinc)
            disSims[[i]][11, j] <- u[[i]][11, j] + RHinc
            cu[[i]][11, j] <- cu[[i]][11, j] + RHinc
            
            ## H
            disSims[[i]][10, j] <- disSims[[i]][10, j] + rtskellam(1, a_dis + b_dis * disSims[[i]][10, j], 
                                                                   a_dis + b_dis * disSims[[i]][10, j], 
                                                                   -disSims[[i]][10, j], u[[i]][10, j] + u[[i]][6, j] - DHinc - RHinc - 
                                                                     disSims[[i]][10, j])  
            Hinc <- disSims[[i]][10, j] - u[[i]][10, j] + DHinc + RHinc
            cu[[i]][10, j] <- cu[[i]][10, j] + Hinc
            
            ## DI given H (MD on incidence)
            DIinc <- disSims[[i]][7, j] - u[[i]][7, j]
            DIinc <- DIinc + rtskellam(1, a_dis + b_dis * DIinc, a_dis + b_dis * DIinc, -DIinc, u[[i]][6, j] - Hinc - DIinc)
            disSims[[i]][7, j] <- u[[i]][7, j] + DIinc
            cu[[i]][7, j] <- cu[[i]][7, j] + DIinc            
            
            ## RI (MD on incidence)
            RIinc <- disSims[[i]][9, j] - u[[i]][9, j]
            RIinc <- RIinc + rtskellam(1, a_dis + b_dis * RIinc, a_dis + b_dis * RIinc, -RIinc, u[[i]][8, j] - RIinc)
            disSims[[i]][9, j] <- u[[i]][9, j] + RIinc
            cu[[i]][9, j] <- cu[[i]][9, j] + RIinc
            
            ## I2 given H, RI and DI
            disSims[[i]][8, j] <- disSims[[i]][8, j] + 
              rtskellam(1, a_dis + b_dis * disSims[[i]][8, j], 
                        a_dis + b_dis * disSims[[i]][8, j], 
                        -disSims[[i]][8, j], u[[i]][8, j] + u[[i]][6, j] - Hinc - DIinc - RIinc - 
                          disSims[[i]][8, j])
            I2inc <- disSims[[i]][8, j] - u[[i]][8, j] + RIinc
            cu[[i]][8, j] <- cu[[i]][8, j] + I2inc
            
            ## I1 given later
            disSims[[i]][6, j] <- disSims[[i]][6, j] + 
              rtskellam(1, a_dis + b_dis * disSims[[i]][6, j], 
                        a_dis + b_dis * disSims[[i]][6, j], -disSims[[i]][6, j], 
                        u[[i]][6, j] + u[[i]][5, j] - I2inc - Hinc - DIinc - disSims[[i]][6, j])
            I1inc <- disSims[[i]][6, j] - u[[i]][6, j] + DIinc + Hinc + I2inc
            cu[[i]][6, j] <- cu[[i]][6, j] + I1inc
            
            ## P given later
            disSims[[i]][5, j] <- disSims[[i]][5, j] + 
              rtskellam(1, a_dis + b_dis * disSims[[i]][5, j], a_dis + b_dis * disSims[[i]][5, j], 
                        -disSims[[i]][5, j], u[[i]][5, j] + u[[i]][2, j] - I1inc - 
                          disSims[[i]][5, j])
            Pinc <- disSims[[i]][5, j] - u[[i]][5, j] + I1inc
            cu[[i]][5, j] <- cu[[i]][5, j] + Pinc
            
            ## RA (MD on incidence)
            RAinc <- disSims[[i]][4, j] - u[[i]][4, j]
            RAinc <- RAinc + rtskellam(1, a_dis + b_dis * RAinc, a_dis + b_dis * RAinc, -RAinc, u[[i]][3, j] - RAinc)
            disSims[[i]][4, j] <- u[[i]][4, j] + RAinc
            cu[[i]][4, j] <- cu[[i]][4, j] + RAinc
            
            ## A given later
            disSims[[i]][3, j] <- disSims[[i]][3, j] + 
              rtskellam(1, a_dis + b_dis * disSims[[i]][3, j], 
                        a_dis + b_dis * disSims[[i]][3, j], -disSims[[i]][3, j], 
                        u[[i]][3, j] + u[[i]][2, j] - Pinc - RAinc - disSims[[i]][3, j])
            Ainc <- disSims[[i]][3, j] - u[[i]][3, j] + RAinc
            cu[[i]][3, j] <- cu[[i]][3, j] + Ainc
            
            ## E
            disSims[[i]][2, j] <- disSims[[i]][2, j] + 
              rtskellam(1, a_dis + b_dis * disSims[[i]][2, j], a_dis + b_dis * disSims[[i]][2, j], 
                        -disSims[[i]][2, j], u[[i]][2, j] + u[[i]][1, j] - Ainc - Pinc - 
                          disSims[[i]][2, j])
            Einc <- disSims[[i]][2, j] - u[[i]][2, j] + Ainc + Pinc
            cu[[i]][2, j] <- cu[[i]][2, j] + Einc
            
            ## S
            disSims[[i]][1, j] <- u[[i]][1, j] - Einc
            cu[[i]][1, j] <- cu[[i]][1, j] - Einc
          } else {
            ## DH
            DHinc <- disSims[[i]][12, j] - u[[i]][12, j]
            cu[[i]][12, j] <- cu[[i]][12, j] + DHinc
            
            ## RH 
            RHinc <- disSims[[i]][11, j] - u[[i]][11, j]
            cu[[i]][11, j] <- cu[[i]][11, j] + RHinc
            
            ## H
            Hinc <- disSims[[i]][10, j] - u[[i]][10, j] + DHinc + RHinc
            cu[[i]][10, j] <- cu[[i]][10, j] + Hinc
            
            ## DI 
            DIinc <- disSims[[i]][7, j] - u[[i]][7, j]
            cu[[i]][7, j] <- cu[[i]][7, j] + DIinc
            
            ## RI
            RIinc <- disSims[[i]][9, j] - u[[i]][9, j]
            cu[[i]][9, j] <- cu[[i]][9, j] + RIinc
            
            ## I2 
            I2inc <- disSims[[i]][8, j] - u[[i]][8, j] + RIinc
            cu[[i]][8, j] <- cu[[i]][8, j] + I2inc
            
            ## I1 
            I1inc <- disSims[[i]][6, j] - u[[i]][6, j] + DIinc + Hinc + I2inc
            cu[[i]][6, j] <- cu[[i]][6, j] + I1inc
            
            ## P given later
            Pinc <- disSims[[i]][5, j] - u[[i]][5, j] + I1inc
            cu[[i]][5, j] <- cu[[i]][5, j] + Pinc
            
            ## RA
            RAinc <- disSims[[i]][4, j] - u[[i]][4, j]
            cu[[i]][4, j] <- cu[[i]][4, j] + RAinc
            
            ## A 
            Ainc <- disSims[[i]][3, j] - u[[i]][3, j] + RAinc
            cu[[i]][3, j] <- cu[[i]][3, j] + Ainc
            
            ## E
            Einc <- disSims[[i]][2, j] - u[[i]][2, j] + Ainc + Pinc
            cu[[i]][2, j] <- cu[[i]][2, j] + Einc
            
            ## S
            disSims[[i]][1, j] <- u[[i]][1, j] - Einc
            cu[[i]][1, j] <- cu[[i]][1, j] - Einc
          }
          # if(!identical(disSims[[i]], u[[i]])) browser()
          
          ## calculate log observation error weights
          #obsDiffs <- c(
          #  DIinc - (pluck(data, paste0("DI", j, "obs"))[t] - ifelse(t > 1, pluck(data, paste0("DI", j, "obs"))[t - 1], 0)),
#            DHinc - (pluck(data, paste0("DH", j, "obs"))[t] - ifelse(t > 1, pluck(data, paste0("DH", j, "obs"))[t - 1], 0))
 #         )  
  #        weights[i] <- weights[i] + sum(dpois(obsDiffs, 1 + obsScale * c(DIinc, DHinc), log = TRUE))
        }
        # ## check counts
        # checkCounts(disSims[[i]], cu[[i]], N)
      }
      tempOut[[t]] <- disSims
      
      ## calculate log-likelihood contribution
      #ll <- ll + log_sum_exp(weights, mn = TRUE)
      
      ## if zero likelihood then return
      #if(!is.finite(ll)) {
      #  return(ll)
      #}
      
      ## normalise weights
      #weights <- exp(weights - log_sum_exp(weights))
      
      ## resample
      #inds <- apply(rmultinom(npart, 1, weights), 2, function(x) which(x == 1))
      #REPLACE RESAMPLING WITH RETURNING ALL PARTICLES A LA RUNNING THE MODEL NPART TIMES
      u <- disSims[1:npart]
      cu <- cu[1:npart]
    }
    disSims <- transpose(tempOut)
    length(disSims)
    length(disSims[[1]])
    
    ## collapse to  correct form
    disSims <- map(disSims, ~{
      map(., ~as.vector(t(.))) %>%
        reduce(rbind) %>%
        {cbind(1:nrow(.), .)}
    })
  disSims
}
