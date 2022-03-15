## run model for each set of design points in pars
## pars: matrix / data frame of parameters in order: 
##       nu, nuA, pE, pEP, pA, pP, pI1, pI1H, pI1D, pI2, pH, pHD
##       each parameter except nu/nuA have entries for
##       each age-group e.g. nu, nuA, pE1, pE2, pE3 etc.
## C:    contact matrix for mixing between age-classes
## data: data frame of data to fit to, with columns
##       DI1, DI2, ..., DI8, DH1, ..., DH8
## u1_moves: matrix of commuter data, with columns:
##            home LAD, work LAD
## u1: 3D array of commuter data, with dimensions:
##            nclasses x nages x nrow(u1_moves)
## ndays: the number of days to fit to
## npart: the number of particles
## MD:    turn model discrepancy process on (TRUE) or off (FALSE)
## obsScale: scaling parameter for Poisson observation process (see code)
## a1, a2, b: parameters for Skellam observation process
## saveAll: a logical specifying whether to return all states (if FALSE then returns just observed states))
## ncores:  the number of cores for OpenMP parallelisation (if NA then defaults to all available cores)

PF <- function(pars, C, data, u1_moves, u1, ndays, npart = 10, MD = TRUE, a1 = 0.01, a2 = 0.2, b = 0.1, 
               a_dis = 0.05, b_dis = 0.5, saveAll = NA, ncores = NA) {
    
    ## convert indicators
    if(is.na(saveAll)) {
        saveAllint <- 0
    } else {
        saveAllint <- ifelse(saveAll, 2, 1)
    }
    MDint <- ifelse(MD, 1, 0)
    
    ## check number of requested cores
    if(is.na(ncores)) {
        ncores <- parallel::detectCores()
    }
    
    ## run particle filter for each set of inputs
    runs <- lapply(1:nrow(pars), function(k, pars, C, u1_moves, u1, npart, ndays, data, MD, a1, a2, b, a_dis, b_dis, saveAll, ncores) {
        
        ## extract observations
        data <- select(data, t, (starts_with("DI") | starts_with("DH")) & ends_with("obs")) %>%
            {rbind(rep(0, ncol(.)), .)} %>%
            mutate(across(!t, ~. - lag(.))) %>%
            slice(-1) %>%
            as.matrix()
        
        ## remove time since not used in code
        data <- data[, -1]
        
        ## set pars
        pars <- unlist(pars[k, ])
        
        ## run particle filter
        ll <- PF_cpp(pars, C, data, 12L, 8L, 339L, u1_moves, u1, ndays,
            npart, MD, a1, a2, b, a_dis, b_dis, saveAll, ncores)
        ll
    }, pars = pars, C = C, u1_moves = u1_moves, u1 = u1, npart = npart, ndays = ndays, data = data, MD = MDint, 
       a1 = a1, a2 = a2, b = b, a_dis = a_dis, b_dis = b_dis, saveAll = saveAllint, ncores = ncores)
    if(!is.na(saveAll)) {
        ll <- map(runs, "ll")
        runs <- map(runs, "particles") %>%
            map(~{
                x <- list()
                for(i in 1:npart) {
                    x[[i]] <- .[(i - 1) * (ndays + 1) + 1:(ndays + 1)]
                }
                x
            })
        ll <- do.call("c", ll)
        return(list(ll = ll, particles = runs))
    } else {
        ll <- do.call("c", runs)
        return(ll)
    }
}
