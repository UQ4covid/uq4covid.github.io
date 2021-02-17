## NGM (order: E, A, P, I1, I2)
NGM <- function(R0 = NA, nu = NA, C, N, nuA, gammaE, pEA, gammaA, gammaP, gammaI1, pI1I2, gammaI2) {
    
    ## check that exactly one of R0 or nu is specified
    stopifnot(!is.na(R0) | !is.na(nu))
    stopifnot(is.na(R0) | is.na(nu))
    
    ## extract lengths
    nage <- length(N)
    stopifnot(all(dim(C) == nage))
    
    ## check inputs
    stopifnot(!missing(C) & !missing(N))
    if(missing(nuA)) nuA <- 0
    if(missing(gammaE)) gammaE <- 0
    if(missing(pEA)) pEA <- 0
    if(missing(gammaA)) gammaA <- 0
    if(missing(gammaP)) gammaP <- 0
    if(missing(gammaI1)) gammaI1 <- 0
    if(missing(pI1I2)) pI1I2 <- 0
    if(missing(gammaI2)) gammaI2 <- 0
    stopifnot(all(c(nuA, gammaE, pEA, gammaA, 
        gammaP, gammaI1, pI1I2, gammaI2) >= 0))
    
    ## check for empty pathways
    stopifnot(gammaE > 0)
    pathA <- ifelse(pEA == 0, FALSE, TRUE)
    if(pathA) stopifnot(gammaA > 0)
    pathI1 <- ifelse(pEA == 1, FALSE, TRUE)
    if(pathI1) stopifnot(gammaP > 0 & gammaI1 > 0)
    pathI2 <- ifelse(pI1I2 == 0 | !pathI1, FALSE, TRUE)
    if(pathI2) stopifnot(gammaI2 > 0)
    
    ## set number of relevant stages
    nstages <- 1
    if(pathA) nstages <- nstages + 1
    if(pathI1) nstages <- nstages + 2
    if(pathI2) nstages <- nstages + 1
        
    ## set up empty F matrix
    F <- matrix(0, nage * nstages, nage * nstages)
    
    stage <- 1
    if(pathA) {
        ## add A to E
        F[1:nage, (nage + stage):((stage + 1) * nage)] <- t(t(nuA * N * C) / ifelse(N == 0, 1, N))
        stage <- stage + 1
    }
    if(pathI1) {
        ## add P to E
        F[1:nage, (stage * nage + 1):((stage + 1) * nage)] <- t(t(N * C) / ifelse(N == 0, 1, N))        
        stage <- stage + 1
        ## add I1 to E
        F[1:nage, (stage * nage + 1):((stage + 1) * nage)] <- t(t(N * C) / ifelse(N == 0, 1, N))
        stage <- stage + 1
    }
    if(pathI2) {
        ## add I2 to E
        F[1:nage, (stage * nage + 1):((stage + 1) * nage)] <- t(t(N * C) / ifelse(N == 0, 1, N))
        stage <- stage + 1
    }
    
    ## set up empty V matrix
    V <- matrix(0, nage * nstages, nage * nstages)
    
    ## add E to E terms
    diag(V)[1:nage] <- gammaE
    
    stage <- 1
    if(pathA) {
        ## add E to A terms
        V[(nage + stage):((stage + 1) * nage), 1:nage] <- diag(nage) * -pEA * gammaE
        ## add A to A terms
        V[(nage + stage):((stage + 1) * nage), (nage + stage):((stage + 1) * nage)] <- diag(nage) * gammaA
        stage <- stage + 1
    }
    if(pathI1) {
        ## add E to P terms
        V[(stage * nage + 1):((stage + 1) * nage), 1:nage] <- diag(nage) * -(1 - pEA) * gammaE
        stagecol <- ifelse(pathA, 2, 1) 
        ## add P to P terms
        V[(stage * nage + 1):((stage + 1) * nage), (stagecol * nage + 1):((stagecol + 1) * nage)] <- diag(nage) * gammaP
        stage <- stage + 1
        ## add P to I1 terms
        V[(stage * nage + 1):((stage + 1) * nage), (stagecol * nage + 1):((stagecol + 1) * nage)] <- diag(nage) * -gammaP
        stagecol <- stagecol + 1
        ## add I1 to I1 terms
        V[(stage * nage + 1):((stage + 1) * nage), (stagecol * nage + 1):((stagecol + 1) * nage)] <- diag(nage) * gammaI1
    }
    if(pathI2) {
        stage <- stage + 1
        ## add I1 to I2 terms
        V[(stage * nage + 1):((stage + 1) * nage), (stagecol * nage + 1):((stagecol + 1) * nage)] <- diag(nage) * -pI1I2 * gammaI1
        stagecol <- stagecol + 1
        ## add I2 to I2 terms
        V[(stage * nage + 1):((stage + 1) * nage), (stagecol * nage + 1):((stagecol + 1) * nage)] <- diag(nage) * gammaI2
    }
    
    if(!is.na(R0)) {
        ## calculate NGM
        K <- F %*% solve(V)
        
        ## calculate and return nu
        nu <- R0 / max(eigen(K)$values)
        return(nu)
    } else {        
        ## calculate NGM
        K <- nu * F %*% solve(V)
        
        ## calculate and return nu
        R0 <- max(eigen(K)$values)
        return(R0)
    }
}

