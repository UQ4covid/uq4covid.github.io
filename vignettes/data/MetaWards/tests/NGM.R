## NGM (order: E, A, P, I1, I2)
NGM <- function(R0 = NA, nu = NA, C, S0, N, nuA, gammaE, pEA, gammaA, gammaP, gammaI1, pI1I2, gammaI2) {
    
    ## check that exactly one of R0 or nu is specified
    stopifnot(!is.na(R0) | !is.na(nu))
    stopifnot(is.na(R0) | is.na(nu))
    
    ## check inputs
    stopifnot(!missing(C) & !missing(N) & !missing(S0))
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
    
    ## extract lengths
    nage <- length(N)
    stopifnot(all(dim(C) == nage))
    stopifnot(length(S0) == nage & all(S0 <= N))
    
    ## check for empty pathways
    stopifnot(gammaE > 0)
    pathA <- ifelse(pEA == 0, FALSE, TRUE)
    if(pathA) stopifnot(gammaA > 0)
    pathI1 <- ifelse(pEA == 1, FALSE, TRUE)
    if(pathI1) stopifnot(gammaP > 0 & gammaI1 > 0)
    pathI2 <- ifelse(pI1I2 == 0 | !pathI1, FALSE, TRUE)
    if(pathI2) stopifnot(gammaI2 > 0)
    
    ## set up empty F matrix
    ## (order: E, A, P, I1, I2)
    F <- matrix(0, nage * 5, nage * 5)
    
    ## add A to E
    F[1:nage, (nage + 1):(2 * nage)] <- t(t(nuA * S0 * C) / ifelse(N == 0, 1, N))
    ## add P to E
    F[1:nage, (2 * nage + 1):(3 * nage)] <- t(t(S0 * C) / ifelse(N == 0, 1, N))
    ## add I1 to E
    F[1:nage, (3 * nage + 1):(4 * nage)] <- t(t(S0 * C) / ifelse(N == 0, 1, N))
    ## add I2 to E
    F[1:nage, (4 * nage + 1):(5 * nage)] <- t(t(S0 * C) / ifelse(N == 0, 1, N))
    
    ## set up empty V matrix
    V <- matrix(0, nage * 5, nage * 5)
    
    ## add E to E terms
    diag(V)[1:nage] <- gammaE
    
    ## add E to A terms
    V[(nage + 1):(2 * nage), 1:nage] <- diag(nage) * -pEA * gammaE
    ## add A to A terms
    V[(nage + 1):(2 * nage), (nage + 1):(2 * nage)] <- diag(nage) * gammaA
    
    ## add E to P terms
    V[(2 * nage + 1):(3 * nage), 1:nage] <- diag(nage) * -(1 - pEA) * gammaE
    ## add P to P terms
    V[(2 * nage + 1):(3 * nage), (2 * nage + 1):(3 * nage)] <- diag(nage) * gammaP
    
    ## add P to I1 terms
    V[(3 * nage + 1):(4 * nage), (2 * nage + 1):(3 * nage)] <- diag(nage) * -gammaP
    ## add I1 to I1 terms
    V[(3 * nage + 1):(4 * nage), (3 * nage + 1):(4 * nage)] <- diag(nage) * gammaI1
    
    ## add I1 to I2 terms
    V[(4 * nage + 1):(5 * nage), (3 * nage + 1):(4 * nage)] <- diag(nage) * -pI1I2 * gammaI1
    ## add I2 to I2 terms
    V[(4 * nage + 1):(5 * nage), (4 * nage + 1):(5 * nage)] <- diag(nage) * gammaI2
    
    ## remove pathways that don't exist
    rem <- NULL
    if(!pathA) rem <- (nage + 1):(2 * nage)
    if(!pathI1) rem <- c(rem, (2 * nage + 1):(5 * nage))
    if(!pathI2) rem <- c(rem, (4 * nage + 1):(5 * nage))
    if(!is.null(rem)) {
        rem <- unique(rem)
        F <- F[-rem, -rem]
        V <- V[-rem, -rem]
    }
    
    if(!is.na(R0)) {
        ## calculate NGM
        K <- F %*% solve(V)
        
        ## calculate and return nu
        nu <- R0 / max(eigen(K)$values)
        return(list(nu = nu, K = nu * K))
    } else {        
        ## calculate NGM
        K <- nu * F %*% solve(V)
        
        ## calculate and return nu
        R0 <- max(eigen(K)$values)
        return(list(R0 = R0, K = K))
    }
}

