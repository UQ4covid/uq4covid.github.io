## NGM (order: E, A, P, I1, I2)
NGM <- function(R0 = NA, nu = NA, S0, N, nuA, gammaE, pEP, gammaA, gammaP, gammaI1, pI1H, pI1D, gammaI2) {
  
    ## check that exactly one of R0 or nu is specified
    stopifnot(length(R0) == 1 & length(nu) == 1 & length(nuA) == 1)
    stopifnot(!is.na(R0) | !is.na(nu))
    stopifnot(is.na(R0) | is.na(nu))
    
    ## check inputs
    stopifnot(!missing(N) & !missing(S0))
    if(missing(nuA)) nuA <- 0
    if(missing(gammaE)) gammaE <- 0
    if(missing(pEP)) pEP <- 0
    if(missing(gammaA)) gammaA <- 0
    if(missing(gammaP)) gammaP <- 0
    if(missing(gammaI1)) gammaI1 <- 0
    if(missing(pI1H)) pI1H <- 0
    if(missing(pI1D)) pI1D <- 0
    if(missing(gammaI2)) gammaI2 <- 0
    stopifnot(all(c(nuA, gammaE, pEP, gammaA, 
                    gammaP, gammaI1, pI1H, pI1D, gammaI2) >= 0))
    
    ## extract lengths
    stopifnot(length(N) == 1 & length(S0) == 1)
    stopifnot(S0 <= N & S0 >= 1)
    
    ## check for empty pathways
    stopifnot(gammaE > 0)
    pathA <- ifelse(pEP == 1, FALSE, TRUE)
    if(pathA) stopifnot(gammaA > 0)
    pathI1 <- ifelse(pEP == 0, FALSE, TRUE)
    if(pathI1) stopifnot(gammaP > 0 & gammaI1 > 0)
    pathI2 <- ifelse((1 - pI1H - pI1D) == 0 | !pathI1, FALSE, TRUE)
    if(pathI2) stopifnot(gammaI2 > 0)
    
    ## set up empty F matrix
    ## (order: E, A, P, I1, I2)
    F <- matrix(0, 5, 5)
    
    ## add A to E
    F[1, 2] <- S0 * nuA / N
    ## add P to E
    F[1, 3] <- S0 / N
    ## add I1 to E
    F[1, 4] <- S0 / N
    ## add I2 to E
    F[1, 5] <- S0 / N
    
    ## set up empty V matrix
    V <- matrix(0, 5, 5)
    
    ## add E to E terms
    V[1, 1] <- gammaE
    
    ## add E to A terms
    V[2, 1] <- -(1 - pEP) * gammaE
    ## add A to A terms
    V[2, 2] <- gammaA
    
    ## add E to P terms
    V[3, 1] <- -pEP * gammaE
    ## add P to P terms
    V[3, 3] <- gammaP
    
    ## add P to I1 terms
    V[4, 3] <- -gammaP
    ## add I1 to I1 terms
    V[4, 4] <- gammaI1
    
    ## add I1 to I2 terms
    V[5, 4] <- -(1 - pI1H - pI1D) * gammaI1
    ## add I2 to I2 terms
    V[5, 5] <- gammaI2
    
    ## remove pathways that don't exist
    rem <- NULL
    if(!pathA) rem <- 2
    if(!pathI1) rem <- c(rem, 3:5)
    if(!pathI2) rem <- c(rem, 5)
    if(!is.null(rem)) {
        rem <- unique(rem)
        F <- F[-rem, -rem]
        V <- V[-rem, -rem]
    }
    
    if(!is.na(R0)) {
        ## calculate NGM
        K <- F %*% solve(V)
        
        ## extract eigenvalues
        eig <- eigen(K)$values
        if(is.complex(eig)) {
            ## extract real eigenvalues
            eig <- as.numeric(eig[Im(eig) == 0])
            stopifnot(length(eig) >= 1)
        }
        
        ## calculate and return nu
        nu <- R0 / max(eig)
        return(list(nu = nu, K = nu * K))
    } else {        
        ## calculate NGM
        K <- nu * F %*% solve(V)
        
        ## extract eigenvalues
        eig <- eigen(K)$values
        if(is.complex(eig)) {
            ## extract real eigenvalues
            eig <- as.numeric(eig[Im(eig) == 0])
            stopifnot(length(eig) >= 1)
        }
        
        ## calculate and return R0
        R0 <- max(eig)
        return(list(R0 = R0, K = K))
    }
}
