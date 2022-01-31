## load libraries
library(skellam)

## truncated Skellam density
dtskellam <- function(x, lambda1, lambda2, LB = -Inf, UB = Inf, log = FALSE) {
    
    ## checks on inputs
    if(missing(x) | missing(lambda1) | missing(lambda2)) stop("'x', 'lambda1' and 'lambda2' must be specified")
    if(length(lambda1) != length(x) & length(lambda1) > 1) stop("'lambda1' must be of length 1 or length(x)")
    if(length(lambda2) != length(x) & length(lambda2) > 1) stop("'lambda2' must be of length 1 or length(x)")
    if(length(LB) != length(x) & length(LB) > 1) stop("'LB' must be of length 1 or length(x)")
    if(length(UB) != length(x) & length(UB) > 1) stop("'UB' must be of length 1 or length(x)")
    
    ## expand inputs
    if(length(lambda1) != length(x)) lambda1 <- rep(lambda1, length(x))
    if(length(lambda2) != length(x)) lambda2 <- rep(lambda2, length(x))
    if(length(LB) != length(x)) LB <- rep(LB, length(x))
    if(length(UB) != length(x)) UB <- rep(UB, length(x))
    
    ## more checks
    if(any(LB > UB)) stop("'LB' can't be > 'UB'")
    if(any(apply(cbind(LB, UB), 1, function(x) identical(x[1], x[2])) & !is.finite(LB))) stop("'LB' and 'UB' misspecified")
    if(any(c(lambda1, lambda2) <= 0)) stop("'lambda1' and 'lambda2' must be > 0")
    
    ## generate non-truncated densities
    ldens <- dskellam(x, lambda1, lambda2, log = TRUE)
    
    ## adjust components for truncation as required
    UBind <- which(!is.finite(LB) & is.finite(UB))
    LBind <- which(is.finite(LB) & !is.finite(UB))
    bothind <- which(!is.finite(LB) & !is.finite(UB))
    if(length(UBind) > 0) {
        ldens[UBind] <- ldens[UBind] - pskellam(UB[UBind], lambda1[UBind], lambda2[UBind], lower.tail = TRUE, log.p = TRUE)
        ldens[UBind] <- ifelse(x[UBind] > UB[UBind], -Inf, ldens[UBind])
    }
    if(length(LBind) > 0) {
        ldens[LBind] <- ldens[LBind] - pskellam(LB[LBind], lambda1[LBind], lambda2[LBind], lower.tail = FALSE, log.p = TRUE)
        ldens[LBind] <- ifelse(x[LBind] < LB[LBind], -Inf, ldens[LBind])
    }
    if(length(bothind) > 0) {
        ldens[bothind] <- ldens[bothind] - log(diff(pskellam(c(LB[bothind], UB[bothind]), lambda1[bothind], lambda2[bothind], lower.tail = TRUE, log.p = FALSE)))
        ldens[bothind] <- ifelse(x[bothind] < LB[bothind] | x[bothind] > UB[bothind], -Inf, ldens[bothind])
    }
    if(log) {
        return(ldens)
    } else {
        return(exp(ldens))
    }
}

## truncated Skellam sampling
rtskellam <- function(lambda1, lambda2, LB = -Inf, UB = Inf, ntries = 1000) {
    
    ## checks on inputs
    if(LB > UB) stop("'LB' can't be > 'UB'")
    if(identical(LB, UB) & !is.finite(LB)) stop("'LB' and 'UB' misspecified")
    if(any(c(lambda1, lambda2) <= 0)) stop("'lambda1' and 'lambda2' must be > 0")
    if(length(lambda1) != 1 | length(lambda2) != 1 | length(LB) != 1 | length(UB) != 1) {
        stop("Inputs must be of length 1 currently")
    } 
    
    if(LB == -Inf & UB == Inf) {
        return(rskellam(1, lambda1, lambda2))
    }
    x <- rskellam(1, lambda1, lambda2)
    k <- 1
    while((x < LB | x > UB) & k < ntries) {
        x <- rskellam(1, lambda1, lambda2)
        k <- k + 1
    }
    if(k == ntries) {
        ## if brute force doesn't work, then try inverse transform sampling
        u <- runif(1, 0, 1)
        if(is.finite(LB) & !is.finite(UB)) {
            stop("Need to check this")
            browser()
            ## normalising constant
            norm <- pskellam(LB, lambda1, lambda2, lower.tail = FALSE, log.p = TRUE)
            UB <- 1000
            x <- cumsum(exp(dskellam(LB:UB, lambda1, lambda2, log = TRUE) - norm))
            ind <- which(x >= u)
            while(length(ind) == 0) {
                LB <- UB
                UB <- UB + 1000
                x <- x[length(x)] + cumsum(exp(dskellam(LB:UB, lambda1, lambda2, log = TRUE) - norm))
                ind <- which(x >= u)
            }
            x <- c(LB:UB)[ind[1]]
        } else {
            if(!is.finite(LB) & is.finite(UB)) {
                stop("Need to check this")
                browser()
                ## normalising constant
                norm <- pskellam(UB, lambda1, lambda2, lower.tail = TRUE, log.p = TRUE)
                LB <- -1000
                x <- cumsum(rev(exp(dskellam(LB:UB, lambda1, lambda2, log=  TRUE) - norm)))
                ind <- which(x >= u)
                while(length(ind) == 0) {
                    UB <- LB
                    LB <- LB - 1000
                    x <- x[length(x)] + cumsum(rev(exp(dskellam(LB:UB, lambda1, lambda2, log.p = TRUE) - norm)))
                    ind <- which(x >= u)
                }
                x <- rev(LB:UB)[ind[1]]
            } else {
                ## normalising constant
                norm <- dskellam(LB:UB, lambda1, lambda2, log = TRUE)
                maxn <- max(norm)
                norm <- maxn + log(sum(exp(norm - maxn)))
                x <- cumsum(exp(dskellam(LB:UB, lambda1, lambda2, log = TRUE) - norm))
                ind <- which(x >= u)
                x <- c(LB:UB)[ind[1]]
            }
        }
#        print(paste("lambda1 =", lambda1, "lambda2 =", lambda2, "LB =", LB, "UB =", UB))
#        stop(paste0("Can't sample from truncated Skellam after ", ntries, " tries"))
    }
    x
}
