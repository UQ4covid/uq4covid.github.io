## load libraries
library(skellam)

## truncated Skellam sampling
trSkellam <- function(p1, p2, LB = -Inf, UB = Inf, ntries = 1000) {
    x <- rskellam(1, p1, p2)
    k <- 1
    while((x < LB | x > UB) & k < ntries) {
        x <- rskellam(1, p1, p2)
        k <- k + 1
    }
    if(k == ntries) {
        ## if brute force doesn't work, then try inverse transform sampling
        u <- runif(1, 0, 1)
        if(is.finite(LB) & !is.finite(UB)) {
            stop("Need to check this")
            browser()
            ## normalising constant
            norm <- pskellam(LB, p1, p2, lower.tail = FALSE, log.p = TRUE)
            UB <- 1000
            x <- cumsum(exp(dskellam(LB:UB, p1, p2, log = TRUE) - norm))
            ind <- which(x >= u)
            while(length(ind) == 0) {
                LB <- UB
                UB <- UB + 1000
                x <- x[length(x)] + cumsum(exp(dskellam(LB:UB, p1, p2, log = TRUE) - norm))
                ind <- which(x >= u)
            }
            x <- c(LB:UB)[ind[1]]
        } else {
            if(!is.finite(LB) & is.finite(UB)) {
                stop("Need to check this")
                browser()
                ## normalising constant
                norm <- pskellam(UB, p1, p2, lower.tail = TRUE, log.p = TRUE)
                LB <- -1000
                x <- cumsum(rev(exp(dskellam(LB:UB, p1, p2, log=  TRUE) - norm)))
                ind <- which(x >= u)
                while(length(ind) == 0) {
                    UB <- LB
                    LB <- LB - 1000
                    x <- x[length(x)] + cumsum(rev(exp(dskellam(LB:UB, p1, p2, log.p = TRUE) - norm)))
                    ind <- which(x >= u)
                }
                x <- rev(LB:UB)[ind[1]]
            } else {
                ## normalising constant
                norm <- dskellam(LB:UB, p1, p2, log = TRUE)
                maxn <- max(norm)
                norm <- maxn + log(sum(exp(norm - maxn)))
                x <- cumsum(exp(dskellam(LB:UB, p1, p2, log = TRUE) - norm))
                ind <- which(x >= u)
                x <- c(LB:UB)[ind[1]]
            }
        }
#        print(paste("p1 =", p1, "p2 =", p2, "LB =", LB, "UB =", UB))
#        stop(paste0("Can't sample from truncated Skellam after ", ntries, " tries"))
    }
    x
}
