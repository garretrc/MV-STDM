function (modFilt) 
{
    eps <- .Machine$double.eps^0.4
    mod <- c(modFilt[match(c("m", "U.C", "D.C", "a", "U.R", "D.R"), 
        names(modFilt))], modFilt$mod[match(c("GG", "W", "JGG", 
        "JW", "X"), names(modFilt$mod))])
    n <- length(mod$U.R)
    p <- NCOL(mod$m)
    mtsp <- tsp(mod$m)
    if (p == 1) {
        dim(mod$m) <- c(n + 1, 1)
        dim(mod$a) <- c(n, 1)
    }
    if (is.null(mod$JGG)) 
        tvGG <- FALSE
    else {
        tvGG <- TRUE
        nz <- mod$JGG != 0
        mod$JGG <- cbind(row(mod$JGG)[nz], col(mod$JGG)[nz], 
            mod$JGG[nz])
    }
    if (is.null(mod$JW)) 
        tvW <- FALSE
    else {
        tvW <- TRUE
        nz <- mod$JW != 0
        mod$JW <- cbind(row(mod$JW)[nz], col(mod$JW)[nz], mod$JW[nz])
    }
    tvGW <- tvGG || tvW
    theta <- matrix(0, n + 1, p)
    if (!tvW) {
        tmp <- La.svd(mod$W, nu = 0)
        Dw <- sqrt(tmp$d)
        Dw <- pmax(Dw, eps)
        Dw.inv <- 1/Dw
        sqrtWinv <- Dw.inv * tmp$vt
        if (!tvGG) 
            tG.Winv <- t(mod$GG) %*% crossprod(sqrtWinv)
    }
    theta[n + 1, ] <- mod$m[n + 1, ] + mod$U.C[[n + 1]] %*% matrix(mod$D.C[n + 
        1, ] * rnorm(p))
    for (i in (n:1)) {
        if (tvGG) 
            mod$GG[mod$JGG[, -3, drop = FALSE]] <- mod$X[i, mod$JGG[, 
                3]]
        if (tvW) {
            mod$W[mod$JW[, -3, drop = FALSE]] <- mod$X[i, mod$JW[, 
                3]]
            tmp <- La.svd(mod$W, nu = 0)
            Dw <- sqrt(tmp$d)
            Dw <- pmax(Dw, eps)
            Dw.inv <- 1/Dw
            sqrtWinv <- Dw.inv * tmp$vt
        }
        if (tvGW) 
            tG.Winv <- t(mod$GG) %*% crossprod(sqrtWinv)
        D.inv <- 1/mod$D.C[i, ]
        D.inv[abs(D.inv) == Inf] <- 0
        tmp <- La.svd(rbind(sqrtWinv %*% mod$GG %*% mod$U.C[[i]], 
            diag(x = D.inv, nrow = length(D.inv))), nu = 0)
        U.H <- mod$U.C[[i]] %*% t(tmp$vt)
        D.H <- 1/tmp$d
        D.H[abs(D.H) == Inf] <- 0
        h <- mod$m[i, ] + crossprod(D.H * t(U.H)) %*% tG.Winv %*% 
            (t(theta[i + 1, , drop = F]) - mod$a[i, ])
        theta[i, ] <- h + U.H %*% matrix(D.H * rnorm(p))
    }
    if (!is.null(mtsp)) {
        theta <- drop(theta)
        tsp(theta) <- mtsp
        class(theta) <- if (p > 1) 
            c("mts", "ts")
        else "ts"
    }
    return(theta = theta)
}