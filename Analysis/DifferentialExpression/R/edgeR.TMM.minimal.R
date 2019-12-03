## extracted code from the edgeR software to explore TMM normalization



## test stability of TMM factors:

## - for counts matrix:
## counts.matrix = read.table("Trinity_trans.counts.matrix", header=T, row.names=1)
## sort(calcNormFactors(counts.matrix[,sample(ncol(counts.matrix))]))

## - for TPM matrix

## tpm.matrix = read.table("Trinity_trans.TPM.not_cross_norm", header=T, row.names=1)
## sort(calcNormFactors(tpm.matrix[,sample(ncol(tpm.matrix))]))


calcNormFactors <-
function (object, lib.size = NULL, method = c("TMM"), refColumn = NULL, logratioTrim = 0.3,
          sumTrim = 0.05, doWeighting = TRUE, Acutoff = -1e+10, p = 0.75,
    ...)
{
    x <- as.matrix(object)
    if (any(is.na(x)))
        stop("NA counts not permitted")
    nsamples <- ncol(x)
    if (is.null(lib.size)) {
        lib.size <- colSums(x)
    }
    else {
        if (anyNA(lib.size))
            stop("NA lib.sizes not permitted")
        if (length(lib.size) != nsamples) {
            if (length(lib.size) > 1L)
                warning("calcNormFactors: length(lib.size) doesn't match number of samples",
                  call. = FALSE)
            lib.size <- rep(lib.size, length = nsamples)
        }
    }

    method <- match.arg(method)
    allzero <- .rowSums(x > 0, nrow(x), nsamples) == 0L
    if (any(allzero))
        x <- x[!allzero, , drop = FALSE]

    f75 <- calcFactorQuantile(data = x, lib.size = lib.size,
                              p = 0.75)
    if (is.null(refColumn)) refColumn <- which.min(abs(f75 -
                                                       mean(f75)))
    if (length(refColumn) == 0L | refColumn < 1 | refColumn >
        nsamples) refColumn <- 1L
    f <- rep(NA, nsamples)
    for (i in 1:nsamples) f[i] <- calcFactorTMM(obs = x[,
                                                        i], ref = x[, refColumn], libsize.obs = lib.size[i],
                                                libsize.ref = lib.size[refColumn], logratioTrim = logratioTrim,
                                                sumTrim = sumTrim, doWeighting = doWeighting, Acutoff = Acutoff)

    f <- f/exp(mean(log(f)))
    names(f) <- colnames(x)
    f
}


calcFactorTMM <-
function (obs, ref, libsize.obs = NULL, libsize.ref = NULL, logratioTrim = 0.3,
    sumTrim = 0.05, doWeighting = TRUE, Acutoff = -1e+10)
{
    obs <- as.numeric(obs)
    ref <- as.numeric(ref)
    if (is.null(libsize.obs))
        nO <- sum(obs)
    else nO <- libsize.obs
    if (is.null(libsize.ref))
        nR <- sum(ref)
    else nR <- libsize.ref
    logR <- log2((obs/nO)/(ref/nR))
    absE <- (log2(obs/nO) + log2(ref/nR))/2
    v <- (nO - obs)/nO/obs + (nR - ref)/nR/ref
    fin <- is.finite(logR) & is.finite(absE) & (absE > Acutoff)
    logR <- logR[fin]
    absE <- absE[fin]
    v <- v[fin]
    if (max(abs(logR)) < 1e-06)
        return(1)
    n <- length(logR)
    loL <- floor(n * logratioTrim) + 1
    hiL <- n + 1 - loL
    loS <- floor(n * sumTrim) + 1
    hiS <- n + 1 - loS
    keep <- (rank(logR) >= loL & rank(logR) <= hiL) & (rank(absE) >=
        loS & rank(absE) <= hiS)
    if (doWeighting)
        f <- sum(logR[keep]/v[keep], na.rm = TRUE)/sum(1/v[keep],
            na.rm = TRUE)
    else f <- mean(logR[keep], na.rm = TRUE)
    if (is.na(f))
        f <- 0
    2^f
}

calcFactorQuantile <-
function (data, lib.size, p = 0.75)
{
    y <- t(t(data)/lib.size)
    f <- apply(y, 2, function(x) quantile(x, p = p))
}


