## Downloaded from: http://bioinfo.cipf.es/noiseq/doku.php?id=downloads on 11-11-2011

## NOISeq paper available here: http://www.ncbi.nlm.nih.gov/pubmed/21903743
##     Differential expression in RNA-seq: A matter of depth. 
##     Tarazona S, García-Alcalde F, Dopazo J, Ferrer A, Conesa A.
##     Genome Res. 2011 Oct 28

####################################################################
###################           NOISeq            ####################
####################################################################

# By Sonia Tarazona
# Last modified: 20-IV-2011




## MAIN FUNCTION


noiseq <- function (datos1, datos2, k = 0.5, norm = "rpkm",
                    long = 1000, q = 0.90, repl = "tech",                  
                    pnr = 0.2, nss = 5, v = 0.02, lc = 1)

# datos1: Matrix containing gene counts and as many columns as samples in
#         group 1.

# datos2: Matrix containing gene counts and as many columns as samples in
#         group 2. Row names must be the same than in datos1.

# k:      When counts = 0, 0 will be changed to k. By default, k = 0.5.
  
# norm:   Normalization method. It can be one of "rpkm" (default), "uqua"
#         (upper quartile), "tmm" (trimmed mean of M) or "n" (no normalization).

# long:   Vector containing genes length, whose names are gene names in datos1
#         and datos2. If long = 1000, no correction by length is applied.

# q:      Threshold to determine differentially expressed genes.
#         By default, q = 0.95.

# pnr:    Percentage of total reads (seq.depth) for each simulated sample.
#         Only needed when noise = "simul". By default, pnr = 1.

# nss:    Number of simulated samples (>= 2). By default, nss = 5.
#         If nss = 0, real samples are used to compute noise.

# v:      Variability in sample total reads used to simulate samples.
#         By default, v = 0.02. Sample total reads is computed as a
#         random value from a uniform distribution in the interval
#         [(pnr-v)*sum(counts), (pnr+v)*sum(counts)]

# lc:     Length correction in done by dividing expression by length^lc.
#         By default, lc = 1. 

# repl:   "tech" if you have technical replicates.
#         "bio" if the replicates are biological (there must be replicates!).

{

  n1 <- ncol(as.matrix(datos1))
  n2 <- ncol(as.matrix(datos2))
  g.sinL <- names(which(is.na(long)))


  if (norm == "n") {      # no normalization
    datos1 <- round(datos1, 100)
    datos2 <- round(datos2, 100)
  }


  if (is.null(k)) {
      m1 <- min(datos1[noceros(datos1, num = FALSE)], na.rm = TRUE)
      m2 <- min(datos2[noceros(datos2, num = FALSE)], na.rm = TRUE)
      mm <- min(m1, m2)
      k <- mm/2    
  } 

  
    # Total counts for each gene:
  suma1 <- rowSums(datos1)
  suma2 <- rowSums(datos2)
  

    # All genes
  todos <- rownames(as.matrix(datos1))

    # Genes with counts in any condition
  concounts <- names(which(suma1+suma2 > 0))

  if (length(long) > 1) {
    long <- long[concounts]
  }


  if (repl == "tech") {  ### technical replicates
    suma1 <- suma1[concounts]
    suma2 <- suma2[concounts]    

    #-------------------------------------------------------------------------#
    # Normalization of counts for each condition (aggregating replicates)

    if (norm == "rpkm") {      # RPKM
      suma1.norm <- rpkm(suma1, long = long, k = k, lc = lc)
      suma2.norm <- rpkm(suma2, long = long, k = k, lc = lc)
    }

    
    if (norm == "uqua") {
      suma.norm <- uqua(cbind(suma1, suma2), long = long, lc = lc, k = k)
      suma1.norm <- suma.norm[ ,1]
      suma2.norm <- suma.norm[ ,2]
    }

    
    if (norm == "tmm") {
      suma.norm <- tmm(as.matrix(cbind(suma1, suma2)), long = long,
                       lc = lc, k = k)      
      suma1.norm <- suma.norm[ ,1]
      suma2.norm <- suma.norm[ ,2]
    }
    
  }


   

    #-------------------------------------------------------------------------#

    ## Noise distribution

  if ((n1+n2)>2) {   # with real samples

    datitos <- cbind(datos1, datos2)
    datitos <- datitos[concounts,]

    gens.sin0 <- setdiff(concounts, g.sinL)

    if (norm == "n") {       # no normalization
      datitos.0 <- sinceros(datitos, k = k)
      datitos.norm <- datitos.0[gens.sin0, ]
    }

    if (norm == "rpkm") {      # RPKM
      datitos.0 <- rpkm(datitos, long = long, k = k, lc = lc)
      datitos.norm <- datitos.0[gens.sin0, ]
    }

    if (norm == "uqua") {      # Upper Quartile
      datitos.0 <- uqua(datitos, long = long, lc = lc, k = k)
      datitos.norm <- datitos.0[gens.sin0, ]
    }

    if (norm == "tmm") {
      datitos.0 <- tmm(datitos, long = long, lc = lc, k = k)
      datitos.norm <- datitos.0[gens.sin0, ]
    }

    datos1.norm <- datitos.norm[ ,1:n1]
    datos2.norm <- datitos.norm[ ,(n1+1):(n1+n2)]

    if (n1 > 1) {
      MD1 <- MD(dat = datos1.norm)
    } else { MD1 <- NULL }

    if (n2 > 1) {
      MD2 <- MD(dat = datos2.norm)
    } else { MD2 <- NULL }

    
  } else {                 # with simulated samples

    if (nss == 0) {
      nss <- 5
    }

    datos.sim <- sim.samples(counts1 = sinceros(suma1, k = k),
                             counts2 = sinceros(suma2, k = k),
                             pnr = pnr, nss = nss, v = v)

    nn <- sapply(datos.sim, ncol)

    dat.sim.norm <- vector("list", length = 2)

    datitos <- cbind(datos.sim[[1]], datos.sim[[2]])

    sumita <- rowSums(datitos)
    g.sin0 <- names(which(sumita > 0))
    gens.sin0 <- setdiff(g.sin0, g.sinL)

    if (norm == "n") {       # no normalization
      datitos.0 <- sinceros(datitos, k = k)
      datitos.norm <- datitos.0[gens.sin0, ]
    }

    if (norm == "rpkm") {      # RPKM
      datitos.0 <- rpkm(datitos, long = long, k = k, lc = lc)
      datitos.norm <- datitos.0[gens.sin0, ]
    }

    if (norm == "uqua") {      # Upper Quartile
      datitos.0 <- uqua(datitos, long = long, lc = lc, k = k)
      datitos.norm <- datitos.0[gens.sin0, ]
    }

    if (norm == "tmm") {
      datitos.0 <- tmm(datitos, long = long, lc = lc, k = k)
      datitos.norm <- datitos.0[gens.sin0, ]
    }

    dat.sim.norm[[1]] <- datitos.norm[ ,1:nn[1]]
    dat.sim.norm[[2]] <- datitos.norm[ ,(nn[1]+1):sum(nn)]

    MD1 <- MD(dat = dat.sim.norm[[1]])
    MD2 <- MD(dat = dat.sim.norm[[2]])

  }

  Mr <- c(as.numeric(MD1$M), as.numeric(MD2$M))
  Dr <- c(as.numeric(MD1$D), as.numeric(MD2$D))

  
  
    #-------------------------------------------------------------------------#

    ## M and D for different experimental conditions

  if (repl == "tech" & norm != "n") {
    
    MDs <- MD(dat = cbind(suma1.norm, suma2.norm))

  } else {

    if (norm == "n" & (n1+n1) == 2) {
      datos1.norm <- sinceros(as.matrix(datos1)[concounts,], k = k)
      datos2.norm <- sinceros(as.matrix(datos2)[concounts,], k = k)
    }

    resum1.norm <- rowMeans(as.matrix(datos1.norm))
    resum2.norm <- rowMeans(as.matrix(datos2.norm))

    MDs <- MD(dat = cbind(resum1.norm, resum2.norm))

  }
    


 
    #-------------------------------------------------------------------------#
    
    ## Probability of differential expression
    
    prob.concounts <- probdeg(MDs$M, MDs$D, Mr, Dr)

    prob <- prob.concounts[todos]
    names(prob) <- todos


    ## Completing M and D
    Ms0 <- as.numeric(MDs$M)
    names(Ms0) <- rownames(MDs$M)
    Ms <- Ms0[todos]
    names(Ms) <- todos

    Ds0 <- as.numeric(MDs$D)
    names(Ds0) <- rownames(MDs$D)
    Ds <- Ds0[todos]
    names(Ds) <- todos

        
    
    ## Differentially expressed genes
    deg <- DEG.q(prob, q = q)


    ## Results
    list("probab" = prob, "deg" = deg, "Ms" = Ms, "Ds" = Ds,
       "Mn" = Mr, "Dn" = Dr)

  }





#*****************************************************************************#

## Reading data

readData <- function (file, header = FALSE, cond1, cond2) {

  mydata <- vector("list", length = 2)

  myfile <- read.delim(file = file, header = header, as.is = TRUE)

  mynames <- as.character(myfile[,1])

  mydata[[1]] <- as.matrix(myfile[,cond1])
  mydata[[2]] <- as.matrix(myfile[,cond2])

  rownames(mydata[[1]]) <- rownames(mydata[[2]]) <- mynames

  mydata

}


readInfo <- function (file, header = FALSE) {

  myfile <- read.delim(file = file, header = header, as.is = TRUE)

  mynames <- as.character(myfile[,1])

  myinfo <- myfile[,2]

  names(myinfo) <- mynames

  myinfo
  
}








#*****************************************************************************#


## Replacing counts=0 with counts=k


sinceros <- function (datos, k) {
  datos0 <- datos

  if (is.null(k)) {

    mini0 <- min(datos[noceros(datos, num = FALSE, k = 0)])

    kc <- mini0/2

    datos0[datos0 == 0] <- kc

  } else {

    datos0[datos0 == 0] <- k

  }
  
  datos0
}



#****************************************************************************#



## RPKM normalization

rpkm <- function (datos, long = 1000, k = 0, lc = 1) {
  
  total <- colSums(as.matrix(datos))

  datos0 <- sinceros(datos, k)                   

  datos.norm <- (t(t(datos0)/total)*10^9)/(long^lc)
  
  na.omit(datos.norm)   
}



#****************************************************************************#


## UQUA: Upper-quartile normalization (Bullard et al. et Sandrine Dudoit, 2010)

uqua <- function (datos, long = 1000, lc = 1, k = 0) {
  
# lc: Length correction. Expression is divided by long^lc. lc can be any real number.

  L <- long^lc

  datos0 <- sinceros(datos, k)

  if (ncol(as.matrix(datos)) > 1) {
    sumatot <- rowSums(datos)
    supertot <- sum(sumatot)
    counts0 <- which(sumatot == 0)

    if (length(counts0) > 0) {
      datitos <- datos[-counts0,]
    } else {
      datitos <- datos
    }

    q3 <- apply(datitos, 2, quantile, probs = 0.75)
    d <- q3*supertot/sum(q3)

    datos.norm <- (t(t(datos0)/d)*10^9)/L

  } else {

    datos.norm <- datos0/L

  }
  
  na.omit(datos.norm)  
  
}



#****************************************************************************#



## TMM: Trimmed Mean of M values normalization (Robinson & Oshlack, 2010)

tmm <- function (datos, long = 1000, lc = 1, k = 0, refColumn = NULL,
                 logratioTrim = .3, sumTrim = 0.05, doWeighting = TRUE,
                 Acutoff = -1e10) {

  # lc: Length correction. Expression is divided by long^lc. lc can be any real number.

  library(edgeR)

  L <- long^lc

  datos0 <- sinceros(datos, k)

  if (ncol(as.matrix(datos)) > 1) {

    fk <- calcNormFactors(as.matrix(datos), method = "TMM", refColumn = refColumn,
                          logratioTrim = logratioTrim, sumTrim = sumTrim,
                          doWeighting = doWeighting, Acutoff = Acutoff)
    
    datos.norm <- (t(t(datos0)*fk)*10^3)/L

  } else {

    datos.norm <- datos0/L

  }
  
  na.omit(datos.norm)  
    
}




#*****************************************************************************#





## Computing M and D

MD <- function (dat = dat, selec = c(1:nrow(dat))) {
  pares <- as.matrix(combn(ncol(dat), 2))

  if (NCOL(pares) > 30) {  
    sub30 <- sample(1:NCOL(pares), size = 30, replace = FALSE)
    pares <- pares[,sub30]
  }
  
  mm <- NULL
  dd <- NULL
  for (i in 1:ncol(pares)) {
    a <- dat[selec,pares[1,i]]
    b <- dat[selec,pares[2,i]]
    mm <- cbind(mm, log(a/b, 2))
    dd <- cbind(dd, abs(a-b))
  }
  list("M" = mm, "D" = dd)
}




#****************************************************************************#




## Computing M and A

MA <- function (dat = dat, selec = c(1:nrow(dat))) {
  pares <- as.matrix(combn(ncol(dat), 2))

  if (NCOL(pares) > 30) {  
    sub20 <- sample(1:NCOL(pares), size = 30, replace = FALSE)
    pares <- pares[,sub20]
  }
  
  mm <- NULL
  aa <- NULL
  for (i in 1:ncol(pares)) {
    a <- dat[selec,pares[1,i]]
    b <- dat[selec,pares[2,i]]
    mm <- cbind(mm, log(a/b, 2))
    aa <- cbind(aa, (log2(a)+log2(b))/2)
  }
  list("M" = mm, "A" = aa)
}







#***************************************************************************#



## Probability for a gene of being differentially expressed

# alternative = c("two.sided", "less", "greater")
# compare =  c("diff", "up", "down")


n.menor <- function (x, S1, S2) {
 
  length(which(S1 <= x[1] &  S2 <= x[2]))

}



#****************************************************************************#


busca <- function (x, S) {
  which(S[,1] == x[1] & S[,2] == x[2])
}



#****************************************************************************#


## Probability for a gene to be differentially expressed 

probdeg <- function (Mg, Dg, Mn, Dn, prec = 2) {

  # Mg, Dg -> signal
  # Mn, Dn -> noise
  # prec = precission (number of digits to round M and D)

  tot <- length(Mn)   # number of points in noise distribution

  Mruido <- abs(round(Mn, prec))
  Druido <- round(Dn, prec)
  Mgen <- abs(round(Mg, prec))
  Dgen <- round(Dg, prec)
  
  MDgen <- unique(cbind(Mgen, Dgen))

  Nres <- apply(MDgen, 1, n.menor, S1 = Mruido, S2 = Druido)

  lugares <- apply(cbind(Mgen,Dgen), 1, busca, S = MDgen)

  Nconj <- Nres[lugares]

  names(Nconj) <- names(lugares)

  Nconj / tot

}




#***************************************************************************#





## To extract DEG

DEG.q <- function (x, q = 1, flag = NULL) {
  nnn <- names(which(x >= q))
  setdiff(nnn, flag)
}


DEG.list <- function(lista, q) {
  lapply(lista, DEG.q, q = q)
}



#****************************************************************************#



## Function to intersect multiple sets

int.mult <- function(lista, todos = NULL) {
  
  if(is.null(todos)) {
    todos <- unlist(lista)
  }

  comunes <- todos

  for(i in 1:length(lista)) {
    comunes <- intersect(comunes, lista[[i]])
  }

  comunes
}



#*****************************************************************************#



## Function of union of multiple sets

uni.mult <- function(lista) {
  
  todos <- lista[[1]]
  
  for(i in 2:length(lista)) {
    todos <- union(todos, lista[[i]])
  }

  todos
}



#*****************************************************************************#




## To calculate number of genes in common for two sets

long.int <- function(x, y) { length(intersect(x, y)) }



#****************************************************************************#


## To transform counts from two samples into a contingency table

tablacont <- function (x, total) {  # x; total: 2 x 1
  m1 <- round(c(x[1], total[1]-x[1]), 0)
  m2 <- round(c(x[2], total[2]-x[2]), 0)
  tt <- rbind(m1, m2)
  colnames(tt) <- c("yes", "no")
  tt
}




#****************************************************************************#

## To generate simulated samples:

# pnr: Percentage of total reads (seq.depth). By default, pnr = 1.
# nss: Number of simulated samples (>= 2). By default, nss = 5.
# v:   Variability in sample total reads.

# Sample total reads is computed as a random value from a uniform distribution
# in the interval [(pnr-v)*sum(counts), (pnr+v)*sum(counts)]


sim.samples <- function(counts1, counts2 = NULL, pnr = 1, nss = 5, v = 0.02) {
  seqdep <- c(sum(counts1), sum(counts2))
  num.reads1 <- (pnr + c(-v,v))*seqdep[1]

  muestras <- vector("list")
  
  muestras$c1 <- NULL
  for (s in 1:nss) {
    tama <- round(runif(1, num.reads1[1], num.reads1[2]), 0)
    muestras$c1 <- cbind(muestras$c1,
                         rmultinom(1, size = tama, prob = counts1))
  }

  if(!is.null(counts2)) {
    num.reads2 <- (pnr + c(-v,v))*seqdep[2]
    muestras$c2 <- NULL
    for (s in 1:nss) {
      tama <- round(runif(1, num.reads2[1], num.reads2[2]), 0)
      muestras$c2 <- cbind(muestras$c2,
                           rmultinom(1, size = tama, prob = counts2))
    }
  }

  muestras
}



#*****************************************************************************#



## Fisher's exact test


fisher.ngs <- function(datos1, datos2, norm = "rpkm", long = NULL,
                       adjP = "fdr", sig.level = 0.05) {

  if (is.null(long)) {
    long <- 1000
  } 

  suma1 <- apply(as.matrix(datos1), 1, sum)
  suma2 <- apply(as.matrix(datos2), 1, sum)
  #names(suma1) <- names(suma2) <- rownames(as.matrix(datos1))
  tot <- c(sum(suma1), sum(suma2))

  if (norm  == "n") {       # no normalization
    dat <- cbind(suma1, suma2)
  }

  if (norm == "rpkm") {      # RPKM
    norm1 <- rpkm(suma1, long = long, k = 0)
    norm2 <- rpkm(suma2, long = long, k = 0)
    dat <- cbind(norm1, norm2)
  }

  if (norm == "uqua") {      # Upper Quartile
    dat <- uqua(cbind(suma1,suma2), long = long, k = 0)
  }

  totfi <- apply(dat, 2, sum)

  fish.pv <- NULL

  for (g in 1:nrow(dat)) {
    ttt <- tablacont(dat[g,], totfi)
    fish.pv <- c(fish.pv, fisher.test(ttt)$p.value)
  }

  names(fish.pv) <- rownames(dat)

  fish.adjP <- p.adjust(fish.pv, method = adjP)

  deg <- DEG.q(1-fish.adjP, q = 1 - sig.level)

  list("p.value" = fish.pv, "adj.p.val" = fish.adjP, "deg" = deg)
}



#***************************************************************************#



## Plot MD noise vs deg

MD.plot <- function (Ms = Ms, Ds = Ds, Mn = Mn, Dn = Dn,
                     xlim = range(Ms), ylim = range(Ds),
                     tit = "") {

  plot(Mn, Dn, pch = ".", main = tit, xlab = "M", ylab = "D",
       xlim = xlim, ylim = ylim)
  points(Ms, Ds, col = 2, pch = 20)

  legend("topright", c("noise", "deg"), col = 1:2, pch = 15,
         bg = "lightgrey")

}



#************************************************************************************#



## Plot to see dependence on length

Lplot <- function (deg, long, ylim = c(0,1), confi = 0)  {
  
# confi:  Confidence level for wilcoxon test. If confi = 0, no test is done.
  
  gens <- names(long)

  if(confi > 0) {
    wilcox.res <- wilcox.test(long[deg], long[setdiff(gens, deg)],
                              alternative = "two.sided", mu = 0,
                              paired = FALSE, exact = NULL, correct = TRUE,
                              conf.int = TRUE, conf.level = confi)
  }


  tabla <- table(long)
  acum <- cumsum(tabla)
  corte <- 0
  N <- 300

  while(N <= length(na.omit(long))) {
    fN <- acum[acum >= N]
    corte <- c(corte, as.numeric(names(fN[1])))
    N <- fN[1] + 300
  }

  corte[length(corte)] <- max(na.omit(long))

  bins.long <- cut(long, breaks = corte)  # divides length into 300 genes bins
  names(bins.long) <- gens

  medianas <- aggregate(long, by = list(bins.long), median)
  colnames(medianas) <- c("bin", "mediana")

  bins <- medianas$bin
  medianas <- medianas[,-1]
  names(medianas) <- bins

  long.deg <- matrix(NA, length(medianas), 2)
  long.deg[,1] <- medianas
  rownames(long.deg) <- bins
  colnames(long.deg) <- c("median.leng", "%deg")

  gen.bin <- na.omit(bins.long)

  for (i in bins) {
    gege <- names(gen.bin)[which(gen.bin == i)]
    num <- length(gege)
    dede <- length(intersect(gege, deg))
    long.deg[i,2] <- dede/num
  }

  plot(long.deg, ylim = ylim, xlab = "median gene length", main = "")
  
  if (confi > 0) {
    legend("top", legend = c(wilcox.res$method,
                     paste("p-value:", format(wilcox.res$p.value,digists = 4)),
                     paste(confi*100, "% confidence interval: [",
                           round(wilcox.res$conf.int[1], 2), "; ",
                           round(wilcox.res$conf.int[2], 2), "]", sep = "")),
           bty = "n")

  }
}




#************************************************************************************#


## Plot to see dependence on length for several methods

Lplot2 <- function (deg, long, ylim = c(0,1), confi = 0, ng = 300, tit = "",
                    plotdata = NULL, x = "bottomright", y = NULL, ss) {

  # deg: List containing several sets of differentially expressed genes.
  # long: Vector containing genes length.
  # ng: Number of genes per length bin.
  # confi:  Confidence level for wilcoxon test. If confi = 0, no test is done.
  # tit: Title for the plot.
  # x,y: The x and y co-ordinates to be used to position the legend.
  # ss: vector containing pch for plot symbols
        
  
  if (is.null(plotdata)) {

    gens <- names(long)

    ## Wilcoxon's test

    if(confi > 0) {

      wilcox.res <- NULL
      
      for(j in 1:length(deg)) {

        resul <- wilcox.test(long[deg[[j]]],
                             long[setdiff(gens, deg[[j]])],
                             alternative = "two.sided", mu = 0,
                             paired = FALSE, exact = NULL,
                             correct = TRUE, conf.int = TRUE,
                             conf.level = confi)

        wilcox.res <- rbind(wilcox.res, unlist(resul))
      }

      rownames(wilcox.res) <- names(deg)
      colnames(wilcox.res) <- names(unlist(resul))
      wilcox.res <- as.data.frame(wilcox.res)
      
    } else { wilcox.res = NULL }


    ## Length bins

    tabla <- table(long)
    acum <- cumsum(tabla)
    corte <- 0
    N <- ng

    while(N <= length(na.omit(long))) {
      fN <- acum[acum >= N]
      corte <- c(corte, as.numeric(names(fN[1])))
      N <- fN[1] + ng
    }

    corte[length(corte)] <- max(na.omit(long))

    bins.long <- cut(long, breaks = corte)  # divides length into ng genes bins
    names(bins.long) <- gens

    medianas <- aggregate(long, by = list(bins.long), median)
    colnames(medianas) <- c("bin", "mediana")

    bins <- medianas$bin
    medianas <- medianas[,-1]
    names(medianas) <- bins

    

    ## %deg in each length bin

    long.deg <- matrix(NA, length(medianas), 1+length(deg))
    long.deg[,1] <- medianas
    rownames(long.deg) <- bins
    colnames(long.deg) <- c("median.leng", names(deg))

    gen.bin <- na.omit(bins.long)

    for (i in bins) {

      gege <- names(gen.bin)[which(gen.bin == i)]
      num <- length(gege)

      for(j in 1:length(deg)) {
        dede <- length(intersect(gege, deg[[j]]))
        long.deg[i,1+j] <- dede/num
      }
    }

    return(list("wilcoxon" = wilcox.res, "Ldeg" = long.deg))
    
  } else {

    long.deg <- plotdata$Ldeg
    wilcox.res <- plotdata$wilcoxon

  }


  ## Plot

  k <- 1.2

  plot(long.deg[,1:2], ylim = ylim, xlab = "median gene length",
       ylab = "%deg", main = tit,
       col = 1, pch = ss[1], cex.axis = k, cex.lab = k, cex.main = k+0.2)

  for (j in 3:ncol(long.deg)) {

    points(long.deg[,c(1,j)], col = j-1, pch = ss[j-1])

  }

  if (confi > 0) {

    intconf <- cbind(rownames(wilcox.res),
                     round(as.numeric(as.character(wilcox.res[,7])),1),
                     round(as.numeric(as.character(wilcox.res[,8])),1))

    leyenda <- apply(intconf, 1, function(x) {
      paste(x[1], ": [", x[2], "; ", x[3], "]", sep = "") } )

    legend(x = x, y = y, legend = leyenda, bty = "n",
           title = paste("Methods & ", confi*100, "% confidence interval
 for length difference between deg and non-deg", sep = ""),
           col = 1:(ncol(long.deg)-1), pch = ss, cex = k)

  } else {

    legend(x = x, y = y, legend = colnames(long.deg)[-1],
           col = 1:(ncol(long.deg)-1), pch = ss, cex = k)

  }

}





#************************************************************************************#




## To count number of non-zero elements

noceros <- function (x, num = TRUE, k = 0) {
  
  nn <- length(which(x > k))
  
  if (num) {
    nn
    
  } else {
    if(nn > 0) { which(x > k) } else { NULL }
  }
}




#************************************************************************************#





## Saturation Plot


satur.plot <- function (datos1, datos2, ylim = NULL, k = 0, tit = "Saturation",
                        cex.main = cex.main, cex.lab = cex.lab, cex.axis = cex.axis,
                        cex = cex,
                        legend = c(deparse(substitute(datos1)),
                                   deparse(substitute(datos2)))) {

  n1 <- ncol(as.matrix(datos1))
  n2 <- ncol(as.matrix(datos2))

  
  if (n1 > 1) {

    muestra1 <- as.list(1:n1)
    for (i in 2:n1) {
      
      combi1 <- combn(n1, i, simplify = FALSE)

      if( length(combi1) > 20 ) {
        sub20 <- sample(1:length(combi1), size = 20, replace = FALSE)
        combi1 <- combi1[sub20]
      }
      
      muestra1 <- append(muestra1, combi1)
    }

    varias1 <- vector("list", length = length(muestra1))
    names(varias1) <- sapply(muestra1, function(x) {
      paste("C1.", x, collapse = "", sep = "")})

    for (i in 1:length(muestra1)) {
      varias1[[i]] <- apply(as.matrix(datos1[,muestra1[[i]]]), 1, sum)
    }

    satura1 <- data.frame("muestra" = names(varias1),
                          "seq.depth" = sapply(varias1, sum),
                          "noceros" = sapply(varias1, noceros, k = k))
  }

  
  if (n2 > 1) {

    if (n1 == n2) {
      muestra2 <- muestra1
      
    } else {
      
      muestra2 <- as.list(1:n2)

      for (i in 2:n2) {

        combi2 <- combn(n2, i, simplify = FALSE)

        if (length(combi2) > 20) {
          sub20 <- sample(1:length(combi2), size = 20, replace = FALSE)
          combi2 <- combi2[sub20]
        }

        muestra2 <- append(muestra2, combi2)
      }
    }

    varias2 <- vector("list", length = length(muestra2))
    names(varias2) <- sapply(muestra2, function(x) {
      paste("C2.", x, collapse = "", sep = "")})

    for (i in 1:length(muestra2)) {
      varias2[[i]] <- apply(as.matrix(datos2[,muestra2[[i]]]), 1, sum)
    }

    satura2 <- data.frame("muestra" = names(varias2),
                          "seq.depth" = sapply(varias2, sum),
                          "noceros" = sapply(varias2, noceros, k = k))
  }


  
  if (n1 == 1) {
    
    total1 <- sum(datos1)
    satura1 <- NULL
    
    for (i in 1:9) {     # 10%, 20%, ..., 90% reads (apart 100% is calculated)     
      muestra1 <- rmultinom(10, size = round(total1*i/10,0), prob = datos1)
      detec1 <- mean(apply(muestra1, 2, noceros, k = k))
      satura1 <- rbind(satura1, c(round(total1*i/10,0), detec1))
    }
    satura1 <- rbind(satura1, c(total1, noceros(datos1, k = k)))
    colnames(satura1) <- c("seq.depth", "noceros")
    satura1 <- as.data.frame(satura1)
  }


  if (n2 == 1) {
    
    total2 <- sum(datos2)
    satura2 <- NULL
    
    for (i in 1:9) {     # 10%, 20%, ..., 90% reads (apart 100% is calculated)     
      muestra2 <- rmultinom(10, size = round(total2*i/10,0), prob = datos2)
      detec2 <- mean(apply(muestra2, 2, noceros, k = k))
      satura2 <- rbind(satura2, c(round(total2*i/10,0), detec2))
    }
    satura2 <- rbind(satura2, c(total2, noceros(datos2, k = k)))
    colnames(satura2) <- c("seq.depth", "noceros")
    satura2 <- as.data.frame(satura2)
  } 
  

  if (is.null(ylim)) {
    ylim <- c(0, NROW(datos1))
  }

  SS1 <- range(satura1$seq.depth)
  SS2 <- range(satura2$seq.depth)

  xM <- max(SS1[2], SS2[2])
  xm <- min(SS1[1], SS2[1])
  

  plot(satura1$seq.depth, satura1$noceros, pch = 16, col = 2, ylim = ylim,
       xlim = c(xm, xM), main = tit, type = "b", xlab = "Sequencing Depth",
       ylab = paste("Number of genes with reads >", k),
       cex.main = cex.main, cex.lab = cex.lab, cex.axis = cex.axis)

  points(satura2$seq.depth, satura2$noceros, pch = 16, col = 4, type = "b")
  
  legend("top", legend = legend, text.col = c(2,4), bty = "n",
         lty = 1, lwd = 2, col = c(2, 4), cex = cex)
  
}





#************************************************************************************#




## Data for saturation plot with different biotypes


saturbio.dat <- function (datos1, datos2, k = 0, infobio, biotypes, nsim = 5) {

  # infobio = vector containing biotype for each gene in "datos"
  # biotypes = list containing groups of biotypes to be studied

  n1 <- NCOL(datos1)
  n2 <- NCOL(datos2)

  biog <- lapply(biotypes, function(x) { which(is.element(infobio, x)) })

  satura1 <- satura2 <- vector("list", length = length(biotypes))
  names(satura1) <- names(satura2) <- names(biotypes)
  

 if (n1 > 1) {  # when replicates are available in condition 2

    varias1 <- vector("list", length = n1) # counts for each sequencing depth
    names(varias1) <- paste(1:n1, "rep", sep = "")

    for(i in 1:n1) {
      
      muestra1 <- combn(n1, i, simplify = FALSE)

      if (length(muestra1) > 20) {
        sub20 <- sample(1:length(muestra1), size = 20, replace = FALSE)
        muestra1 <- muestra1[sub20]
      }

      sumrepli <- NULL

      for (com in 1:length(muestra1)) {
        sumrepli <- cbind(sumrepli, rowSums(as.matrix(datos1[,muestra1[[com]]])))
      }

      varias1[[i]] <- as.matrix(sumrepli)
      
    }
  }


  if (n1 == 1) {   # simulating replicates for datos1

    varias1 <- vector("list", length = nsim) # counts for each sequencing depth
    names(varias1) <- paste("dp", 1:nsim, sep = "")
    
    total1 <- sum(datos1)    

    for (i in 1:(nsim-1)) {    # total1*i/nbox reads (100% is calculated apart)

      muestra1 <- rmultinom(10, size = round(total1*i/nsim,0), prob = datos1)
      
      varias1[[i]] <- muestra1
    }

    varias1[[nsim]] <- as.matrix(datos1)
  }

  seq.depth1 <- sapply(varias1, function(x) { mean(colSums(x)) })

  # computing saturation for each biotype (datos1)
  for (j in 1:length(satura1))  {

    conbio1 <- lapply(varias1, function(x) { as.matrix(x[biog[[j]],]) })

    satura1[[j]] <- sapply(conbio1, function(x) { mean(apply(x, 2, noceros, k = k)) })
  }

  
  
  if (n2 > 1) {   # when replicates are available in condition 2

    varias2 <- vector("list", length = n2) # counts for each sequencing depth
    names(varias2) <- paste(1:n2, "rep", sep = "")

    for(i in 1:n2) {
      
      muestra2 <- combn(n2, i, simplify = FALSE)
      
      if (length(muestra2) > 20) {
        sub20 <- sample(1:length(muestra2), size = 20, replace = FALSE)
        muestra2 <- muestra2[sub20]
      }

      sumrepli <- NULL

      for (com in 1:length(muestra2)) {
        sumrepli <- cbind(sumrepli, rowSums(as.matrix(datos2[,muestra2[[com]]])))
      }

      varias2[[i]] <- as.matrix(sumrepli)
      
    }
  }

  
  if (n2 == 1) {  # replicates have to be simulated

    varias2 <- vector("list", length = nsim) # counts for each sequencing depth
    names(varias2) <- paste("dp", 1:nsim, sep = "")
    
    total2 <- sum(datos2)    

    for (i in 1:(nsim-1)) {    # total1*i/nbox reads (100% is calculated apart)

      muestra2 <- rmultinom(10, size = round(total2*i/nsim,0), prob = datos2)
      
      varias2[[i]] <- muestra2
    }

    varias2[[nsim]] <- as.matrix(datos2)
  }  

  seq.depth2 <- sapply(varias2, function(x) { mean(colSums(x)) })
  
  # computing saturation for each biotype (datos2)
  for (j in 1:length(satura2))  {

    conbio2 <- lapply(varias2, function(x) { as.matrix(x[biog[[j]]]) })
    
    satura2[[j]] <- sapply(conbio2, function(x) { mean(apply(x, 2, noceros, k = k)) })
  }

  

  ## computing detection increasing per million reads

  newdet1 <- newdet2 <- vector("list", length = length(biotypes))
  names(newdet1) <- names(newdet2) <- names(biotypes)
    
  # condition 1
  for (j in 1:length(newdet1))  {

    puntos1 <- data.frame("x" = seq.depth1, "y" = satura1[[j]])

    pendi <- NULL

    for(i in 2:nrow(puntos1)) {
      pendi <- c(pendi, (puntos1$y[i]-puntos1$y[i-1])/(puntos1$x[i]-puntos1$x[i-1]))
    }

    pendimil1 <- c(NA, pendi)*1000000

    newdet1[[j]] <- pendimil1

  }

  # condition 2
  for (j in 1:length(newdet2))  {

    puntos2 <- data.frame("x" = seq.depth2, "y" = satura2[[j]])

    pendi <- NULL

    for(i in 2:nrow(puntos2)) {
      pendi <- c(pendi, (puntos2$y[i]-puntos2$y[i-1])/(puntos2$x[i]-puntos2$x[i-1]))
    }

    pendimil2 <- c(NA, pendi)*1000000

    newdet2[[j]] <- pendimil2

  }


  ### Results

  satura <- list("cond1" = satura1, "cond2" = satura2,
                 "bionum" = sapply(biog, length),
                 "depth1" = seq.depth1, "depth2" = seq.depth2,
                 "newdet1" = newdet1, "newdet2" = newdet2)
  satura

}





#************************************************************************************#






## Saturation plot with different biotypes

saturbio.plot <- function(depth1, depth2, sat1, sat2, newdet1 = NULL, newdet2 = NULL,
                          xlim = NULL, yleftlim = NULL, yrightlim = NULL,
                          main = "Saturation", lwdL = 2, lwdR = 10,
                          xlab = "Sequencing depth",
                          legend = c("sample1", "sample2"), bionum = NULL,
                          ylabL = "Number of detected features",
                          ylabR = "New detections per million reads",
                          cex.main = 1, cex.lab = 1, cex.axis = 1, cex = 1) {
  
  # yleftlim for plot and plot.y2
  if (is.null(yleftlim)) {
    yleftlim <- c(min(c(sat1, sat2)), max(c(sat1, sat2)))
  }

  # xlim for plot
  if (is.null(xlim)) {
    SS1 <- range(depth1)
    SS2 <- range(depth2)

    xM <- max(SS1[2], SS2[2])
    xm <- min(SS1[1], SS2[1])

    xlim <- c(xm, xM)
  }

  percen1 <- round(100*max(sat1)/bionum, 1)
  percen2 <- round(100*max(sat2)/bionum, 1)


  if(is.null(newdet1)) {

    # PLOT for detections
    plot(depth1, sat1, pch = 16, col = 2, ylim = yleftlim, lwd = lwdL,
         xlim = xlim, main = main, type = "b", xlab = xlab,
         ylab = ylabL, cex.main = cex.main, cex.lab = cex.lab, cex.axis = cex.axis)

    points(depth2, sat2, pch = 16, col = 4, type = "b")

    legend("top",
           legend = paste(legend, ": ", c(percen1, percen2), "% detected", sep = ""),
           text.col = c(2,4), bty = "n", lty = 1, lwd = lwdL, col = c(2, 4), cex = cex)
    
  } else {

    # yrightlim for plot.y2
    if (is.null(yrightlim)) {
      yrightlim <- c(0, max(10,max(na.omit(c(newdet1, newdet2)))))
    }

    # PLOT with 2 axis
    plot.y2(x = depth1, yright = newdet1, yleft = sat1, type = c("h", "b"),
            lwd = c(lwdR, lwdL), xlab = xlab, xlim = xlim,
            yrightlim = yrightlim, yleftlim = yleftlim,
            yylab = c(ylabR, ylabL), pch = c(1,19), col = c("pink",2),
            main = main, x2 = depth2, yright2 = newdet2, yleft2 = sat2,
            col2 = c("lightblue1",4),
            cex.main = cex.main, cex.lab = cex.lab, cex.axis = cex.axis, cex = cex)

    rect(xlim[1]+diff(range(xlim))*0.3, max(yleftlim)*1.5,
         max(xlim)*1.5, yleftlim[1]+diff(range(yleftlim))*0.85,
         col = "grey90", border = "grey90") 

    text(x = xlim[1]+diff(range(xlim))*0.55, y = max(yleftlim), "Left axis", font = 3)
    text(x = xlim[1]+diff(range(xlim))*0.75, y = max(yleftlim), "Right axis", font = 3)
    text(x = xlim[1]+diff(range(xlim))*0.95, y = max(yleftlim), "%detected", font = 3)

    text(x = xlim[1]+diff(range(xlim))*0.4, y = yleftlim[1]+diff(range(yleftlim))*0.95,
         legend[1], font = 2)
    points(x = xlim[1]+diff(range(xlim))*0.55,
           y = yleftlim[1]+diff(range(yleftlim))*0.95, lty = 1, pch = 16, col = 2)
    points(x = xlim[1]+diff(range(xlim))*0.75,
           y = yleftlim[1]+diff(range(yleftlim))*0.95, pch = 15, col = "pink")
    text(x = xlim[1]+diff(range(xlim))*0.95,
         y = yleftlim[1]+diff(range(yleftlim))*0.95, percen1)

    text(x = xlim[1]+diff(range(xlim))*0.4,
         y = yleftlim[1]+diff(range(yleftlim))*0.9, legend[2], font = 2)
    points(x = xlim[1]+diff(range(xlim))*0.55,
           y = yleftlim[1]+diff(range(yleftlim))*0.9, lty = 1, pch = 16, col = 4)
    points(x = xlim[1]+diff(range(xlim))*0.75,
           y = yleftlim[1]+diff(range(yleftlim))*0.9, pch = 15, col = "lightblue1")
    text(x = xlim[1]+diff(range(xlim))*0.95,
         y = yleftlim[1]+diff(range(yleftlim))*0.9, percen2)


    
    #legend(x = sum(range(xlim))*0.4, y = max(yleftlim), legend = legend,
    #       lty = 1, pch = 16, col = c(2,4), title = "Left axis", lwd = 2)
    
    #legend(x = sum(range(xlim))*0.7, y = max(yleftlim), legend = legend,
    #       pch = 15, col = c("pink","lightblue1"), title = "Right axis")

  }
}






#************************************************************************************#



## Plot to compare counts distributions for two experimental conditions


cd.plot <- function (cond1, cond2, legend = c("group1", "group2"),
                     tit = "Counts distributions", xlim = c(0,100), ylim = c(0,100)) {

  if ( max(ncol(as.matrix(cond1)), ncol(as.matrix(cond2))) == 1) {

    suma <- cond1 + cond2

    detect <- which(suma > 0)

    suma1.0 <- cond1[detect]
    suma2.0 <- cond2[detect]

    qq <- (1:100)

    cum1 <- cumsum(sort(suma1.0, decreasing = TRUE))/sum(suma1.0)
    cum2 <- cumsum(sort(suma2.0, decreasing = TRUE))/sum(suma2.0)

    nu1 <- length(suma1.0)
    nu2 <- length(suma2.0)

    yy1 <- cum1[round(nu1*qq/100, 0)]*100
    yy2 <- cum2[round(nu2*qq/100, 0)]*100

    plot(qq, yy1, xlab = "% detected genes", ylab = "% cumulative reads",
         type = "b", col = 2, main = tit, xlim = xlim, ylim = ylim,
         cex.main = 1.7, cex.lab = 1.5, cex.axis = 1.3)

    points(qq, yy2, type = "b", col = 4)

#    ks.resul <- ks.test(suma1.0, suma2.0)

#    text(50, 70, "Kolmogorov-Smirnov test", cex = 1.5)

#    text(50, 65, paste("p-value =", ks.resul$p.value), cex = 1.5)

    legend("bottom", legend = legend, text.col = c(2,4), bty = "n",
           lty = 1, lwd = 2, col = c(2,4), cex = 1.5)
    

  } else {


    suma <- apply(cbind(cond1, cond2), 1, sum)

    detect <- which(suma > 0)

    cond1.0 <- as.matrix(cond1)[detect,]

    cond2.0 <- as.matrix(cond2)[detect,]

    qq <- (1:100)

    cum1 <- apply(cond1.0, 2, function(x) {
                      cumsum(sort(x, decreasing = TRUE)) / sum(x) } )
    cum2 <- apply(cond2.0, 2, function(x) {
                      cumsum(sort(x, decreasing = TRUE)) / sum(x) } )

    nu1 <- nrow(cond1.0)
    nu2 <- nrow(cond2.0)

    yy1 <- cum1[round(nu1*qq/100, 0),]*100
    yy2 <- cum2[round(nu2*qq/100, 0),]*100

    plot(qq, yy1[,1], xlab = "% detected genes", ylab = "% cumulative reads",
         type = "l", col = 2, main = tit, ylim = c(0,100), lwd = 1)

    for (i in 2:ncol(cond1.0)) {
      lines(qq, yy1[,i], col = 2, lwd = 1)
    }

    for (i in 1:ncol(cond2.0)) {
      lines(qq, yy2[,i], col = 4, lwd = 1)
    }

    legend("bottom", legend = legend, text.col = c(2,4), bty = "n",
           lty = 1, lwd = 2, col = c(2,4))
  }
}





#************************************************************************************#





## Mean length for detected genes Plot   REVISAR!!!!!!!!!!!


DL.plot <- function (datos1, datos2, long, k = 0, ylim = range(na.omit(long)),
                     legend = c(deparse(substitute(datos1)),
                       deparse(substitute(datos2))), tit = "") {

  n1 <- ncol(as.matrix(datos1))

  muestra1 <- as.list(1:n1)
  for (i in 2:n1) {
    muestra1 <- append(muestra1, combn(n1, i, simplify = FALSE))
  }

  varias1 <- vector("list", length = length(muestra1))
  names(varias1) <- sapply(muestra1, function(x) {
    paste("C1", x, collapse = ".", sep = "")})

  for (i in 1:length(muestra1)) {
    varias1[[i]] <- apply(as.matrix(datos1[,muestra1[[i]]]), 1, sum)
  }

  no0.1 <- lapply(varias1, noceros, k = k, num = FALSE)

  satura1 <- data.frame("seq.depth" = sapply(varias1, sum),
                        "lengths" = sapply(no0.1,
                          function (x) { median(na.omit(long[x])) }))

  n2 <- ncol(as.matrix(datos2))

  if (n1 == n2) {
    muestra2 <- muestra1
  } else {
    muestra2 <- as.list(1:n2)
    for (i in 2:n2) {
      muestra2 <- append(muestra2, combn(n2, i, simplify = FALSE))
    }
  }

  varias2 <- vector("list", length = length(muestra2))
  names(varias2) <- sapply(muestra2, function(x) {
    paste("C2", x, collapse = ".", sep = "")})

  for (i in 1:length(muestra2)) {
    varias2[[i]] <- apply(as.matrix(datos2[,muestra2[[i]]]), 1, sum)
  }

  no0.2 <- lapply(varias2, noceros, k = k, num = FALSE)

  satura2 <- data.frame("seq.depth" = sapply(varias2, sum),
                        "lengths" = sapply(no0.2,
                          function (x) { median(na.omit(long[x])) }))

  SS1 <- range(satura1$seq.depth)
  SS2 <- range(satura2$seq.depth)

  xM <- max(SS1[2], SS2[2])
  xm <- min(SS1[1], SS2[1])
  
  plot(satura1, pch = 16, col = 2, ylim = ylim, xlim = c(xm, xM),
       main = tit, type = "b", xlab = "Sequencing Depth",
       ylab = "Median length of genes with reads > 0")
  
  points(satura2, pch = 16, col = 4, type = "b")

  abline(h = median(long, na.rm = TRUE), lty = 2)
  
  legend("top", legend = legend, text.col = c(2,4), bty = "n",
         lty = 1, lwd = 2, col = c(2, 4))
  
}






#************************************************************************************#




## Mean length for detected genes Plot according to BIOTYPES


DLbio.dat <- function (datos1, datos2, long, k = 0, infobio, biotypes, nsim = 5)  {

  # infobio = vector containing biotype for each gene in "datos"
  # biotypes = list containing groups of biotypes to be studied
  # nsim = number of replicates to be simulated

  n1 <- NCOL(datos1)
  n2 <- NCOL(datos2)

  biog <- lapply(biotypes, function(x) { which(is.element(infobio, x)) })

  satura1 <- satura2 <- vector("list", length = length(biotypes))
  names(satura1) <- names(satura2) <- names(biotypes)

  newdet1 <- newdet2 <- NULL
  
  
  ## when replicates are available in condition 1
  if (n1 > 1) { 

    varias1 <- vector("list", length = n1) # counts for each sequencing depth
    names(varias1) <- paste(1:n1, "rep", sep = "")

    for(i in 1:n1) {
      
      muestra1 <- combn(n1, i, simplify = FALSE)

      if (length(muestra1) > 20) {
        sub20 <- sample(1:length(muestra1), size = 20, replace = FALSE)
        muestra1 <- muestra1[sub20]
      }

      sumrepli <- NULL

      for (com in 1:length(muestra1)) {
        sumrepli <- cbind(sumrepli, rowSums(as.matrix(datos1[,muestra1[[com]]])))
      }

      varias1[[i]] <- as.matrix(sumrepli)
      
    }
  }
  

  ## replicates simulated for condition 1
  if (n1 == 1) {  # replicates have to be simulated

    varias1 <- vector("list", length = nsim) # counts for each sequencing depth
    names(varias1) <- paste("dp", 1:nsim, sep = "")
    
    total1 <- sum(datos1)    

    for (i in 1:(nsim-1)) {    # total1*i/nbox reads (100% is calculated apart)

      muestra1 <- rmultinom(10, size = round(total1*i/nsim,0), prob = datos1)
      
      varias1[[i]] <- as.matrix(muestra1)
    }

    varias1[[nsim]] <- as.matrix(datos1)
  }


  ## sequencing depth for each new sample
  seq.depth1 <- sapply(varias1, function(x) { mean(colSums(x)) })


  ## computing length for each biotype (datos1)
    for (j in 1:length(satura1))  {

      conbio1 <- lapply(varias1, function(x) { as.matrix(x[biog[[j]]]) })
      long1 <- long[biog[[j]]]
      conbio1.0 <- lapply(conbio1,
                          function(x) {
                            apply(x, 2, function(y) {
                              median(long1[noceros(y, k = k, num = FALSE)], na.rm = TRUE)
                            })})

      dtc1 <- vector("list", length = length(conbio1))  # per each seq. depth

      for (ss in 1:length(conbio1)) {
        dtc1[[ss]] <- unique(unlist(apply(conbio1[[ss]], 2,
                                          noceros, k = k, num = FALSE)))
      }

      nous1 <- NA

      for (nn in 1:(length(dtc1)-1)) {
        nuevos <- setdiff(dtc1[[nn+1]], dtc1[[nn]])
        nous1 <- c(nous1, median(na.omit(long1[nuevos])))
      }                                   

      satura1[[j]] <- sapply(conbio1.0, mean, na.rm = TRUE)

      newdet1 <- rbind(newdet1, nous1)
    }

  rownames(newdet1) <- names(satura1)
  colnames(newdet1) <- paste("rep", 1:ncol(newdet1))



  ## when replicates are available in condition 2
  if (n2 > 1) { 

    varias2 <- vector("list", length = n2) # counts for each sequencing depth
    names(varias2) <- paste(1:n2, "rep", sep = "")

    for(i in 1:n2) {
      
      muestra2 <- combn(n2, i, simplify = FALSE)
      
      if (length(muestra2) > 20) {
        sub20 <- sample(1:length(muestra2), size = 20, replace = FALSE)
        muestra2 <- muestra2[sub20]
      }

      sumrepli <- NULL

      for (com in 1:length(muestra2)) {
        sumrepli <- cbind(sumrepli, rowSums(as.matrix(datos2[,muestra2[[com]]])))
      }

      varias2[[i]] <- as.matrix(sumrepli)
      
    }
  }

  

  ## replicates simulated for condition 2
  if (n2 == 1) {  # replicates have to be simulated

    varias2 <- vector("list", length = nsim) # counts for each sequencing depth
    names(varias2) <- paste("dp", 1:nsim, sep = "")
    
    total2 <- sum(datos2)    

    for (i in 1:(nsim-1)) {    # total1*i/nbox reads (100% is calculated apart)

      muestra2 <- rmultinom(10, size = round(total2*i/nsim,0), prob = datos2)
      
      varias2[[i]] <- as.matrix(muestra2)
    }

    varias2[[nsim]] <- as.matrix(datos2)
  }


  ## sequencing depth for each new sample (datos2)
  seq.depth2 <- sapply(varias2, function(x) { mean(colSums(x)) })


  ## computing saturation for each biotype (datos2)
  for (j in 1:length(satura2))  {

    conbio2 <- lapply(varias2, function(x) { as.matrix(x[biog[[j]]]) })
    long2 <- long[biog[[j]]]
    conbio2.0 <- lapply(conbio2,
                        function(x) {
                          apply(x, 2, function(y) {
                            median(long2[noceros(y, k = k, num = FALSE)], na.rm = TRUE)
                          })})

    dtc2 <- vector("list", length = length(conbio2))  # per each seq. depth

    for (ss in 1:length(conbio2)) {
      dtc2[[ss]] <- unique(unlist(apply(conbio2[[ss]], 2,
                                        noceros, k = k, num = FALSE)))
    }

    nous2 <- NA

    for (nn in 1:(length(dtc2)-1)) {
      nuevos <- setdiff(dtc2[[nn+1]], dtc2[[nn]])
      nous2 <- c(nous2, median(na.omit(long2[nuevos])))
    }

    satura2[[j]] <- sapply(conbio2.0, mean, na.rm = TRUE)

    newdet2 <- rbind(newdet2, nous2)
  }

  rownames(newdet2) <- names(satura2)
  colnames(newdet2) <- paste("rep", 1:ncol(newdet2))

  

  ## computing global length
  #length.biotype <- sapply(biog, function (x) { median(na.omit(long[x])) })
  length.biotype <- NULL
  
  total12 <- rowSums(cbind(datos1, datos2))
  totno0 <- noceros(total12, num = FALSE, k = 0)

  long.mayor150 <- which(long > 150)
  long.menor150 <- which(long <= 150)

  for (i in 1:length(biog)) {
    det.menor150 <- int.mult(list(biog[[i]], totno0, long.menor150))
    mayor150 <- intersect(long.mayor150, biog[[i]])
    estos <- union(det.menor150, mayor150)
    longestos <- na.omit(long[estos])
    length.biotype <- c(length.biotype, median(longestos))
  }
  
  satura <- list("cond1" = satura1, "cond2" = satura2,
                 "bionum" = sapply(biog, length),
                 "depth1" = seq.depth1, "depth2" = seq.depth2,
                 "newdet1" = newdet1, "newdet2" = newdet2,
                 "biolength" = length.biotype)
  satura
}





#************************************************************************************#





## PLOT: Mean length for detected genes Plot according to BIOTYPES

DLbio.plot <- function (depth1, depth2, sat1, sat2, biolong, k = 0,
                        xlim = NULL, ylim = NULL,
                        legend = c("sample1", "sample2"), main = "",
                        ylab = "Median length of detected genes",
                        cex.main = 1, cex.lab = 1, cex.axis = 1, cex = 1,
                        xlab = "Sequencing depth")  {

  # ylim for plot
  if (is.null(ylim)) {

    ylim <- c(min(c(na.omit(sat1), na.omit(sat2), biolong)),
              max(c(na.omit(sat1), na.omit(sat2), biolong)))

    ylim <- ylim + 0.1 * diff(range(ylim)) * c(-1,1)    
  }  

  # xlim for plot
  if (is.null(xlim)) {
    SS1 <- range(depth1)
    SS2 <- range(depth2)

    xM <- max(SS1[2], SS2[2])
    xm <- min(SS1[1], SS2[1])

    xlim <- c(xm, xM)
  }

  # PLOT
  if(is.na(biolong)) {
    
    plot(1:5, 1:5, type = "n", axes = FALSE, main = main, xlab = "", ylab = "",
         cex.main = cex.main)
    text(3, 4, "Biotype not found in the dataset", adj = 0.5, cex = cex.main, font = 2)
    
  } else {

    if( diff(range(ylim)) < 100 ) {   # correcting ylim
      ylim <- mean(ylim) + 50*c(-1,1)
      ylim[1] <- max(ylim[1], 0)
    }      

    plot(depth1, sat1, pch = 16, col = 2, ylim = ylim, xlim = xlim,
         main = main, type = "o", xlab = xlab, ylab = ylab,
         cex.main = cex.main, cex.lab = cex.lab, cex.axis = cex.axis)

    points(depth2, sat2, pch = 16, col = 4, type = "o")

    abline(h = biolong, lty = 2, col = "grey")

    text(mean(xlim), biolong + 0.02 * diff(range(ylim)),
         "median global length", col = "grey", cex = cex)

    legend("top", legend = legend, text.col = c(2,4), bty = "n",
           lty = 1, lwd = 2, col = c(2, 4), cex = cex, horiz = TRUE)

  }
  
}




#************************************************************************************#





## Detecting outliers

# lim: Values outside the interval lim = c(a,b) will be tagged as outliers.
#      "bp" = Limits are computed as in boxplots. q -/+ k*IR
#

# k:   Coefficient to compute "bp" limits. By default, k = 1.5.




outliers <- function (x, lim = "bp", k = 1.5) {
  
  if ( lim == "bp" ) {

    Q <- quantile(na.omit(x), c(0.25,0.75))
    IR <- diff(Q)

    lim <- c(Q[1] - k*IR, Q[2] + k*IR)

  }

  outl.left <- which(x < lim[1])
  outl.right <- which(x > lim[2])

  outl <- c(outl.left, outl.right)

  list("left" = outl.left, "right" = outl.right, "both" = outl)

}

  

    

######################################################################################



## NOISeq-sim (to study different performance aspects)   POR ARREGLAR!!!!!!!!!!!


noiseqsim <- function (datos1, datos2, k = 0.5, norm = "rpkm", long = 1000,
                       q = 0.90, pnr = c(1,1), nss = 5, v = 0.02, lc = 1,
                       expcond = 0)

  # expcond: 0 (if both samples are used to estimate noise distribution)
  #          1 (if sample 1 is used to estimate noise distribution)
  #          2 (if sample 2 is used to estimate noise distribution)
  
  {

    g.sinL <- names(which(is.na(long)))
      
    # Normalization
    if (norm == "n") {       # no normalization
      datos1.norm <- sinceros(datos1, k = k)
      datos2.norm <- sinceros(datos2, k = k)
    }
    if (norm == "rpkm") {      # RPKM
      datos1.norm <- rpkm(datos1, long = long, k = k, lc = lc)
      datos2.norm <- rpkm(datos2, long = long, k = k, lc = lc)
    }

    if (norm == "uqua") {
      datos.norm <- uqua(cbind(datos1, datos2), long = long, lc = lc, k = k)
      datos1.norm <- datos.norm[ ,1]
      datos2.norm <- datos.norm[ ,2]
    }

    
    # M and D for different experimental conditions
    MDs <- MD(dat = cbind(datos1.norm, datos2.norm))

    
    # Noise distribution with simulated samples

    if (nss == 0) {
        nss <- 5
      }

    if (expcond == 0) {
      datos.sim <- sim.samples2(counts1 = datos1, counts2 = datos2,
                                pnr = pnr, nss = nss, v = v)
      datitos <- cbind(datos.sim[[1]], datos.sim[[2]])
      dat.sim.norm <- vector("list", length = 2)
    }

    if (expcond == 1) {
      datos.sim <- sim.samples2(counts1 = datos1, counts2 = NULL,
                                pnr = pnr, nss = nss, v = v)
      datitos <- datos.sim[[1]]
      dat.sim.norm <- vector("list", length = 1)
    }

    if (expcond == 2) {
      datos.sim <- sim.samples2(counts1 = datos2, counts2 = NULL,
                                pnr = pnr, nss = nss, v = v)
      datitos <- datos.sim[[1]]
      dat.sim.norm <- vector("list", length = 1)
    }

    nn <- sapply(datos.sim, ncol)

    sumita <- apply(datitos, 1, sum)
    g.sin0 <- names(which(sumita > 0))
    gens.sin0 <- setdiff(g.sin0, g.sinL)

    if (norm == "n") {       # no normalization
      datitos.0 <- sinceros(datitos, k = k)
      datitos.norm <- datitos.0[gens.sin0, ]
    }

    if (norm == "rpkm") {      # RPKM
      datitos.0 <- rpkm(datitos, long = long, k = k, lc = lc)
      datitos.norm <- datitos.0[gens.sin0, ]
    }

    if (norm == "uqua") {      # Upper Quartile
      datitos.0 <- uqua(datitos, long = long, lc = lc, k = k)
      datitos.norm <- datitos.0[gens.sin0, ]
    }

    dat.sim.norm[[1]] <- datitos.norm[ ,1:nn[1]]
    MD1 <- MD(dat = dat.sim.norm[[1]])

    if (expcond == 0) {
      dat.sim.norm[[2]] <- datitos.norm[ ,(nn[1]+1):sum(nn)]
      MD2 <- MD(dat = dat.sim.norm[[2]])
    } else { MD2 <- NULL }

    Mr <- c(as.numeric(MD1$M), as.numeric(MD2$M))
    Dr <- c(as.numeric(MD1$D), as.numeric(MD2$D))
   

    # Differential expression

    prob <- probdeg(MDs$M, MDs$D, Mr, Dr)
    deg <- DEG.q(prob, q = q)

    list("probab" = prob, "deg" = deg, "Ms" = MDs$M, "Ds" = MDs$D,
       "Mn" = Mr, "Dn" = Dr)

  }





#************************************************************************************#




## To generate simulated samples (version 2):

# pnr: Percentage of total reads (seq.depth). By default, pnr = 1.
# nss: Number of simulated samples (>= 2). By default, nss = 5.
# v:   Variability in sample total reads.

# Sample total reads is computed as a random value from a uniform distribution
# in the interval [(pnr-v)*sum(counts), (pnr+v)*sum(counts)]


sim.samples2 <- function(counts1, counts2 = NULL, pnr = 1, nss = 5, v = 0.02) {

  seqdep <- c(sum(counts1), sum(counts2))

  num.reads1 <- (pnr + c(-v,v))*seqdep[1]

  muestras <- vector("list")
  
  muestras$c1 <- NULL
  for (s in 1:nss) {
    tama <- round(runif(1, num.reads1[1], num.reads1[2]), 0)
    muestras$c1 <- cbind(muestras$c1,
                         rmultinom(1, size = tama, prob = counts1))
  }

  if(!is.null(counts2)) {
    num.reads2 <- (pnr + c(-v,v))*seqdep[2]
    muestras$c2 <- NULL
    for (s in 1:nss) {
      tama <- round(runif(1, num.reads2[1], num.reads2[2]), 0)
      muestras$c2 <- cbind(muestras$c2,
                           rmultinom(1, size = tama, prob = counts2))
    }
  }

  muestras
}








###################################################################################


#### Results for ROC curve

rocdata <- function(comparing, method.type, features, positives, qroc = seq(1,0,-0.005),
                    aroc = seq(0,1, 0.005))  {

  # positives = Vector of features presenting real differential expression.
  # qroc = Threshold probabilities to compute points for ROC curve for those methods
  #        with differential expression probabilities output.
  # aroc = Threhold significances to compute points for ROC curve for those methods
  #        whose output are p-values.
  # comparing = List with length the number of methods to be compared. Each element of
  #             the list is a vector cointaining the probabilities or p-values for
  #             each gene or biological feature under study.
  # method.type = Vector of the same length of "comparing" with "q" or "a" in each
  #               element, according to the type of method in "comparing".
  # features = List with length the number of methods to be compared. Each element of
  #            the list is a vector cointaining the names of the biological features
  #            (genes or whatever) associated to the probabilities or p-values in
  #            the list "comparing".

  P <- length(positives)
  N <- length(comparing[[1]]) - P

  roc <- vector("list", length =  length(comparing))
  names(roc) <- names(comparing)

  for (k in 1:length(comparing)) {

    taula <- NULL
    gg <- as.character(features[[k]])

    if (method.type[k] == "q") {

      for (i in 1:length(qroc)) {
        deg <- gg[which(comparing[[k]] >= qroc[i])]
        TP <- long.int(deg, positives)
        TPR <- TP / P
        FPR <- (length(deg)-TP) / N
        taula <- rbind(taula, c(FPR,TPR))
      }
    }

    if (method.type[k] == "a") {

      for (i in 1:length(aroc)) {
        deg <- gg[which(comparing[[k]] <= aroc[i])]
        TP <- long.int(deg, positives)
        TPR <- TP / P
        FPR <- (length(deg)-TP) / N
        taula <- rbind(taula, c(FPR,TPR))
      }
    }

    roc[[k]] <- taula

  }

  roc

}




#************************************************************************************#



#### ROC curve plot

rocplot <- function(data, col, lin, xlim = c(0,1), ylim = c(0,1),
                    xlab = "FPR", ylab = "TPR", main = "") {

  # data = List coming from the "rocdata" function. Each element of the list is a two
  #        columns matrix with (x,y)-coordinates for the plot.
  # col =  Vector of colors for lines with the same length as ROC data list.
  # lin = Vector with type of lines, with the same length as ROC data list.

  plot(data[[1]], type = "l", col = col[1], lty = lin[1], lwd = 2,
       xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, main = main)

  if (length(data) > 1) {

    for (i in 2:length(data)) {
      lines(data[[i]], col = col[i], lwd = 2, lty = lin[i])
    }
  }

  abline(a = 0, b = 1, lty = 3)

  legend("bottomright", names(data), lwd = 3, col = col, lty = lin)

}





#************************************************************************************#


#### Results for FDR curve

fdrdata <- function(comparing, method.type, num, positives, features) {

  # positives = Vector of features presenting real differential expression.
  # num = Vector containing the number of selected features to be drawn in X-axis.
  # comparing = List with length the number of methods to be compared. Each element of
  #             the list is a vector cointaining the probabilities or p-values for
  #             each gene or biological feature under study.
  # method.type = Vector of the same length of "comparing" with "q" or "a" in each
  #               element, according to the type of method in "comparing".
  # features = List with length the number of methods to be compared. Each element of
  #            the list is a vector cointaining the names of the biological features
  #            (genes or whatever) associated to the probabilities or p-values in
  #            the list "comparing".

  P <- length(positives)
  N <- length(comparing[[1]]) - P

  fdr <- vector("list", length =  length(comparing))
  names(fdr) <- names(comparing)

  for (k in 1:length(comparing)) {

    taula <- NULL
    gg <- as.character(features[[k]])

    if (method.type[k] == "q") {
      ordre <- rank(1-comparing[[k]])
    }

    if (method.type[k] == "a") {
      ordre <- rank(comparing[[k]])
    }

    taula <- c(0,0)

    for (i in 2:length(num)) {
      deg <- gg[which(ordre <= num[i])]
      TP <- long.int(deg, positives)
      FDR <- 1- (TP / length(deg))
      taula <- rbind(taula, c(length(deg),FDR))
    }

    fdr[[k]] <- taula

  }

  fdr
  
}



#************************************************************************************#



#### FDR plot

fdrplot <- function(data, col, lin, xlim = c(0,1000), ylim = c(0,1),
                    xlab = "Selected features", ylab = "FDR", main = "") {

  # data = List coming from the "fdrdata" function. Each element of the list is a two
  #        columns matrix with (x,y)-coordinates for the plot.
  # col =  Vector of colors for lines with the same length as ROC data list.
  # lin = Vector with type of lines, with the same length as ROC data list.

  plot(data[[1]], type = "l", col = col[1], lty = lin[1], lwd = 2,
       xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, main = main)

  if (length(data) > 1) {

    for (i in 2:length(data)) {
      lines(data[[i]], col = col[i], lwd = 2, lty = lin[i])
    }
  }

  legend("bottomright", names(data), lwd = 3, col = col, lty = lin)

}








#************************************************************************************#





## Counts for detected genes Plot according to BIOTYPES


countsbio.dat <- function (datos1, datos2, k = 0, infobio, biotypes, nbox = 5)  {

  # infobio = vector containing biotype for each gene in "datos"
  # biotypes = list containing groups of biotypes to be studied
  # nbox = number of different depths to be plotted (only used for simulated data)

  n1 <- NCOL(datos1)
  n2 <- NCOL(datos2)

  biog <- lapply(biotypes, function(x) { which(is.element(infobio, x)) })
  # which genes belong to each biotype

  satura1 <- satura2 <- vector("list", length = length(biotypes))
  names(satura1) <- names(satura2) <- names(biotypes)


  ## replicates available for condition 1
  if (n1 > 1) {

    varias1 <- vector("list", length = n1) # counts for each sequencing depth
    names(varias1) <- paste(1:n1, "rep", sep = "")
    
    for(i in 1:n1) {

      muestra1 <- combn(n1, i, simplify = FALSE)

      if (length(muestra1) > 20) {
        sub20 <- sample(1:length(muestra1), size = 20, replace = FALSE)
        muestra1 <- muestra1[sub20]
      }

      sumrepli <- NULL

      for (com in 1:length(muestra1)) {
        sumrepli <- cbind(sumrepli, rowSums(as.matrix(datos1[,muestra1[[com]]])))
      }

      varias1[[i]] <- rowMeans(sumrepli)
    }
  }


  ## replicates simulated for condition 1
  if (n1 == 1) {  # replicates have to be simulated

    varias1 <- vector("list", length = nbox) # counts for each sequencing depth
    names(varias1) <- paste("dp", 1:nbox, sep = "")
    
    total1 <- sum(datos1)    

    for (i in 1:(nbox-1)) {    # total1*i/nbox reads (100% is calculated apart)

      muestra1 <- rmultinom(10, size = round(total1*i/nbox,0), prob = datos1)
      
      varias1[[i]] <- rowMeans(muestra1)
    }

    varias1[[nbox]] <- datos1
  }

  seq.depth1 <- sapply(varias1, sum)  # sequencing depth for each new sample


  
  
  ## selecting detected genes for each biotype (datos1)
  for (j in 1:length(satura1))  {

    satura1[[j]] <- vector("list", length = length(varias1))
    names(satura1[[j]]) <- names(varias1)

    # selecting genes in bioclass j for each sample
    conbio1 <- lapply(varias1, function(x) { x[biog[[j]]] })  

    for (i in 1:length(conbio1)) {
            
      # selecting the genes with counts > k
      noK <- noceros(conbio1[[i]], k = k, num = FALSE)

      if (is.null(noK)) {
        satura1[[j]][[i]] <- NA
      } else {
        satura1[[j]][[i]] <- conbio1[[i]][noK]
      }
    }
  }
  

  ## replicates available for condition 2
  if (n2 > 1) {

    varias2 <- vector("list", length = n2)
    names(varias2) <- paste(1:n2, "rep", sep = "")
    
    for (i in 1:n2) {

      muestra2 <- combn(n2, i, simplify = FALSE)
      
      if (length(muestra2) > 20) {
        sub20 <- sample(1:length(muestra2), size = 20, replace = FALSE)
        muestra2 <- muestra2[sub20]
      }

      sumrepli2 <- NULL

      for (com in 1:length(muestra2)) {
        sumrepli2 <- cbind(sumrepli2, rowSums(as.matrix(datos2[,muestra2[[com]]])))
      }

      varias2[[i]] <- rowMeans(sumrepli2)
    }
  }

    
  ## replicates simulated for condition 2
  if (n2 == 1) {  # replicates have to be simulated

    varias2 <- vector("list", length = nbox) # counts for each sequencing depth
    names(varias2) <- paste("dp", 1:nbox, sep = "")
    
    total2 <- sum(datos2)    

    for (i in 1:(nbox-1)) {    # total1*i/nbox reads (100% is calculated apart)

      muestra2 <- rmultinom(10, size = round(total2*i/nbox,0), prob = datos2)
      
      varias2[[i]] <- rowMeans(muestra2)
    }

    varias2[[nbox]] <- datos2
  }

  seq.depth2 <- sapply(varias2, sum)  # sequencing depth for each new sample

  
  ## selecting detected genes for each biotype (datos2)
  for (j in 1:length(satura2))  {

    satura2[[j]] <- vector("list", length = length(varias2))
    names(satura2[[j]]) <- names(varias2)

    # selecting genes in bioclass j for each sample
    conbio2 <- lapply(varias2, function(x) { x[biog[[j]]] })  

    for (i in 1:length(conbio2)) {
            
      # selecting the genes with counts > k
      noK <- noceros(conbio2[[i]], k = k, num = FALSE)

      if (is.null(noK)) {
        satura2[[j]][[i]] <- NA
      } else {
        satura2[[j]][[i]] <- conbio2[[i]][noK]
      }
    }
  }
  

  ## results
  satura <- list("cond1" = satura1, "cond2" = satura2,
                 "bionum" = sapply(biog, length),
                 "depth1" = seq.depth1, "depth2" = seq.depth2)                 
  satura
}




#************************************************************************************#





## PLOT: Mean length for detected genes Plot according to BIOTYPES

countbio.plot <- function (depth1, depth2, sat1, sat2, bionum, 
                           legend = c("sample1", "sample2"), main = "",
                           ylab = "# counts of detected features",
                           xlab = "Sequencing Depth", ylim = NULL, las = 0,
                           cex.main = 1, cex.lab = 1, cex.axis = 1, cex = 1)  {

  if (bionum == 0) {

    plot(1:5, 1:5, type = "n", main = main, cex.main = cex.main, axes = FALSE,
         xlab = "", ylab = "")
    text(3, 4, "Biotype not found in the dataset", font = 2, adj = 0.5,
         cex = cex.main)
    
  } else {

    # ylim for plot
    if (is.null(ylim)) {
      ylim <- c(min(c(na.omit(unlist(sat1)), na.omit(unlist(sat2)))),
                max(c(na.omit(unlist(sat1)), na.omit(unlist(sat2)))))
      ylim <- ylim + 0.05 * diff(range(ylim)) * c(-1,1)
    }

    # boxplots data

    depth <- c(depth1, depth2)
    d.sort <- sort(depth)
    d.order <- order(depth)

    n1 <- length(depth1)
    n2 <- length(depth2)

    databox <- vector("list", length = n1+n2)
    names(databox) <- d.sort

    colo <- NULL

    for (i in 1:(n1+n2)) {

      if ( d.order[i] > n1 ) {

        dd <- d.order[i] - n1
        databox[[i]] <- sat2[[dd]]
        colo <- c(colo, 4)

      } else {

        dd <- d.order[i]
        databox[[i]] <- sat1[[dd]]
        colo <- c(colo, 2)
      }
    }

    # BOXPLOT
    boxplot(databox, col = colo, ylim = ylim, main = main,
            type = "b", xlab = xlab, ylab = ylab,
            cex.main = cex.main, cex.lab = cex.lab, cex.axis = cex.axis)

    legend("top", legend = legend, text.col = c(2,4), bty = "o", bg = "white",
           fill = c(2, 4), cex = cex, horiz = TRUE)

    cuantos <- sapply(databox, length)
    cuantos1 <- which(cuantos == 1)
    sumant <- sapply(databox, sum, na.rm = TRUE)
    sumant0 <- which(sumant == 0)
    cuantosNA <- intersect(cuantos1, sumant0)

    cuantos[cuantosNA] <- 0

    mtext(cuantos, 3, at = 1:length(databox), cex = 0.6*cex, las = las)

  }
}










###################################################################################


#### Results for PRECISION-RECALL (PR) curve

prdata <- function(comparing, method.type, features, positives, qroc = seq(1,0,-0.005),
                   aroc = seq(0,1, 0.005))  {

  # positives = Vector of features presenting real differential expression.
  # qroc = Threshold probabilities to compute points for ROC curve for those methods
  #        with differential expression probabilities output.
  # aroc = Threhold significances to compute points for ROC curve for those methods
  #        whose output are p-values.
  # comparing = List with length the number of methods to be compared. Each element of
  #             the list is a vector cointaining the probabilities or p-values for
  #             each gene or biological feature under study.
  # method.type = Vector of the same length of "comparing" with "q" or "a" in each
  #               element, according to the type of method in "comparing".
  # features = List with length the number of methods to be compared. Each element of
  #            the list is a vector cointaining the names of the biological features
  #            (genes or whatever) associated to the probabilities or p-values in
  #            the list "comparing".

  P <- length(positives)
  N <- length(comparing[[1]]) - P

  pr <- vector("list", length =  length(comparing))
  names(pr) <- names(comparing)

  for (k in 1:length(comparing)) {

    taula <- NULL
    gg <- as.character(features[[k]])
    n <- length(gg)

    if (method.type[k] == "q") {

      for (i in 1:length(qroc)) {
        deg <- gg[which(comparing[[k]] >= qroc[i])]
        TP <- long.int(deg, positives)
        Pd <- length(deg)
        FP <- Pd-TP
        TN <- N-FP
        TPR <- TP / P # recall
        preci <- TP / Pd
        F1 <- 2*TPR*preci/(TPR+preci)
        nume <- TP*TN-FP*(P-TP)
        deno <- sqrt(Pd)*sqrt(P)*sqrt(N)*sqrt(n-Pd)        
        MCC <- nume/deno
        taula <- rbind(taula, c(TPR, preci, 1-qroc[i], F1, MCC))
      }
      colnames(taula) <- c("recall", "precision", "1-q", "F1", "MCC")
    }

    if (method.type[k] == "a") {

      for (i in 1:length(aroc)) {
        deg <- gg[which(comparing[[k]] <= aroc[i])]
        TP <- long.int(deg, positives)
        Pd <- length(deg)
        FP <- Pd-TP
        TN <- N-FP
        TPR <- TP / P
        preci <- TP / length(deg)
        F1 <- 2*TPR*preci/(TPR+preci)
        nume <- TP*TN-FP*(P-TP)
        deno <- sqrt(Pd)*sqrt(P)*sqrt(N)*sqrt(n-Pd)
        MCC <- nume/deno        
        taula <- rbind(taula, c(TPR, preci, aroc[i], F1, MCC))
      }
      colnames(taula) <- c("recall", "precision", "alpha", "F1", "MCC")
    }

    
    pr[[k]] <- taula

  }

  pr

}




#************************************************************************************#



#### PRECISION-RECALL curve

prplot <- function(data, col, lin, xlim = c(0,1), ylim = c(0,1), lty = 1,
                    xlab = "Recall", ylab = "Precision", main = "", plottype = "PR") {

  # data = List coming from the "prdata" function. Each element of the list is a two
  #        columns matrix with (x,y)-coordinates for the plot.
  # col =  Vector of colors for lines with the same length as ROC data list.
  # lin = Vector with type of lines, with the same length as ROC data list.
  # plottype = If PR, Precision-Recall curve is plotted.
  #            If F1, q/a values are plotted against F1-score.
  #            If FDR, q/a values are plotted against False Discovery Rate.
  #            If MCC, q/a values are plotted against Matthew's Correlation Coefficient.

  if (length(lty) == 1) {
    lty <- rep(lty, length(col))
  }

  if (plottype == "PR") {

    plot(data[[1]][,1:2], type = "l", col = col[1], lwd = 2, lty = lty[1],
         xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, main = main)

    if (length(data) > 1) {

      for (i in 2:length(data)) {
        lines(data[[i]][,1:2], col = col[i], lwd = 2, lty = lty[i])
      }
    }

    legend("bottomleft", names(data), lwd = 3, col = col, lty = lty, bty = "n")

  }

  if (plottype == "F1") {

    plot(data[[1]][,3:4], type = "l", col = col[1], lty = lty[1],
         lwd = 2, xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, main = main)

    if (length(data) > 1) {

      for (i in 2:length(data)) {
        lines(data[[i]][,3:4], col = col[i], lwd = 2, lty = lty[i])
      }      
    }
    legend("bottom", names(data), lwd = 3, col = col, ncol = 2, lty = lty, bty = "n")
  }

  if (plottype == "FDR") {

    plot(data[[1]][,3], 1-data[[1]][,2], type = "l", col = col[1], lty = lty[1],
         lwd = 2, xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, main = main)

    if (length(data) > 1) {

      for (i in 2:length(data)) {
        lines(data[[i]][,3], 1-data[[i]][,2], col = col[i], lwd = 2, lty = lty[i])
      }      
    }
    legend("topleft", names(data), lwd = 3, col = col, ncol = 2, bty = "n", lty = lty)
  }

  if (plottype == "MCC") {

    plot(data[[1]][,3], data[[1]][,5], type = "l", col = col[1], lty = lty[1], 
         lwd = 2, xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, main = main)

    if (length(data) > 1) {

      for (i in 2:length(data)) {
        lines(data[[i]][,3], data[[i]][,5], col = col[i], lwd = 2, lty = lty[i])
      }      
    }
    legend("bottom", names(data), lwd = 3, col = col, ncol = 2, lty = lty, bty = "n")
  }  

}




#************************************************************************************#



#########################################################################################
#################   Plot with 2 different Y axis (left and right)   #####################
#########################################################################################


# By Ajay Shah (taken from [R] Plot 2 time series with different y axes (left and right),
# in https://stat.ethz.ch/pipermail/r-help/2004-March/047775.html) 

# Modified by: Sonia Tarazona

### PARAMETERS (default):
# x: data to be drawn on X-axis
# yright: data to be drawn on Y right axis
# yleft: data to be drawn on Y left axis
# yrightlim (range(yright, na.rm = TRUE)): ylim for rigth Y-axis
# yleftlim (range(yleft, na.rm = TRUE)): ylim for left Y-axis
# xlab (NULL): Label for X-axis
# yylab (c("","")): Labels for right and left Y-axis
# pch (c(1,2)): Type of symbol for rigth and left data
# col (c(1,2)): Color for rigth and left data
# linky (TRUE): If TRUE, points are connected by lines.
# smooth (0): Friedman's super smoothing
# lwds (1): Line width for smoothed line
# length (10): Number of tick-marks to be drawn on axis
# ...: Other graphical parameters to be added by user (such as main, font, etc.)
###



plot.y2 <- function(x, yright, yleft, yrightlim = range(yright, na.rm = TRUE),
                    yleftlim = range(yleft, na.rm = TRUE), xlim = range(x, na.rm = TRUE),
                    xlab = NULL, yylab = c("",""), lwd = c(2,2),
                    pch = c(1,2), col = c(1,2), type = c("b","b"),
                    linky = TRUE, smooth = 0,
                    lwds = 1, length = 10, ...,
                    x2 = NULL, yright2 = NULL, yleft2 = NULL, col2 = c(3,4))
{
  #par(mar = c(5,2,4,2), oma = c(0,3,0,3))

  ## Plotting RIGHT axis data

  plot(x, yright, ylim = yrightlim, axes = FALSE, ylab = "", xlab = xlab, xlim = xlim,
       pch = pch[1], type = type[1], lwd = lwd[1], col = col[1], ...)
  
  axis(4, pretty(yrightlim, length), col = 1, col.axis = 1)

  if (is.null(yright2) == FALSE) {
    points(x2, yright2, type = type[1], pch = pch[1], lwd = lwd[1], col = col2[1], ...)
  }
  
  #if (linky) lines(x, yright, col = col[1], ...)
  
  if (smooth != 0) lines(supsmu(x, yright, span = smooth), col = col[1], lwd = lwds, ...)
  
  if(yylab[1]=="") {
    mtext(deparse(substitute(yright)), side = 4, outer = FALSE, line = 2,
          col = 1,...)
  } else {
    mtext(yylab[1], side = 4, outer = FALSE, line = 2, col = 1, ...)
  }
  

  par(new = T)

  ## Plotting LEFT axis data
  
  plot(x, yleft, ylim = yleftlim, axes = FALSE, ylab = "" , xlab = xlab, xlim = xlim,
       pch = pch[2], type = type[2], lwd = lwd[2], col = col[2], ...)
  
  box()
  
  axis(2, pretty(yleftlim, length), col = 1, col.axis = 1)

  if (is.null(yleft2) == FALSE) {
    points(x2, yleft2, type = type[2], pch = pch[2], lwd = lwd[2], col = col2[2], ...)
  }
  

  #if (linky) lines(x, yleft, col = col[2], ...)
  
  if (smooth != 0) lines(supsmu(x, yleft, span = smooth), col = col[2], lwd=lwds, ...)
  
  if(yylab[2] == "") {
    mtext(deparse(substitute(yleft)), side = 2, outer = FALSE, line = 2, col = 1, ...)
  } else {
    mtext(yylab[2], side = 2, outer = FALSE, line = 2, col = 1, ...)
  }
  
  
  ## X-axis
  axis(1, at = pretty(xlim, length))
  
   
}









###############################################################################
###############################################################################





## Global Saturation Plot including new detection per million reads


satur.plot2 <- function (datos, yleftlim = NULL, yrightlim = NULL, k = 5,
                         tit = "Saturation", cex.main = cex.main, scale = 10^6,
                         xlab = "Sequencing depth (million reads)",
                         cex.lab = cex.lab, cex.axis = cex.axis, cex = cex,
                         legend = c(deparse(substitute(datos1)),
                                    deparse(substitute(datos2)))) {


# For datos1

  bioseq1 <- NULL

  for (i in 1:length(datos$cond1)) {
    bioseq1 <- rbind(bioseq1, datos$cond1[[i]])
  }

  tot1 <- colSums(bioseq1, na.rm = TRUE)

  satura1 <- data.frame("seq.depth" = datos$depth1,
                        "detections" = tot1)



  # new detections (sample 1)

  puntos1 <- data.frame("x" = satura1$seq.depth, "y" = satura1$detections)

  pendi <- NULL

  for (i in 2:nrow(puntos1)) {

    nuevo <- max(0, (puntos1$y[i]-puntos1$y[i-1])/
                   (puntos1$x[i]-puntos1$x[i-1]))

    pendi <- c(pendi, nuevo)
  }

  newdet1 <- c(NA, pendi)*1000000





# For datos2

  bioseq2 <- NULL

  for (i in 1:length(datos$cond2)) {
    bioseq2 <- rbind(bioseq2, datos$cond2[[i]])
  }

  tot2 <- colSums(bioseq2, na.rm = TRUE)

  satura2 <- data.frame("seq.depth" = datos$depth2,
                        "detections" = tot2)

 

    # new detections (sample 2)

  puntos2 <- data.frame("x" = satura2$seq.depth, "y" = satura2$detections)

  pendi <- NULL

  for (i in 2:nrow(puntos2)) {

    nuevo <- max(0, (puntos2$y[i]-puntos2$y[i-1])/
                    (puntos2$x[i]-puntos2$x[i-1]))

    pendi <- c(pendi, nuevo)
  }

  newdet2 <- c(NA, pendi)*1000000


  


    ## PLOT LIMITS

    # yleftlim for plot.y2
    if (is.null(yleftlim)) { 
      yleftlim <- c(0, max(c(tot1, tot2),na.rm = TRUE))
    }


    # xlim
    SS1 <- range(satura1$seq.depth/scale)
    SS2 <- range(satura2$seq.depth/scale)

    xM <- max(SS1[2], SS2[2])
    xm <- min(SS1[1], SS2[1])


    # yrightlim for plot.y2
    if (is.null(yrightlim)) {
      maxi <- max(na.omit(c(newdet1, newdet2)))
      if (maxi < 10) {
        yrightlim <- c(0, max(10,maxi))
      } else {
        yrightlim <- c(0, maxi*1.1)
      }
    }

    

    ## PLOT for SAMPLES 1 and 2

    #par(xpd=TRUE, mar=par()$mar+c(0,0,2,0))

    plot.y2(x = satura1$seq.depth/scale, yright = newdet1, yleft = satura1$detections,
            type = c("h", "o"), lwd = c(10,2), pch = c(1, 19), 
            yrightlim = yrightlim, yleftlim = yleftlim,
            xlim = c(xm, xM), main = tit, 
            xlab = xlab, col = c("pink", 2),
            yylab = c("New detection per million reads",
                      paste("Number of genes with reads >", k)), cex = cex,
            cex.main = cex.main, cex.lab = cex.lab, cex.axis = cex.axis,
            x2 = satura2$seq.depth/scale, yright2 = newdet2,
            yleft2 = satura2$detections, col2 = c("lightblue1", 4))

    legend.text <- c(paste(legend[1], "(left axis)"),
                     paste(legend[2], "(left axis)"),
                     paste(legend[1], "(right axis)"),
                     paste(legend[2], "(right axis)"))

    legend.pch <- c(16, 16, 15, 15)

    legend.col <- c(2, 4, "pink", "lightblue1")

    legend("top",
           legend = legend.text, pch = legend.pch,
           col = legend.col, ncol = 2, bty = "n")
 
    par(mar=c(5, 4, 4, 2) + 0.1)

  list("sample1" = satura1, "sample2" = satura2)  ## results

}





#***************************************************************************#





####  FEATURES DETECTION per BIOTYPE


biodetection <- function(misdatos1, misdatos2 = NULL, infobio, k = 0,
                         leg = c("Sample 1", "Sample 2"), main = NULL) {

  detect1 <- misdatos1 > k

  genome <- 100*table(infobio)/sum(table(infobio))

  ordre <- order(genome, decreasing = TRUE)

  perdet1 <- genome*table(infobio, detect1)[names(genome),2]/
                    table(infobio)[names(genome)]

  perdet2 <- 100*table(infobio, detect1)[names(genome),2] /
                 sum(table(infobio, detect1)[,2])

  ceros <- rep(0, length(genome))

  biotable <- as.matrix(rbind(genome[ordre], perdet1[ordre],
                              perdet2[ordre], ceros))
  rownames(biotable) <- c("genome", "detectionVSgenome", "detectionVSsample",
                          "ceros")

 #ymax1 <- round(max(biotable[,1:3])*1.2,0)
  ymax1 <- ceiling(max(biotable[,1:3], na.rm = TRUE))
  ymax1sin <- max(biotable[,-c(1:3)], na.rm = TRUE)

 
  if (!is.null(misdatos2)) {


    detect2 <- misdatos2 > k

    perdet3 <- genome*table(infobio, detect2)[names(genome),2]/
                      table(infobio)[names(genome)]

    perdet4 <- 100*table(infobio, detect2)[names(genome),2] /
                   sum(table(infobio, detect2)[,2])

    biotable2 <- as.matrix(rbind(genome[ordre], perdet3[ordre],
                                 perdet4[ordre], ceros))
    rownames(biotable2) <- c("genome", "detectionVSgenome", "detectionVSsample",
                             "ceros")

 
  #ymax3 <- round(max(biotable2[,1:3])*1.2,0)
  ymax2 <- ceiling(max(biotable2[,1:3], na.rm = TRUE))
  ymax2sin <- max(biotable2[,-c(1:3)], na.rm = TRUE)

  ymaxL <- max(ymax1, ymax2)
  ymaxR <- ceiling(max(ymax1sin, ymax2sin))

  # scaling data on the right (datos1)
  biotable2b <- biotable2
  biotable2b[,-c(1:3)] <- biotable2b[,-c(1:3)]*ymaxL/ymaxR

} else {
 
  ymaxL <- ymax1
  ymaxR <- ceiling(ymax1sin)

}

  # scaling data on the right (datos1)
  biotable1 <- biotable
  biotable1[,-c(1:3)] <- biotable1[,-c(1:3)]*ymaxL/ymaxR



  if (is.null(misdatos2)) {

    # Plot (1 sample)
    par(mar = c(10, 4, 2, 2))

    barplot(biotable1[c(1,3),], main = main,
            xlab = NULL, ylab = "%features", axis.lty = 1, legend = FALSE,
            beside = TRUE, col = c("grey", 2), las = 2,
            ylim = c(0, ymaxL), border = c("grey", 2))

    barplot(biotable1[c(2,4),], main = main,
            xlab = NULL, ylab = "%features", axis.lty = 1, legend = FALSE,
            beside = TRUE, col = c(2, 1), las = 2, density = 30,
            ylim = c(0, ymaxL), border = 2, add = TRUE)

    axis(side=4, at = pretty(c(0,ymaxL), n = 5), labels = round(pretty(c(0,ymaxL), n = 5)*ymaxR/ymaxL, 1))

    abline(v = 9.5, col = 3, lwd = 2, lty = 2)

    text(x = 9, y = ymaxL*1.1, "Left axis", col = 3, font = 3, adj = 1)
    text(x = 10, y = ymaxL*1.1, "Right axis", col = 3, font = 3, adj = 0)

    legend(x = "topright", bty = "n", horiz = FALSE,
           fill = c("grey", 2, 2), density = c(NA,30,NA),
           border = c("grey", 2, 2),
           legend = c("% in genome", "detected", "% in sample"))

    #dev.off()

  } else {

    # Plot (2 samples)

    par(mar = c(10, 4, 2, 2))

 
    # Datos1
    barplot(biotable1[c(1,3),],
            main = paste(main, leg[1]),
            xlab = NULL, ylab = "%features", axis.lty = 1, legend = FALSE,
            beside = TRUE, col = c("grey", 2), las = 2,
            ylim = c(0, ymaxL), border = c("grey", 2))

    barplot(biotable1[c(2,4),],
            main = paste(main, leg[1]),
            xlab = NULL, ylab = "%features", axis.lty = 1, legend = FALSE,
            beside = TRUE, col = c(2, 1), las = 2, density = 30,
            ylim = c(0, ymaxL), border = 2, add = TRUE)

    axis(side=4, at = pretty(c(0,ymaxL), n = 5), labels = round(pretty(c(0,ymaxL), n = 5)*ymaxR/ymaxL, 1))

    abline(v = 9.5, col = 3, lwd = 2, lty = 2)

    text(x = 9, y = ymaxL*1.1, "Left axis", col = 3, font = 3, adj = 1)
    text(x = 10, y = ymaxL*1.1, "Right axis", col = 3, font = 3, adj = 0)

    legend(x = "topright", bty = "n", horiz = FALSE,
           fill = c("grey", 2, 2), density = c(NA,30,NA), border = c("grey", 2, 2),
           legend = c("% in genome", "detected", "% in sample"))


    # Datos2
    barplot(biotable2b[c(1,3),],
            main = paste(main, leg[2]),
            xlab = NULL, ylab = "%features", axis.lty = 1, legend = FALSE,
            beside = TRUE, col = c("grey", 4), las = 2,
            ylim = c(0, ymaxL), border = c("grey", 4))

    barplot(biotable2b[c(2,4),],
            main = paste(main, leg[2]),
            xlab = NULL, ylab = "%features", axis.lty = 1, legend = FALSE,
            beside = TRUE, col = c(4, 1), las = 2, density = 30,
            ylim = c(0, ymaxL), border = 4, add = TRUE)

    axis(side=4, at = pretty(c(0,ymaxL), n = 5), labels = round(pretty(c(0,ymaxL), n = 5)*ymaxR/ymaxL, 1))

    text(x = 9, y = ymaxL*1.1, "Left axis", col = 3, font = 3, adj = 1)
    text(x = 10, y = ymaxL*1.1, "Right axis", col = 3, font = 3, adj = 0)

    abline(v = 9.5, col = 3, lwd = 2, lty = 2)

    legend(x = "topright", bty = "n", horiz = FALSE,
           fill = c("grey", 4, 4), density = c(NA,30,NA),
           border = c("grey", 4, 4),
           legend = c("% in genome", "detected", "% in sample"))
  }
}


###########################################################################
###########################################################################



## DLbio.plot with second Y-axis (median length of new detections)

DLbio2.plot <- function(depth1, depth2, sat1, sat2, newdet1 = NULL, newdet2 = NULL,
                        xlim = NULL, ylim = NULL, biolong, k = 0,
                        main = NULL, xlab = "Sequencing depth",
                        samples = c("sample1", "sample2"), 
                        ylab = "Median length", cex.main = 1, cex.lab = 1,
                        cex.axis = 1, cex = 1) {  

  # ylim for plot
  if (is.null(ylim)) {

    ylim <- range(c(sat1, sat2, biolong, newdet1, newdet2), na.rm = TRUE)

    ylim <- ylim + 0.2 * diff(ylim) * c(-1,1)
    ylim[1] <- max(0, ylim[1])
  }  

  # xlim for plot
  if (is.null(xlim)) {
    SS1 <- range(depth1)
    SS2 <- range(depth2)

    xM <- max(SS1[2], SS2[2])
    xm <- min(SS1[1], SS2[1])

    xlim <- c(xm, xM)
  }

  # PLOT
  if(is.na(biolong)) {
    
    plot(1:5, 1:5, type = "n", axes = FALSE, main = main, xlab = "", ylab = "",
         cex.main = cex.main)
    text(3, 4, "Biotype not found in the dataset", adj = 0.5, cex = cex.main, font = 2)
    
  } else {

    if( diff(ylim) < 100 ) {   # correcting ylim
      ylim <- mean(ylim) + 50*c(-1,1)
      ylim[1] <- max(ylim[1], 0)
    }

    plot(depth1, newdet1, type = "h", lwd = 5, col = "pink", xlim = xlim, ylim = ylim,
         main = main, xlab = xlab, ylab = ylab, cex.main = cex.main,
         cex.lab = cex.lab, cex.axis = cex.axis)

    lines(depth2, newdet2, type = "h", lwd = 5, col = "lightblue1")

    points(depth1, sat1, pch = 16, col = 2, type = "o")

    points(depth2, sat2, pch = 16, col = 4, type = "o")

    abline(h = biolong, lty = 2, col = "grey")

    text(mean(xlim), biolong + 0.02 * diff(ylim),
         "median global length", col = "grey", cex = cex)

    legend("top", legend = c("All detected", samples, "New detections", samples),
           fill = c("white", 2,4, "white", "pink", "lightblue1"),
           border = "white",
           bty = "n", cex = cex, ncol = 2)    

  }
  
}




##************************************************************************************##


## RANKING STATISTIC FOR Gene Set Analysis

ranking <- function(results) {

  M <- results$Ms

  D <- results$Ds

  prob <- results$probab

  ## Changing NA by 0
  M[is.na(M)] <- 0
  D[is.na(D)] <- 0
  prob[is.na(prob)] <- 0


  ## Ranking

#   ranking1 <- M*prob
# 
#   ranking2 <- sign(M)*prob
# 
#   ranking3 <- M*D
# 
#   ranking4 <- M*D*prob

  ranking5 <- sqrt(M*M + D*D)*sign(M)


  ## Ranking results
  #list(ranking1, ranking2, ranking3, ranking4, ranking5)
  theranking <- data.frame(names(ranking5), ranking5)
  rownames(theranking) <- NULL
  colnames(theranking) <- c("ID", "statistic")
  
  theranking

}





