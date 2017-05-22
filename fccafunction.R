library(CCA)
##### functions start ###################
exonNet <-
  function (dat, tor = 0.3, cutMethod = c("value", "number"), cutValue = 0.05, 
            cutNum = NULL) 
  {   
    indexlist <- split(1:length(dat[, 1]), dat[, 1])
    #index <- combn(1:length(indexlist), 2)
    
    nlen<-length(dat[,1])
    clen<-choose(nlen,2)
    t=1
    index<-matrix(NA,2,clen)
    for (i in 1:(nlen-1))
    {
      
      for (j in (i+1):nlen)
      {
        index[1,t]<-dat[i,1]
        index[2,t]<-dat[j,1]
        t=t+1
      }
    }
    
    net = Reduce("rbind", apply(index, 2, pCCAExon, dat, indexlist, 
                                tor = tor))
    if (cutMethod == "value") {
      net_cutoff <- cutoff(net, cutValue)
    }
    if (cutMethod == "number") {
      num = cutNum
      net_cutoff <- cutoff(net, num = num, opt = cutMethod)
    }
    return(net_cutoff)
  }



pCCAExon <-
  function (index, dat, indexlist, tor = 0.3) 
  {
    pos1 = indexlist[[index[1]]]
    pos2 = indexlist[[index[2]]]
    name1 = names(indexlist)[index[1]]
    name2 = names(indexlist)[index[2]]
    X = t(as.matrix(dat[pos1, -c(1, 2)]))
    Y = t(as.matrix(dat[pos2, -c(1, 2)]))
    if (dim(X)[1] != dim(Y)[1]) {
      stop("The sample size of X and Y should be the same!!!")
    }
    n = dim(X)[1]
    p = dim(X)[2]
    q = dim(Y)[2]
    if (sum(X == 0) > p * n/2) {
      return(data.frame(gene1 = name1, gene2 = name2, pcor = 0))
    }
    if (sum(Y == 0) > q * n/2) {
      return(data.frame(gene1 = name1, gene2 = name2, pcor = 0))
    }
    if (p > q) {
      t <- p
      p <- q
      q <- t
    }
    X <- scale(X)
    Y <- scale(Y)
    mcor <- matcor(X, Y)
    if (all(abs(mcor$Xcor) <= tor)) {
      warning("The Xcor is nearly null correlation!!!")
    }
    if (all(abs(mcor$Ycor) <= tor)) {
      warning("The Ycor is nearly null correlation!!!")
    }
    if (all(abs(mcor$XYcor) <= tor)) {
      warning("The XYcor is nearly null correlation!!!")
      return(data.frame(gene1 = name1, gene2 = name2, pcor = 0))
    }
    if (n < q) {
      para <- tryCatch(suppressWarnings(estim.regul(X, Y, plt = TRUE, 
                                                    grid1 = seq(1e-04, 0.2, l = 5), grid2 = seq(0, 0.2, 
                                                                                                l = 5))), error = function(e) NULL)
      if (is.null(para)) {
        return(data.frame(gene1 = name1, gene2 = name2, pcor = 0))
      }
      else {
        ccor <- rcc(X, Y, para$lambda1, para$lambda2)$cor
      }
    }
    else {
      ccor <- tryCatch(cc(X, Y), error = function(e) NULL)
      if (is.null(ccor)) {
        return(data.frame(gene1 = name1, gene2 = name2, pcor = 0))
      }
      else {
        ccor <- ccor$cor
        
      }
    }
    pvec <- dimenReduc(ccor, p, q, n)[, "p"]
    pvec[is.na(pvec)] = 0.9999
    if (pvec[1] >= 0.05) 
      return(data.frame(gene1 = name1, gene2 = name2, pcor = 0))
    idxP <- rep(0, length(pvec))
    pos <- pvec > 0
    idxP[pos] = -log(pvec[pos])
    wccor <- sum(ccor * idxP)/sum(idxP)
    return(data.frame(gene1 = name1, gene2 = name2, pcor = wccor))
  }


cutoff <-
  function (mat, val = 0.05, num = 0, colnum = 3, opt = "value") 
  {
    if (dim(mat)[2] == 5) {
      mat <- cutoffCore(mat, val, num, colnum = 5, opt)
    }
    cutoffCore(mat, val, num, colnum, opt)
  }


cutoffCore <-
  function (mat, val = 0.05, num = 0, colnum = 3, opt = "value") 
  {
    vec <- mat[, colnum]
    vec[is.na(vec)] <- 0
    mat[, colnum] <- vec
    if (opt == "value") {
      if (val > 0.5) {
        val <- 1 - val
      }
      if (any(vec < 0)) {
        a <- quantile(vec, 1 - val/2)
        b <- quantile(vec, val/2)
        pos <- (vec >= a | vec <= b)
      }
      else {
        pos <- (vec > quantile(vec, 1 - val))
      }
    }
    else if (opt == "number") {
      pos <- rank(abs(vec)) > length(vec) - num
    }
    mat[!pos, colnum] <- 0
    mat
  }


drawNet<-function (net) 
{
  library(Rgraphviz)
  net = net[net[, 3] != 0, ]
  df <- data.frame(from = net[, 1], to = net[, 2], weight = net[, 
                                                                3])
  g <- graphBAM(df)
  plot(g)
}

posNet <-function (dat, tor = 0.3, cutMethod = c("value", "number"), cutValue = 0.05, 
                   cutNum = NULL) 
{
  indexlist <- split(1:length(dat[, 1]), dat[, 1])
  index <- combn(1:length(indexlist), 2)
  net = Reduce("rbind", apply(index, 2, pCCAallel, dat, indexlist, 
                              tor = tor))
  if (cutMethod == "value") {
    net_cutoff <- cutoff(net, cutValue)
  }
  if (cutMethod == "number") {
    num = cutNum
    net_cutoff <- cutoff(net, num = num, opt = cutMethod)
  }
  return(net_cutoff)
}


dimenReduc <-
  function (ccor, p, q, n) 
  {
    ev <- ccor^2
    ev2 <- 1 - ev
    m <- n - 3/2 - (p + q)/2
    w <- cbind(NULL)
    for (i in 1:p) {
      w <- cbind(w, prod(ev2[i:p]))
    }
    d1 <- cbind(NULL)
    d2 <- cbind(NULL)
    f <- cbind(NULL)
    for (i in 1:p) {
      s <- sqrt((p^2 * q^2 - 4)/(p^2 + q^2 - 5))
      si <- 1/s
      df1 <- p * q
      d1 <- cbind(d1, df1)
      df2 <- m * s - p * q/2 + 1
      d2 <- cbind(d2, df2)
      r <- (1 - w[i]^si)/w[i]^si
      f <- cbind(f, r * df2/df1)
      p <- p - 1
      q <- q - 1
    }
    pv = tryCatch(suppressWarnings(pf(f, d1, d2, lower.tail = FALSE)))
    dmat <- cbind(t(w), t(f), t(d1), t(d2), t(pv))
    colnames(dmat) <- c("WilksL", "F", "df1", "df2", "p")
    rownames(dmat) <- c(seq(1:length(w)))
    return(dmat)
  }


pCCAallel <-function (index, dat, indexlist, tor = 0.3) 
{
  pCCAExon(index, dat, indexlist, tor = tor)
}

###############functions end###########################################
















