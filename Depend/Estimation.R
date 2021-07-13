# Estimation functions ----
# --
# Fast computation of pair correlation function, based on the K function.
pcf.c = function(X,fast = T, correction = "iso", method = "c",...) {
  if(fast) {
    p = pcf.fv(Kest(X,correction = correction,...), method = method)
    p$iso = p$pcf
  } else {
    p = pcf(X,...)
  }
  return(p)
}

# PCF estimate with informative range of r-values.
trimmed_pcf = function(X, STDS, ...) {
  dnb = nndist(X); dnb = dnb[dnb > 0]
  
  rlow = mean(STDS)/4
  rmax = mean(STDS)*10
  rinc = (rmax-rlow)/5e2 # 500 equally spaced query points. Should be enough for most cases.
  
  pcf = pcf.c(unmark(X), r = seq(0,rmax, by = rinc), ...)
  range = which(pcf$r > rlow)
  
  return(list("PCF" = pcf, "Range" = range))
}

# Markcorrelation function based on mark correlation integral.
Kmark.c = function(X,f, method = "c",...) {
  p = pcf.fv(Kmark(X,f, ...), method = method)
  p$iso = p$pcf
  return(p)
}

g2.estim = function(M,alpha, R = max(M), plt = F, nsim = 1e5) {
  MM = sort(M)
  uM = unique(MM)
  oecdf = ecdf(MM)(uM)
  
  uM  = c(0,uM); oecdf = c(0,oecdf)
  
  necdf = uM/R
  g2ecdf = (oecdf-(1-alpha)*necdf)/alpha
  
  fq.o = approxfun(y = uM, x = oecdf)
  fq = approxfun(y = uM, x = g2ecdf)
  
  ss.g2o = abs(fq.o(runif(nsim,0,1))-fq.o(runif(nsim,0,1)))
  og.g2o = sapply(uM, function(u) mean(ss.g2o <= u))
  
  ss.g2 = abs(fq(runif(nsim,0,1))-fq(runif(nsim,0,1)))
  og.g2 = sapply(uM, function(u) mean(ss.g2 <= u))
  
  og.ne = function(x) 1-(1-x/R)^2
  
  ss.mix = abs(fq(runif(nsim,0,1))-runif(nsim, 0, R))
  og.mix = sapply(uM, function(u) mean(ss.mix <= u))
  
  fo.g2o = approxfun(x = uM, y = og.g2o, rule = 2)
  fo.g2  = approxfun(x = uM, y = og.g2, rule = 2)
  fo.ne  = og.ne
  fo.mix = approxfun(x = uM, y = og.mix, rule = 2)
  
  if(plt) {
    foo = abs(sample(M,1e5, T)-sample(M,1e5, T))
    s0 = sapply(uM, function(u) mean(foo <= u))
    plot(fo(uM)~uM, type = "l", col = "green", lwd = 1.5, xlab = "u", ylab = "Fu")
    points(s0~uM, type = "l", col = "red", lwd = 1.5)
    legend(x = "bottomright", legend = c("gamma2", "gamma2o"), col = c("green", "red"), lty = 1)
  }
  return(list(g2o = fo.g2o, g2 = fo.g2, g2e = fo.ne, g2ez = fo.mix)) 
}

# Experimental function to estimate bias in localization uncertainties.
sigma.bias.estim = function(summat, sigma, r, r.range, alpha, into, sseq = seq(0.25, 2.75, by = 0.25/10)) {
  bst.s = 1
  bst.er = Inf
  for(sc in sseq) {
    At = get.A(r[r.range], sigma*sc, 1e3)
    W = apply(summat, 1, function(s) funsplit(s, At*alpha/into, 0, fixbc = T))
    er = mean((summat-W%*%t(At)*alpha/into)^2)
    if(er < bst.er) {
      bst.er = er
      bst.s = sc
    }
  }
  return(bst.s)
}

correct.inrate = function(r0,b) {
  f = function(l) (exp(l*b)-l*b-1)/(l*(exp(l*b)-1))-1/r0
  xd = uniroot(f, c(r0/10, r0))
  return(xd$root)
}

# Estimate model parameters.
estimate.IBCpp = function(X, framerate, sig = attributes(X)$s, eta = 1, spatial_lag = 1, nopt = 50, ...) {
  if(!all(is.integer(X$marks))) {catprint("Please provide point pattern with frame-marks."); stop("Non-integer marks detected.")}

  X$marks = X$marks/framerate
  delt = 1/framerate
  
  breakl(); catprint("Estimating summary functions.")
  #######
  largeR = mean(sig)*spatial_lag
  closep = closepairs(X, largeR)

  ii = closep$i
  ij = closep$j
  ds = closep$d

  tlag = abs(X$marks[ii]-X$marks[ij])
  rs = seq(0,largeR, length.out = 500)[-1]
  NR = length(rs)
  us = seq(0,max(X$marks), by = delt)
  NU = length(us)

  m = matrix(0, nrow = NR, ncol = NU)
  i = 0
  j = 0

  for(rr in rs) {
    overprint(paste(NR-i, "       ", sep = ""))
    i = i+1

    tsub = tlag[ds <= rr]
    if(length(tsub) == 0) {
      ee = Vectorize(function(x) return(0))
    } else {
      ee = ecdf(tsub)
    }
    m[i,] = ee(us)*length(tsub)/intensity(X)/(X$n-1)
  }
  
  g2s = g2.estim(X$marks, eta)
  A = get.A(rs, sig, 1e4)
  g22 = g2s$g2(us)
  g2o = g2s$g2o(us)
  KM = m[,length(us)]
  
  # Note: may want to add some smoothing to ss in future.
  ss = apply(m, 2, function(x) diff(x)/rs[-1]/(2*pi)*1/(diff(rs)[1]))
  
  FM = ss[,length(us)]
  for(i in 1:ncol(m)) {
    ss[,i] = ss[,i]-(FM-1)*g22[i]-g2o[i]
  }
  ok = ss
  
  AH = A[-1]
  ws = apply(ok,2, function(x) sum(x*AH)/sum(AH^2))/eta
  
  
  ##
  catprint("")
  breakl(); catprint("Estimating rate-parameters.")

  d1 = ws*intensity(X)

  g22[g22 == 0] = min(g22[g22 != 0])
  weights = clamp((ws/(g22)), 0, Inf)^2
  weights = weights/sum(weights)
  
  ff = function(par) {
    nc = moment.approx(par, delt)[3]
    gam1guess = dist1approx(par,delt)(us)
    guess = (gam1guess-g22)*nc
    o = sum(weights*(guess-d1)^2)
    return(o)
  }

  f = function(par) ff(ncparam(par,delt))
  yes = boptim.try(c(3,3), f, nopt, plt = F)
  opar = ncparam(abs(yes$par), delt)
  lr = 1/((mean(X$marks)-(1-eta)*sum(range(X$marks))/2)/eta-gam2d.AB(opar, delt))

  ratepar = c(lr, opar)
  cat("-")
  overprint(ratepar)
  cat("-")
  names(ratepar) = c("r_F", "r_B", "r_D", "r_R")
  return(list(par = ratepar))
}
