# Helper functions ----
# --

# Nice printing.
overprint = function(s) {
  cat("\r",s)
}

catprint = function(s) {
  cat(s, "\n")
}

estprint = function(e,n,d = 3) {
  e = round(e, d)
  s = paste(n, "<-", e, sep = " ")
  for(si in s) {catprint(si)}
}

breakl = function() {
  catprint("______________________________________")
}

# Write line to csv file
write.to.file = function(l,file, main = NULL) {
  if(!file.exists(file)) {
    file.create(file)
    if(!is.null(main)) {
      write.table(main,file = file, append = T, col.names = F, row.names = F, quote = F)
    }
  } else {
    write.table(l,file = file, append = T, col.names = F, row.names = F, quote = F)
  }
}

# Multi-start optimization with performance graph.
boptim.try = function(starter.v, fun, n = 1e3, plt = T, silent = F, prt = F, use.sof = T, reset = 100, mean = 0, ...) { 
  if(length(starter.v) == 2) {starter = function() rnorm(starter.v[1],mean, starter.v[2]) } else {starter = starter.v}
  b = Inf
  avs = c()
  sof = rep(0,length(starter()))
  for(i in 1:n) {
    copt = tryCatch(optim(sof+starter(), fun, ...), error = function(e) "ER")
    if(copt[1] == "ER") {
      print(copt)
      n = n+1
      sof = rep(0, length(starter()))
      if(!silent) {
        print(paste("Function returned an error on iteration ",i, sep = ""))
      }
    } else {
      cval = copt$value
      avs[i] = cval
      if(cval < b) {
        params = copt$par
        b = cval
        if(prt) {print(c(round(params,3), log(b) ))}
        if(use.sof) {sof = params; if(any(sof > reset)) {sof = rep(reset, length(params))}}
      }
    }
    
  }
  if(plt) {
    hist(log(avs), breaks = round(min(c(n,100))))
    abline(v = log(b), col = "red", lty = 2, lwd = 2)
  }
  return(list(par = params, val = b))
}

boptim = function(starter, fun, n = 1e3, plt = T, ...) { 
  b = Inf
  avs = numeric(n)
  for(i in 1:n) {
    copt = optim(starter(), fun, ...)
    cval = copt$value
    avs[i] = cval
    if(cval < b) {
      params = copt$par
      b = cval
    }
  }
  if(plt) {
    hist(log(avs), breaks = min(c(n/3,100)))
    abline(v = log(b), col = "red", lty = 2, lwd = 2)
  }
  return(list(par = params, val = b))
}

# Clamps a value
clamp = function(v, l, u) {
  pmin(pmax(v, l), u)
}

# Splits a function into an affine combination of 2 given functions
funsplit = function(O,A,B, fixbc = F) {
  if(!fixbc) {
    M = cbind(A,B,1)
    S = t(M)%*%M
    out = solve(S)%*%t(M)%*%O; out = as.numeric(out)
    names(out) = letters[1:3]
  } else {
    M = cbind(A)
    S = t(M)%*%M
    out = solve(S)%*%t(M)%*%(O-B); out = as.numeric(out)
    names(out) = "a"
    
  }
  return(out)
}
funsplitx = function(O,A,B,C = NULL) {
  if(!is.null(C)) {
    M = cbind(A,B,C)
    S = t(M)%*%M
    out = solve(S)%*%t(M)%*%O; out = as.numeric(out)
    names(out) = letters[1:3]
    return(out)
  } else {
    M = cbind(A,B)
    S = t(M)%*%M
    out = solve(S)%*%t(M)%*%O; out = as.numeric(out)
    names(out) = letters[1:2]
    return(out)
  }
}

# PSF autoconvolution via simulation from collection of uncertainties.
AN = function(r,s1,s2) {
  con = 2*pi*s1*s2*sqrt(1/s1^2+1/s2^2)*sqrt(s1^2+s2^2)
  rter = exp(-r^2/(2*(s1^2+s2^2)))
  return(1/con*rter)
}
get.A = function(r, sigdst, nsim = 1e3) {
  fx = function() {sv = sample(sigdst,2,T); o = AN(r, sv[1],sv[2]); return(o)}
  m = replicate(nsim, fx(), simplify = T)
  o = rowMeans(m)
  return(o)
}
APS = function(r,sigs,ns = 1e3) {
  s1 = sample(sigs, ns, T)
  s2 = sample(sigs, ns, T)
  
  S = cbind(s1,s2)
  
  o = apply(S, 1, function(v) AN(r,v[1],v[2]))
  return(rowMeans(o))
}
good_sample = function(v,size,replace,m) {
  if(length(v) == 0) {return(NA)}
  if(length(v) == 1) {return(v)}
  
  return(sample(v,size, replace))
}
splitf.A = function(A, r, B, sm) {
  o = t(apply( sm,1, function(k) funsplitx(k, A, B , 1) ))
  return(o)
}

splitf.A.noint = function(A, r, B, sm) {
  o = t(apply( sm,1, function(k) funsplitx(k, A, B) ))
  return(o)
}

PCPALM = function(dt,un, noptim = 5) { # PC-PALM method of Sengupta et. al.
  sig = mean(un)
  
  pcf = trimmed_pcf(dt, un)
  rs = pcf$Range
  pcv = pcf$PCF$iso[rs]
  
  r = pcf$PCF$r[rs]
  psf = 1/(4*pi*sig^2)*exp(-r^2/(4*sig^2))
  
  gmodel = function(par, plt = F) {
    par = abs(par)
    a = par[1]
    b = par[2]
    rho = (par[3]/intensity(dt))
    gx = function(r) {return(a*exp(-r/b)+1)}
    
    intt = function(x) integrate(function(r) exp(-r^2/(4*sig^2))*2*pi*besselI(x*r/(2*sig^2),0)*gx(r)*r, 0, max(r)*2)$value
    
    conv = sapply(r, intt)
    
    o = psf*(rho+conv)
    
    if(plt) {
      plot(pcv~r, type = "l", col = "green", lwd = 2)
      points(o~r, type = "l", col = "orange", lwd = 2, lty = 2)
      points(psf*rho+1~r, type = "l", col = "red", lwd = 1.5, lty = 2)
      points(conv*psf~r, type = "l", lwd = 1.5, lty = 2)
    }
    return(o)
  }
  
  f = function(par) {
    par = abs(par)
    v = pcv-gmodel(par)
    return(sum(v^2))
  }
  
  a = boptim.try(c(3,5), f, noptim, plt = F)
  gmodel(a$par[1:3], T)
  nc = abs(a$par[3])
  return(list("Poisson_G" = nc, "Geometric_G" = nc/2+1))
}

dt_to_pp = function(dt) {
  pp = ppp(dt$x, dt$y, marks = dt$frame, window = owin(range(dt$x), range(dt$y)))
  return(pp)
}

BlinkingCSRTest = function(pp, sd, params, framerate, eta = 1, r = seq(0,300,length.out = 150), nsim = 500) {
  delt = 1/framerate
  nblnk = moment.approx(params[2:4], delt)[1]
  
  N = round(pp$n*eta/nblnk)
  
  Lobsi = Lest(pp, r = r, correction = "border")
  Lobsi = Lobsi$border-Lobsi$r
  
  simf = function() sim_blinking_pp(function() runifpoint(N, win = pp$window), 
                                    pars.to.simulator(params[2:4], 0, framerate, T), 
                                    function() rexp(1, params[1]), framerate, noise = 1-eta, s = sd)
  
  fi = function(i) {
    sim.X <- simf()$O
    ff = Lest(sim.X, r = r, correction = "border")
    return(ff$border-ff$r)
  }
  
  #library(future.apply)
  #plan(multisession, workers = 6)
  #oo = future_sapply(1:500, fi, future.seed = T)
  oo = sapply(1:nsim, fi)
  
  cseti = create_curve_set(list(r=r, obs=Lobsi, sim_m=oo))
  resi = global_envelope_test(cseti, type="erl")
  return(resi)
}
