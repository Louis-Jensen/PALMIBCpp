# Obtain simulator for continuous markov transitions in blinking model.
mark.sim.fun.complete = function(v, prt = F) {
  nd = (length(v)-1)/2
  
  rFB = v[1]
  v = v[-1]
  
  out.rates = v[seq(1,nd*2, by = 2)]
  back.rates = v[seq(2, nd*2, by = 2)]
  
  states = c("F",paste("D", 1:nd, sep = ""), "B")
  
  rmat = matrix(0, nrow = nd+2, ncol = nd+2)
  rmat[1,] = c(0, out.rates,rFB)
  
  rmat[2:(nd+1),1]   = back.rates
  diag(rmat) = -rowSums(rmat)
  
  colnames(rmat) = states
  rownames(rmat) = states
  if(prt) {
    print("Embedded chain transition matrix:")
    print(markovchain::generatorToTransitionMatrix(rmat))
    print("- - -")
    print("Off-diagonal generator matrix:")
    print(rmat-diag(diag(rmat)))
  }
  
  sim.f = new("ctmc", states = states,
              byrow = TRUE, generator = rmat,
              name = "Temporal State Model")
  
  return(sim.f)
}

# Simulate continuous time signal until the bleaching event occurs.
coolsim = function (ctmc, T = 0, include.T0 = TRUE) {
  state <- "F"
  trans <- generatorToTransitionMatrix(ctmc@generator)
  states <- c()
  time <- c()
  if (include.T0 == TRUE) {
    states <- c(states, state)
    time <- c(time, 0)
    wait = c()
  }
  t <- 0
  i <- 1
  while (i <= Inf) {
    idx <- which(ctmc@states == state)
    wt = rexp(1, -ctmc@generator[idx, idx])
    t <- t + wt
    state <- ctmc@states[sample(1:dim(ctmc), 1, prob = trans[idx,])]
    states <- c(states, state)
    time <- c(time, t)
    wait = c(wait,wt)
    if(state == "B") {
      break
    }
    i <- i + 1
  }
  wait = c(wait,0)
  df <- data.frame(states,time,wait)
  names(df) <- c("states", "time","wait")
  return(df)
}

# Simulate discrete, observed temporal signal from a blinking cluster.
discretesim = function(f,tau,frate, frames = T, plt = F ) {
  flen = 1/frate
  
  starttime = runif(1,0,flen)
  
  contsim = coolsim(f)
  
  contsim$time = contsim$time + starttime
  if(plt) {plot(c(F,contsim$states == "F", F)~c(0,contsim$time,1e7), type = "s", xlim = c(0,contsim$time[nrow(contsim)]+1.5*1.5*flen), lwd = 2, col = "black", ylim = c(-0.1,1.2), yaxt = "n", ylab = "", xaxt = "n", xlab = "")}
  
  entrytimes = contsim[which(contsim$states == "F"),]$time
  exittimes  = contsim[which(contsim$states == "F")+1,]$time
  time.spent = exittimes-entrytimes
  
  entryframes = ceiling(entrytimes*frate)
  exitframes  = ceiling(exittimes*frate)
  
  potential.frames = (entryframes[1]:exitframes[length(exitframes)])*flen
  
  overlap  = sapply(potential.frames, function(x) sum(x-entrytimes[entrytimes < x & entrytimes > x-flen]))
  underlap = sapply(potential.frames, function(x) sum(exittimes[exittimes < x & exittimes > x-flen]-x))

  prolap = (overlap+underlap)/flen
  
  wr = which(prolap == 0)
  prolap[wr] = sign(prolap[wr+1]) == -1
  wr1 = which(prolap == 1)
  for(i in wr1) {
    s = T
    j = 0
    while(1) {
      j = j+1
      s = prolap[i-j] == 0
      if(!s) {break}
      prolap[i-j] = 1
    }
  }
  
  prolap[prolap < 0] = 1+prolap[prolap < 0]
  prolap[prolap == 0] = -1
  wf = which(prolap >= tau)
  
  fr.out = potential.frames[wf]
  
  if(plt) {xm = contsim$time[nrow(contsim)]+2*1.5*flen; vs = seq(0,xm, by = flen); nvs = setdiff(vs, fr.out)}
  if(plt) {abline(v = vs, lty = 2, col = rgb(0,0,0,0.25))}
  if(plt) {points(rep(0,length(nvs))~nvs, pch = 19, cex = 1.25, col = "red")}
  if(plt) {points(rep(1,length(fr.out))~c(fr.out), pch = 19, cex = 1.25, col = "red", type = "h", lwd = 2)}
  if(plt) {points(rep(1,length(fr.out))~c(fr.out), pch = 19, cex = 1.25, col = "red", type = "p")}
  #if(plt) {points(rep(1+0.1,length(fr.out))~c(fr.out), pch = as.character(c(1:length(fr.out))), col = "red")}
  
  if(frames) {
    return(fr.out/flen)
  } else {
    return(fr.out)
  }
  
}

# nice plots from discrete blinking cluster simulations.
plotdsim = function(simframes, add = F) {
  simframes = simframes-min(simframes)+1
  mf = max(simframes)
  all.frames = 1:mf
  black.frames = setdiff(all.frames,simframes)
  
  ysignal = rep(1, length(simframes))
  yblack = rep(0, length(black.frames))
  
  yo = c(ysignal,yblack)
  xo = c(simframes,black.frames)
  
  if(add == F) {
    plot(yo~xo, type = "p", ylim = c(-0.1,1.1), pch = 19)
  } else {
    points(yo~xo, type = "p", ylim = c(-0.1,1.1), pch = 19)
  }
}

# Parameters to simulator ----
# ---
pars.to.simulator = function(par, tau, framerate, frames = F) {
  f = mark.sim.fun.complete(par, prt = F)
  gen = function() as.numeric(discretesim(f,tau,framerate, frames = frames))
  return(gen)
}

### Monte-Carlo estimation of gamma1(f) for various functions, f.
# MGF ----
# --
gam1venum = function(simul,v) {
  nu = length(v)
  g = length(simul)
  if(g <= 1) {
    return(c(rep(0,nu),g))
  }
  
  fout = rep(0, nu)
  for(j1 in 1:(g-1)) {
    for(j2 in (j1+1):g) {
      fout = fout + 2*(exp(-v*abs(simul[j1]-simul[j2])) )
    }
  }
  
  return(c(fout,g))
}
est.gam1v = function(v,nsim,simulator) {
  o = replicate(nsim, gam1venum(simulator(),v))
  nu = length(v)
  
  enums = o[1:nu,]
  gnums = o[nu+1,]
  
  return(list("E" = enums,"G" = gnums))
}

# CDF ----
# --
gam1uenum = function(simul,u) {
  nu = length(u)
  g = length(simul)
  if(g <= 1) {
    return(c(rep(0,nu),g))
  }
  
  fout = rep(0, nu)
  for(j1 in 1:(g-1)) {
    for(j2 in (j1+1):g) {
      fout = fout + 2*(abs(simul[j1]-simul[j2]) <= u)
    }
  }
  
  return(c(fout,g))
}

est.gam1u = function(u,nsim,simulator) {
  o = replicate(nsim, gam1uenum(simulator(),u))
  nu = length(u)
  
  enums = o[1:nu,]
  gnums = o[nu+1,]
  
  return(list("E" = enums,"G" = gnums))
}

# Mean
# --
gam1momenum = function(simul) {
  nu = 1
  g = length(simul)
  if(g <= 1) {
    return(c(0,g))
  }
  
  fout = 0
  for(j1 in 1:(g-1)) {
    for(j2 in (j1+1):g) {
      fout = fout + 2*(abs(simul[j1]-simul[j2]))
    }
  }
  return(c(fout,g))
}

est.gam1mom = function(nsim,simulator) {
  o = replicate(nsim, gam1momenum(simulator()))
  nu = 1
  
  enums = o[1:nu,]
  gnums = o[nu+1,]
  
  return(list("E" = enums,"G" = gnums))
}

# General f ----
# --
gam1gen = function(simul,f) {
  nu = 1
  g = length(simul)
  if(g <= 1) {
    return(c(0,g))
  }
  
  fout = 0
  for(j1 in 1:(g-1)) {
    for(j2 in (j1+1):g) {
      fout = fout + f(simul[j1],simul[j2])+f(simul[j2],simul[j1])
    }
  }
  return(c(fout,g))
}

est.gam1gen = function(nsim,simulator,f) {
  o = replicate(nsim, gam1gen(simulator(),f))
  nu = 1
  
  enums = o[1:nu,]
  gnums = o[nu+1,]
  
  return(list("E" = enums,"G" = gnums))
}

# Estimation of reappearance number distribution, G ----
# --
est.gdist = function(nsim,simulator) {
  o = replicate(nsim, gam1uenum(simulator(),1))
  
  gnums = o[2,]
  
  return(gnums)
}

approx.gdist = function(nsim, par, dlt) {
  pt = pars.to.terms(par)
  p = pt$p
  rout = 1/pt$mf
  
  f1s = function() {
    nb = rgeom(1, p)+1
    
    t1 = nb+sum(rexp(nb, rout)/dlt)
    if(nb > 1) {
      t2 = sapply(pt$back.rates, function(r) rexp(nb-1, r))
      if(nb > 2) {t2 = rowSums(t2)}
      
      t2 = sum((t2/dlt <= 1)*(1-t2/dlt))
    } else {
      t2 = 0
    }
    
    return(t1-t2)
  }
  
  o = replicate(nsim, f1s())
  return(o)
}

# Generic summarizing function ----
# --
estfunc = function(par,u,nsim,g2,tau,framerate,combinefun = NULL, type = "u", normf = T) {
  if(length(g2) == 1 & (type == "u" | type == "v")) {g2 = rep(g2,length(u))}
  f1   = mark.sim.fun.complete(par, 1, prt = F)
  
  if(normf) {
    gen1 = function() discretesim(f1,tau,framerate, frames = F)
  } else {
    gen1 = function() discretesim(f1,tau,framerate)
  }
  
  if(type == "u") {
    g1list = est.gam1u(u,nsim,gen1)
  } else if(type == "v"){
    g1list = est.gam1v(u,nsim,gen1)
  } else {
    g1list = est.gam1mom(nsim,gen1)
    g1list$E = matrix(g1list$E, nrow = 1)
  }
  
  g1u = rowMeans(g1list$E)
  gs = g1list$G
  
  if(any(is.null(combinefun))) {
    return(g1list)
  } else {
    out = combinefun(g1u,gs,g2)
    return(out)
  }
}



