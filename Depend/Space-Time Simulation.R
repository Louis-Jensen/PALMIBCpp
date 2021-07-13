# Given a temporal dynamics generator, adds realistic blinking clusters to any point process.
add.cluster.layer = function(X,generator, s= 0,m1parent,fr) { 
  s = s[[1]]
  if(length(s) == 1) {s = c(s,s)}
  xl = X$x
  yl = X$y
  w = X$window
  
  xo = numeric(0)
  yo = numeric(0)
  wp = numeric(0)
  mks = numeric(0)
  for(i in 1:X$n) {
    np = 0
    while(np < 1) {
      marksim = generator()
      np = length(marksim)
    }
    
    parentmark = ceiling(m1parent()*fr)
    nx = xl[i]+rnorm(np,0, sample(s,np,T) )
    ny = yl[i]+rnorm(np,0, sample(s,np,T) )
    xo = c(xo,nx)
    yo = c(yo,ny)
    wp = c(wp,rep(i, np))
    mks = c(mks, parentmark+marksim)
  }
  Xo = ppp(x = xo, y = yo, marks = mks, window = owin(range(xo), range(yo)))
  attributes(Xo)$parentid = wp
  
  return(Xo)
}

noise_add = function(X, sigs, frames, k = 1, cutr = 0, use.frames = F, keepall = F) {
  if(cutr > 0) {
    X$marks = X$marks+cutr[1]
    sigs = sigs[frames > cutr[1]]
    frames = frames[frames > cutr[1]]
  }
  
  if(use.frames) {
    monka = dbscan::frNN(matrix(frames), k, X$marks)
    
    o = unlist(lapply(monka$id, function(l) good_sample(l,1,T, m = mean(sigs))))
    o[is.na(o)] = sample(length(frames), sum(is.na(o)), T)
    
    s = sigs[o]
    
    dfx = data.frame(X$x, X$y, X$marks)
    dfx$s = s
  } else {
    if(length(sigs) == X$n) {
      dfx = data.frame(X$x, X$y, X$marks, s = sigs)
    } else {
      dfx = data.frame(X$x, X$y, X$marks, s = good_sample(sigs, X$n, T))
    }
    
  }
  
  dfx$sh1 = rnorm(nrow(dfx), 0, dfx$s)
  dfx$sh2 = rnorm(nrow(dfx), 0, dfx$s)
  dfx$X.x = dfx$X.x+dfx$sh1
  dfx$X.y = dfx$X.y+dfx$sh2
  
  if(keepall) {
    sreport = dfx$s
    OX = ppp(x = dfx$X.x, y = dfx$X.y, marks = dfx$X.marks, window = owin(range(dfx$X.x), range(dfx$X.y)))
    OX$marks = OX$marks-cutr[1]
    
    attributes(OX)$of = X$marks
    attributes(OX)$s = sreport
    retl = list(D = OX, kept = "ALL")
  } else {
    kept = which(inside.owin(dfx$X.x, dfx$X.y, X$window))
    
    sreport = dfx$s[kept]
    OX = ppp(x = dfx$X.x, y = dfx$X.y, marks = dfx$X.marks, window = X$window)
    OX$marks = OX$marks-cutr[1]
    
    attributes(OX)$of = X$marks
    attributes(OX)$s = sreport
    retl = list(D = OX, kept = kept)
  }
  
  return(retl)
}

# Simulates blinking dataset.
sim_blinking_pp = function(X_generator, T_generator, ACT_T, framerate, noise = 0, s = 0, frames = 1) {
  X = X_generator()
  W = X$window
  
  Z = add.cluster.layer(X, T_generator, 0, ACT_T, framerate)
  
  cp = attributes(Z)$parent
  if(!(all(s == 0) )) {
    if(length(s) > 1) {
      Z = noise_add(Z, s, frames, use.frames = length(frames) > 1)
    } else {
      Z = noise_add(Z, s, frames, use.frames = F)
    }
  }
  cp = cp[Z$kept]
  Z = Z$D
  
  NZ = Z$n
  NE = round(1/(1/noise - 1)*NZ)
  
  E = runifpoint(NE, W)
  E$marks = runif(NE, 0, max(Z$marks))
  attributes(E)$s = good_sample(s, E$n, T)
  
  O = superimpose(Z, E)
  
  Z$marks = Z$marks/framerate
  E$marks = E$marks/framerate
  O$marks = O$marks/framerate
  
  return(list(X = X, Z = Z, E = E, O = O, Labels = cp))
}


#sim = sim_blinking_pp(function() rpoispp(1, win = square(20)), 
#                pars.to.simulator(c(3,8,0.8), 0, 25, T), 
#                function() rexp(1, 0.01), 25, noise = 0.01, s = 0.05)

# Simulate from fitted model
# fitsim = function() {
#   nr = length(fit$params$rates)
#   X = fit$data
#   cutr = fit$extras$cut
#   fr = fit$extras$framerate
#   
#   if(type == "p") {
#     #parproc = runifpoint(round(fit$params$intensities[4]*area(X$window)), win = X$window)
#     parproc = rpoispp(fit$params$intensities[4], win = X$window)
#   } else if(type == "dimer"){
#     parproc = rThomas(fit$params$intensities[4]/mu, scale = s, mu = mu, win = fit$data$window)
#   }
#   else {
#     parproc = rthin(X, area(X$window)*fit$params$intensities[4]/X$n)
#   }
#   print(parproc$n)
# 
#   np = rpoispp(fit$params$intensities[3], win = X$window)
#   np$marks = runif(np$n, 0, max(X$marks))
#   nproc.o = add.cluster.layer(parproc,pars.to.simulator(fit$params$rates[-nr], tau, fr, T),  0, function() rexp(1, fit$params$rates[nr]),fr )
#   nproc = noise_add(X = nproc.o, sigs = sigs, frames, k = 1, cutr[1], T)
#   onadd = attributes(nproc)$s
#   nproc$marks = nproc$marks/fr
#   nproc = superimpose(nproc, np)
#   attributes(nproc)$s = onadd
# 
#   return(nproc)
# }

