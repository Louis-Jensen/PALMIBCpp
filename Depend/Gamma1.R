bleach_cdf = function(par) {
  r = abs(par)
  
  S = matrix(c(-sum(r[1:2]),r[2],r[3],-r[3]), nrow = 2, ncol = 2, byrow = T)
  S0 = -rowSums(S)
  fx = function(x) as.numeric(1-c(1,0)%*%expm(S*x)%*%c(1,1))
  fxv = Vectorize(fx)
  return(fxv)
}

bleach_qf = function(par, fineness = 1/25) {
  u = ubget(par)*2
  yy = seq(0,u,by = fineness)
  xx = bleach_cdf(par)(yy)
  
  return(approxfun(x = xx, y = yy, rule = 2, yleft = 0))
}

phase_type_interval = function(par, alpha = 0.025) {
  lb = par[1]
  ld = par[2]
  lr = par[3]
  m = (lr+ld)/((lb+ld)*lr-ld*lr)
  v1 = 1/(lb^2*lr)*(ld+lr+(lb*ld+ld^2+ld*lr)/lr)
  v = 2*v1-m^2

  r = m+sqrt(v/alpha)
  return(r)
}

mgf = function(v, lamb) {
  lamb/(lamb-v)
}

mgn = function(v, p) {
  p*exp(v)/(1-(1-p)*exp(v))
}

gam2d.AB = function(par, delt) {
  pt = pars.to.terms(par)
  GM = pt$mf/delt+1/2
  A = 1/GM*(pt$mf2/(2*delt)+ pt$mf+3*delt/8)
  B = 1/(pt$mninf*GM)*( pt$mf/delt+1/2 )*(1/2*(pt$mninf2-pt$mninf)*(pt$mf+pt$mr)+pt$mninf*delt*1/2)
  return(A+B)
}

gam1.asymp.x = function(v,mgff,mgfr,mninf,mninf2,mf,mf2, spec) {
  mgFR = mgfr(-v)*mgff(-v)
  enum = mninf*(mgff(-v)+v*mf-1)+((mgff(-v)-1)/(mgFR-1))^2*mgfr(-v)*(spec(v)-mgFR*mninf+mninf-1  )
  denom = v^2*(mninf*mf2 + mf^2*(mninf2-mninf))
  o = 2*(enum)/denom
  return(o)
}

gam1.asymp = function(v, par) {
  terms = pars.to.terms(abs(par))
  o = gam1.asymp.x(v, terms$mgff, terms$mgfr, terms$mninf, terms$mninf2, terms$mf, terms$mf2, terms$spec)
  return(o)
}

gam1.asymp.d1 = function(v,par) {
  par = abs(par)
  lB = par[1]
  lF = par[2]
  lR = par[3]
  lO = lF+lB
  out = (lB*(lR+v))/(lB * lR+(lB+lF+lR)*v+v^2)
  
  return(out)
}

gam1approx.x = function(v,mgff,mgfr,mninf,mninf2,mf,mf2, spec, delt) {
  eps = 1; teps = 1
  mgFR = mgfr(-v)*mgff(-v)
  a = 2*mninf*(exp(delt*v)*(mgff(-v)*exp(-delt*v*eps)-1) + (exp(delt*v)-1)*(mf/delt+eps))*1/(1-exp(delt*v))^2
  b = mgfr(-v)*(spec(v)-mgFR*mninf+mninf-1)/(1-mgFR)^2
  c = 2*(exp(delt*v*(1+teps)))*(1-mgff(-v)*exp(-delt*v*eps))^2*1/(1-exp(delt*v))^2
  d = mninf2*(mf/delt+eps)^2+mninf*((mf2-mf^2)/delt^2-mf/delt-eps)
  return((a+b*c)/d)
}

gam1approx = function(v, par, delt) {
  terms = pars.to.terms(abs(par))
  o = gam1approx.x(v, terms$mgff,terms$mgfr,terms$mninf,terms$mninf2,terms$mf,terms$mf2, terms$spec, delt)
  o[is.nan(o)] = 1
  return(o)
}
### New FFT-based PCF and CDF functions for gamma 1. Very fast!
ubget = function(par) {
  lb = par[1]
  ld = par[2]
  lr = par[3]
  m = (lr+ld)/((lb+ld)*lr-ld*lr)
  v1 = 1/(lb^2*lr)*(ld+lr+(lb*ld+ld^2+ld*lr)/lr)
  v = 2*v1-m^2
  
  r = m+6.5*sqrt(v)
  return(r)
}

dist1approx = function(par, delt) {
  par = abs(par)
  a = 0
  b = ubget(par)
  NM = (b-a)/delt
  BM = floor(log(NM, 2))
  N = 2^BM
  
  j = 0:(N-1)
  k = j
  dy = (b-a)/N
  y = a + k*dy
  
  phi = (-1)^(-2*a/(b-a)*j)*gam1approx(1i*2*pi/(b-a)*(j-N/2), par, delt)
  
  C = 1/(b-a)*(-1)^as.integer(((a/(b-a)+j/N)*N))
  dv = Re(fft(phi, inverse = T)*C)
  
  dltseq = seq(0,b,by = delt)
  fu = approxfun(x = y, y = cumsum(dv*dy), rule = 2, yleft = 0, yright = 1)
  yo = approxfun(y = fu(dltseq), x = dltseq, rule = 2, yleft = 0, yright = 1)
}

pcf1approx = function(par, delt) {
  par = abs(par)
  a = 0
  b = ubget(par)
  NM = (b-a)/delt
  BM = floor(log(NM, 2))
  N = 2^BM
  
  j = 0:(N-1)
  k = j
  dy = (b-a)/N
  y = a + k*dy
  
  phi = (-1)^(-2*a/(b-a)*j)*gam1approx(1i*2*pi/(b-a)*(j-N/2), par, delt)
  
  C = 1/(b-a)*(-1)^as.integer(((a/(b-a)+j/N)*N))
  dv = Re(fft(phi, inverse = T)*C)
  
  dltseq = seq(0,b,by = delt)
  fu = approxfun(x = y, y = dv*dy, rule = 2, yleft = 0, yright = 0)
  yo = approxfun(y = fu(dltseq)/sum(fu(dltseq)), x = dltseq, rule = 2, yleft = 0, yright = 0)
  
  return(yo)
}

nctol1.helper = function(nc,l2,l3,d) {# parametrize in terms of nc, ld, and lr
  m1 = (l3*d+exp(-l3*d)-1)/(l3*d)
  m2 = (2-2*exp(-l3*d)-2*l3*d+(l3*d)^2)/(l3*d)^2
  
  n = nc
  o = -(1/(2*d^2*n))*(-d-2*d^2*l2+3*d^2*l2*m1-d^2*l2*m2+d*n+ 
                        d^2* l2* n - 
                        d^2* l2* m1* n - sqrt(4 *d^2 *(2 + 4 *d *l2 + 2 *d^2* l2^2 - 
                                                         4* d* l2* m1 - 4* d^2* l2^2* m1 + 2* d^2* l2^2* m1^2)*n + (d + 
                                                                                                                      2*d^2 *l2 - 3 *d^2 *l2 *m1 + d^2* l2 *m2 - d*n - d^2*l2*n + 
                                                                                                                      d^2*l2*m1*n)^2))
  return(o)
}

ncparam = function(par, delt) {
  par = abs(par)
  nc = par[1]
  l2 = par[2]
  l3 = par[3]
  d = delt
  lb = nctol1.helper(nc,l2,l3,d)
  return(c(lb,l2,l3))
}

### New FFT based PCF and CDF functions for gamma 1. Very fast!
meanapprox = function(par, delt, s = 0.01) {
  return(-diff(gam1approx(c(-s,s), par, delt))/(2*s))
}

moment.approx = function(par, delt) {
  pt = pars.to.terms(par)
  
  l3 = abs(par[3])
  
  mur1 = (l3*delt+exp(-l3*delt)-1)/(l3*delt)
  mur2 = (2-2*exp(-l3*delt)-2*l3*delt+(l3*delt)^2)/(l3*delt)^2
  
  m1 = pt$mninf*(pt$mf/delt+1)-(pt$mninf-1)*mur1
  
  m2 = 
    pt$mninf2*(pt$mf/delt+1)^2+pt$mninf*(pt$mf2-pt$mf^2)/delt^2+
    (pt$mninf2+1-2*pt$mninf)*mur1^2+(pt$mninf-1)*(mur2-mur1^2)-
    2*(pt$mf/delt+1)*mur1*(pt$mninf2-pt$mninf)
  
  ov = c(m1,m2,m2/m1-1)
  names(ov) = c("M1", "M2", "nc")
  return(ov)
}

pars.to.terms = function(par) {
  par = abs(par)
  nd = (length(par)-1)/2
  
  rFB = par[1]
  par = par[-1]
  
  out.rates  = par[seq(1,nd*2, by = 2)]
  back.rates = par[seq(2, nd*2, by = 2)]
  
  rOut = sum(out.rates)+rFB
  p = rFB/rOut
  
  mix = out.rates/sum(out.rates)
  
  mninf = 1/p
  mninf2 = (1-p)/p^2+1/p^2
  
  mf = 1/rOut
  mf2 = 2/rOut^2
  
  mr = 0; for(i in 1:nd) {mr = mr+mix[i]*1/back.rates[i]}
  
  mgff = function(v) mgf(v,rOut) 
  mgfr = function(v) {
    o = 0
    for(i in 1:nd) {
      o = o+mix[i]*mgf(v,back.rates[i])
    }
    return(o)
  }
  CDFr = function(v) {
    o = 0
    for(i in 1:nd) {
      o = o+mix[i]*pexp(v,back.rates[i])
    }
    return(o)
  }
  PDFr = function(v) {
    o = 0
    for(i in 1:nd) {
      o = o+mix[i]*dexp(v,back.rates[i])
    }
    return(o)
  }
  spec = function(v) mgn(log(mgff(-v)*mgfr(-v)),p)
  
  return(list(rOut = rOut,mf = mf, mf2 = mf2, mr = mr, mninf = mninf, mninf2 = mninf2, mgff = mgff, mgfr = mgfr, pdfr = PDFr, cdfr = CDFr, spec = spec, p = p, mix = mix, ndark = nd, out.rates = out.rates, back.rates = back.rates))
}
