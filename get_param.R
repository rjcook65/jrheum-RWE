# ====================================================
# Derive Simulation Parameters
# 
# 
# ====================================================

# Dataset of combinations of X1, Z, B, A, X1_o
# The order is how they appear in the rows or columns of the qmatrix (col or row 1 => state 1)
states_dataset.f <- function() {
  state_comb <- data.frame(X1 = rep(c(rep(0,2), rep(1,2)),5), 
                           Z = c(rep(0, 8), rep(1,8), rep(0, 4)),
                           B = rep(c(0,1), 10),
                           A = c(rep(0, 16), rep(1, 4)), 
                           X1_o = c(rep(c(0, 0, 1, 1, 1, 1, 0, 0), 2), rep(0, 2), rep(1, 2)))
  state_comb$states <- c(1:nrow(state_comb))
  return(state_comb)
}

# Build an 8 x 8 matrix for a multistate model for X1, B, and X1o
q_X1BXo.f <- function(X2, 
                    alpha0, gamma01, gamma02, # X2 1 -> 0
                    alpha1, gamma11, gamma12, # X2 0 -> 1
                    alpha3, gamma31, gamma32) {# B
  
  qmat <- matrix(0, nrow=8, ncol=8)
  
  qmat[1, 2] <- alpha3*exp(X2*gamma32)
  qmat[3, 4] <- alpha3*exp(gamma31 + X2*gamma32)
  
  qmat[1, 7] <- alpha1*exp(X2*gamma12)
  qmat[2, 8] <- alpha1*exp(gamma11 + X2*gamma12)
  qmat[3, 5] <- alpha0*exp(X2*gamma02)
  qmat[4, 6] <- alpha0*exp(gamma01 + X2*gamma02)
  
  qmat[5, 6] <- alpha3*exp(X2*gamma32)
  qmat[7, 8] <- alpha3*exp(gamma31 + X2*gamma32)
  
  qmat[5, 3] <- alpha1*exp(X2*gamma12)
  qmat[6, 4] <- alpha1*exp(gamma11 + X2*gamma12)
  qmat[7, 1] <- alpha0*exp(X2*gamma02)
  qmat[8, 2] <- alpha0*exp(gamma01 + X2*gamma02)
  
  diag(qmat) <- 0
  diag(qmat) <- (-1)*apply(qmat, 1, sum)
  
  return(qmat)
}

# Build an 16 x 16 matrix for a multistate model for X1, Z, B, and X1o
q_X1ZBXo.f <- function(X2, 
                     alpha0, gamma01, gamma02, # X2 1 -> 0
                     alpha1, gamma11, gamma12, # X2 0 -> 1
                     alpha2, gamma21, gamma22, gamma23, # Z
                     alpha3, gamma31, gamma32) {# B
  
  qmat <- matrix(0, nrow=16, ncol=16)
  
  qmat[1:8, 1:8] <- q_X1BXo.f(X2, 
                            alpha0, gamma01, gamma02, 
                            alpha1, gamma11, gamma12, 
                            alpha3, gamma31, gamma32)
  qmat[1, c(9)] <- alpha2*exp(X2*gamma22) 
  qmat[2, c(10)] <- alpha2*exp(X2*gamma22 + gamma23) 
  qmat[3, c(11)] <- alpha2*exp(gamma21 + X2*gamma22) 
  qmat[4, c(12)] <- alpha2*exp(gamma21 + X2*gamma22 + gamma23)
  
  qmat[5, c(13)] <- alpha2*exp(X2*gamma22) 
  qmat[6, c(14)] <- alpha2*exp(X2*gamma22 + gamma23) 
  qmat[7, c(15)] <- alpha2*exp(gamma21 + X2*gamma22)
  qmat[8, c(16)] <- alpha2*exp(gamma21 + X2*gamma22 + gamma23) 
  diag(qmat) <- 0
  
  diag(qmat) <- (-1)*apply(qmat, 1, sum)
  
  return(qmat)
}


# Build a matrix for a multistate model for X1, Z, B and A, X1o
q_X1ZBA.f <- function(maxvisits,
                    X2, 
                    alpha0, gamma01, gamma02, # X2 1 -> 0
                    alpha1, gamma11, gamma12, # X2 0 -> 1
                    alpha2, gamma21, gamma22, gamma23, # Z
                    alpha3, gamma31, gamma32, #B
                    alpha4, gamma41, gamma42, gamma43) {# A
  
  maxstates <- 16
  
  qq <- matrix(0, nrow=16, ncol=20)
  qq[1:16, 1:16] <- q_X1ZBXo.f(X2, 
                           alpha0, gamma01, gamma02, # X2 1 -> 0
                           alpha1, gamma11, gamma12, # X2 0 -> 1
                           alpha2, gamma21, gamma22, gamma23, # Z
                           alpha3, gamma31, gamma32)
  qq[1, c(17)] <- alpha4*exp(X2*gamma42) 
  qq[2, c(18)] <- alpha4*exp(X2*gamma42) 
  qq[3, c(19)] <- alpha4*exp(gamma41 + X2*gamma42) 
  qq[4, c(20)] <- alpha4*exp(gamma41 + X2*gamma42) 
  
  qq[5, c(17)] <- alpha4*exp(X2*gamma42) 
  qq[6, c(18)] <- alpha4*exp(X2*gamma42) 
  qq[7, c(19)] <- alpha4*exp(gamma41 + X2*gamma42) 
  qq[8, c(20)] <- alpha4*exp(gamma41 + X2*gamma42) 
  
  diag(qq) <- 0
  
  qmat <- diag((maxvisits+1)*maxstates)
  for (k in 0:(maxvisits-1)) {
    row.range <- (maxstates*k) + c(1:16)
    col.range <- (maxstates*k) + c(1:20)
    qmat[row.range, col.range] <- qq
  }
  
  row.range <- (maxstates*maxvisits) + c(1:16)
  col.range <- (maxstates*maxvisits) + c(1:maxstates)
  qmat[row.range, col.range] <- qq[c(1:16), c(1:maxstates)]
  
  diag(qmat) <- (-1)*apply(qmat, 1, sum)
  return(qmat)
}

Etau.f <- function(alpha0, gamma01, gamma02, # X2 1 -> 0
                   alpha1, gamma11, gamma12, # X2 0 -> 1
                   alpha2, gamma21, gamma22, gamma23, # Z
                   alpha3, gamma31, gamma32, #B
                   alpha4, gamma41, gamma42, gamma43,
                   tau, PX2eq1, maxvisits) {

  q0 <- q_X1ZBA.f(maxvisits=maxvisits,
                X2 = 0, 
                alpha0=alpha0, gamma01=gamma01, gamma02=gamma02, # X2 1 -> 0
                alpha1=alpha1, gamma11=gamma11, gamma12=gamma12, # X2 0 -> 1
                alpha2=alpha2, gamma21=gamma21, gamma22=gamma22, gamma23=gamma23, # Z
                alpha3=alpha3, gamma31=gamma31, gamma32=gamma32, #B
                alpha4=alpha4, gamma41=gamma41, gamma42=gamma42, gamma43=gamma43)
  q1 <- q_X1ZBA.f(maxvisits=maxvisits,
                X2 = 1, 
                alpha0=alpha0, gamma01=gamma01, gamma02=gamma02, # X2 1 -> 0
                alpha1=alpha1, gamma11=gamma11, gamma12=gamma12, # X2 0 -> 1
                alpha2=alpha2, gamma21=gamma21, gamma22=gamma22, gamma23=gamma23, # Z
                alpha3=alpha3, gamma31=gamma31, gamma32=gamma32, #B
                alpha4=alpha4, gamma41=gamma41, gamma42=gamma42, gamma43=gamma43)
  
  P0 <- MatrixExp(q0, t=tau)
  P1 <- MatrixExp(q1, t=tau)
  vals <- mclapply(0:maxvisits, function(k, inmaxvisits, inP0, inP1, inPX2eq1) {
    maxstates <- 16
    
    pickA.col <- (maxstates*k) + c(1:maxstates)
    pickD.col <- (maxstates*k) + c(9:16)
    pickC.col <- (maxstates*k) + seq(2,16,by=2)
    
    pD0 <- sum(inP0[1, pickD.col])
    pD1 <- sum(inP1[1, pickD.col])
    pD  <- pD0*(1 - inPX2eq1) + pD1*inPX2eq1
    
    pC0 <- sum(inP0[1, pickC.col])
    pC1 <- sum(inP1[1, pickC.col])
    pC  <- pC0*(1 - inPX2eq1) + pC1*inPX2eq1
    
    
    pE0 <- sum(inP0[1, pickA.col])
    pE1 <- sum(inP1[1, pickA.col])
    pE  <- pE0*(1 - inPX2eq1) + pE1*inPX2eq1
    
    return(data.frame(k, pE, pD, pC))
  }, inmaxvisits=maxvisits, inP0=P0, inP1=P1, inPX2eq1=PX2eq1, mc.cores=1)
  vals <- do.call("rbind", vals)
  
  outdata <- NULL
  outdata$vals  <- vals
  outdata$EAtau <- sum(vals$k*vals$pE)
  outdata$PD    <- sum(vals$pD)
  outdata$PC    <- sum(vals$pC)
  return(outdata) 
}

getparam.f <- function(r10, sum10,
                       EAtau, PBio, PZ, 
                       gamma01, gamma02, # X1 1 -> 0
                       gamma11, gamma12, # X1 0 -> 1
                       gamma21, gamma22, gamma23, # Z
                       gamma31, gamma32, # B
                       gamma41, gamma42, gamma43, # A
                       tau, PX2eq1, maxvisits) {
  
  alpha0 <- (1 + 1/r10)/sum10
  alpha1 <- r10*alpha0
  
  funcZB.f <- function(x, 
                      alpha0, gamma01, gamma02, # X2 1 -> 0
                      alpha1, gamma11, gamma12, # X2 0 -> 1
                      gamma21, gamma22, gamma23,
                      gamma31, gamma32, #B
                      PX2eq1, tau, constZ, constB) {
    y <- numeric(2)
    
    
    alpha2 <- exp(x[1])
    alpha3 <- exp(x[2])
    
    
    q0 <- q_X1ZBXo.f(X2 = 0, 
                  alpha0=alpha0, gamma01=gamma01, gamma02=gamma02, # X2 1 -> 0
                  alpha1=alpha1, gamma11=gamma11, gamma12=gamma12, # X2 0 -> 1
                  alpha2=alpha2, gamma21=gamma21, gamma22=gamma22, gamma23=gamma23, # Z
                  alpha3=alpha3, gamma31=gamma31, gamma32=gamma32) # B
    q1 <- q_X1ZBXo.f(X2 = 1, 
                  alpha0=alpha0, gamma01=gamma01, gamma02=gamma02, # X2 1 -> 0
                  alpha1=alpha1, gamma11=gamma11, gamma12=gamma12, # X2 0 -> 1
                  alpha2=alpha2, gamma21=gamma21, gamma22=gamma22, gamma23=gamma23, # Z
                  alpha3=alpha3, gamma31=gamma31, gamma32=gamma32) # B
    
    P0 <- MatrixExp(q0, t=tau)
    P1 <- MatrixExp(q1, t=tau)
    
    eq0Z <- sum(P0[1,c(9:16)])
    eq1Z <- sum(P1[1,c(9:16)])
    
    eqZ <- ( (eq0Z*(1-PX2eq1)) + (eq1Z*PX2eq1) ) - constZ
    eq0B <- sum(P0[1,seq(2,16,by=2)])
    eq1B <- sum(P1[1,seq(2,16,by=2)])
    
    eqB <- ( (eq0B*(1-PX2eq1)) + (eq1B*PX2eq1) ) - constB
    y[1] <- eqZ
    y[2] <- eqB
    return(y)
  }
  est <- nleqslv(c(0.1, 0.1),funcZB.f, 
             alpha0=alpha0, gamma01=gamma01, gamma02=gamma02,
             alpha1=alpha1, gamma11=gamma11, gamma12=gamma12,
             gamma21=gamma21, gamma22=gamma22, gamma23=gamma23,
             gamma31=gamma31, gamma32=gamma32,
             PX2eq1=PX2eq1, tau=tau, constZ=PZ, constB=PBio, control=list(ftol=10^(-100)))$x
  
  alpha2 <- exp(est[1])
  alpha3 <- exp(est[2])
  
  funcA.f <- function(x, 
                      alpha0, gamma01, gamma02,
                      alpha1, gamma11, gamma12,
                      alpha2, gamma21, gamma22, gamma23, # Z
                      alpha3, gamma31, gamma32,
                      gamma41, gamma42, gamma43,
                      PX2eq1, tau, maxvisits, const) {
    alpha4 <- x
    
    vals <- Etau.f(alpha0=alpha0, gamma01=gamma01, gamma02=gamma02, # X2 1 -> 0
                   alpha1=alpha1, gamma11=gamma11, gamma12=gamma12, # X2 0 -> 1
                   alpha2=alpha2, gamma21=gamma21, gamma22=gamma22, gamma23=gamma23, # Z
                   alpha3=alpha3, gamma31=gamma31, gamma32=gamma32, #B
                   alpha4=alpha4, gamma41=gamma41, gamma42=gamma42, gamma43=gamma43,
                   tau=tau, PX2eq1=PX2eq1, maxvisits=maxvisits)
    
    eq <- vals$EAtau - const
    return(eq)
  }

  alpha4 <- uniroot(funcA.f, interval=c(1e-08,100), extendInt = "yes",
                    alpha0=alpha0, gamma01=gamma01, gamma02=gamma02,
                    alpha1=alpha1, gamma11=gamma11, gamma12=gamma12,
                    alpha2=alpha2, gamma21=gamma21, gamma22=gamma22, gamma23=gamma23, # Z
                    alpha3=alpha3, gamma31=gamma31, gamma32=gamma32,
                    gamma41=gamma41, gamma42=gamma42, gamma43=gamma43,
                    PX2eq1=PX2eq1, tau=tau, maxvisits=maxvisits, const=EAtau)$root
  
  
  outstats <- data.frame(EAtau, PBio, PZ, 
                         alpha0, gamma01, gamma02, # X1 1 -> 0
                         alpha1, gamma11, gamma12, # X1 0 -> 1
                         alpha2, gamma21, gamma22, gamma23, # Z
                         alpha3, gamma31, gamma32, # B
                         alpha4, gamma41, gamma42, gamma43, # A
                         tau, PX2eq1, maxvisits)
  return(outstats)
}
