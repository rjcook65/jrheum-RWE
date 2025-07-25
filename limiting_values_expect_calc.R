# =====================================================
# Calculate the limiting values
#
# =====================================================

limiting.val.calc.f <- function(maxvisits, 
                                V2o, # V^o_2(t)
                                alpha0, gamma01, gamma02, # X1 1 -> 0
                                alpha1, gamma11, gamma12, # X1 0 -> 1
                                alpha2, gamma21, gamma22, gamma23, # Z
                                alpha3, gamma31, gamma32, # B
                                alpha4, gamma41, gamma42, gamma43, # A
                                tau, PX2eq1, defclass, ncores) {
  
  # V2o could be a scalar or vector with 2 to 3 entries depending on the working model (A, B, or C)
  
  UZ.f <- function(x, V2o, 
                   alpha2, gamma21, gamma22, gamma23, # Z
                   inwt, inP0mat, inP1mat, inmat, PX2eq1) {
    
    # Establishing the parameters to solve for and initializing the expected observed score functions
    # UZ_beta23, UZ_beta22, UZ_beta23 are the portion of the observed score U2 within the integral
    # Since V^o_2(t) is a vector, U2 is a vector, _beta1 etc., indicate which entry in the vector 
    # UZ_beta1 etc. are associated with
    
    beta23 <- x[1] #B(t) term for working models
    UZ_beta23 <- 0 
    if (V2o >= 2){
      beta22 <- x[2] #X2(t) term for working models
      UZ_beta22 <- 0
      if (V2o == 3){
        beta21 <- x[3] #X1o(t) term for working models
        UZ_beta21 <- 0
      }
    }
    
    # Note that integration will be approximated using Gauss Legendre Integration which 
    # requires summing the full function in the integratal at various chosen time points with corresponding weights
    # In other words, all UZ functions will be a summation at different times. This will be done using a for loop. 
    # Therefore, UZ_beta1 etc. are initialized at 0. 
    
    for (k in 1:length(inwt)) { # length of for loop corresponds to number of weights/nodes (number depends on precision vs computation time)

      if (V2o == 1){ # MODEL A
        # 1
        # print("for loop")
        # print(beta23)
        # print(inP1mat[k, 1:10])
        term0 <- sum( inmat$B*exp(inmat$X1*gamma21 + 0*gamma22 + inmat$B*gamma23)*inP0mat[k, inmat$col] )
        term1 <- sum( inmat$B*exp(inmat$X1*gamma21 + 1*gamma22 + inmat$B*gamma23)*inP1mat[k, inmat$col] )
        ex1 <- alpha2*( (term0*(1 - PX2eq1)) + (term1*PX2eq1) )
        
        # 2
        term0 <- sum( exp(inmat$X1*gamma21 + 0*gamma22 + inmat$B*gamma23)*inP0mat[k, inmat$col] )
        term1 <- sum( exp(inmat$X1*gamma21 + 1*gamma22 + inmat$B*gamma23)*inP1mat[k, inmat$col] )
        ex2 <- alpha2*( (term0*(1 - PX2eq1)) + (term1*PX2eq1) )
        
        # 3
        term0 <- sum( exp(inmat$B*beta23)*inP0mat[k, inmat$col] )
        term1 <- sum( exp(inmat$B*beta23)*inP1mat[k, inmat$col] )
        ex3 <- term0*(1 - PX2eq1) + term1*PX2eq1
        
        # 4
        term0 <- sum( inmat$B*exp(inmat$B*beta23)*inP0mat[k, inmat$col] )
        term1 <- sum( inmat$B*exp(inmat$B*beta23)*inP1mat[k, inmat$col] )
        ex4 <- term0*(1 - PX2eq1) + term1*PX2eq1
        
        U1_beta23 <- ex1 - ((ex2/ex3)*ex4)
        #print(U1_beta23)
        UZ_beta23 <- UZ_beta23 + (U1_beta23*inwt[k])
        
      } else if (V2o == 2){ # MODEL B

        # 1
        term0 <- sum( 0*exp(inmat$X1*gamma21 + 0*gamma22 + inmat$B*gamma23)*inP0mat[k, inmat$col] )
        term1 <- sum( 1*exp(inmat$X1*gamma21 + 1*gamma22 + inmat$B*gamma23)*inP1mat[k, inmat$col] )
        ex1.beta22 <- alpha2*( (term0*(1 - PX2eq1)) + (term1*PX2eq1) ) 
        
        term0 <- sum( inmat$B*exp(inmat$X1*gamma21 + 0*gamma22 + inmat$B*gamma23)*inP0mat[k, inmat$col] )
        term1 <- sum( inmat$B*exp(inmat$X1*gamma21 + 1*gamma22 + inmat$B*gamma23)*inP1mat[k, inmat$col] )
        ex1.beta23 <- alpha2*( (term0*(1 - PX2eq1)) + (term1*PX2eq1) )
        
        
        # 2
        term0 <- sum( exp(inmat$X1*gamma21 + 0*gamma22 + inmat$B*gamma23)*inP0mat[k, inmat$col] )
        term1 <- sum( exp(inmat$X1*gamma21 + 1*gamma22 + inmat$B*gamma23)*inP1mat[k, inmat$col] )
        ex2 <- alpha2*( (term0*(1 - PX2eq1)) + (term1*PX2eq1) )
        
        # 3
        term0 <- sum( exp(0*beta22 + inmat$B*beta23)*inP0mat[k, inmat$col] )
        term1 <- sum( exp(1*beta22 + inmat$B*beta23)*inP1mat[k, inmat$col] )
        ex3 <- term0*(1 - PX2eq1) + term1*PX2eq1
        
        # 4
        term0 <- sum( 0*exp(0*beta22 + inmat$B*beta23)*inP0mat[k, inmat$col] )
        term1 <- sum( 1*exp(1*beta22 + inmat$B*beta23)*inP1mat[k, inmat$col] )
        ex4.beta22 <- term0*(1 - PX2eq1) + term1*PX2eq1
        
        term0 <- sum( inmat$B*exp(0*beta22 + inmat$B*beta23)*inP0mat[k, inmat$col] )
        term1 <- sum( inmat$B*exp(1*beta22 + inmat$B*beta23)*inP1mat[k, inmat$col] )
        ex4.beta23 <- term0*(1 - PX2eq1) + term1*PX2eq1
        
        
        U1_beta22 <- ex1.beta22 - ((ex2/ex3)*ex4.beta22)
        U1_beta23 <- ex1.beta23 - ((ex2/ex3)*ex4.beta23)
        
        UZ_beta22 <- UZ_beta22 + (U1_beta22*inwt[k])
        UZ_beta23 <- UZ_beta23 + (U1_beta23*inwt[k])
        
      }else if (V2o == 3){ # MODEL C
        # 1
        term0 <- sum( inmat$X1_o*exp(inmat$X1*gamma21 + 0*gamma22 + inmat$B*gamma23)*inP0mat[k, inmat$col] )
        term1 <- sum( inmat$X1_o*exp(inmat$X1*gamma21 + 1*gamma22 + inmat$B*gamma23)*inP1mat[k, inmat$col] )
        ex1.beta21 <- alpha2*( (term0*(1 - PX2eq1)) + (term1*PX2eq1) ) #
        
        #term0 <- sum( 0*exp(inmat$X1*gamma21 + 0*gamma22 + inmat$B*gamma23)*inP0mat[k, inmat$col] )
        term1 <- sum( 1*exp(inmat$X1*gamma21 + 1*gamma22 + inmat$B*gamma23)*inP1mat[k, inmat$col] )
        ex1.beta22 <- alpha2*( (term1*PX2eq1) ) #(term0*(1 - PX2eq1)) + 
        
        term0 <- sum( inmat$B*exp(inmat$X1*gamma21 + 0*gamma22 + inmat$B*gamma23)*inP0mat[k, inmat$col] )
        term1 <- sum( inmat$B*exp(inmat$X1*gamma21 + 1*gamma22 + inmat$B*gamma23)*inP1mat[k, inmat$col] )
        ex1.beta23 <- alpha2*( (term0*(1 - PX2eq1)) + (term1*PX2eq1) )
        
        # 2
        term0 <- sum( exp(inmat$X1*gamma21 + 0*gamma22 + inmat$B*gamma23)*inP0mat[k, inmat$col] )
        term1 <- sum( exp(inmat$X1*gamma21 + 1*gamma22 + inmat$B*gamma23)*inP1mat[k, inmat$col] )
        ex2 <- alpha2*( (term0*(1 - PX2eq1)) + (term1*PX2eq1) )
        
        # 3
        term0 <- sum( exp(inmat$X1_o*beta21 + 0*beta22 + inmat$B*beta23)*inP0mat[k, inmat$col] )
        term1 <- sum( exp(inmat$X1_o*beta21 + 1*beta22 + inmat$B*beta23)*inP1mat[k, inmat$col] )
        ex3 <- term0*(1 - PX2eq1) + term1*PX2eq1
        
        # 4
        term0 <- sum( inmat$X1_o*exp(inmat$X1_o*beta21 + 0*beta22 + inmat$B*beta23)*inP0mat[k, inmat$col] )
        term1 <- sum( inmat$X1_o*exp(inmat$X1_o*beta21 + 1*beta22 + inmat$B*beta23)*inP1mat[k, inmat$col] )
        ex4.beta21 <- term1*PX2eq1 + term0*(1 - PX2eq1) 
        
        #term0 <- sum( 0*exp(0*beta22 + inmat$B*beta23)*inP0mat[k, inmat$col] )
        term1 <- sum( 1*exp(inmat$X1_o*beta21 + 1*beta22 + inmat$B*beta23)*inP1mat[k, inmat$col] )
        ex4.beta22 <- term1*PX2eq1 #+ term0*(1 - PX2eq1) 
        
        term0 <- sum( inmat$B*exp(inmat$X1_o*beta21 + 0*beta22 + inmat$B*beta23)*inP0mat[k, inmat$col] )
        term1 <- sum( inmat$B*exp(inmat$X1_o*beta21 + 1*beta22 + inmat$B*beta23)*inP1mat[k, inmat$col] )
        ex4.beta23 <- term0*(1 - PX2eq1) + term1*PX2eq1
        
        
        U1_beta21 <- ex1.beta21 - ((ex2/ex3)*ex4.beta21)
        U1_beta22 <- ex1.beta22 - ((ex2/ex3)*ex4.beta22)
        U1_beta23 <- ex1.beta23 - ((ex2/ex3)*ex4.beta23)
        
        UZ_beta21 <- UZ_beta21 + (U1_beta21*inwt[k])
        UZ_beta22 <- UZ_beta22 + (U1_beta22*inwt[k])
        UZ_beta23 <- UZ_beta23 + (U1_beta23*inwt[k])
      }
      
    }      
    
    if (V2o == 1){
      eq <- UZ_beta23
      print(c(beta23, eq))
      
    }else if (V2o == 2){
      eq <- c(UZ_beta23, UZ_beta22)
      print(c(beta22, beta23, eq))
      
    }else if (V2o == 3){
      eq <- c(UZ_beta23, UZ_beta22, UZ_beta21)
      print(c(beta21, beta22, beta23, eq))
    }
  
    return(eq)
  }
  
  pick <- defclass[defclass$Z == 0 & defclass$A == 0,]
  
  mat <- NULL
  for (k in 0:maxvisits) {
    mat <- rbind(mat, data.frame(col = ((16*k) + pick$states),
                                 B = pick$B,
                                 X1 = pick$X1,
                                 X1_o = pick$X1_o))
  }
  
  q0 <- q_X1ZBA.f(maxvisits,
                  X2=0,
                  stats$alpha0, stats$gamma01, stats$gamma02,
                  stats$alpha1, stats$gamma11, stats$gamma12,
                  stats$alpha2, stats$gamma21, stats$gamma22, stats$gamma23,
                  stats$alpha3, stats$gamma31, stats$gamma32,
                  stats$alpha4, stats$gamma41, stats$gamma42, stats$gamma43)
  
  q1 <- q_X1ZBA.f(maxvisits,
                  X2=1,
                  stats$alpha0, stats$gamma01, stats$gamma02,
                  stats$alpha1, stats$gamma11, stats$gamma12,
                  stats$alpha2, stats$gamma21, stats$gamma22, stats$gamma23,
                  stats$alpha3, stats$gamma31, stats$gamma32,
                  stats$alpha4, stats$gamma41, stats$gamma42, stats$gamma43)
  
  ngauss <- 100 # number of nodes
  gauss <- gauleg.f(ngauss, 0, tau) # num of nodes, x1 (lower limit of integration)
  # x2 (higher limit of integration)
  
  Pmat <- mclapply(1:ngauss, function(kth, inq0, inq1, ingauss) {  
    
    P0 <- MatrixExp(inq0, t=ingauss$nodes[kth])
    P1 <- MatrixExp(inq1, t=ingauss$nodes[kth])
    # print(ingauss$nodes[kth])
    # print(P1[1,1:10])
    return( c(ingauss$nodes[kth], ingauss$weights[kth], P0[1,], P1[1,]) )
    
  }, inq0=q0, inq1=q1, ingauss=gauss, mc.cores=1)
  Pmat <- do.call("rbind", Pmat)
  Pmat <- Pmat[order(Pmat[,1]),]
  NCOL <- ncol(q0)
  
  tt <- Pmat[,1]
  wt <- Pmat[,2]
  P0mat <- Pmat[,(c(1:NCOL)+2)]
  P1mat <- Pmat[,(c(1:NCOL)+2+NCOL)]
  #print(length(Pmat[1,]))
  # print("P1mat")
  # print(P1mat[1,1:10])
  
  if (V2o == 1){ # MODEL A Starting Values
    start = 1
  }else if(V2o == 2){ # MODEL B Starting Values
    start = c(1,1)
  }else if(V2o == 3){ # MODEL C Starting Values
    start = c(1,1,1)
  }
  
  est <- nleqslv(x=start,
                 fn=UZ.f, global="dbldog", xscalm="fixed",
                 control=list(trace=0, btol=0.01, maxit=5000, xtol=1e-08, ftol=1e-08),
                 V2o = V2o,
                 alpha2=alpha2, gamma21=gamma21, gamma22=gamma22, gamma23=gamma23,
                 inwt=wt, inP0mat=P0mat, inP1mat=P1mat, inmat=mat,
                 PX2eq1=PX2eq1)
  print(est)
  beta23 <- est$x[1]
  beta22 <- est$x[2]
  beta21 <- est$x[3]
  
  print( c(beta21, beta22, beta23) )
  
  return( data.frame(beta21, beta22, beta23) )
}