library(plotly)
library(beepr)

library(nleqslv)
library(doParallel)
library(foreach)
library(latex2exp)
n_cores <- detectCores()
cluster <- makeCluster(n_cores-1)
registerDoParallel(cluster)

source('get_param.r')
source('func.gauleg.r')
source('limiting_values_expect_calc.r')

# Other fixed values
maxvisits <- 20
r10 <- 1/3
sum10 <- 0.25
EAtau <- 8 #### Setting 1 = 8, Setting 2 = 4
PBio <- 0.8 
PZ <- 0.4
tau <- 1 
PX2eq1 <- 0.5

# Fixed marker values
gamma02 <- log(0.8)
gamma12 <- log(1.2)

# Fixed failure values
gamma22 <- log(1.5)
gamma23 <- log(0.5)

# Fixed bio values 
gamma32 <- 0

# Fixed visit values 
gamma41 <- 0 #log0.5 Setting 2
gamma42 <- 0
gamma43 <- 0

r1 <- seq(0.5, 1.8, by = 0.2) 
r3 <- seq(0.5, 2, by = 0.2) 

r2 <- 2
resA1 <- foreach(ir3=r3, .combine='rbind') %do% { #row
  
  foreach(ir1=r1, .combine='cbind') %dopar% { #col
    
    library(nleqslv)
    library(msm)
    library(doParallel)
    
    stats <- getparam.f(r10=r10, sum10=sum10,
                        EAtau=EAtau, PBio=PBio, PZ=PZ,
                        gamma01=log(ir1), gamma02=gamma02, # X1 1 -> 0
                        gamma11=-log(ir1), gamma12=gamma12, # X1 0 -> 1
                        gamma21=log(r2), gamma22=gamma22, gamma23=gamma23, # Z
                        gamma31=log(ir3), gamma32=gamma32, # B
                        gamma41=gamma41, gamma42=gamma42, gamma43=gamma43, # A
                        tau = tau, PX2eq1=PX2eq1, maxvisits=maxvisits)
    
    
    est <- limiting.val.calc.f(stats$maxvisits,
                               V2o = 1, # V^o_2(t)
                               stats$alpha0, stats$gamma01, stats$gamma02, # X1 1 -> 0
                               stasts$alpha1, stats$gamma11, stats$gamma12, # X1 0 -> 1
                               stats$alpha2, stats$gamma21, stats$gamma22, stats$gamma23, # Z
                               stats$alpha3, stats$gamma31, stats$gamma32, # B
                               stats$alpha4, stats$gamma41, stats$gamma42, stats$gamma43, # A
                               stats$tau, stats$PX2eq1, defclass=states_dataset.f(), 1)$beta23
    (est - gamma23)/gamma23*100
  }
}
stopCluster(cl = cluster)

###################### B
r2 <- 2
resB1 <- foreach(ir3=r3, .combine='rbind') %do% { #row
  
  foreach(ir1=r1, .combine='cbind') %dopar% { #col
    
    library(nleqslv)
    library(msm)
    library(doParallel)
    
    stats <- getparam.f(r10=r10, sum10=sum10,
                        EAtau=EAtau, PBio=PBio, PZ=PZ,
                        gamma01=log(ir1), gamma02=gamma02, # X1 1 -> 0
                        gamma11=-log(ir1), gamma12=gamma12, # X1 0 -> 1
                        gamma21=log(r2), gamma22=gamma22, gamma23=gamma23, # Z
                        gamma31=log(ir3), gamma32=gamma32, # B
                        gamma41=gamma41, gamma42=gamma42, gamma43=gamma43, # A
                        tau = tau, PX2eq1=PX2eq1, maxvisits=maxvisits)
    
    
    est <- limiting.val.calc.f(stats$maxvisits,
                               V2o = 2, # V^o_2(t)
                               stats$alpha0, stats$gamma01, stats$gamma02, # X1 1 -> 0
                               stasts$alpha1, stats$gamma11, stats$gamma12, # X1 0 -> 1
                               stats$alpha2, stats$gamma21, stats$gamma22, stats$gamma23, # Z
                               stats$alpha3, stats$gamma31, stats$gamma32, # B
                               stats$alpha4, stats$gamma41, stats$gamma42, stats$gamma43, # A
                               stats$tau, stats$PX2eq1, defclass=states_dataset.f(), 1)$beta23
    (est - gamma23)/gamma23*100
  }
}
stopCluster(cl = cluster)


###################### C
r2 <- 2
resC1 <- foreach(ir3=r3, .combine='rbind') %do% { #row
  
  foreach(ir1=r1, .combine='cbind') %dopar% { #col
    
    library(nleqslv)
    library(msm)
    library(doParallel)
    
    stats <- getparam.f(r10=r10, sum10=sum10,
                        EAtau=EAtau, PBio=PBio, PZ=PZ,
                        gamma01=log(ir1), gamma02=gamma02, # X1 1 -> 0
                        gamma11=-log(ir1), gamma12=gamma12, # X1 0 -> 1
                        gamma21=log(r2), gamma22=gamma22, gamma23=gamma23, # Z
                        gamma31=log(ir3), gamma32=gamma32, # B
                        gamma41=gamma41, gamma42=gamma42, gamma43=gamma43, # A
                        tau = tau, PX2eq1=PX2eq1, maxvisits=maxvisits)
    
    
    est <- limiting.val.calc.f(stats$maxvisits,
                               V2o = 3, # V^o_2(t)
                               stats$alpha0, stats$gamma01, stats$gamma02, # X1 1 -> 0
                               stasts$alpha1, stats$gamma11, stats$gamma12, # X1 0 -> 1
                               stats$alpha2, stats$gamma21, stats$gamma22, stats$gamma23, # Z
                               stats$alpha3, stats$gamma31, stats$gamma32, # B
                               stats$alpha4, stats$gamma41, stats$gamma42, stats$gamma43, # A
                               stats$tau, stats$PX2eq1, defclass=states_dataset.f(), 1)$beta23
    (est - gamma23)/gamma23*100
  }
}
