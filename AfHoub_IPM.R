#### African houbara IPM ####

library(R2jags)



### Read in the data ####
load("AfHoubIPM_data.Rdata")


### IPM script ####
sink(paste0(extension,"houbara_IPM.jags"))
cat("
model{

#------------------------------------------------------------------------
# Integrated population model for African houbara
# - Age-structured model with 5 age classes
#     - Juvenile (birds in their first year)
#     - 2nd year adults (birds in their second year)
#     - 3rd year adults
#     - 4th year adults
#     - 5th year adults
# - Post-breeding census, no separation of the sexes

#------------------------------------------------------------------------
##### Some parameters
# N.EWA = median estimates for abundance from HDSM
# N.tot = IPM abundance estimate
# N.c = number of captive-bred birds
# N.w = number of wild birds
# N.juv.c = number of captive-bred juveniles
# N.ad2.c = number of captive-bred 2 year olds (etc)
# beta.age.c = age- and year-specific survival of captive-bred birds
# beta.age.w = age- and year-specific survival of wild birds
# R.ad2.c = productivity for 2-year old captive bred birds, product of breeding parameters (defined within code)

##### Indices
# t in 1:nYears (7 years)
# a in 1:nAgeCats (5 age categories)
# o in 1:2 (2 origin categories, wild/captive bred)
# i in 1:nind (number of individuals, may be nind.c or nind.w for captive-bred/wild)


#------------------------------------------------------------------------
# 1. Specify priors 
#------------------------------------------------------------------------

# Specify ratios for initial population size
# -------------------------------------------------
ratio.origin <- 0.67     # 2:1 cbr:wild
ratio.sex <- 0.5         # 1:1 male:female

# Age ratios from SSD
propAge1_w <- 0.12       # proportion of wild birds that are 1 year old
propAge2_w <- 0.03       # proportion of wild birds that are 2 years old
propAge3_w <- 0.03       # proportion of wild birds that are 3 years old
propAge4_w <- 0.01       # proportion of wild birds that are 4 years old
propAge5_w <- 0.81       # proportion of wild birds that are 5+ years old
propAge1_c <- 0.15       # proportion of captive-bred birds that are 1 year old
propAge2_c <- 0.11       # proportion of captive-bred birds that are 2 years old
propAge3_c <- 0.09       # proportion of captive-bred birds that are 3 years old
propAge4_c <- 0.06       # proportion of captive-bred birds that are 4 years old
propAge5_c <- 0.59       # proportion of captive-bred birds that are 5 years old


# Initial population sizes
# -------------------------------------------------
N.est ~ dnorm(N.EWA[1], tau.y)                             # N.est = letting year 1 abundance vary each iteration

N.juv.c[1] ~ dpois(mean.juv.c[1])                          # CBR juvs
mean.juv.c[1] <- N.est * ratio.origin * propAge1_c
N.ad2.c[1] ~ dpois(mean.ad2.c[1])                          # CBR 2yos
mean.ad2.c[1] <- N.est * ratio.origin * propAge2_c
N.ad3.c[1] ~ dpois(mean.ad3.c[1])                          # CBR 3yos
mean.ad3.c[1] <- N.est * ratio.origin * propAge3_c
N.ad4.c[1] ~ dpois(mean.ad4.c[1])                          # CBR 4yos
mean.ad4.c[1] <- N.est * ratio.origin * propAge4_c
N.ad5.c[1] ~ dpois(mean.ad5.c[1])                          # CBR 5+yos
mean.ad5.c[1] <- N.est * ratio.origin * propAge5_c

N.juv.w[1] ~ dpois(mean.juv.w[1])                          # wild juvs
mean.juv.w[1] <- N.est * (1-ratio.origin) * propAge1_c
N.ad2.w[1] ~ dpois(mean.ad2.w[1])                          # wild 2yos
mean.ad2.w[1] <- N.est * (1-ratio.origin) * propAge2_c
N.ad3.w[1] ~ dpois(mean.ad3.w[1])                          # wild 3yos
mean.ad3.w[1] <- N.est * (1-ratio.origin) * propAge3_c
N.ad4.w[1] ~ dpois(mean.ad4.w[1])                          # wild 4yos
mean.ad4.w[1] <- N.est * (1-ratio.origin) * propAge4_c
N.ad5.w[1] ~ dpois(mean.ad5.w[1])                          # wild 5+yos
mean.ad5.w[1] <- N.est * (1-ratio.origin) * propAge5_c


# Priors for observation error
# ------------------------------------------------
tau.y <- pow(sigma.y, -2)
sigma.y ~ dunif(0, 10000)
sigma2.y <- pow(sigma.y, 2)


# Priors for breeding parameters
# -------------------------------------------------
for(a in 1:nAgeCats){   # age categories
  Br1Pr[a] ~ dunif(0, 1)              # Female breeding probability: captive-bred only
}

for(a in 1:nAgeCats){
  for(y in 1:nYears){
    N.BrAt[a,y] ~ dnorm(1, 0.001)     # Number of breeding attempts
    CS[a,y] ~ dnorm(1, 0.001)         # Clutch size
  }
}

for(o in 1:2){    # origin: wild or captive-bred
  for(a in 1:nAgeCats){   # age categories
    HR[o,a] ~ dunif(0, 1)            # Hatch rate
    JS[o,a] ~ dbeta(7, 4)            # Survival of hatchlings to count
  }
}
   
alpha0 ~ dnorm(0,.001)   # Intercept prior for nest survival model
    
for (i in 1:8){ # age/origin categories
    beta.age.origin[i] ~ dnorm(0, 0.001)  # effect of age/origin on nest survival
} #i


# CBR juvenile survival
for(t in 1:nYears){
  CBR.j.s[t] ~ dbeta(5, 5)   # Prior for annual survival of released birds between release and count
}



# Priors for captive-bred survival parameters
# -------------------------------------------------
f0.c ~ dunif(0.7, 1)                                        # Prior for site fidelity
r0.c ~ dunif(0, 1)                                          # Prior for dead recovery
p0.c ~ dunif(0, 1)                                          # Prior for detection

for(a in 1:max.age){
  for(t in 1:(nYears)){
    beta.age.c[a,t] <- mu.s.c[a] + epsilon.s.c[a,t]         # fixed age and random year effect
    beta.age.c_ps[a,t] <- exp(beta.age.c[a,t])/(1+exp(beta.age.c[a,t]))   # age and year-specific est. on probability scale
    epsilon.s.c[a,t] ~ dnorm(0, tau.s.c[a])
  } #t
  mean.s.c[a] ~ dunif(0, 1)                                 # Prior for mean age-spec. survival
  mu.s.c[a] <- log(mean.s.c[a]/(1-mean.s.c[a]))
  sigma.s.c[a] ~ dunif(0, 20)                               # Prior for age-spec. sd
  tau.s.c[a] <- pow(sigma.s.c[a], -2)
  sigma2.s.c[a] <- pow(sigma.s.c[a], 2)
} #a



# Priors for wild survival parameters
# -------------------------------------------------
f0.w ~ dunif(0.7, 1)                                        # Prior for site fidelity
r0.w ~ dunif(0, 1)                                          # Prior for dead recovery
p0.w ~ dunif(0, 1)                                          # Prior for detection

for(a in 1:max.age){
  for(t in 1:(nYears)){
    beta.age.w[a,t] <- mu.s.w[a] + epsilon.s.w[a,t]         # fixed age and random year effect
    beta.age.w_ps[a,t] <- exp(beta.age.w[a,t])/(1+exp(beta.age.w[a,t]))   # age and year-specific est. on probability scale
    epsilon.s.w[a,t] ~ dnorm(0, tau.s.w[a])
  } #t
  mean.s.w[a] ~ dunif(0, 1)                                 # Prior for mean age-spec. survival
  mu.s.w[a] <- log(mean.s.w[a]/(1-mean.s.w[a]))
  sigma.s.w[a] ~ dunif(0, 10)                               # Prior for age-spec. sd
  tau.s.w[a] <- pow(sigma.s.w[a], -2)
  sigma2.s.w[a] <- pow(sigma.s.w[a], 2)
} #a




#------------------------------------------------------------------------
# 2. Specify likelihoods
#------------------------------------------------------------------------

# 2.2.1 # Likelihood for the survival model: captive bred 
#------------------------------------------------------

for (i in 1:nind.c) {
   for (t in f.surv.c[i]:(nYears)){
     logit(s.c[i,t]) <- beta.age.c[age.c[i,t],t]     # survival with fixed effect of age, random year
     r.c[i,t] <- r0.c                                # dead recovery
     pr.c[i,t] <- p0.c                               # detection probability
     F.c[i,t] <- f0.c                                # site fidelity
   } #t
} #ind


# -------------------------------------------------
# Define state-transition and observation matrices  
for (i in 1:nind.c){
   # Define probabilities of state S(t+1) given S(t)
   for (t in f.surv.c[i]:(nYears)){
      ps.c[1,i,t,1] <- s.c[i,t]*F.c[i,t]
      ps.c[1,i,t,2] <- s.c[i,t]*(1-F.c[i,t])
      ps.c[1,i,t,3] <- (1-s.c[i,t])*r.c[i,t]
      ps.c[1,i,t,4] <- (1-s.c[i,t])*(1-r.c[i,t])
      
      ps.c[2,i,t,1] <- 0
      ps.c[2,i,t,2] <- s.c[i,t]
      ps.c[2,i,t,3] <- (1-s.c[i,t])*r.c[i,t]
      ps.c[2,i,t,4] <- (1-s.c[i,t])*(1-r.c[i,t])
      
      ps.c[3,i,t,1] <- 0
      ps.c[3,i,t,2] <- 0
      ps.c[3,i,t,3] <- 0
      ps.c[3,i,t,4] <- 1
      
      ps.c[4,i,t,1] <- 0
      ps.c[4,i,t,2] <- 0
      ps.c[4,i,t,3] <- 0
      ps.c[4,i,t,4] <- 1
# -------------------------------------------------
      # Define probabilities of O(t) given S(t)
      po.c[1,i,t,1] <- pr.c[i,t]
      po.c[1,i,t,2] <- 0
      po.c[1,i,t,3] <- 1-pr.c[i,t]
      po.c[2,i,t,1] <- 0
      po.c[2,i,t,2] <- 0
      po.c[2,i,t,3] <- 1
      po.c[3,i,t,1] <- 0
      po.c[3,i,t,2] <- 1
      po.c[3,i,t,3] <- 0
      po.c[4,i,t,1] <- 0
      po.c[4,i,t,2] <- 0
      po.c[4,i,t,3] <- 1
      } #t
   } #i
# -------------------------------------------------
# Likelihood 
for (i in 1:nind.c){
   # Define latent state at first capture
   z.surv.c[i,f.surv.c[i]] <- y.surv.c[i,f.surv.c[i]]
   for (t in (f.surv.c[i]+1):(nYears+1)){
      # State process: draw S(t) given S(t-1)
      z.surv.c[i,t] ~ dcat(ps.c[z.surv.c[i,t-1], i, t-1,])
      # Observation process: draw O(t) given S(t)
      y.surv.c[i,t] ~ dcat(po.c[z.surv.c[i,t], i, t-1,])
      } #t
   } #i


# 2.2.2 # Likelihood for the survival model: wild
#------------------------------------------------------

for (i in 1:nind.w) {
   for (t in f.surv.w[i]:(nYears)){
     logit(s.w[i,t]) <- beta.age.w[age.w[i,t],t]     # survival with fixed effect of age, random year
     r.w[i,t] <- r0.w                                # dead recovery
     pr.w[i,t] <- p0.w                               # detection probability
     F.w[i,t] <- f0.w                                # site fidelity
   } #t
} #ind


# -------------------------------------------------
# Define state-transition and observation matrices  
for (i in 1:nind.w){
   # Define probabilities of state S(t+1) given S(t)
   for (t in f.surv.w[i]:(nYears)){
      ps.w[1,i,t,1] <- s.w[i,t]*F.w[i,t]
      ps.w[1,i,t,2] <- s.w[i,t]*(1-F.w[i,t])
      ps.w[1,i,t,3] <- (1-s.w[i,t])*r.w[i,t]
      ps.w[1,i,t,4] <- (1-s.w[i,t])*(1-r.w[i,t])
      
      ps.w[2,i,t,1] <- 0
      ps.w[2,i,t,2] <- s.w[i,t]
      ps.w[2,i,t,3] <- (1-s.w[i,t])*r.w[i,t]
      ps.w[2,i,t,4] <- (1-s.w[i,t])*(1-r.w[i,t])
      
      ps.w[3,i,t,1] <- 0
      ps.w[3,i,t,2] <- 0
      ps.w[3,i,t,3] <- 0
      ps.w[3,i,t,4] <- 1
      
      ps.w[4,i,t,1] <- 0
      ps.w[4,i,t,2] <- 0
      ps.w[4,i,t,3] <- 0
      ps.w[4,i,t,4] <- 1
# -------------------------------------------------
      # Define probabilities of O(t) given S(t)
      po.w[1,i,t,1] <- pr.w[i,t]
      po.w[1,i,t,2] <- 0
      po.w[1,i,t,3] <- 1-pr.w[i,t]
      po.w[2,i,t,1] <- 0
      po.w[2,i,t,2] <- 0
      po.w[2,i,t,3] <- 1
      po.w[3,i,t,1] <- 0
      po.w[3,i,t,2] <- 1
      po.w[3,i,t,3] <- 0
      po.w[4,i,t,1] <- 0
      po.w[4,i,t,2] <- 0
      po.w[4,i,t,3] <- 1
      } #t
   } #i
# -------------------------------------------------
# Likelihood 
for (i in 1:nind.w){
   # Define latent state at first capture
   z.surv.w[i,f.surv.w[i]] <- y.surv.w[i,f.surv.w[i]]
   for (t in (f.surv.w[i]+1):(nYears+1)){
      # State process: draw S(t) given S(t-1)
      z.surv.w[i,t] ~ dcat(ps.w[z.surv.w[i,t-1], i, t-1,])
      # Observation process: draw O(t) given S(t)
      y.surv.w[i,t] ~ dcat(po.w[z.surv.w[i,t], i, t-1,])
      } #t
   } #i



# 2.3 # Likelihood for the breeding model
#--------------------------------------------------
   ### 1) Probability of females breeding ####
  for(a in 1:nAgeCats){   # age categories
          Br1Tr[a] ~ dbin(Br1Pr[a], TrF[a])
  } 
   
   ### 2) Number of attempts per female ####
   for(a in 1:nAgeCats){   # age categories
    for(y in 1:nYears){
      NestAttempts[a,y] ~ dpois(N.BrAt[a,y]*NAfemales[a,y])
      }
   }
   
   ### 3) Clutch size  ####
   for(a in 1:nAgeCats){   # age categories
       for(y in 1:nYears){
          clutchAgg[a,y] ~ dpois(CS[a,y]*clutchNo[a,y])
      }
   }

  ### 4) Hatching rate ####
   for(o in 1:2){    # origin: wild or captive-bred
      for(a in 1:nAgeCats){   # age categories
          hatchlings[o,a] ~ dbin(HR[o,a], eggs[o,a])
      }
   }
   
   ### 5) Logistic expoure model for daily nest survival ####
    for(d in 1:length(ns.succ)){
        eta[d] <- alpha0 + beta.age.origin[ns.age.origin[d]]      # Fixed age/origin effect, random year effect
        expeta[d] <- exp(eta[d])
        s[d] <- expeta[d] / (1+ expeta[d])
        theta[d] <- pow(s[d], ns.interval[d])
        ns.succ[d] ~ dbern(theta[d])
    } #d
    
   for(o in 1:8){    # age/origin category, 8 = 4 age * 2 origin
      daily.ns[o] <- (exp(alpha0+beta.age.origin[o])/(1+exp(alpha0+beta.age.origin[o])))
      annual.ns[o] <- pow((exp(alpha0+beta.age.origin[o])/(1+exp(alpha0+beta.age.origin[o]))), 23)
   }
   
   
   
   ### 7) Final productivity estimation ####
   ## Multiply each breeding parameter for overall productivity estimation
   for(t in 1:nYears){
      R.ad2.c[t] <- Br1Pr[1] * N.BrAt[1,t] * CS[1,t] * HR[1,1] * annual.ns[1] * JS[1,1]
      R.ad3.c[t] <- Br1Pr[2] * N.BrAt[2,t] * CS[2,t] * HR[1,2] * annual.ns[2] * JS[1,2]
      R.ad4.c[t] <- Br1Pr[3] * N.BrAt[3,t] * CS[3,t] * HR[1,3] * annual.ns[3] * JS[1,3]
      R.ad5.c[t] <- Br1Pr[4] * N.BrAt[4,t] * CS[4,t] * HR[1,4] * annual.ns[4] * JS[1,4]
      
      R.ad2.w[t] <- Br1Pr[1] * N.BrAt[1,t] * CS[1,t] * HR[2,1] * annual.ns[1] * JS[2,1]
      R.ad3.w[t] <- Br1Pr[2] * N.BrAt[2,t] * CS[2,t] * HR[2,2] * annual.ns[2] * JS[2,2]
      R.ad4.w[t] <- Br1Pr[3] * N.BrAt[3,t] * CS[3,t] * HR[2,3] * annual.ns[3] * JS[2,3]
      R.ad5.w[t] <- Br1Pr[4] * N.BrAt[4,t] * CS[4,t] * HR[2,4] * annual.ns[4] * JS[2,4]
   }
   


# 2.4 # System process
#--------------------------------------------------

for(t in 1:nYears){
  N.c[t] <- N.juv.c[t] + N.ad2.c[t] + N.ad3.c[t] + N.ad4.c[t] + N.ad5.c[t]    # Number of captive bred birds as sum of each age group
  N.w[t] <- N.juv.w[t] + N.ad2.w[t] + N.ad3.w[t] + N.ad4.w[t] + N.ad5.w[t]    # Number of wild birds as sum of each age group
  N.tot[t] <- N.c[t] + N.w[t]                                                 # Total number of birds as sum of captive-bred and wild groups
}

for(t in 2:nYears){
  mean.juv.c[t] <- Ncbr[t] * CBR.j.s[t]                                                  # captive bred juvs from releases only
  mean.ad2.c[t] <- N.juv.c[t-1] * beta.age.c_ps[1,t]                                     # captive bred 2yos
  mean.ad3.c[t] <- N.ad2.c[t-1] * beta.age.c_ps[2,t]                                     # captive bred 3yos
  mean.ad4.c[t] <- N.ad3.c[t-1] * beta.age.c_ps[3,t]                                     # captive bred 4yos
  mean.ad5.c[t] <- (N.ad4.c[t-1] * beta.age.c_ps[4,t])+                                  # captive bred 5+yos
                   (N.ad5.c[t-1] * beta.age.c_ps[5,t])
  mean.juv.w[t] <- ((((N.juv.w[t-1] * beta.age.w_ps[1,t]) * ratio.sex) * R.ad2.w[t])+    # wild juvs from last yr's wild juvs
                    (((N.ad2.w[t-1] * beta.age.w_ps[2,t]) * ratio.sex) * R.ad3.w[t])+    # wild juvs from last yr's wild 2yos
                    (((N.ad3.w[t-1] * beta.age.w_ps[3,t]) * ratio.sex) * R.ad4.w[t])+    # wild juvs from last yr's wild 3+yos
                    (((N.ad4.w[t-1] * beta.age.w_ps[4,t]) * ratio.sex) * R.ad5.w[t])+    # wild juvs from last yr's wild 2yos
                    (((N.ad5.w[t-1] * beta.age.w_ps[5,t]) * ratio.sex) * R.ad5.w[t])+    # wild juvs from last yr's wild 3+yos
                    (((N.juv.c[t-1] * beta.age.c_ps[1,t]) * ratio.sex) * R.ad2.c[t])+    # wild juvs from last yr's CBR juvs
                    (((N.ad2.c[t-1] * beta.age.c_ps[2,t]) * ratio.sex) * R.ad3.c[t])+    # wild juvs from last yr's CBR 2yos
                    (((N.ad3.c[t-1] * beta.age.c_ps[3,t]) * ratio.sex) * R.ad4.c[t])+    # wild juvs from last yr's CBR 3yos
                    (((N.ad4.c[t-1] * beta.age.c_ps[4,t]) * ratio.sex) * R.ad5.c[t])+    # wild juvs from last yr's CBR 4yos
                    (((N.ad5.c[t-1] * beta.age.c_ps[5,t]) * ratio.sex) * R.ad5.c[t]))    # wild juvs from last yr's CBR 5+yos
  mean.ad2.w[t] <- N.juv.w[t-1] * beta.age.w_ps[1,t]                                     # wild 2yos
  mean.ad3.w[t] <- N.ad2.w[t-1] * beta.age.w_ps[2,t]                                     # wild 3yos
  mean.ad4.w[t] <- N.ad3.w[t-1] * beta.age.w_ps[3,t]                                     # wild 4yos
  mean.ad5.w[t] <- (N.ad4.w[t-1] * beta.age.w_ps[4,t])+                                  # wild 5+yos
                   (N.ad5.w[t-1] * beta.age.w_ps[5,t])
  
  N.juv.c[t] ~ dpois(mean.juv.c[t])T(0,)
  N.ad2.c[t] ~ dpois(mean.ad2.c[t])T(0,)
  N.ad3.c[t] ~ dpois(mean.ad3.c[t])T(0,)
  N.ad4.c[t] ~ dpois(mean.ad4.c[t])T(0,)
  N.ad5.c[t] ~ dpois(mean.ad5.c[t])T(0,)
  N.juv.w[t] ~ dpois(mean.juv.w[t])T(0,)
  N.ad2.w[t] ~ dpois(mean.ad2.w[t])T(0,)
  N.ad3.w[t] ~ dpois(mean.ad3.w[t])T(0,)
  N.ad4.w[t] ~ dpois(mean.ad4.w[t])T(0,)
  N.ad5.w[t] ~ dpois(mean.ad5.w[t])T(0,)
}


#------------------------------------------------------------------------  
# Observation model for the count data (linking counts estimated from the
# HDSM (N.EWA) to N.tot)
#--------------------------------------------------

for (t in 2:nYears){
  N.EWA[t] ~ dnorm(N.tot[t], tau.y)
}

#------------------------------------------------------------------------
# 3. Derived parameters
#------------------------------------------------------------------------

# Population growth rate: total:
for (t in 1:(nYears-1)){
  lambda[t] <- N.tot[t+1]/N.tot[t]
  logla[t] <- log(lambda[t])
}
mlam <- exp((1/(nYears-1))*sum(logla[1:(nYears-1)]))

# Population growth rate: captive bred:
for (t in 1:(nYears-1)){
  lambda.c[t] <- N.c[t+1]/N.c[t]
  logla.c[t] <- log(lambda.c[t])
}
mlam.c <- exp((1/(nYears-1))*sum(logla.c[1:(nYears-1)]))

# Population growth rate: wild:
for (t in 1:(nYears-1)){
  lambda.w[t] <- N.w[t+1]/N.w[t]
  logla.w[t] <- log(lambda.w[t])
}
mlam.w <- exp((1/(nYears-1))*sum(logla.w[1:(nYears-1)]))



}
",fill=TRUE)
sink()


# Initial values
inits <- function(){list(p0.c = runif(1,0,1), f0.c = runif(1,0.9,1), r0.c = runif(1,0,1), z.surv.c = ld.init(rCH.c, f.c),
                         p0.w = runif(1,0,1), f0.w = runif(1,0.9,1), r0.w = runif(1,0,1), z.surv.w = ld.init(rCH.w, f.w))}


params <- c("N.tot","N.est","tau.y","sigma.y","sigma2.y",
            "N.c","N.juv.c","N.ad2.c","N.ad3.c","N.ad4.c","N.ad5.c",
            "N.w","N.juv.w","N.ad2.w","N.ad3.w","N.ad4.w","N.ad5.w",
            "mean.juv.c","mean.ad2.c","mean.ad3.c","mean.ad4.c","mean.ad5.c",
            "mean.juv.w","mean.ad2.w","mean.ad3.w","mean.ad4.w","mean.ad5.w",
            "R.ad2.w","R.ad3.w","R.ad4.w","R.ad5.w","R.ad2.c","R.ad3.c","R.ad4.c","R.ad5.c",
            "Br1Pr","N.BrAt","CS","HR","annual.ns","JS","CBR.j.s",
            "s0.c","beta.age.c_ps","mean.s.c","mu.s.c","sigma2.s.c","epsilon.s.c",
            "s0.w","beta.age.w_ps","mean.s.w","mu.s.w","sigma2.s.w","epsilon.s.w",
            "lambda","logla","mlam","lambda.c","logla.c","mlam.c","lambda.w","logla.w","mlam.w")


ni <- 1000000   ;   nb <- 800000   ;   nt <- 100   ;   nc <- 2
outIPM <- jags(jags.data, inits, params, paste0(extension,"houbara_IPM.jags"),
               n.thin=nt, n.chains=nc, n.burnin=nb,n.iter=ni)

print(outIPM)


