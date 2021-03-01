#### African houbara IPM ####

library(R2jags)
library(tidyr)
library(dplyr)
library(lubridate)

extension <- "../"    ## if running on mac
# extension <- ""       ## if running on sbatch
setwd("~/Google Drive/Houbara/Models/")  # remove if running on sbatch
source(paste0(extension,'Scripts/IPM functions.R'))



###------------------------------------------------------------------------#
## PREPARE COUNT DATA ####
count.data <- read.csv(paste0(extension,'Data/Count data/sco2001_2019.csv'), as.is = T)

### Subset to relevant years
count.data <- count.data %>% filter(Year > 2010 & Year < 2018)
count.data$ObsPt_Date <- as.Date(count.data$ObsPt_Date, format = "%d/%m/%Y")
count.data$Month <- month(count.data$ObsPt_Date)
count.data <- count.data %>% filter(Month >=8)  ## counts from end Aug-Dec only

# ## Turn NAs into zeros
# count.data$Nb[is.na(count.data$Nb)] <- 0  # specify zero counts as such
# 
# ## Calculate distance to observations:
# for (i in 1:nrow(count.data)) {
#   count.data$Distance[i] <- (trackDistance(count.data$MS_LON[i], count.data$MS_LAT[i],
#                                            count.data$Obj_LON[i], count.data$Obj_LAT[i]))
# }
# 
# 
# ## Truncate distance at 1500m
# hist(count.data$Distance, breaks=100)  # truncate at 1.5km
# count.data$Nb[which(count.data$Distance >1.5)] <- 0  # set observations to 0 if they were out of range
# count.data$Distance[which(count.data$Distance >1.5)] <- NA  # set distance to NA if out of range
# count.data <- count.data[which(is.na(count.data$Distance) | count.data$Distance <= 1.5),]
# hist(count.data$Distance, breaks=100)
# 
# ## Take all observations during the first survey of each site
# count.data <- count.data[order(count.data$Year, count.data$MS_Name, count.data$ObsPt_Date), ]
# first.visits <-
#   count.data %>%   ## get ObsPt_Sta_IDs of first surveys of each site
#   group_by(MS_ID, Year) %>%
#   slice(1)
# first.visits$ObsPt_Sta_ID <- factor(first.visits$ObsPt_Sta_ID)
# 
# ## Subset count data to that from these IDs
# count.data <- count.data[which(count.data$ObsPt_Sta_ID %in% first.visits$ObsPt_Sta_ID),]
# hist(count.data$Distance, breaks=100)
# 
# 
# ## From here just work on the non-zero data ####
# ncap <- pull(count.data %>% group_by(Year) %>%            # number of sites with 0s each year
#                filter(Nb==0) %>% summarise(ncap = n()) %>%
#                ungroup())
# count.data <- count.data %>% filter(Nb > 0)
# 
# ## Add new site ID starting from 1 for each year
# new_ids <- data.frame(MS_ID = count.data$MS_ID, Year = count.data$Year)
# new_ids <- distinct(new_ids)
# new_ids <- new_ids %>% group_by(Year) %>%
#   mutate(new_site_id = sequence(n()))
# count.data <- left_join(count.data, new_ids, by = c("Year", "MS_ID"))  # add new ID back to count data
# 
# 
# ### Prep the data to go in the model ####
# ## organise into yearly observations:
# table(count.data$Year)  ## max = 188, augment by 300
# nobs <- as.vector(table(count.data$Year))                 # number of groups each year
# maxObs <- max(nobs)                                       # max number of groups in any year
# years <- unique(count.data$Year)                          # years of study
# nYears <- length(years)                                   # numbers of years
# 
# ## Augment the data
# nz <- 300                                                 # number of "pseudogroups" added to each year
# y <- matrix(data = NA, ncol = nYears, nrow = maxObs+nz)   # y is matrix of length max groups in a year and width n years
# for(i in 1:nYears) {                                      # fill matrix with obs from each year
#   y.year <- c(count.data$Nb[which(count.data$Year == years[i])], rep(0, nz))
#   y[,i] <- c(y.year, rep(NA, (maxObs+nz)-length(y.year)))  # y is the observations, augmented obs, then NAs for the rest
# }
# groupsize <- y                                            # group size is y
# groupsize[groupsize == 0] <- NA                           # change augmented 0s to NAs
# groupsize <- groupsize-1                                  # input groupsize-1 as data so it fits a Poisson dist.
# y[y>1] <- 1                                               # y is 1s or 0s (or NAs, after augmented data)
# # nind = nobs, not needed in this model without 0s
# d <- y
# for(i in 1:nYears) {                                      # fill matrix with distances
#   d.year <- c(count.data$Distance[which(count.data$Year == years[i])], rep(NA, nz))
#   d[,i] <- c(d.year,rep(NA, (maxObs+nz)-length(d.year)))  # distances, NA for augmented obs, then NAs for the rest
# }
# zst <- y                                                  # Starting values for data augmentation variable
# 
# 
# ### bin distance data
# B <- 1.5                                           # Radius of circle
# delta <- 0.25                                      # width of distance bins
# midpt <- seq(delta/2, B, delta)                    # make mid-points and chop up data
# dclass <- d %/% delta + 1                          # convert distances to cat. distances
# nD <- length(midpt)                                # number of distance intervals
# 
# 
# ### site data
# site <- y
# for(i in 1:nYears) {                               # fill matrix with site obs belong to
#   site.year <- c(count.data$new_site_id[which(count.data$Year == years[i])], rep(NA, nz))
#   site[,i] <- c(site.year,rep(NA, (maxObs+nz)-length(site.year)))  # sites, NA for augmented obs, then NAs for the rest
# }
# nsites <- pull(count.data %>%                      # get n. sites surveyed each year
#                  filter(!is.na(new_site_id)) %>%
#                  group_by(Year) %>%
#                  summarise(n = n_distinct(new_site_id)))
# 
# maxsites <- max(nsites)
# maxobs <- max(nobs)


###------------------------------------------------------------------------#
## PREPARE ABUNDANCE DATA ####
N.EWA <- c(35269, 20146, 13194, 15226, 23437, 17735, 13900)    # median estimates from HDS model



###------------------------------------------------------------------------#
## PREPARE RELEASE DATA ####
Ncbr <- c(10895, 7982, 11131, 13230, 11895, 11427, 6586)


###------------------------------------------------------------------------#
## PREPARE SURVIVAL DATA ####
pivot <- read.csv(paste0(extension,'Data/Pivot files/AfHoubPivot_annual.csv'))
firstLocs <- read.csv(paste0(extension,"Data/Tracking/Houbara_firstOccasions.csv"))
firstLocs$relQuarter <- as.factor(quarter(as.Date(firstLocs$Date)))

## subset to since 2000 (second column is first year 1997, so 16th column is 2011)
pivot <- cbind(pivot[,1],pivot[,16:ncol(pivot)])
colnames(pivot) <- c("Ind_ID",seq(1:8))
pivot <- pivot[rowSums(pivot[,2:(ncol(pivot)-1)]) != 0, ] # remove rows for birds outside this window
pivot <- pivot[which(pivot[,2] != 2),] # remove birds dead at first occasion
pivot <- pivot[rowSums(pivot[,2:ncol(pivot)]==1)!=0,]  # remove birds not seen alive during the study period

# set occasions after death to be 0
g <-vector()
for (i in 1:nrow(pivot)) {
  g[i] <- ifelse(2 %in% pivot[i,2:ncol(pivot)],
                 (which(pivot[i,] == 2))+1, NA)
}

for(i in row_number(pivot[which(2 %in% pivot[1,2:ncol(pivot)])])) {
  if(!is.na(g[i])) {
    for(t in g[i]:ncol(pivot)) {
      pivot[i,t] <- 0
    }
  }
}

# remove the last column, newly added:
pivot <- pivot[,1:(ncol(pivot)-1)]

## add info to pivot and subset
pivot$sex <- firstLocs$Sex[match(pivot$Ind_ID, firstLocs$Ind_ID)]
pivot$age.at.rel <- firstLocs$age_years[match(pivot$Ind_ID, firstLocs$Ind_ID)]
pivot$origin <- firstLocs$Origin[match(pivot$Ind_ID, firstLocs$Ind_ID)]
pivot$relArea <- firstLocs$area[match(pivot$Ind_ID, firstLocs$Ind_ID)]
pivot$relQtr <- firstLocs$relQuarter[match(pivot$Ind_ID, firstLocs$Ind_ID)]



#### CAPTIVE BRED SUBSET ####
pivot.c <- pivot[which(pivot$origin == "captive bred") ,]  # captive bred only
pivot.c <- pivot.c[which(pivot.c$age.at.rel == 0), ]
pivot.c <- pivot.c[which(pivot.c$relArea != "outside"), ]

pivot.c$Ind_ID <- factor(pivot.c$Ind_ID)


## extract variables for covars:
relPer.c <- pivot.c$relQtr
relArea.c <- ifelse(pivot.c$relArea == "hunting",1,2)


### Format pivot for the models ####
head(pivot.c)
CH.c <- as.matrix(pivot.c[,2:(ncol(pivot.c)-5)])
colnames(CH.c) <- 1:ncol(CH.c)

f.c <- apply(CH.c, 1, get.first)  # this is temporary, to be used to subsample. Will be replaced line 201.

## subset to smaller number of individuals
pivot.c$firstYear <- f.c
pivot.c <- pivot.c %>% group_by(firstYear) %>%
  sample_n(80)
pivot.c <- pivot.c[,1:(ncol(pivot.c)-1)]

CH.c <- as.matrix(pivot.c[,2:(ncol(pivot.c)-5)])
colnames(CH.c) <- 1:ncol(CH.c)

f.c <- apply(CH.c, 1, get.first)

rCH.c <- CH.c  # Recoded CH
rCH.c[rCH.c==0] <- 3

# age
age.at.rel.c <- pivot.c$age.at.rel
age.at.rel.c <- age.at.rel.c+1  ## get rid of 0 ages
age.c <- matrix(NA, ncol = dim(CH.c)[2], nrow = dim(CH.c)[1])
for (i in 1:nrow(CH.c)){   ## for each ind
  age.c[i,f.c[i]] <- age.at.rel.c[i]
  for (t in (f.c[i]+1):ncol(CH.c)){  ## for each year since the year after capture
    age.c[i,t] <- (age.at.rel.c[i])+(t-f.c[i])
  }  # t
}  # i
max.age <- 5
age.c[age.c>max.age] <- max.age



#### WILD SUBSET ####
pivot.w <- pivot[which(pivot$origin == "wild") ,]
pivot.w <- pivot.w[which(!is.na(pivot.w$age.at.rel)), ]
pivot.w <- pivot.w[which(pivot.w$relArea != "outside"), ]

pivot.w$Ind_ID <- factor(pivot.w$Ind_ID)


## extract variables for covars:
relPer.w <- pivot.w$relQtr
relArea.w <- ifelse(pivot.w$relArea == "hunting",1,2)


### Format pivot for the models ####
head(pivot.w)
CH.w <- as.matrix(pivot.w[,2:(ncol(pivot.w)-5)])
colnames(CH.w) <- 1:ncol(CH.w)

rCH.w <- CH.w  # Recoded CH
rCH.w[rCH.w==0] <- 3
f.w <- apply(CH.w, 1, get.first)


# age
age.at.rel.w <- pivot.w$age.at.rel
age.at.rel.w <- age.at.rel.w+1  ## get rid of 0 ages
age.w <- matrix(NA, ncol = dim(CH.w)[2], nrow = dim(CH.w)[1])
for (i in 1:nrow(CH.w)){   ## for each ind
  age.w[i,f.w[i]] <- age.at.rel.w[i]
  for (t in (f.w[i]+1):ncol(CH.w)){  ## for each year since the year after capture
    age.w[i,t] <- (age.at.rel.w[i])+(t-f.w[i])
  }  # t
}  # i
max.age <- 5
age.w[age.w>max.age] <- max.age



###------------------------------------------------------------------------#
## PREPARE BREEDING DATA ####
nAgeCats <- max.age-1

### Breeding data: ####
## 1) Number of tracked females in the breeding season each year by each age category ####
houbara_TrF <- read.csv(paste0(extension,"Data/Breeding data/Processed parameters/N_females_avail.csv"))
houbara_TrF <- houbara_TrF %>% filter(Year > 2010 & Year < 2018) %>%
  filter(Origin == "captive bred")
table(houbara_TrF$Age_seasons)

# TrF <- array(data = as.vector(table(houbara_TrF$Origin, houbara_TrF$Age_seasons, houbara_TrF$Year)),
#              c(2,nAgeCats,7))
TrF <- array(data = as.vector(table(houbara_TrF$Age_seasons)), nAgeCats)


### 2) Number of breeding attempts by tracked females each year by each age category ####
houbara_Br1Tr <- houbara_TrF %>% filter(Breed_status == 1) %>%  # subset only to females that bred
  filter(Year > 2010 & Year < 2018) %>%
  filter(Origin == "captive bred")


Br1Tr <- array(data = as.vector(table(houbara_Br1Tr$Age_seasons)), nAgeCats)
# varies by age


## 3) Number of nesting attempts by known females (whether tracked or not) ####
nestAttempts <- read.csv(paste0(extension,"Data/Breeding data/Processed parameters/nests_AllAttempts.csv"))
nestAttempts <- nestAttempts %>% filter(Year > 2010 & Year < 2018)
table(nestAttempts$Age_seasons, nestAttempts$Year)

NestAttempts <- array(data = as.vector(table(nestAttempts$Age_seasons, nestAttempts$Year)),
                      c(nAgeCats, 7))

NestAttemptFemales <- nestAttempts %>% 
  group_by(Year_ID) %>%
  slice(1) %>%
  ungroup()
table(NestAttemptFemales$Age_seasons, NestAttemptFemales$Year)
NAfemales <- array(data = as.vector(table(NestAttemptFemales$Age_seasons, NestAttemptFemales$Year)),
                   c(nAgeCats, 7))
# varies by age and year


## 4) Total number of eggs in clutches each year by each age and origin category ####
eggsInClutches <- read.csv(paste0(extension,"Data/Breeding data/Processed parameters/eggsInClutches.csv"))
eggsInClutches <- eggsInClutches %>% filter(Year > 2010 & Year < 2018)

clutchSize <- eggsInClutches %>% 
  filter(Year > 2010 & Year < 2018) %>% 
  group_by(Year, Age_seasons) %>%
  summarise(eggs = sum(Nb_Item, na.rm=T)) %>%
  ungroup()

clutchAgg <- array(data = clutchSize$eggs,
                   c(nAgeCats, 7))


## 5) Total number of clutches each year by each age and origin category ####
clutchNo <- eggsInClutches %>% 
  filter(Year > 2010 & Year < 2018) %>% 
  group_by(Year, Age_seasons) %>%
  summarise(clutches = n()) %>%
  ungroup()

clutchNo <- array(data = clutchNo$clutches,
                  c(nAgeCats, 7))
# varies by age and year


## 6) Hatch rate ####
hatchRate <- read.csv(paste0(extension,"Data/Breeding data/Processed parameters/hatchRate_nestVisitsOnly.csv"))
hatchRate$Origin[which(hatchRate$Origin == "captive bred")] <- "Captive bred"
hatchRate$Origin[which(hatchRate$Origin == "wild")] <- "Wild"
head(hatchRate)

hatchRate2 <- hatchRate %>% 
  filter(Year > 2002 & Year < 2020) %>% 
  group_by(Age_seasons, Origin) %>%
  summarise(MaxEggs = sum(MaxEggs, na.rm=T), MaxHatch = sum(MaxHatched, na.rm=T)) %>%
  ungroup()

hatchlings <- array(data = hatchRate2$MaxHatch,
                    c(2,nAgeCats))

eggs <- array(data = hatchRate2$MaxEggs,
              c(2,nAgeCats))
# varies by age and origin


## 7) Nest survival data ####
nestSurvival <- read.csv(paste0(extension,"Data/Breeding data/Processed parameters/nestVisitsWithEggs.csv"))
nestSurvival <- nestSurvival %>% filter(Year > 2002 & Year < 2020) %>%
  filter(!is.na(Age_seasons) & !is.na(Origin))

# for calculating sample sizes:
ns_nests <- nestSurvival %>%
  group_by(MS_Name) %>%
  slice(1)
table(ns_nests$Age_seasons, ns_nests$Origin)

nestSurvival$ageOrigin <- case_when(
  nestSurvival$Age_seasons == 1 & nestSurvival$Origin == "Captive bred" ~ 1,
  nestSurvival$Age_seasons == 2 & nestSurvival$Origin == "Captive bred" ~ 2,
  nestSurvival$Age_seasons == 3 & nestSurvival$Origin == "Captive bred" ~ 3,
  nestSurvival$Age_seasons == 4 & nestSurvival$Origin == "Captive bred" ~ 4,
  nestSurvival$Age_seasons == 1 & nestSurvival$Origin == "Wild" ~ 5,
  nestSurvival$Age_seasons == 2 & nestSurvival$Origin == "Wild" ~ 6,
  nestSurvival$Age_seasons == 3 & nestSurvival$Origin == "Wild" ~ 7,
  nestSurvival$Age_seasons == 4 & nestSurvival$Origin == "Wild" ~ 8
  
)

ns.succ <- nestSurvival$success
ns.interval <- nestSurvival$interval
ns.age.origin <- nestSurvival$ageOrigin
# varies by age and origin



### IPM script ####
sink(paste0(extension,"Models/houbara_IPM.jags"))
cat("
model{

#------------------------------------------------------------------------
# Integrated population model for African houbara 
# - Age-structured model with 3 age classes
#     - Juvenile (birds in their first year)
#     - 2nd year adults (birds in their second year)
#     - 3rd year adults (birds in their third year)
# - Post-breeding census, no separation of the sexes
# - Captive-bred demographic rates only for now
#------------------------------------------------------------------------

#------------------------------------------------------------------------
# 1. Specify priors 
#------------------------------------------------------------------------

# Ratios for initial population size
# -------------------------------------------------
ratio.origin <- 0.67     # 2:1 cbr:wild
ratio.sex <- 0.5         # 1:1 male:female
propAge1_w <- 0.12       # age ratios from SSD
propAge2_w <- 0.03
propAge3_w <- 0.03
propAge4_w <- 0.01
propAge5_w <- 0.81
propAge1_c <- 0.15       # age ratios from SSD
propAge2_c <- 0.11
propAge3_c <- 0.09
propAge4_c <- 0.06
propAge5_c <- 0.59


# Initial population sizes
# -------------------------------------------------
N.est ~ dnorm(N.EWA[1], tau.y)

N.juv.c[1] ~ dpois(mean.juv.c[1])                          # CBR juvs
mean.juv.c[1] <- N.est * ratio.origin * propAge1_c
N.ad2.c[1] ~ dpois(mean.ad2.c[1])                          # CBR 2yos
mean.ad2.c[1] <- N.est * ratio.origin * propAge2_c
N.ad3.c[1] ~ dpois(mean.ad3.c[1])                          # CBR 3+yos
mean.ad3.c[1] <- N.est * ratio.origin * propAge3_c
N.ad4.c[1] ~ dpois(mean.ad4.c[1])                          # CBR 4yos
mean.ad4.c[1] <- N.est * ratio.origin * propAge4_c
N.ad5.c[1] ~ dpois(mean.ad5.c[1])                          # CBR 5+yos
mean.ad5.c[1] <- N.est * ratio.origin * propAge5_c

N.juv.w[1] ~ dpois(mean.juv.w[1])                          # wild juvs
mean.juv.w[1] <- N.est * (1-ratio.origin) * propAge1_c
N.ad2.w[1] ~ dpois(mean.ad2.w[1])                          # wild 2yos
mean.ad2.w[1] <- N.est * (1-ratio.origin) * propAge2_c
N.ad3.w[1] ~ dpois(mean.ad3.w[1])                          # wild 3+yos
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
  Br1Pr[a] ~ dunif(0, 1)           # Female breeding probability: captive-bred only
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
    JS[o,a] ~ dbeta(7, 4)
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
  N.c[t] <- N.juv.c[t] + N.ad2.c[t] + N.ad3.c[t] + N.ad4.c[t] + N.ad5.c[t]
  N.w[t] <- N.juv.w[t] + N.ad2.w[t] + N.ad3.w[t] + N.ad4.w[t] + N.ad5.w[t]
  N.tot[t] <- N.c[t] + N.w[t]
}

for(t in 2:nYears){
  mean.juv.c[t] <- Ncbr[t] * CBR.j.s[t]                                               # captive bred juvs from releases only
  mean.ad2.c[t] <- N.juv.c[t-1] * beta.age.c_ps[1,t]                                  # captive bred 2yos
  mean.ad3.c[t] <- N.ad2.c[t-1] * beta.age.c_ps[2,t]                                  # captive bred 3yos
  mean.ad4.c[t] <- N.ad3.c[t-1] * beta.age.c_ps[3,t]                                  # captive bred 4yos
  mean.ad5.c[t] <- (N.ad4.c[t-1] * beta.age.c_ps[4,t])+                                 # captive bred 5+yos
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
  mean.ad2.w[t] <- N.juv.w[t-1] * beta.age.w_ps[1,t]                                  # wild 2yos
  mean.ad3.w[t] <- N.ad2.w[t-1] * beta.age.w_ps[2,t]                                  # wild 3yos
  mean.ad4.w[t] <- N.ad3.w[t-1] * beta.age.w_ps[3,t]                                  # wild 4yos
  mean.ad5.w[t] <- (N.ad4.w[t-1] * beta.age.w_ps[4,t])+                               # wild 5+yos
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


## Bundle data and produce summary
str(bugs.data <- list (nYears=length(Ncbr), max.age=max.age, Ncbr=Ncbr, y.surv.c=rCH.c, f.surv.c=f.c,
                       nind.c=dim(rCH.c)[1], age.c=age.c, y.surv.w=rCH.w, f.surv.w=f.w, nind.w=dim(rCH.w)[1], 
                       age.w=age.w, nAgeCats=nAgeCats, Br1Tr=Br1Tr, TrF=TrF, NAfemales=NAfemales, eggs=eggs,
                       NestAttempts=NestAttempts, clutchAgg=clutchAgg, clutchNo=clutchNo, hatchlings=hatchlings,
                       ns.succ=ns.succ, ns.interval=ns.interval, ns.age.origin=ns.age.origin, N.EWA=N.EWA))

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

# ni <- 1000   ;   nb <- 500   ;   nt <- 2   ;   nc <- 2
ni <- 1000000   ;   nb <- 9000000   ;   nt <- 100   ;   nc <- 2
outIPM <- jags(bugs.data, inits, params, paste0(extension,"Models/houbara_IPM.jags"),
               n.thin=nt, n.chains=nc, n.burnin=nb,n.iter=ni)

print(outIPM)
# traceplot(outIPM)
filename <- paste0(extension,"Rdata/IPM/AfHoubIPM_CORRECTED_withoutHDS_survivalAnnual_breedingAnnual_",Sys.Date(),".Rdata")
save(outIPM, bugs.data, file=filename)


