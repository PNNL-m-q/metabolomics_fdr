# This script compares the capabilities of two empirical bayes (EB) methods,
# one simple (aka the GMM) and the other complex (aka the hierarchical model), in their ability to estimate the FDR. 
# The complex EB model is based on the work of Jeong et. al (2011). We use an adapted version of
# their original code (https://github.com/jjs3098/CNU-Bioinformatics-Lab)

# Required packages: 
# require(truncnorm)
# require(tidyverse)
# require(rlist)
# require(gridExtra)
# require(parallel)
# require(ggpubr)

# Assign the misc directory path to misc_folder_path
misc_folder_path <- "/misc/folder/path/"

### Function(s) -------------------------------------------------------------

# Jeong et al functions (modified)
source(paste0(misc_folder_path, "MetPC_modified.R"))

# Data simulator (truncated normal / truncated normal)
datsim <- function(refdat, true_params, num_match){
  
  # Extract true parameters
  rho     <- true_params$rho  # probability that a given ref metabolite is within the sample
  eta0    <- true_params$eta0  # probability a given reference metabolite is matched (assuming bj == 0)
  beta0   <- true_params$beta0    # logisitic regression parameter to determine prob a reference metabolite is matched (assuming bj != 0)
  beta1   <- true_params$beta1   # logisitic regression parameter to determine prob a reference metabolite is matched (assuming bj != 0)
  beta2   <- true_params$beta2  # logisitic regression parameter to determine prob a reference metabolite is matched (assuming bj != 0)
  eta1    <- true_params$eta1  # probability a given reference metabolite is matched (assuming bjs == 1)
  alpha0  <- true_params$alpha0  # logisitic regression parameter to determine prob a reference metabolite is matched (assuming bjs != 1)
  alpha1  <- true_params$alpha1 # logisitic regression parameter to determine prob a reference metabolite is matched (assuming bjs != 1)
  alpha2  <- true_params$alpha2  # logisitic regression parameter to determine prob a reference metabolite is matched (assuming bjs != 1)
  tau     <- true_params$tau  # probability that a given match is correct (i.e. true positive)
  muT     <- true_params$muT    # mean of true positive distribution
  varT    <- true_params$varT  # variance of true positive distribution
  muF     <- true_params$muF   # mean of false positive distribution
  varF    <- true_params$varF # variance of false positive distribution
  
  nRef <- nrow(refdat)
  simY <- rbinom(nRef, 1, rho)
  
  refdat <- refdat %>% mutate(simY = simY)
  
  # Next, we simulate the propensity of each reference metabolite to match given the simulated Y
  # and assumed competition scores (computed based on real reference data)
  simZ <- vector("numeric", length = nrow(refdat))
  for(i in 1:nrow(refdat)){
    if(refdat$simY[i] == 0){
      if(refdat$bj[i] == 0){
        simZ[i] <- rbinom(1, 1, eta0)
        
      } else if(refdat$bj[i] != 0){
        gmma <- 1 - 1/(1 + exp(beta0 + beta1*refdat$bj[i] + beta2*refdat$bj[i]^2))
        simZ[i] <- rbinom(1, 1, gmma)
        
      }
    } else if(refdat$simY[i] == 1){
      if(refdat$bjs[i] == 1){
        simZ[i] <- rbinom(1, 1, eta1)
        
      } else if(refdat$bjs[i] != 1){
        lmda <- 1 - 1/(1 + exp(alpha0 + alpha1*refdat$bjs[i] + alpha2*refdat$bjs[i]^2))
        simZ[i] <- rbinom(1, 1, lmda)
      }
    }
  }
  rm(i, gmma, lmda)
  
  refdat <- refdat %>% mutate(simZ = simZ)
  
  # Next, for each refmet with simZ == 1, we assume there were five 
  # corresponding matches. A constant probability of the correctness of
  # match (tau) is assumed across all matches where simY == 1. When 
  # simY == 0 and simZ == 1, the probability of the correct match is 0.
  # Corresponding similarity scores are also simulated from a mixture of
  # truncated normal distributions
  hitdat <- vector("list", length = sum(refdat$simZ == 1))
  names(hitdat) <- refdat$Ref[refdat$simZ == 1]
  for(ref in names(hitdat)){
    if(refdat$simY[refdat$Ref == ref] == 0){
      simW <- rbinom(num_match, 1, 0)
      simS <- vector("numeric", length = num_match)
      for(j in 1:length(simW)){
        # Simulates exclusively from false positive distribution since
        # simW is always 0 in this case
        simS[j] <- simW[j]*rtruncnorm(1, a = 0, b = 90, mean = muT, sd = sqrt(varT)) + 
          (1-simW[j])*rtruncnorm(1, a = 0, b = 90, mean = muF, sd = sqrt(varF))
      }
      hitdat[[ref]] <- data.frame(Ref  = rep(ref, num_match),
                                  simW = simW,
                                  simS = simS)
    } else if(refdat$simY[refdat$Ref == ref] == 1){
      simW <- rbinom(num_match, 1, tau)
      simS <- vector("numeric", length = num_match)
      for(j in 1:length(simW)){
        simS[j] <- simW[j]*rtruncnorm(1, a = 0, b = 90, mean = muT, sd = sqrt(varT)) + 
          (1-simW[j])*rtruncnorm(1, a = 0, b = 90, mean = muF, sd = sqrt(varF))
      }
      hitdat[[ref]] <- data.frame(Ref  = rep(ref, num_match),
                                  simW = simW,
                                  simS = simS)
    }
  }
  rm(ref, simW, simS)
  
  hitdat_all <- Reduce("rbind.data.frame", hitdat)
  
  return(list(refdat      = refdat,
              hitdat      = hitdat_all,
              true_params = true_params))
}

# Function to draw from a right-truncated gamma distribution
rgammat <- function(n, rtrunc, shape, scale = 1) {
  
  F.b <- pgamma(rtrunc, shape = shape, scale = scale)
  
  u <- runif(n, min = 0, max = F.b)
  
  qgamma(u, shape = shape, scale = scale)
  
}

# Data simulator (gamma / mirrored gamma)
datsim_nonorm <- function(refdat, true_params, num_match){
  
  # Extract true parameters
  rho     <- true_params$rho  # probability that a given ref metabolite is within the sample
  eta0    <- true_params$eta0  # probability a given reference metabolite is matched (assuming bj == 0)
  beta0   <- true_params$beta0    # logisitic regression parameter to determine prob a reference metabolite is matched (assuming bj != 0)
  beta1   <- true_params$beta1   # logisitic regression parameter to determine prob a reference metabolite is matched (assuming bj != 0)
  beta2   <- true_params$beta2  # logisitic regression parameter to determine prob a reference metabolite is matched (assuming bj != 0)
  eta1    <- true_params$eta1  # probability a given reference metabolite is matched (assuming bjs == 1)
  alpha0  <- true_params$alpha0  # logisitic regression parameter to determine prob a reference metabolite is matched (assuming bjs != 1)
  alpha1  <- true_params$alpha1 # logisitic regression parameter to determine prob a reference metabolite is matched (assuming bjs != 1)
  alpha2  <- true_params$alpha2  # logisitic regression parameter to determine prob a reference metabolite is matched (assuming bjs != 1)
  tau     <- true_params$tau  # probability that a given match is correct (i.e. true positive)
  shapeT  <- true_params$shapeT  # shape parameter of true positive distribution
  scaleT  <- true_params$scaleT  # scale parameter of true positive distribution
  shapeF  <- true_params$shapeF  # shape parameter of false positive distribution
  scaleF  <- true_params$scaleF  # scale parameter of false positive distribution
  
  nRef <- nrow(refdat)
  simY <- rbinom(nRef, 1, rho)
  
  refdat <- refdat %>% mutate(simY = simY)
  
  # Next, we simulate the propensity of each reference metabolite to match given the simulated Y
  # and assumed competition scores (computed based on real reference data)
  simZ <- vector("numeric", length = nrow(refdat))
  for(i in 1:nrow(refdat)){
    if(refdat$simY[i] == 0){
      if(refdat$bj[i] == 0){
        simZ[i] <- rbinom(1, 1, eta0)
        
      } else if(refdat$bj[i] != 0){
        gmma <- 1 - 1/(1 + exp(beta0 + beta1*refdat$bj[i] + beta2*refdat$bj[i]^2))
        simZ[i] <- rbinom(1, 1, gmma)
        
      }
    } else if(refdat$simY[i] == 1){
      if(refdat$bjs[i] == 1){
        simZ[i] <- rbinom(1, 1, eta1)
        
      } else if(refdat$bjs[i] != 1){
        lmda <- 1 - 1/(1 + exp(alpha0 + alpha1*refdat$bjs[i] + alpha2*refdat$bjs[i]^2))
        simZ[i] <- rbinom(1, 1, lmda)
      }
    }
  }
  rm(i, gmma, lmda)
  
  refdat <- refdat %>% mutate(simZ = simZ)
  
  # Next, for each refmet with simZ == 1, we assume there were five 
  # corresponding matches. A constant probability of the correctness of
  # match (tau) is assumed across all matches where simY == 1. When 
  # simY == 0 and simZ == 1, the probability of the correct match is 0.
  # Corresponding similarity scores are simulated from a mixture of
  # a gamma and mirrored gamma distribution. 
  hitdat <- vector("list", length = sum(refdat$simZ == 1))
  names(hitdat) <- refdat$Ref[refdat$simZ == 1]
  for(ref in names(hitdat)){
    if(refdat$simY[refdat$Ref == ref] == 0){
      simW <- rbinom(num_match, 1, 0)
      simS <- vector("numeric", length = num_match)
      for(j in 1:length(simW)){
        # Simulates exclusively from false positive distribution since
        # simW is always 0 in this case
        # Note that we model the false positive distribution as mirrored gamma, shifted by the max score (90)
        # We also truncate each gamma distribution at 90
        simS[j] <- simW[j]*rgammat(1, rtrunc = 90, shape = shapeT, scale = scaleT) + 
          (1-simW[j])*(90-rgammat(1, rtrunc = 90, shape = shapeF, scale = scaleF))
      }
      hitdat[[ref]] <- data.frame(Ref  = rep(ref, num_match),
                                  simW = simW,
                                  simS = simS)
    } else if(refdat$simY[refdat$Ref == ref] == 1){
      simW <- rbinom(num_match, 1, tau)
      simS <- vector("numeric", length = num_match)
      for(j in 1:length(simW)){
        simS[j] <- simW[j]*rgammat(1, rtrunc = 90, shape = shapeT, scale = scaleT) + 
          (1-simW[j])*(90-rgammat(1, rtrunc = 90, shape = shapeF, scale = scaleF))
      }
      hitdat[[ref]] <- data.frame(Ref  = rep(ref, num_match),
                                  simW = simW,
                                  simS = simS)
    }
  }
  rm(ref, simW, simS)
  
  hitdat_all <- Reduce("rbind.data.frame", hitdat)
  
  return(list(refdat      = refdat,
              hitdat      = hitdat_all,
              true_params = true_params))
}

# SimpleEB: Define E(Z|X), where Z is the latent variable representative of whether observation is tp or fp
ez <- function(data, paramsk){
  
  pi_k   <- paramsk[1]
  mufp_k <- paramsk[2] 
  sdfp_k <- sqrt(paramsk[3])
  mutp_k <- paramsk[4]
  sdtp_k <- sqrt(paramsk[5])
  
  edens_tp <- dnorm(data, mean = mutp_k, sd = sdtp_k)
  edens_fp <- dnorm(data, mean = mufp_k, sd = sdfp_k)
  
  ezgx <- (pi_k*edens_tp)/(pi_k*edens_tp + (1-pi_k)*edens_fp) # expected value of latent Z, given observed X
  
  return(ezgx)
}

# SimpleEB: Define M step, i.e. how params are updated 
m_update <- function(data, paramsk){
  
  pi_k    <- paramsk[1]
  mufp_k  <- paramsk[2] 
  varfp_k <- paramsk[3]
  mutp_k  <- paramsk[4]
  vartp_k <- paramsk[5]
  
  ezgx <- ez(data   = data, 
             paramsk = paramsk)
  
  # pi update
  pi <- mean(ezgx)
  
  # mutp update
  mutp <- sum(data*ezgx)/sum(ezgx)
  
  # mufp update
  mufp <- sum(data*(1-ezgx))/sum(1-ezgx)
  
  #vartp update
  vartp <- sum(ezgx*(data - mutp)^2)/sum(ezgx)
  
  #varfp update
  varfp <- sum((1-ezgx)*(data - mufp)^2)/sum(1-ezgx)
  
  return(c(pi, mufp, varfp, mutp, vartp))
  
}

# SimpleEB: Define observed loglikelihood:
obsll <- function(data, params){
  
  pi   <- params[1]
  mufp <- params[2] 
  sdfp <- sqrt(params[3])
  mutp <- params[4]
  sdtp <- sqrt(params[5])
  
  dens_tp <- dnorm(data, mean = mutp, sd = sdtp, log = FALSE)
  dens_fp <- dnorm(data, mean = mufp, sd = sdfp, log = FALSE)
  
  ll <- sum(log(pi*dens_tp + (1-pi)*dens_fp))
  
  return(ll)
}

# SimpleEB: Value initializer (based on k means)
init_valgen <- function(data){
  
  clustlabs <- kmeans(data,2)$cluster
  
  if(mean(data[clustlabs == 1]) < mean(data[clustlabs == 2])) {
    
    imutp <- mean(data[clustlabs == 1])
    imufp <- mean(data[clustlabs == 2])
    
    ivartp <- var(data[clustlabs == 1])
    ivarfp <- var(data[clustlabs == 2])
    
    ipi <- mean(clustlabs == 1)
    
    return(c(ipi, imufp, ivarfp, imutp, ivartp))
    
  } else{
    
    imutp <- mean(data[clustlabs == 2])
    imufp <- mean(data[clustlabs == 1])
    
    ivartp <- var(data[clustlabs == 2])
    ivarfp <- var(data[clustlabs == 1])
    
    ipi <- mean(clustlabs == 2)
    
    return(c(ipi, imufp, ivarfp, imutp, ivartp))
  }
  
}

# SimpleEB: EM algorithm
emfunc <- function(data, eps, max.iters = NULL){
  
  # Initialize parameters
  iterparams <- matrix(init_valgen(data), ncol = 5)
  
  # Initialize counter
  k <- 1
  
  # Determine value of obsll based on initial params
  objval <- obsll(data, params = iterparams[k,])  
  
  while(TRUE){
    
    # update parameters
    iterparams <- rbind(iterparams, m_update(data = data, paramsk = iterparams[k,]))
    
    # recompute obsll
    objval[k+1] <- obsll(data, params = iterparams[k+1,])
    
    # Check if converged
    if(abs(objval[k+1] - objval[k]) < eps){
      break
    }
    
    if(!is.null(max.iters)){
      if(k == max.iters){
        break
      }
    }
    
    k <- k+1
  }
  
  return(list(objval = objval,
              final.params = iterparams[nrow(iterparams),],
              iter.params = iterparams))
}

# SimpleEB: PEP computation
# To compute the FDR, we first need to compute the "posterior error probability",
# which is defined as P(Z = 0 | X), where Z is an indicator denoting between true and false positives.
seb_pepfunc <- function(data, fparams){
  
  pi_t   <- fparams[1]
  mufp_f <- fparams[2] 
  sdfp_f <- sqrt(fparams[3])
  mutp_f <- fparams[4]
  sdtp_f <- sqrt(fparams[5])
  
  fdens_tp <- dnorm(data, mean = mutp_f, sd = sdtp_f)
  fdens_fp <- dnorm(data, mean = mufp_f, sd = sdfp_f)
  
  pep <- ((1-pi_t)*fdens_fp)/((1-pi_t)*fdens_fp + pi_t*fdens_tp)
  return(pep)
}

# Simple EB: FDR estimation
seb_fdrfunc <- function(dataset, thresh){
  tempdat <- dataset[dataset$simS <= thresh,]
  fdr <- mean(tempdat$PEP)
  return(fdr)
}

# Simple EB: True FDR
seb_actualfdr <- function(dataset, thresh){
  tempdat <- dataset[dataset$simS <= thresh,]
  fdr <- mean(tempdat$simW == 0)
  return(fdr)
}

# Complex EB: Value initializer 
ceb_init_valgen <- function(simdata){
  
  clustlabs <- kmeans(simdata$hitdat$simS,2)$cluster
  
  if(mean(simdata$hitdat$simS[clustlabs == 1]) < mean(simdata$hitdat$simS[clustlabs == 2])) {
    
    imutp <- mean(simdata$hitdat$simS[clustlabs == 1])
    imufp <- mean(simdata$hitdat$simS[clustlabs == 2])
    
    ivartp <- var(simdata$hitdat$simS[clustlabs == 1])
    ivarfp <- var(simdata$hitdat$simS[clustlabs == 2])
    
    itau <- mean(clustlabs == 1)
    
  } else{
    
    imutp <- mean(simdata$hitdat$simS[clustlabs == 2])
    imufp <- mean(simdata$hitdat$simS[clustlabs == 1])
    
    ivartp <- var(simdata$hitdat$simS[clustlabs == 2])
    ivarfp <- var(simdata$hitdat$simS[clustlabs == 1])
    
    itau <- mean(clustlabs == 2)
    
  }
  
  ibeta  <- glm(simZ ~ bj + I(bj^2), data = simdata$refdat[simdata$refdat$bj != 0,], family = binomial(link = "logit"))$coefficients
  ialpha <- glm(simZ ~ bjs + I(bjs^2), data = simdata$refdat[simdata$refdat$bjs != 1,], family = binomial(link = "logit"))$coefficients
  ieta0  <- mean(simdata$refdat$simZ[simdata$refdat$bj == 0])
  ieta1  <- mean(simdata$refdat$simZ[simdata$refdat$bjs == 1])
  
  iterparams <- c(0.5, ialpha, ieta1, ibeta, ieta0, itau, imutp, imufp, ivartp, ivarfp)
  return(iterparams)
}

# Complex EB: model fit
ceb_func <- function(simdata, iter){
  Sdat <- simdata$hitdat %>% mutate(i = 1:nrow(simdata$hitdat),
                                    Score = simS,
                                    j = 1)
  for(i in 1:nrow(Sdat)){
    Sdat$j[i] <- which(simdata$refdat$Ref == Sdat$Ref[i])
  }
  Sdat <- Sdat %>% select(j,i,Score)
  
  obdata <- list(Z=simdata$refdat$simZ,
                 S=Sdat,
                 b=simdata$refdat$bj,
                 b_star=simdata$refdat$bjs)
  
  ivs <- ceb_init_valgen(simdata = simdata)
  
  # Note that ivs[1:9] are to initialize alpha through eta0 below.
  # However, model sometimes does not fit with those initializations.
  # So we only initialize tau through sigmaF
  result <- MetIDv2(iter   = iter, 
                    l_list = simdata$refdat %>% rename(Name = Ref), 
                    cutoff = 0, 
                    ob     = obdata, 
                    rho    = 0.5,
                    alpha  = c(1,1,1),
                    eta1   = 0.5,
                    beta   = c(1,1,1),
                    eta0   = 0.5,
                    tau    = ivs[10],
                    muT    = ivs[11], 
                    muF    = ivs[12], 
                    sigmaT = ivs[13], 
                    sigmaF = ivs[14])
  
  iterparams <- cbind(result$parameters$rho, 
                      result$parameters$alpha,
                      result$parameters$eta1,
                      result$parameters$beta,
                      result$parameters$eta0,
                      result$parameters$tau,
                      result$parameters$muT,
                      result$parameters$sigmaT,
                      result$parameters$muF,
                      result$parameters$sigmaF)
  
  return(list(allres = result,
              obdata = obdata,
              final.params = iterparams[nrow(iterparams),],
              iter.params = iterparams))
}

# Complex EB: FDR estimation
ceb_fdrfunc <- function(dataset, thresh){
  tempdat <- dataset[dataset$PY >= thresh,]
  lfdr    <- 1-tempdat$PY
  fdr <- mean(lfdr)
  return(fdr)
}

# Complex EB: True FDR
ceb_actualfdr <- function(dataset, thresh){
  present_names <- dataset$Ref[dataset$PY >= thresh]
  fdr <- 1-mean(dataset$simY[dataset$Ref %in% present_names])
  return(fdr)
}

# Mega function to compare FDR estimation capabilities of the simple and complex EB methods jointly.
# This is essentially a wrapper for all of the code in sections presented below (from Data Simulation onwards)
bigsim <- function(rep, refdat, true_params, num_match, normal = TRUE, seb.eps = 1e-6, seb.maxiters = NULL, ceb.iters = 200){
  
  set.seed(rep)
  
  if(normal){
    simdat <- datsim(refdat = refdat,
                     true_params = true_params,
                     num_match = num_match)
  } else{
    simdat <- datsim_nonorm(refdat = refdat,
                            true_params = true_params,
                            num_match = num_match)
  }
  
  ## Simple EB fit
  sebfit <- emfunc(data = simdat$hitdat$simS, eps = seb.eps, max.iters = seb.maxiters)
  
  FDRdat_seb <- simdat$hitdat %>% 
    mutate(PEP  = seb_pepfunc(data = simdat$hitdat$simS, fparams = sebfit$final.params))
  FDRdat_seb$FDR <- sapply(FDRdat_seb$simS, seb_fdrfunc, dataset = FDRdat_seb)
  FDRdat_seb$tFDR <- sapply(FDRdat_seb$simS, seb_actualfdr, dataset = FDRdat_seb) 
  FDRdat_seb <- FDRdat_seb %>% arrange(simS)
  
  # Computes the average absolute difference between the estimated and actual FDR (based on SEB model)
  seb_FDRerr <- sum(abs(FDRdat_seb$FDR-FDRdat_seb$tFDR))/nrow(FDRdat_seb)
  
  ## Complex EB fit
  cebfit <- ceb_func(simdata = simdat, iter = ceb.iters)
  
  FDRdat_ceb <- cebfit$allres$idtable %>% 
    rename(PY = Confidence,
           Ref = Name) %>% 
    arrange(PY) %>%
    left_join(simdat$refdat)
  
  FDRdat_ceb$FDR <- sapply(FDRdat_ceb$PY, ceb_fdrfunc, dataset = FDRdat_ceb)
  FDRdat_ceb$tFDR <- sapply(FDRdat_ceb$PY, ceb_actualfdr, dataset = FDRdat_ceb)
  
  # Computes the average absolute difference between the estimated and actual FDR (based on CEB model)
  ceb_FDRerr <- sum(abs(FDRdat_ceb$FDR-FDRdat_ceb$tFDR))/nrow(FDRdat_ceb)
  
  return(list(simdat     = simdat,
              sebfit     = sebfit,
              FDRdat_seb = FDRdat_seb,
              cebfit     = cebfit,
              FDRdat_ceb = FDRdat_ceb,
              FDRerr = c(seb_FDRerr, ceb_FDRerr)))
}

# -------------------------------------------------------------------------

### Pre-Sim Data Prep -------------------------------------------------------

# Need to define reference library and competition scores on which to base simulated data on
# We use the reference data provided by CoreMS as the basis for our simulated data

# Library-library cosine correlation scores
ll_scoremat <- read.csv(paste0(misc_folder_path, "CosineCorrelationLowResGCMSDataBase.csv"), row.names = 1, header = TRUE)
rownames(ll_scoremat) <- trimws(rownames(ll_scoremat), which = "both") # Important! Trims leading and trailing spaces
ll_scoremat[lower.tri(ll_scoremat)] <- t(ll_scoremat)[lower.tri(t(ll_scoremat))]

# For compatability with the software written by Jeong et al. (2011), the cosine scores in ll_scoremat must be transformed
# according to the following: New Score = 180/pi*acos(Old Score). Note that the acos() function returns a NaN value if the
# input is even epsilon larger than 1. We therefore cap each value in the matrix to 1 prior to transforming
ll_scoremat[ll_scoremat > 1] <- 1
xfll_scoremat <- 180/pi*acos(ll_scoremat)
rm(ll_scoremat)

# For the sake of simplifying the simulation (i.e. working with a smaller dataset), retain only the first 
# 200 metabolites
xfll_scoremat <- xfll_scoremat[1:200, 1:200]

# Compute competition scores using OAF functions, cutoff of 30 is the default specified by
# the original authors
ak     <- makeak(xfll_scoremat,cutoff = 30)
lslist <- makelslist(xfll_scoremat,cutoff = 30)
bj     <- makebj(lslist,ak)
bjs    <- makebjs(lslist,ak)
names(bj) <- rownames(xfll_scoremat)
names(bjs) <- rownames(xfll_scoremat)
rm(ak, lslist)

refdat <- data.frame(Ref = rownames(xfll_scoremat),
                     bj  = bj,
                     bjs = bjs)

### Simulation Study ----------------------------------------------------------

set.seed(1123)
n.cores <- detectCores()
clust <- makeCluster(n.cores-2)
allfuncs <- setdiff(ls(), c("bigsim", "bj", "bjs", "refdat", "xfll_scoremat", "clust", "n.cores"))
clusterExport(clust, allfuncs)
clusterEvalQ(clust, {
  library(truncnorm)
  library(tidyverse)
})
replist <- sapply(1:100, list)

# Normal (n), clear separation (cs) of TP and FP distributions 
results_ncs <- parLapply(clust, replist, bigsim, 
                         refdat = refdat,
                         true_params = list(rho    = 0.1,
                                            eta0   = 0.5,  
                                            beta0  = 3,    
                                            beta1  = 2,   
                                            beta2  = 1,  
                                            eta1   = 0.7, 
                                            alpha0 = 6,  
                                            alpha1 = 4, 
                                            alpha2 = 2,  
                                            tau    = 0.8,  
                                            muT    = 5,    
                                            varT   = 3^2, 
                                            muF    = 40,   
                                            varF   = 10^2),
                         num_match = 5,
                         normal = TRUE,
                         seb.eps = 1e-6, seb.maxiters = 500, ceb.iters = 500)

stopCluster(clust)

### Figure S1 ---------------------------------------------------------------

# Visualizing true vs actual FDR (GMM)
templist <- list.ungroup(list.select(results_ncs, FDRdat_seb))
seb_mae <- vector("numeric", length = length(results_ncs))
for(i in 1:length(templist)){
  templist[[i]] <- templist[[i]] %>% mutate(Seedno = i)
  seb_mae[i] <- median(abs(templist[[i]]$FDR-templist[[i]]$tFDR))
}
templist <- list.rbind(templist)
tempdat <- list.rbind(list.ungroup(list.select(results_ncs, FDRerr)))

p1 <- ggplot(data = templist, aes(x = simS, y = FDR, group = Seedno)) + geom_line(alpha = 0.3, aes(colour = "Estimated")) + 
  geom_line(aes(y = tFDR, colour = "Truth"), alpha = 0.3) + theme_bw() + ggtitle("Gaussian Mixture Model") + xlab("Similarity Score") +
  annotate("text", x=60, y=0.25, label= paste0("Median MAE: ", round(median(seb_mae),3)), size = 5) + 
  theme(axis.title = element_text(size = 20),
        title = element_text(size = 20),
        legend.title = element_text(size = 12),
        axis.text = element_text(size = 12)) +
  scale_color_manual(name = "Legend", values = c("Estimated" = "black", "Truth" = "red"))



# Visualizing true vs actual FDR (Hierarchical Empirical Bayes model)
templist <- list.ungroup(list.select(results_ncs, FDRdat_ceb))
ceb_mae <- vector("numeric", length = length(results_ncs))
for(i in 1:length(templist)){
  templist[[i]] <- templist[[i]] %>% mutate(Seedno = i)
  ceb_mae[i] <- median(abs(templist[[i]]$FDR-templist[[i]]$tFDR))
}
templist <- list.rbind(templist)
tempdat <- list.rbind(list.ungroup(list.select(results_ncs, FDRerr)))

p2 <- ggplot(data = templist, aes(x = PY, y = FDR, group = Seedno)) + geom_line(alpha = 0.3, aes(colour = "Estimated")) + 
  geom_line(aes(y = tFDR, colour = "Truth"), alpha = 0.3) + theme_bw() + ggtitle("Hierarchical Empirical Bayes Model") + xlab("Pr(Y)") +
  annotate("text", x=0.25, y=0.1, label= paste0("Median MAE: ", round(median(ceb_mae),3)), size = 5) + 
  theme(axis.title = element_text(size = 20),
        title = element_text(size = 15),
        legend.title = element_text(size = 12),
        axis.text = element_text(size = 12)) +
  scale_color_manual(name = "Legend", values = c("Estimated" = "black", "Truth" = "red"))

ggpubr::ggarrange(p1, p2, nrow = 1, common.legend = TRUE, legend = "bottom")
