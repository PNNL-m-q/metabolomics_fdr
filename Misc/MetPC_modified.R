##R package: MetPC##
##Sub function##
#Split Spectra data
# Modified by Jaesik Jeong (Jan. 2 2019)

makelist <- function(data){
  list <- data.frame(do.call(rbind,strsplit(strsplit(as.character(data),split = " ")[[1]],split = ":",fixed = TRUE)))
  list$X1 <- as.numeric(as.character(list$X1))
  list$X2 <- as.numeric(as.character(list$X2))
  list <- list[order(list$X1),]
  spectralist <- list[list$X2!=0,]
  return(spectralist)
}

makelists <- function(data){
  spectralist <- lapply(data,makelist)
  return(spectralist)
}

#Pre-identification : calculate similarity score

cosinescore <- function(list1,list2){
  score_table <- matrix(NA,ncol = length(list1),nrow = length(list2))
  for(i in 1:length(list1)){
    for(j in 1:length(list2)){
      lengthi <- sqrt(sum(list1[[i]]$X2^2))
      lengthj <- sqrt(sum(list2[[j]]$X2^2))
      spectra_table <- merge(list1[[i]],list2[[j]],by="X1")
      cosine <- sum(spectra_table$X2.x*spectra_table$X2.y)/(lengthi*lengthj)
      c_score <- acos(cosine)*180/pi
      if(is.nan(c_score)){c_score <- 0}
      score_table[j,i] <- c_score
    }
  }
  return(score_table)
}

minscore <- function(stable){
  minimumscore <- matrix(0,nrow = ncol(stable), ncol = 3)
  for(i in 1:ncol(stable)){
    minimumscore[i,] <- c(which.min(stable[,i]),i,min(stable[,i]))
  }
  minimumscore <- minimumscore[order(minimumscore[,1]),]
  minimumscore <- data.frame(minimumscore)
  colnames(minimumscore) <- c("j","i","Score")
  return(minimumscore)
}

makeZ <- function(s,jlength){
  z <- rep(0,jlength)
  for(i in 1:jlength){
    zi <- length(s[s[,1]==i,1])
    z[i] <- zi
  }
  z[which(z>0)]=1
  return(z)
}

#Pre-identification : calculate competition score

libraryscore <- function(list){
  score_table <- matrix(0,ncol = length(list),nrow = length(list))
  for(i in 1:(length(list)-1)){
    for(j in (i+1):length(list)){
      lengthi <- sqrt(sum(list[[i]]$X2^2))
      lengthj <- sqrt(sum(list[[j]]$X2^2))
      spectra_table <- merge(list[[i]],list[[j]],by="X1")
      cosine <- sum(spectra_table$X2.x*spectra_table$X2.y)/(lengthi*lengthj)
      c_score <- acos(cosine)*180/pi
      score_table[i,j] <- c_score
      score_table[j,i] <- score_table[i,j]
    }
  }
  return(score_table)
}

makeak <- function(ltable,cutoff){
  ak <- rep(0,nrow(ltable))
  for(i in 1:nrow(ltable)){
    ak[i] <- length(ltable[i,ltable[i,]<cutoff])
  }
  return(ak)
}

makelslist <- function(ltable,cutoff){
  lslist <- list()
  for(i in 1:nrow(ltable)){
    lslist[[i]] <- which(ltable[i,]<cutoff)
  }
  return(lslist)
}

makebj <- function(list,ak){
  inv_ak <- 1/ak
  bj <- rep(0,length(ak))
  for(i in 1:length(ak)){
    bj[i] <- sum(inv_ak[list[[i]][list[[i]]!=i]])
  }
  return(bj)
}

makebjs <- function(list,ak){
  inv_ak <- 1/ak
  bjs <- rep(0,length(ak))
  for(i in 1:length(ak)){
    bjs[i] <- sum(inv_ak[list[[i]]])
  }
  return(bjs)
}

#E-step

logit <- function(b,x){
  1-1/(1+exp(b[1]+b[2]*x+b[3]*x^2))
}

pz_y0 <- function(eta0,bj,beta,z){
  prob <- eta0^(bj==0)*logit(beta,bj)^(bj>0)
  result <- (prob^z)*((1-prob)^(1-z))
  return(result)
}

pz_y1 <- function(eta1,bjs,alpha,z){
  prob <- eta1^(bjs==1)*logit(alpha,bjs)^(bjs!=1)
  result <- (prob^z)*((1-prob)^(1-z))
  return(result)
}

ps_w <- function(s,mu,sigma){
  exp(-((s-mu)^2)/(2*sigma))/(sqrt(2*sigma*pi))
}

pms_w <- function(s,mu1,mu2,sigma1,sigma2,prop){
  prop*ps_w(s,mu1,sigma1)+(1-prop)*ps_w(s,mu2,sigma2)
}

estY <- function(rho,eta0,eta1,bj,bjs,alpha,beta,z,tau,s,muT,muF,sigmaT,sigmaF){
  esty <- z
  for(i in 1:length(z)){
    sv <- prod(pms_w(s[s$j==i,3],muT,muF,sigmaT,sigmaF,tau))
    fv <- prod(ps_w(s[s$j==i,3],muF,sigmaF))
    esty[i] <- (rho*pz_y1(eta1,bjs[i],alpha,z[i])*sv)/((1-rho)*pz_y0(eta0,bj[i],beta,z[i])*fv+rho*pz_y1(eta1,bjs[i],alpha,z[i])*sv)
  }
  return(esty)
}

estW <- function(rho,eta0,eta1,bj,bjs,alpha,beta,z,tau,s,muT,muF,sigmaT,sigmaF){
  estw <- s
  colnames(estw) <- c("j","i","W")
  zj <- z[s$j]
  bjj <- bj[s$j]
  bjsj <- bjs[s$j]
  sv <- ps_w(s$Score,muT,sigmaT)*tau*pz_y1(eta1,bjsj,alpha,zj)*rho
  fv <- ps_w(s$Score,muF,sigmaF)*(pz_y0(eta0,bjj,beta,zj)*(1-rho)+(1-tau)*pz_y1(eta1,bjsj,alpha,zj)*rho)
  estw$W <- sv/(sv+fv)
  return(estw)
}

estYW <- function(rho,eta0,eta1,bj,bjs,alpha,beta,z,tau,s,muT,muF,sigmaT,sigmaF){
  estyw <- s
  colnames(estyw) <- c("j","i","YW")
  for(i in 1:nrow(s)){
    zj <- z[s$j[i]]
    bjj <- bj[s$j[i]]
    bjsj <- bjs[s$j[i]]
    si <- s[-i,]
    v1 <- rho*pz_y1(eta1,bjsj,alpha,zj)*tau*ps_w(s[i,3],muT,sigmaT)*prod(pms_w(si[si$j==s$j[i],3],muT,muF,sigmaT,sigmaF,tau))
    v2 <- (1-rho)*pz_y0(eta0,bjj,beta,zj)*prod(ps_w(s[s$j==s$j[i],3],muF,sigmaF))
    v3 <- rho*pz_y1(eta1,bjsj,alpha,zj)*prod(pms_w(s[s$j==s$j[i],3],muT,muF,sigmaT,sigmaF,tau))
    estyw[i,3] <- v1/(v2+v3)
  }
  return(estyw)
}

estmtY <- function(rho,eta0,eta1,bj,bjs,alpha,beta,z,tau,s,muT1,muT2,muF,sigmaT1,sigmaT2,sigmaF,prop){
  esty <- z
  for(i in 1:length(z)){
    sv <- prod(tau*pms_w(s[s$j==i,3],muT1,muT2,sigmaT1,sigmaT2,prop)+(1-tau)*ps_w(s[s$j==i,3],muF,sigmaF))
    fv <- prod(ps_w(s[s$j==i,3],muF,sigmaF))
    esty[i] <- (rho*pz_y1(eta1,bjs[i],alpha,z[i])*sv)/((1-rho)*pz_y0(eta0,bj[i],beta,z[i])*fv+rho*pz_y1(eta1,bjs[i],alpha,z[i])*sv)
  }
  return(esty)
}

estmtW <- function(rho,eta0,eta1,bj,bjs,alpha,beta,z,tau,s,muT1,muT2,muF,sigmaT1,sigmaT2,sigmaF,prop){
  estw <- s
  colnames(estw) <- c("j","i","W")
  zj <- z[s$j]
  bjj <- bj[s$j]
  bjsj <- bjs[s$j]
  sv <- pms_w(s$Score,muT1,muT2,sigmaT1,sigmaT2,prop)*tau*pz_y1(eta1,bjsj,alpha,zj)*rho
  fv <- ps_w(s$Score,muF,sigmaF)*(pz_y0(eta0,bjj,beta,zj)*(1-rho)+(1-tau)*pz_y1(eta1,bjsj,alpha,zj)*rho)
  estw$W <- sv/(sv+fv)
  return(estw)
}

estmtYW <- function(rho,eta0,eta1,bj,bjs,alpha,beta,z,tau,s,muT1,muT2,muF,sigmaT1,sigmaT2,sigmaF,prop){
  estyw <- s
  colnames(estyw) <- c("j","i","YW")
  for(i in 1:nrow(s)){
    zj <- z[s$j[i]]
    bjj <- bj[s$j[i]]
    bjsj <- bjs[s$j[i]]
    si <- s[-i,]
    v1 <- rho*pz_y1(eta1,bjsj,alpha,zj)*tau*pms_w(s[i,3],muT1,muT2,sigmaT1,sigmaT2,prop)*prod(tau*pms_w(si[si$j==s$j[i],3],muT1,muT2,sigmaT1,sigmaT2,prop)+(1-tau)*ps_w(si[si$j==s$j[i],3],muF,sigmaF))
    v2 <- (1-rho)*pz_y0(eta0,bjj,beta,zj)*prod(ps_w(s[s$j==s$j[i],3],muF,sigmaF))
    v3 <- rho*pz_y1(eta1,bjsj,alpha,zj)*prod(tau*pms_w(s[s$j==s$j[i],3],muT1,muT2,sigmaT1,sigmaT2,prop)+(1-tau)*ps_w(s[s$j==s$j[i],3],muF,sigmaF))
    estyw[i,3] <- v1/(v2+v3)
  }
  return(estyw)
}

estmfY <- function(rho,eta0,eta1,bj,bjs,alpha,beta,z,tau,s,muT,muF1,muF2,sigmaT,sigmaF1,sigmaF2,prop){
  esty <- z
  for(i in 1:length(z)){
    sv <- prod(tau*ps_w(s[s$j==i,3],muT,sigmaT)+(1-tau)*pms_w(s[s$j==i,3],muF1,muF2,sigmaF1,sigmaF2,prop))
    fv <- prod(pms_w(s[s$j==i,3],muF1,muF2,sigmaF1,sigmaF2,prop))
    esty[i] <- (rho*pz_y1(eta1,bjs[i],alpha,z[i])*sv)/((1-rho)*pz_y0(eta0,bj[i],beta,z[i])*fv+rho*pz_y1(eta1,bjs[i],alpha,z[i])*sv)
  }
  return(esty)
}

estmfW <- function(rho,eta0,eta1,bj,bjs,alpha,beta,z,tau,s,muT,muF1,muF2,sigmaT,sigmaF1,sigmaF2,prop){
  estw <- s
  colnames(estw) <- c("j","i","W")
  zj <- z[s$j]
  bjj <- bj[s$j]
  bjsj <- bjs[s$j]
  sv <- ps_w(s$Score,muT,sigmaT)*tau*pz_y1(eta1,bjsj,alpha,zj)*rho
  fv <- pms_w(s$Score,muF1,muF2,sigmaF1,sigmaF2,prop)*(pz_y0(eta0,bjj,beta,zj)*(1-rho)+(1-tau)*pz_y1(eta1,bjsj,alpha,zj)*rho)
  estw$W <- sv/(sv+fv)
  return(estw)
}

estmfYW <- function(rho,eta0,eta1,bj,bjs,alpha,beta,z,tau,s,muT,muF1,muF2,sigmaT,sigmaF1,sigmaF2,prop){
  estyw <- s
  colnames(estyw) <- c("j","i","YW")
  for(i in 1:nrow(s)){
    zj <- z[s$j[i]]
    bjj <- bj[s$j[i]]
    bjsj <- bjs[s$j[i]]
    si <- s[-i,]
    v1 <- rho*pz_y1(eta1,bjsj,alpha,zj)*tau*ps_w(s[i,3],muT,sigmaT)*prod(tau*ps_w(si[si$j==s$j[i],3],muT,sigmaT)+(1-tau)*pms_w(si[si$j==s$j[i],3],muF1,muF2,sigmaF1,sigmaF2,prop))
    v2 <- (1-rho)*pz_y0(eta0,bjj,beta,zj)*prod(pms_w(s[s$j==s$j[i],3],muF1,muF2,sigmaF1,sigmaF2,prop))
    v3 <- rho*pz_y1(eta1,bjsj,alpha,zj)*prod(tau*ps_w(s[s$j==s$j[i],3],muT,sigmaT)+(1-tau)*pms_w(s[s$j==s$j[i],3],muF1,muF2,sigmaF1,sigmaF2,prop))
    estyw[i,3] <- v1/(v2+v3)
  }
  return(estyw)
}

#M-step

E_rho <- function(y){
  sum(y)/length(y)
}

E_eta0 <- function(y,z){
  sum((1-y)*z)/sum(1-y)
}

E_eta1 <- function(y,z){
  sum(y*z)/sum(y)
}

E_tau <- function(y,z,yw){
  n <- nrow(yw)
  cal_table <- cbind(yw,Y=y[yw$j],Z=z[yw$j])
  esttau <- sum(cal_table$YW*cal_table$Z)/sum(cal_table$Y*cal_table$Z)
  return(esttau)
}

estP <- function(piT,muT,muF,sigmaT,sigmaF,s){
  sv <- ps_w(s[,3],muT,sigmaT)*piT
  fv <- ps_w(s[,3],muF,sigmaF)*(1-piT)
  estp <- sv/(sv+fv)
  return(estp)
}

E_prop <- function(p){
  sum(p)/length(p)
}

E_muT <- function(p,s){
  sum(p*s[,3])/sum(p)
}

E_muF <- function(p,s){
  sum((1-p)*s[,3])/sum(1-p)
}

E_sigmaT <- function(p,s,mu){
  sum(p*(s[,3]-mu)^2)/sum(p)
}

E_sigmaF <- function(p,s,mu){
  sum((1-p)*(s[,3]-mu)^2)/sum(1-p)
}

estPs <- function(pis,mu1,mu2,mu3,sigma1,sigma2,sigma3,s){
  v1 <- ps_w(s[,3],mu1,sigma1)*pis[1]
  v2 <- ps_w(s[,3],mu2,sigma2)*pis[2]
  v3 <- ps_w(s[,3],mu3,sigma3)*(1-pis[1]-pis[2])
  estp1 <- v1/(v1+v2+v3)
  estp2 <- v2/(v1+v2+v3)
  estp <- cbind(estp1,estp2)
  return(estp)
}

E_props <- function(ps){
  props1 <- sum(ps[,1])/nrow(ps)
  props2 <- sum(ps[,2])/nrow(ps)
  props <- c(props1,props2)
}

E_mu1 <- function(ps,s){
  sum(ps[,1]*s[,3])/sum(ps[,1])
}

E_mu2 <- function(ps,s){
  sum(ps[,2]*s[,3])/sum(ps[,2])
}

E_mu3 <- function(ps,s){
  sum((1-ps[,1]-ps[,2])*s[,3])/sum(1-ps[,1]-ps[,2])
}

E_sigma1 <- function(ps,s,mu){
  sum(ps[,1]*(s[,3]-mu)^2)/sum(ps[,1])
}

E_sigma2 <- function(ps,s,mu){
  sum(ps[,2]*(s[,3]-mu)^2)/sum(ps[,2])
}

E_sigma3 <- function(ps,s,mu){
  sum((1-ps[,1]-ps[,2])*(s[,3]-mu)^2)/sum(1-ps[,1]-ps[,2])
}

##Main function##
#Peak merging

pmerge <- function(data){
  uniquedata <- NULL
  for(i in 1:nrow(data)){
    uniquedata <- rbind(uniquedata,data[data$Name==data$Name[i]&data$Area==max(data$Area[data$Name==data$Name[i]]),])
  }
  uniquedata <- unique(uniquedata[,colnames(uniquedata)])
  return(uniquedata)
}


#Calculating observed data

cal_ob <- function(sampledata,librarydata,cutoff=30){
  samplelist <- makelists(sampledata)
  librarylist <- makelists(librarydata)
  score <- cosinescore(samplelist,librarylist)
  s <- minscore(score)
  z <- makeZ(s,nrow(score))
  lscore <- libraryscore(librarylist)
  ak <- makeak(lscore,cutoff)
  lslist <- makelslist(lscore,cutoff)
  bj <- makebj(lslist,ak)
  bjs <- makebjs(lslist,ak)
  result <- list(Z=z,S=s,b=bj,b_star=bjs)
  return(result)
}

cal_ob_new <- function(ls_scoremat, ll_scoremat, cutoff=30){
  
  s <- minscore(ls_scoremat)
  z <- makeZ(s,nrow(ls_scoremat))
  
  ak <- makeak(ll_scoremat,cutoff)
  lslist <- makelslist(ll_scoremat,cutoff)
  bj <- makebj(lslist,ak)
  bjs <- makebjs(lslist,ak)
  result <- list(Z=z,S=s,b=bj,b_star=bjs)
  
  return(result)
}


#Estimating parameters

estpar_tf <- function(iter,ob,rho=0.5,eta0=0.5,eta1=0.5,alpha=c(1,1,1),beta=c(1,1,1),tau=0.5,piT=0.5,muT,muF,sigmaT,sigmaF, optalg){
  y <- 0
  z <- ob$Z
  w <- 0
  yw <- 0
  p <- 0
  s <- ob$S
  bj <- ob$b
  bjs <- ob$b_star
  rhos <- rep(0,iter)
  eta0s <- rep(0,iter)
  eta1s <- rep(0,iter)
  alphas <- matrix(0,nrow = iter,ncol = 3)
  betas <- matrix(0,nrow = iter,ncol = 3)
  taus <- rep(0,iter)
  piTs <- rep(0,iter)
  muTs <- rep(0,iter)
  muFs <- rep(0,iter)
  sigmaTs <- rep(0,iter)
  sigmaFs <- rep(0,iter)
  Qa <- function(alpha){
    sum(y*(z*log(logit(alpha,bjs))+(1-z)*log(1-logit(alpha,bjs))))
  }
  Qb <- function(beta){
    sum((1-y)*(z*log(logit(beta,bj))+(1-z)*log(1-logit(beta,bj))))
  }
  for(i in 1:iter){
    y <- estY(rho,eta0,eta1,bj,bjs,alpha,beta,z,tau,s,muT,muF,sigmaT,sigmaF)
    w <- estW(rho,eta0,eta1,bj,bjs,alpha,beta,z,tau,s,muT,muF,sigmaT,sigmaF)
    yw <- estYW(rho,eta0,eta1,bj,bjs,alpha,beta,z,tau,s,muT,muF,sigmaT,sigmaF)
    p <- estP(piT,muT,muF,sigmaT,sigmaF,s)
    rhos[i] <- E_rho(y)
    eta0s[i] <- E_eta0(y,z)
    eta1s[i] <- E_eta1(y,z)
    alphas[i,] <- optim(alpha, Qa, method = optalg, control = list(fnscale = -1))$par
    betas[i,] <- optim(beta, Qb, method = optalg, control = list(fnscale = -1))$par
    taus[i] <- E_tau(y,z,yw)
    piTs[i] <- E_prop(p)
    muTs[i] <- E_muT(p,s)
    muFs[i] <- E_muF(p,s)
    sigmaTs[i] <- E_sigmaT(p,s,muTs[i])
    sigmaFs[i] <- E_sigmaF(p,s,muFs[i])
    rho <- rhos[i]
    eta0 <- eta0s[i]
    eta1 <- eta1s[i]
    alpha <- alphas[i,]
    beta <- betas[i,]
    tau <- taus[i]
    piT <- piTs[i]
    muT <- muTs[i]
    muF <- muFs[i]
    sigmaT <- sigmaTs[i]
    sigmaF <- sigmaFs[i]
  }
  alphas <- data.frame(alphas)
  colnames(alphas) <- c("alpha0","alpha1","alpha2")
  betas <- data.frame(betas)
  colnames(betas) <- c("beta0","beta1","beta2")
  parameter <- list(rho=rhos,eta0=eta0s,eta1=eta1s,alpha=alphas,beta=betas,tau=taus,piT=piTs,muT=muTs,muF=muFs,sigmaT=sigmaTs,sigmaF=sigmaFs)
  return(parameter)
}

estpar_ttf <- function(iter,ob,rho=0.5,eta0=0.5,eta1=0.5,alpha=c(1,1,1),beta=c(1,1,1),tau=0.5,pis=c(0.3,0.3),muT1,muT2,muF,sigmaT1,sigmaT2,sigmaF, optalg){
  y <- 0
  z <- ob$Z
  w <- 0
  yw <- 0
  ps <- 0
  s <- ob$S
  bj <- ob$b
  bjs <- ob$b_star
  rhos <- rep(0,iter)
  eta0s <- rep(0,iter)
  eta1s <- rep(0,iter)
  alphas <- matrix(0,nrow = iter,ncol = 3)
  betas <- matrix(0,nrow = iter,ncol = 3)
  taus <- rep(0,iter)
  piss <- matrix(0,nrow=iter,ncol=2)
  muT1s <- rep(0,iter)
  muT2s <- rep(0,iter)
  muFs <- rep(0,iter)
  sigmaT1s <- rep(0,iter)
  sigmaT2s <- rep(0,iter)
  sigmaFs <- rep(0,iter)
  Qa <- function(alpha){
    sum(y*(z*log(logit(alpha,bjs))+(1-z)*log(1-logit(alpha,bjs))))
  }
  Qb <- function(beta){
    sum((1-y)*(z*log(logit(beta,bj))+(1-z)*log(1-logit(beta,bj))))
  }
  for(i in 1:iter){
    prop <- pis[2]/(1-pis[1])
    y <- estmtY(rho,eta0,eta1,bj,bjs,alpha,beta,z,tau,s,muT1,muT2,muF,sigmaT1,sigmaT2,sigmaF,prop)
    w <- estmtW(rho,eta0,eta1,bj,bjs,alpha,beta,z,tau,s,muT1,muT2,muF,sigmaT1,sigmaT2,sigmaF,prop)
    yw <- estmtYW(rho,eta0,eta1,bj,bjs,alpha,beta,z,tau,s,muT1,muT2,muF,sigmaT1,sigmaT2,sigmaF,prop)
    ps <- estPs(pis,muT1,muT2,muF,sigmaT1,sigmaT2,sigmaF,s)
    rhos[i] <- E_rho(y)
    eta0s[i] <- E_eta0(y,z)
    eta1s[i] <- E_eta1(y,z)
    alphas[i,] <- optim(alpha, Qa, method = optalg, control = list(fnscale = -1))$par
    betas[i,] <- optim(beta, Qb, method = optalg, control = list(fnscale = -1))$par
    taus[i] <- E_tau(y,z,yw)
    piss[i,] <- E_props(ps)
    muT1s[i] <- E_mu1(ps,s)
    muT2s[i] <- E_mu2(ps,s)
    muFs[i] <- E_mu3(ps,s)
    sigmaT1s[i] <- E_sigma1(ps,s,muT1s[i])
    sigmaT2s[i] <- E_sigma2(ps,s,muT2s[i])
    sigmaFs[i] <- E_sigma3(ps,s,muFs[i])
    rho <- rhos[i]
    eta0 <- eta0s[i]
    eta1 <- eta1s[i]
    alpha <- alphas[i,]
    beta <- betas[i,]
    tau <- taus[i]
    pis <- piss[i,]
    muT1 <- muT1s[i]
    muT2 <- muT2s[i]
    muF <- muFs[i]
    sigmaT1 <- sigmaT1s[i]
    sigmaT2 <- sigmaT2s[i]
    sigmaF <- sigmaFs[i]
  }
  alphas <- data.frame(alphas)
  colnames(alphas) <- c("alpha0","alpha1","alpha2")
  betas <- data.frame(betas)
  colnames(betas) <- c("beta0","beta1","beta2")
  parameter <- list(rho=rhos,eta0=eta0s,eta1=eta1s,alpha=alphas,beta=betas,tau=taus,pis=piss,muT1=muT1s,muT2=muT2s,muF=muFs,sigmaT1=sigmaT1s,sigmaT2=sigmaT2s,sigmaF=sigmaFs)
  return(parameter)
}

estpar_tff <- function(iter,ob,rho=0.5,eta0=0.5,eta1=0.5,alpha=c(1,1,1),beta=c(1,1,1),tau=0.5,pis=c(0.3,0.3),muT,muF1,muF2,sigmaT,sigmaF1,sigmaF2, optalg){
  y <- 0
  z <- ob$Z
  w <- 0
  yw <- 0
  ps <- 0
  s <- ob$S
  bj <- ob$b
  bjs <- ob$b_star
  rhos <- rep(0,iter)
  eta0s <- rep(0,iter)
  eta1s <- rep(0,iter)
  alphas <- matrix(0,nrow = iter,ncol = 3)
  betas <- matrix(0,nrow = iter,ncol = 3)
  taus <- rep(0,iter)
  piss <- matrix(0,nrow=iter,ncol=2)
  muTs <- rep(0,iter)
  muF1s <- rep(0,iter)
  muF2s <- rep(0,iter)
  sigmaTs <- rep(0,iter)
  sigmaF1s <- rep(0,iter)
  sigmaF2s <- rep(0,iter)
  Qa <- function(alpha){
    sum(y*(z*log(logit(alpha,bjs))+(1-z)*log(1-logit(alpha,bjs))))
  }
  Qb <- function(beta){
    sum((1-y)*(z*log(logit(beta,bj))+(1-z)*log(1-logit(beta,bj))))
  }
  for(i in 1:iter){
    prop <- pis[2]/(1-pis[1])
    y <- estmfY(rho,eta0,eta1,bj,bjs,alpha,beta,z,tau,s,muT,muF1,muF2,sigmaT,sigmaF1,sigmaF2,prop)
    w <- estmfW(rho,eta0,eta1,bj,bjs,alpha,beta,z,tau,s,muT,muF1,muF2,sigmaT,sigmaF1,sigmaF2,prop)
    yw <- estmfYW(rho,eta0,eta1,bj,bjs,alpha,beta,z,tau,s,muT,muF1,muF2,sigmaT,sigmaF1,sigmaF2,prop)
    ps <- estPs(pis,muT,muF1,muF2,sigmaT,sigmaF1,sigmaF2,s)
    rhos[i] <- E_rho(y)
    eta0s[i] <- E_eta0(y,z)
    eta1s[i] <- E_eta1(y,z)
    alphas[i,] <- optim(alpha, Qa, method = optalg, control = list(fnscale = -1))$par
    betas[i,] <- optim(beta, Qb, method = optalg, control = list(fnscale = -1))$par
    taus[i] <- E_tau(y,z,yw)
    piss[i,] <- E_props(ps)
    muTs[i] <- E_mu1(ps,s)
    muF1s[i] <- E_mu2(ps,s)
    muF2s[i] <- E_mu3(ps,s)
    sigmaTs[i] <- E_sigma1(ps,s,muTs[i])
    sigmaF1s[i] <- E_sigma2(ps,s,muF1s[i])
    sigmaF2s[i] <- E_sigma3(ps,s,muF2s[i])
    rho <- rhos[i]
    eta0 <- eta0s[i]
    eta1 <- eta1s[i]
    alpha <- alphas[i,]
    beta <- betas[i,]
    tau <- taus[i]
    pis <- piss[i,]
    muT <- muTs[i]
    muF1 <- muF1s[i]
    muF2 <- muF2s[i]
    sigmaT <- sigmaTs[i]
    sigmaF1 <- sigmaF1s[i]
    sigmaF2 <- sigmaF2s[i]
  }
  alphas <- data.frame(alphas)
  colnames(alphas) <- c("alpha0","alpha1","alpha2")
  betas <- data.frame(betas)
  colnames(betas) <- c("beta0","beta1","beta2")
  parameter <- list(rho=rhos,eta0=eta0s,eta1=eta1s,alpha=alphas,beta=betas,tau=taus,pis=piss,muT=muTs,muF1=muF1s,muF2=muF2s,sigmaT=sigmaTs,sigmaF1=sigmaF1s,sigmaF2=sigmaF2s)
  return(parameter)
}

#Making the identification table

id_table <- function(iter,ob,rho=0.5,eta0=0.5,eta1=0.5,alpha=c(1,1,1),beta=c(1,1,1),tau=0.5,pis=c(0.3,0.3),muT,muT2=NA,muF,muF2=NA,sigmaT,sigmaT2=NA,sigmaF,sigmaF2=NA,l_list,cutoff=0.8){
  z <- ob$Z
  s <- ob$S
  bj <- ob$b
  bjs <- ob$b_star
  prop <- pis[2]/(1-pis[1])
  if(!(is.na(muT2)|is.na(muF2)|is.na(sigmaT2)|is.na(sigmaF2))){
    stop("Too many values have been entered")
  }else if(is.na(muT2)&is.na(muF2)&is.na(sigmaT2)&is.na(sigmaF2)){
    Py <- estY(rho,eta0,eta1,bj,bjs,alpha,beta,z,tau,s,muT,muF,sigmaT,sigmaF)
  }else if(!(is.na(muF2)|is.na(sigmaF2))){
    Py <- estmfY(rho,eta0,eta1,bj,bjs,alpha,beta,z,tau,s,muT,muF,muF2,sigmaT,sigmaF,sigmaF2,prop)
  }else if(!(is.na(muT2)|is.na(sigmaT2))){
    Py <- estmtY(rho,eta0,eta1,bj,bjs,alpha,beta,z,tau,s,muT,muT2,muF,sigmaT,sigmaT2,sigmaF,prop)
  }else{
    stop("Invalid variable input")
  }
  cut <- which(Py>=cutoff)
  Name <- l_list$Name[cut]
  Confidence <- Py[cut]
  idtable <- data.frame(Name = Name[order(Confidence, decreasing = T)], Confidence = Confidence[order(Confidence, decreasing = T)])
  return(idtable)
}

#Combination function

MetIDv2 <- function(iter,l_list,optalg = "Nelder-Mead",cutoff=0.8,ob,rho=0.5,eta0=0.5,eta1=0.5,alpha=c(1,1,1),beta=c(1,1,1),tau=0.5,piT=0.5,pis=c(0.3,0.3),muT,muT2=NA,muF,muF2=NA,sigmaT,sigmaT2=NA,sigmaF,sigmaF2=NA){
  z <- ob$Z
  s <- ob$S
  bj <- ob$b
  bjs <- ob$b_star
  if(!(is.na(muT2)|is.na(muF2)|is.na(sigmaT2)|is.na(sigmaF2))){
    stop("Too many values have been entered")
  }else if(is.na(muT2)&is.na(muF2)&is.na(sigmaT2)&is.na(sigmaF2)){
    parameters <- estpar_tf(iter,ob,rho,eta0,eta1,alpha,beta,tau,piT,muT,muF,sigmaT,sigmaF,optalg)
    Py <- estY(rho,eta0,eta1,bj,bjs,alpha,beta,z,tau,s,muT,muF,sigmaT,sigmaF)
  }else if(!(is.na(muF2)|is.na(sigmaF2))){
    parameters <- estpar_tff(iter,ob,rho,eta0,eta1,alpha,beta,tau,pis,muT,muF,muF2,sigmaT,sigmaF,sigmaF2,optalg)
    prop <- pis[2]/(1-pis[1])
    Py <- estmfY(rho,eta0,eta1,bj,bjs,alpha,beta,z,tau,s,muT,muF,muF2,sigmaT,sigmaF,sigmaF2,prop)
  }else if(!(is.na(muT2)|is.na(sigmaT2))){
    parameters <- estpar_ttf(iter,ob,rho,eta0,eta1,alpha,beta,tau,pis,muT,muT2,muF,sigmaT,sigmaT2,sigmaF,optalg)
    prop <- pis[2]/(1-pis[1])
    Py <- estmtY(rho,eta0,eta1,bj,bjs,alpha,beta,z,tau,s,muT,muT2,muF,sigmaT,sigmaT2,sigmaF,prop)
  }else{
    stop("Invalid variable input")
  }
  cut <- which(Py>=cutoff)
  Name <- l_list$Name[cut]
  Confidence <- Py[cut]
  idtable <- data.frame(Name = Name[order(Confidence, decreasing = T)],Confidence = Confidence[order(Confidence, decreasing = T)])
  return(list(idtable = idtable,
              parameters = parameters,
              Py = Py))
}

MetID <- function(iter, optalg, l_names,cutoff=0.8,ob,rho=0.5,eta0=0.5,eta1=0.5,alpha=c(1,1,1),beta=c(1,1,1),tau=0.5,piT=0.5,pis=c(0.3,0.3),muT,muT2=NA,muF,muF2=NA,sigmaT,sigmaT2=NA,sigmaF,sigmaF2=NA){
  z <- ob$Z
  s <- ob$S
  bj <- ob$b
  bjs <- ob$b_star
  if(!(is.na(muT2)|is.na(muF2)|is.na(sigmaT2)|is.na(sigmaF2))){
    stop("Too many values have been entered")
  }else if(is.na(muT2)&is.na(muF2)&is.na(sigmaT2)&is.na(sigmaF2)){
    parameters <- estpar_tf(iter,ob,rho,eta0,eta1,alpha,beta,tau,piT,muT,muF,sigmaT,sigmaF,optalg)
    Py <- estY(rho,eta0,eta1,bj,bjs,alpha,beta,z,tau,s,muT,muF,sigmaT,sigmaF)
  }else if(!(is.na(muF2)|is.na(sigmaF2))){
    parameters <- estpar_tff(iter,ob,rho,eta0,eta1,alpha,beta,tau,pis,muT,muF,muF2,sigmaT,sigmaF,sigmaF2,optalg)
    prop <- pis[2]/(1-pis[1])
    Py <- estmfY(rho,eta0,eta1,bj,bjs,alpha,beta,z,tau,s,muT,muF,muF2,sigmaT,sigmaF,sigmaF2,prop)
  }else if(!(is.na(muT2)|is.na(sigmaT2))){
    parameters <- estpar_ttf(iter,ob,rho,eta0,eta1,alpha,beta,tau,pis,muT,muT2,muF,sigmaT,sigmaT2,sigmaF,optalg)
    prop <- pis[2]/(1-pis[1])
    Py <- estmtY(rho,eta0,eta1,bj,bjs,alpha,beta,z,tau,s,muT,muT2,muF,sigmaT,sigmaT2,sigmaF,prop)
  }else{
    stop("Invalid variable input")
  }
  cut <- which(Py>=cutoff)
  # Name <- l_list$Name[cut]
  Name <- l_names[cut]
  Confidence <- Py[cut]
  idtable <- data.frame(Name = Name[order(Confidence, decreasing = T)],Confidence = Confidence[order(Confidence, decreasing = T)])
  return(list(idtable = idtable,
              parameters = parameters,
              Py = Py))
}

#Dimension reduction
#Duplicate metabolite screening

comp_list <- function(samplelists,Num){
  if(Num<2){
    stop("More than one sample is required.")
  }
  complist <- samplelists[[1]]
  for(i in 2:Num){
    complist <- merge(complist[,1:2],samplelists[[i]],by="Name")
  }
  return(complist$Name)
}

#Numbering of metabolite in sample

givenum <- function(m_library,complist){
  num <- rep(0,length(complist))
  for(i in 1:length(complist)){
    num[i] <- which(as.character(m_library$Name)==as.character(complist[i]))
  }
  n_complist <- data.frame(Name=complist,num)
  return(n_complist)
}

# Output reduced values

comp_vec <- function(s,complist,sampledata,cutoff,method){
  s <- s[s$Score<cutoff,]
  compvec <- rep(0,length(complist$num))
  if(missing(method)){
    stop("Enter the ''TIC'' or ''BPC'' method.")
  }
  if(method=="TIC"){
    rpvalue <- function(data,num){sum(data)/num}
  }else if(method=="BPC"){
    rpvalue <- function(data,num){max(data)}
  }else{stop("Enter the ''TIC'' or ''BPC'' method.")}
  for(i in 1:length(complist$num)){
    n <- s[s$j==complist$num[i],2]
    metms <- makelists(sampledata$Spectra[n])
    intensity_data <- metms[[1]]$X2
    if(length(n)>=2){
      for(j in 2:length(n)){
        intensity_data <- c(intensity_data,metms[[j]]$X2)
      }
    }
    compvec[i] <- rpvalue(intensity_data,length(n))
  }
  return(compvec)
}

# Generate trace plot of 4 parameter estimates
plot_pars <- function(pars, iter){
	plot1 <- ggplot(data.frame(Iteration=1:iter, rho=pars$rho), aes(x=Iteration, y=rho))+geom_line()
	plot2 <- ggplot(data.frame(Iteration=1:iter, tau=pars$tau), aes(x=Iteration, y=tau))+geom_line()
	plot3 <- ggplot(data.frame(Iteration=1:iter, muT=pars$muT), aes(x=Iteration, y=muT))+geom_line()
	plot4 <- ggplot(data.frame(Iteration=1:iter, sigmaT=pars$sigmaT), aes(x=Iteration, y=sigmaT))+geom_line()
	grid.arrange(plot1, plot2, plot3, plot4, nrow=2, ncol=2)
}

