###   Basic functions for Real Data Analysis   Runchao Jiang


library(compiler)
enableJIT(3)
library(rgenoud)
library(survival)










### Fps: random design for observed treatment
Fps <- function(DataList)
{
  X <- DataList$X
  A <- DataList$A
  
  ps.fit <- glm(A ~ 1, family=binomial(link="logit"))
  piA <- ps.fit$fitted.values
  theta <- ps.fit$coefficients
  
  return(list(piA=piA, theta=theta))
}





### Freg: fit the cox model for the augmented part 
Freg <- function(DataList)
{
  X <- DataList$X
  A <- DataList$A
  Nsize <- length(A)
  Ttilde <- DataList$Ttilde
  delta <- DataList$delta
  
  fit.cox <- coxph(Surv(time=Ttilde,event=delta) ~ X+A+I(A*X), ties='breslow')
  beta <- fit.cox$coefficients
  
  exp1 <- c(exp(cbind(X, 1, X) %*% beta))
  exp0 <- c(exp(X %*% beta[1:ncol(X)]))
  
  expA <- exp1 * A + exp0 * (1 - A)
  den <- sum(expA) - cumsum(expA) + expA
  Lambda0 <- cumsum(delta / den)
  dLambda0 <- c(Lambda0[1], diff(Lambda0))
  
  return(list(beta=beta, Lambda0=Lambda0, dLambda0=dLambda0, exp1=exp1, exp0=exp0))
}







### FPC: fit the Kaplan-Meier estimate for C
FPC <- function(DataList)
{
  fit <- survfit(Surv(DataList$Ttilde, 1 - DataList$delta) ~ 1)
  PC.unique <- fit$surv
  PC <- rep(PC.unique, times = (fit$n.event + fit$n.censor))
  return(PC) 
}






### IPWE: the un-smoothed IPWE
IPWE.nonsmooth <- function(eta, t0, DataList, ps)
{
  A <- DataList$A
  piA <- ps$piA
  
  g <- as.numeric((cbind(1,DataList$X) %*% eta) >= 0)
  
  w <- (A*g+(1-A)*(1-g))/(piA*A+(1-piA)*(1-A))
  
  num <- w * DataList$delta
  den <- sum(w)-cumsum(w)+w
  
  Nevent <- sum(DataList$Ttilde<=t0)    
  result <- prod(1-num[1:Nevent]/den[1:Nevent])
  
  return(result)
}







### smooth.IPWE: smoothed IPWE
IPWE.smooth <- function(eta, t0, DataList, ps)
{
  A <- DataList$A
  piA <- ps$piA
  
  sd.etaX <- sd(c(cbind(1,DataList$X) %*% eta))
  if (!is.finite(sd.etaX)) return(-1000)
  if (sd.etaX > 0) eta <- eta/sd.etaX else eta <- c(ifelse(eta[1] >= 0, 1, -1), rep(0, ncol(DataList$X)))
  g <- pnorm(c(cbind(1,DataList$X) %*% eta)/((length(A)/4)^(-1/3)))
  
  w <- (A*g+(1-A)*(1-g))/(piA*A+(1-piA)*(1-A))
  
  num <- w * DataList$delta
  den <- sum(w)-cumsum(w)+w
  
  Nevent <- sum(DataList$Ttilde<=t0)
  result <- prod(1-num[1:Nevent]/den[1:Nevent])
  
  return(result)
}







### Genetic.IPWE: optimization for IPWE
Genetic.IPWE <- function(t0, DataList, ps, Smooth, Neta)
{  
  if (!Smooth) fn <- IPWE.nonsmooth else fn <- IPWE.smooth
  temp <- genoud(fn=fn, t0=t0, DataList=DataList, ps=ps,
                 nvars=Neta, 
                 Domains=cbind(rep(-1,Neta),rep(1,Neta)),
                 starting.values=rep(0,Neta),
                 max=TRUE,
                 
                 print.level=0,
                 BFGS=FALSE, 
                 optim.method="Nelder-Mead",
                 P9=0, 
                 
                 unif.seed=1107,
                 int.seed=0130)
  
  eta.est <- temp$par
  etahat <- eta.est/sqrt(sum(eta.est^2))
  Shat.etahat <- temp$value
  
  return(c(etahat, Shat.etahat))
}




###   pre-compute eta-free terms in AIPWE
Fprep <- function(DataList, t0, ps, reg, PC)
{
  A <- DataList$A
  piA <- ps$piA
  w.den <- piA*A+(1-piA)*(1-A)
  w1 <- (A)/w.den
  w0 <- (1-A)/w.den
  
  Nevent <- sum(DataList$Ttilde<=t0)
  idx <- 1:Nevent
  
  PT1 <- exp(- reg$Lambda0[idx] %o% reg$exp1)  # row: time points; col: patients
  PT0 <- exp(- reg$Lambda0[idx] %o% reg$exp0)  # row: time points; col: patients
  num.temp <- PC[idx]*reg$dLambda0[idx]
  num.PT1 <- t(t(PT1) * ((1-w1)*reg$exp1)) * num.temp
  num.PT0 <- t(t(PT0) * ((1-w0)*reg$exp0)) * num.temp
  den.temp <- PC[idx]
  den.PT1 <- t(t(PT1) * (1-w1)) * den.temp
  den.PT0 <- t(t(PT0) * (1-w0)) * den.temp
  
  prep <- list(w.den=w.den, idx=idx, 
               num.PT1=num.PT1, num.PT0=num.PT0,
               den.PT1=den.PT1, den.PT0=den.PT0)
  return(prep)
}






### AIPWE: the non-smoothed AIPWE
AIPWE.nonsmooth <- function(eta, t0, DataList, prep)
{
  A <- DataList$A
  delta <- DataList$delta
  idx <- prep$idx
  
  g <- as.numeric(cbind(1,DataList$X) %*% eta >= 0)
  
  w <- (A*g + (1-A)*(1-g))/prep$w.den
  
  num <- w[idx]*delta[idx]+as.double(prep$num.PT1%*%g+prep$num.PT0%*%(1-g))
  
  den <- (sum(w)-cumsum(w)+w)[idx] + as.double(prep$den.PT1%*%g+prep$den.PT0%*%(1-g))
  
  result <- prod(1-num/den)
  
  return(result)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
  
  
  return(result)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
}






AIPWE.smooth <- function(eta, t0, DataList, prep)
{
  A <- DataList$A
  delta <- DataList$delta
  idx <- prep$idx
  
  sd.etaX <- sd(c(cbind(1,DataList$X)%*%eta))
  if (!is.finite(sd.etaX)) return(-1000)
  if (sd.etaX > 0) eta <- eta/sd.etaX else eta <- c(ifelse(eta[1] >= 0, 1, -1), rep(0, ncol(DataList$X)))
  g <- pnorm(c(cbind(1,DataList$X)%*%eta)/((length(A)/4)^(-1/3)))
  
  w <- (A*g + (1-A)*(1-g))/prep$w.den
  
  num <- w[idx]*delta[idx]+as.double(prep$num.PT1%*%g+prep$num.PT0%*%(1-g))
  
  den <- (sum(w)-cumsum(w)+w)[idx] + as.double(prep$den.PT1%*%g+prep$den.PT0%*%(1-g))
  
  result <- prod(1-num/den)
  
  return(result)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
}





Genetic.AIPWE <- function(t0, DataList, prep, Smooth, Neta)
{  
  if (!Smooth) fn <- AIPWE.nonsmooth else fn <- AIPWE.smooth
  temp <- genoud(fn=fn, t0=t0, DataList=DataList, prep=prep,
                 nvars=Neta, 
                 Domains=cbind(rep(-1,Neta),rep(1,Neta)),
                 starting.values=rep(0,Neta),
                 max=TRUE,
                 
                 print.level=0,
                 BFGS=FALSE, 
                 optim.method="Nelder-Mead",
                 P9=0, 
                 
                 unif.seed=1107,
                 int.seed=0130)
    
  eta.est <- temp$par
  etahat <- eta.est/sqrt(sum(eta.est^2))
  Shat.etahat <- temp$value
  
  return(c(etahat, Shat.etahat))
}







### IPWE.se: se(IPWE)
IPWE.se <- function(eta, t0, DataList, ps)
{  
  A <- DataList$A
  Nsize <- length(A)
  X <- DataList$X
  delta <- DataList$delta
  piA <- ps$piA
  theta <- ps$theta
  
  Z <- cbind(rep(1, Nsize))
  dimZ <- ncol(Z)
  
  g <- as.numeric((cbind(1,X) %*% eta) >= 0)
  
  w <- (A*g+(1-A)*(1-g))/(piA*A+(1-piA)*(1-A)) 
  Nevent <- sum(DataList$Ttilde<=t0)    
  
  lambdagu <- ((w*delta)[1:Nevent])/((sum(w)-cumsum(w)+w)[1:Nevent])
  Lambdagu <- cumsum(lambdagu)
  Sgu <- cumprod(1-lambdagu)
  dSgu <- diff(c(1,Sgu))
  
  h <- num1 <- num2 <- matrix(0, nrow=Nsize, ncol=Nevent)
  
  for (i in 1:Nsize)
  {
    for (u in 1:Nevent)
    {
      temp.dNiu <- ifelse(i==u, delta[i], 0)
      temp.Yiu <- ifelse(i>=u, 1, 0)
      h[i,u] <- Sgu[u]*temp.dNiu + temp.Yiu*dSgu[u]
      num1[i,u] <- h[i,u]*w[i]
    }
  }
  
  num2.u.temp <- array(data=0, dim=c(Nsize, Nevent, dimZ))
  for (i in 1:Nsize)
  {
    for (u in 1:Nevent)
    {
      num2.u.temp[i,u,] <- h[i,u]*w[i]^2*(1-2*g[i])*piA[i]/(1+exp(sum(Z[i,]*theta)))*Z[i,]
    }
  }
  num2.u <- apply(num2.u.temp, c(2,3), mean)
  num2.fixed <- solve(t(piA/(1+exp(c(Z %*% theta)))*Z) %*% Z/Nsize)
  num2.i <- (A-piA)*Z
  
  for (i in 1:Nsize)
  {
    for (u in 1:Nevent)
    {
      num2[i,u] <- rbind(num2.u[u,]) %*% num2.fixed %*% cbind(num2.i[i,])
    }
  }
  
  SguSCu <- (sum(w)-cumsum(w)+w)[1:Nevent]/Nsize
  den <- matrix(rep(Sgu*SguSCu, i), nrow=Nsize, byrow=T)
  
  Avar.Lambdag <- mean((rowSums((num1+num2)/den))^2)  
  Shat <- Sgu[Nevent]
  Avar.Sg <- (Shat)^2*Avar.Lambdag
  se.Sg <- sqrt(Avar.Sg/Nsize)
  
  return(se.Sg)
}







### AIPWE.se: se(AIPWE)
AIPWE.se <- function(eta, t0, DataList, ps, reg, PC)
{  
  A <- DataList$A
  Nsize <- length(A)
  X <- DataList$X
  delta <- DataList$delta
  Ttilde <- DataList$Ttilde
  
  piA <- ps$piA
  theta <- ps$theta
  Z <- cbind(rep(1, Nsize))
  dimZ <- ncol(Z)
  
  g <- as.numeric((cbind(1,X) %*% eta) >= 0)  
  
  w <- (A*g+(1-A)*(1-g))/(piA*A+(1-piA)*(1-A))
  Nevent <- sum(Ttilde<=t0)
  
  beta <- reg$beta
  Lambda0 <- reg$Lambda0
  dLambda0 <- reg$dLambda0
  vg <- cbind(X, g, g*X)
  expg <- reg$exp1 * g + reg$exp0 *(1-g)   # expg[i]=exp(beta'vg[i,]) 
  vA <- cbind(X, A, A*X)
  dimv <- ncol(vA)
  expA <- reg$exp1 * A + reg$exp0 *(1-A)   # expg[i]=exp(beta'vA[i,]) 
  
  # part 0: Y, dN, dM
  LambdaC <- cumsum((1 - delta) / (Nsize : 1))
  dLambdaC <- diff(c(0,LambdaC))
  Y <- dN <- dM <- dMC <- matrix(NA, nrow=Nsize, ncol=Nsize)
  for (i in 1:Nsize)
  {
    for (s in 1:Nsize)
    {
      dN[i,s] <- ifelse(i==s, delta[i], 0)
      Y[i,s] <- ifelse(i>=s, 1, 0)
      dM[i,s] <- dN[i,s]-Y[i,s]*expA[i]*dLambda0[s]
      dMC[i,s] <- ifelse(i==s, 1-delta[i], 0) - Y[i,s]*dLambdaC[s]
    }
  }
  
  # part 1: Eden
  den <- num <-  matrix(NA, nrow=Nsize, ncol=Nevent)
  PT <- t(exp(-Lambda0[1:Nevent] %o% expg)) # PT[i,u]=P(Ti>u|Xi, g(Xi)) by cox model
  for (i in 1:Nsize)
  {
    for (u in 1:Nevent)
    {
      num[i,u] <- w[i]*dN[i,u]+(1-w[i])*PT[i,u]*PC[u]*expg[i]*dLambda0[u]
      den[i,u] <- w[i]*Y[i,u]+(1-w[i])*PT[i,u]*PC[u]
    }
  }
  Eden <- colMeans(den)
  
  # part 2: Sgu=S(u|X, g(X))
  Enum <- colMeans(num)
  Sgu <- cumprod(1-Enum/Eden)
  
  # part 3: hi(u)
  h <- matrix(NA, nrow=Nsize, ncol=Nevent)
  dSgu <- diff(c(1, Sgu))
  for (i in 1:Nsize)
  {
    for (u in 1:Nevent)
    {
      h[i,u] <- num[i,u]*Sgu[u]+den[i,u]*dSgu[u]
    }
  }
  phih <- colSums(t(h)/(Eden*Sgu))
  
  # part 4: theta, partial derivative, Nevent by 3
  ptheta <- array(data=NA, dim=c(Nsize, Nevent, dimZ))
  expexp2 <- piA/(1+exp(c(Z%*%theta)))  # expexp2[i]=exp(theta'Z[i,])/(1+exp(\theta'Z[i,]))^2 
  for (i in 1:Nsize)
  {
    for (u in 1:Nevent)
    {
      temp41 <- PT[i,u]*PC[u]*expg[i]*dLambda0[u]
      temp42 <- PT[i,u]*PC[u]
      ptheta[i,u, ] <- ((dN[i,u]-temp41)*Sgu[u]+(Y[i,u]-temp42)*dSgu[u])*(w[i]^2)*(1-2*g[i])*expexp2[i]*Z[i,]
    }
  }
  Eptheta <- apply(ptheta, c(2,3), mean) 
  
  # part 5: theta, IF, Nsize by 3
  IFtheta <- matrix(NA, nrow=Nsize, ncol=dimZ)
  temp5 <- solve(t(expexp2*Z)%*%Z/Nsize)
  for (i in 1:Nsize)
  {
    IFtheta[i,] <- (A[i]-piA[i])* temp5 %*% Z[i,]
  }
  phitheta <- colSums(Eptheta %*% t(IFtheta)/(Eden*Sgu)) 
  
  # part 6: beta partial derivative, Nevent by 5
  pbeta <- array(data=NA, dim=c(Nsize, Nevent, dimv))
  for (i in 1:Nsize)
  {
    for (u in 1:Nevent)
    {
      pbeta[i,u, ] <- (1-w[i])*PC[u]*PT[i,u]*expg[i]*(dLambda0[u]*Sgu[u]*(-Lambda0[u]*expg[i]+1)-dSgu[u]*Lambda0[u])*vg[i,]
    }
  }
  Epbeta <- apply(pbeta, c(2,3), mean)
  
  # part 7: beta IF, Nsize by 5
  Ecoefvv <- sapply(1:Nsize, function(s) t(Y[,s]*expA*vA)%*%vA/Nsize, simplify="array")   # 5 by 5 by Nsize
  Ecoef <- sapply(1:Nsize, function(s) mean(Y[,s]*expA))   # Nsize
  Ecoefv <- sapply(1:Nsize, function(s) colMeans(Y[,s]*expA*vA))   # 5 by Nsize
  Efrac <- sapply(1:Nsize, function(s) (Ecoefvv[,,s]*Ecoef[s]-Ecoefv[,s]%o%Ecoefv[,s])/(Ecoef[s]^2), simplify="array") 
  EE <- solve(apply(Efrac[,,as.logical(delta)], c(1,2), sum) / Nsize)
  
  IFbeta.temp <- array(NA, dim=c(Nsize, Nsize, dimv))
  for (i in 1:Nsize)
  {
    for (s in 1:Nsize)
    {
      IFbeta.temp[i,s,] <- (vA[i,]-Ecoefv[,s]/Ecoef[s])*dM[i,s]
    }
  }
  IFbeta <- t(EE %*% t(apply(IFbeta.temp, c(1,3), sum)))   # Nsize by 5
  
  phibeta <- colSums(Epbeta %*% t(IFbeta)/(Eden*Sgu))
  
  # part 8: Lambda0, partial derivative, length=Nevent
  pLambda0 <- matrix(NA, nrow=Nsize, ncol=Nevent)
  for (i in 1:Nsize)
  {
    for (u in 1:Nevent)
    {
      pLambda0[i,u] <- (1-w[i])*PT[i,u]*(-expg[i])*PC[u]*(expg[i]*dLambda0[u]*Sgu[u]+dSgu[u])
    }
  }
  EpLambda0 <- colMeans(pLambda0)
  
  # part 9: Lambda0, IF, Nsize by Nevent
  temp9.I <- t(t(dM)/Ecoef)
  IFLambda0.I <- t(apply(temp9.I[, 1:Nevent], 1, cumsum))   
  
  temp9.II <- t(t(Ecoefv)/(Ecoef^2)*delta)/Nsize
  temp9.II.fix <- apply(temp9.II[,1:Nevent], 1, cumsum)   # Nevent by 5
  
  IFLambda0.II <- -t(temp9.II.fix %*% t(IFbeta))
  IFLambda0 <- IFLambda0.I + IFLambda0.II   # Nsize by Nevent
  
  phiLambda0 <- colSums(t(IFLambda0)*EpLambda0/(Eden*Sgu))
  
  # part 10 & 11: dLambda0
  pdLambda0 <- IFdLambda0 <- matrix(NA, nrow=Nsize, ncol=Nevent)
  for (i in 1:Nsize)
  {
    for (u in 1:Nevent)
    {
      pdLambda0[i,u] <- (1-w[i])*PT[i,u]*PC[u]*expg[i]*Sgu[u]
    }
  }
  IFdLambda0 <- t(apply(IFLambda0, 1, function(x) diff(c(0, x)))) # surprisingly correct
  EpdLambda0 <- colMeans(pdLambda0)
  phidLambda0 <- colSums(t(IFdLambda0)*EpdLambda0/(Eden*Sgu))
  
  # part 12 & 13: LambdaC
  pLambdaC <- matrix(NA, nrow=Nsize, ncol=Nevent)
  for (i in 1:Nsize)
  {
    for (u in 1:Nevent)
    {
      pLambdaC[i,u] <- (1-w[i])*PT[i,u]*PC[u]*(-1)*(expg[i]*dLambda0[u]*Sgu[u]+dSgu[u])
    }
  }
  EpLambdaC <- colMeans(pLambdaC)
  
  temp13 <- t(t(dMC)/colMeans(Y))
  IFLambdaC <- (t(apply(temp13[,1:Nevent], 1, cumsum)))
  
  phiLambdaC <- colSums(t(IFLambdaC)*EpLambdaC/(Eden*Sgu))
  
  
  # part 14: combine
  phi <- phih + phitheta + phibeta + phiLambda0 + phidLambda0 + phiLambdaC
  
  
  Avar.Lambdahat <- mean(phi^2)
  Shat <- Sgu[Nevent]
  Avar.Shat <- (Shat)^2*Avar.Lambdahat
  se.Shat <- sqrt(Avar.Shat/Nsize)
  
  return(se.Shat)
}






### diff.se: se(IPWE-trt0/trt1)
IPWE.diff.se <- function(eta, t0, DataList, ps)
{  
  A <- DataList$A
  Nsize <- length(A)
  X <- DataList$X
  delta <- DataList$delta
  piA <- ps$piA
  theta <- ps$theta
  Ttilde <- DataList$Ttilde
  
  Z <- cbind(rep(1, Nsize))
  dimZ <- ncol(Z)
  
  g <- as.numeric((cbind(1,X) %*% eta) >= 0)
  
  w <- (A*g+(1-A)*(1-g))/(piA*A+(1-piA)*(1-A)) 
  Nevent <- sum(DataList$Ttilde<=t0)    
  
  lambdagu <- ((w*delta)[1:Nevent])/((sum(w)-cumsum(w)+w)[1:Nevent])
  Lambdagu <- cumsum(lambdagu)
  Sgu <- cumprod(1-lambdagu)
  dSgu <- diff(c(1,Sgu))
  
  h <- num1 <- num2 <- matrix(0, nrow=Nsize, ncol=Nevent)
  
  for (i in 1:Nsize)
  {
    for (u in 1:Nevent)
    {
      temp.dNiu <- ifelse(i==u, delta[i], 0)
      temp.Yiu <- ifelse(i>=u, 1, 0)
      h[i,u] <- Sgu[u]*temp.dNiu + temp.Yiu*dSgu[u]
      num1[i,u] <- h[i,u]*w[i]
    }
  }
  
  num2.u.temp <- array(data=0, dim=c(Nsize, Nevent, dimZ))
  for (i in 1:Nsize)
  {
    for (u in 1:Nevent)
    {
      num2.u.temp[i,u,] <- h[i,u]*w[i]^2*(1-2*g[i])*piA[i]/(1+exp(sum(Z[i,]*theta)))*Z[i,]
    }
  }
  num2.u <- apply(num2.u.temp, c(2,3), mean)
  num2.fixed <- solve(t(piA/(1+exp(c(Z %*% theta)))*Z) %*% Z/Nsize)
  num2.i <- (A-piA)*Z
  
  for (i in 1:Nsize)
  {
    for (u in 1:Nevent)
    {
      num2[i,u] <- rbind(num2.u[u,]) %*% num2.fixed %*% cbind(num2.i[i,])
    }
  }
  
  SguSCu <- (sum(w)-cumsum(w)+w)[1:Nevent]/Nsize
  den <- matrix(rep(Sgu*SguSCu, i), nrow=Nsize, byrow=T)
  
  phi <- rowSums((num1+num2)/den)
  Shat <- Sgu[Nevent]
  surv.phi <- -Shat * phi
  
  
  # diff with trt 0: for subpopulation with A==0
  sub <- A == 0
  Nsubsize <- sum(sub)
  
  fit.sub <- survfit(Surv(Ttilde[sub], delta[sub]) ~ 1)
  Shat.sub <- fit.sub$surv[sum(fit.sub$time <= t0)]
  sub.Lambda <- -log(rep(fit.sub$surv, times = (fit.sub$n.event + fit.sub$n.censor)))
  sub.dLambda <- diff(c(0, sub.Lambda))
  expand.dLambda <- rep(0, Nsize)
  expand.dLambda[sub] <- sub.dLambda
  
  sub.dN <- sub.Y <- sub.num <- sub.den <- matrix(NA, nrow=Nsize, ncol=Nevent)
  for (i in 1:Nsize)
  {
    for (s in 1:Nevent)
    {
      sub.dN[i,s] <- ifelse(i==s, delta[i], 0)
      sub.Y[i,s] <- ifelse(i>=s, 1, 0)
      sub.num[i,s] <- sub[i] * (sub.dN[i,s]-sub.Y[i,s]*expand.dLambda[s])
      sub.den[i,s] <- sub[i] * sub.Y[i,s]
    }
  }
  E.sub.den <- colMeans(sub.den)
  psi <- rowSums(t(t(sub.num)/E.sub.den)) 
  
  surv.psi <- -Shat.sub * psi
  
  diff.se.surv.0 <- sqrt(mean((surv.phi - surv.psi)^2)/Nsize)
  diff.se.Lambda.0 <- sqrt(mean((phi-psi)^2)/Nsize)
  
  
  
  # diff with trt 1: for subpopulation with A==1
  sub <- A == 1
  Nsubsize <- sum(sub)
  
  fit.sub <- survfit(Surv(Ttilde[sub], delta[sub]) ~ 1)
  Shat.sub <- fit.sub$surv[sum(fit.sub$time <= t0)]
  sub.Lambda <- -log(rep(fit.sub$surv, times = (fit.sub$n.event + fit.sub$n.censor)))
  sub.dLambda <- diff(c(0, sub.Lambda))
  expand.dLambda <- rep(0, Nsize)
  expand.dLambda[sub] <- sub.dLambda
  
  sub.dN <- sub.Y <- sub.num <- sub.den <- matrix(NA, nrow=Nsize, ncol=Nevent)
  for (i in 1:Nsize)
  {
    for (s in 1:Nevent)
    {
      sub.dN[i,s] <- ifelse(i==s, delta[i], 0)
      sub.Y[i,s] <- ifelse(i>=s, 1, 0)
      sub.num[i,s] <- sub[i] * (sub.dN[i,s]-sub.Y[i,s]*expand.dLambda[s])
      sub.den[i,s] <- sub[i] * sub.Y[i,s]
    }
  }
  E.sub.den <- colMeans(sub.den)
  psi <- rowSums(t(t(sub.num)/E.sub.den)) 
  
  surv.psi <- -Shat.sub * psi
  
  diff.se.surv.1 <- sqrt(mean((surv.phi - surv.psi)^2)/Nsize)
  diff.se.Lambda.1 <- sqrt(mean((phi-psi)^2)/Nsize)
  
  
  return(c(diff.se.surv.0=diff.se.surv.0, diff.se.Lambda.0=diff.se.Lambda.0, 
           diff.se.surv.1=diff.se.surv.1, diff.se.Lambda.1=diff.se.Lambda.1))
}









### AIPWE.diff.se: se(AIPWE-trt0/trt1)
AIPWE.diff.se <- function(eta, t0, DataList, ps, reg, PC)
{  
  A <- DataList$A
  Nsize <- length(A)
  X <- DataList$X
  delta <- DataList$delta
  Ttilde <- DataList$Ttilde
  
  piA <- ps$piA
  theta <- ps$theta
  Z <- cbind(rep(1, Nsize))
  dimZ <- ncol(Z)
  
  g <- as.numeric((cbind(1,X) %*% eta) >= 0)  
  
  w <- (A*g+(1-A)*(1-g))/(piA*A+(1-piA)*(1-A))
  Nevent <- sum(Ttilde<=t0)
  
  beta <- reg$beta
  Lambda0 <- reg$Lambda0
  dLambda0 <- reg$dLambda0
  vg <- cbind(X, g, g*X)
  expg <- reg$exp1 * g + reg$exp0 *(1-g)   # expg[i]=exp(beta'vg[i,]) 
  vA <- cbind(X, A, A*X)
  dimv <- ncol(vA)
  expA <- reg$exp1 * A + reg$exp0 *(1-A)   # expg[i]=exp(beta'vA[i,]) 
  
  # part 0: Y, dN, dM
  LambdaC <- cumsum((1 - delta) / (Nsize : 1))
  dLambdaC <- diff(c(0,LambdaC))
  Y <- dN <- dM <- dMC <- matrix(NA, nrow=Nsize, ncol=Nsize)
  for (i in 1:Nsize)
  {
    for (s in 1:Nsize)
    {
      dN[i,s] <- ifelse(i==s, delta[i], 0)
      Y[i,s] <- ifelse(i>=s, 1, 0)
      dM[i,s] <- dN[i,s]-Y[i,s]*expA[i]*dLambda0[s]
      dMC[i,s] <- ifelse(i==s, 1-delta[i], 0) - Y[i,s]*dLambdaC[s]
    }
  }
  
  # part 1: Eden
  den <- num <-  matrix(NA, nrow=Nsize, ncol=Nevent)
  PT <- t(exp(-Lambda0[1:Nevent] %o% expg)) # PT[i,u]=P(Ti>u|Xi, g(Xi)) by cox model
  for (i in 1:Nsize)
  {
    for (u in 1:Nevent)
    {
      num[i,u] <- w[i]*dN[i,u]+(1-w[i])*PT[i,u]*PC[u]*expg[i]*dLambda0[u]
      den[i,u] <- w[i]*Y[i,u]+(1-w[i])*PT[i,u]*PC[u]
    }
  }
  Eden <- colMeans(den)
  
  # part 2: Sgu=S(u|X, g(X))
  Enum <- colMeans(num)
  Sgu <- cumprod(1-Enum/Eden)
  
  # part 3: hi(u)
  h <- matrix(NA, nrow=Nsize, ncol=Nevent)
  dSgu <- diff(c(1, Sgu))
  for (i in 1:Nsize)
  {
    for (u in 1:Nevent)
    {
      h[i,u] <- num[i,u]*Sgu[u]+den[i,u]*dSgu[u]
    }
  }
  phih <- colSums(t(h)/(Eden*Sgu))
  
  # part 4: theta, partial derivative, Nevent by 3
  ptheta <- array(data=NA, dim=c(Nsize, Nevent, dimZ))
  expexp2 <- piA/(1+exp(c(Z%*%theta)))  # expexp2[i]=exp(theta'Z[i,])/(1+exp(\theta'Z[i,]))^2 
  for (i in 1:Nsize)
  {
    for (u in 1:Nevent)
    {
      temp41 <- PT[i,u]*PC[u]*expg[i]*dLambda0[u]
      temp42 <- PT[i,u]*PC[u]
      ptheta[i,u, ] <- ((dN[i,u]-temp41)*Sgu[u]+(Y[i,u]-temp42)*dSgu[u])*(w[i]^2)*(1-2*g[i])*expexp2[i]*Z[i,]
    }
  }
  Eptheta <- apply(ptheta, c(2,3), mean) 
  
  # part 5: theta, IF, Nsize by 3
  IFtheta <- matrix(NA, nrow=Nsize, ncol=dimZ)
  temp5 <- solve(t(expexp2*Z)%*%Z/Nsize)
  for (i in 1:Nsize)
  {
    IFtheta[i,] <- (A[i]-piA[i])* temp5 %*% Z[i,]
  }
  phitheta <- colSums(Eptheta %*% t(IFtheta)/(Eden*Sgu)) 
  
  # part 6: beta partial derivative, Nevent by 5
  pbeta <- array(data=NA, dim=c(Nsize, Nevent, dimv))
  for (i in 1:Nsize)
  {
    for (u in 1:Nevent)
    {
      pbeta[i,u, ] <- (1-w[i])*PC[u]*PT[i,u]*expg[i]*(dLambda0[u]*Sgu[u]*(-Lambda0[u]*expg[i]+1)-dSgu[u]*Lambda0[u])*vg[i,]
    }
  }
  Epbeta <- apply(pbeta, c(2,3), mean)
  
  # part 7: beta IF, Nsize by 5
  Ecoefvv <- sapply(1:Nsize, function(s) t(Y[,s]*expA*vA)%*%vA/Nsize, simplify="array")   # 5 by 5 by Nsize
  Ecoef <- sapply(1:Nsize, function(s) mean(Y[,s]*expA))   # Nsize
  Ecoefv <- sapply(1:Nsize, function(s) colMeans(Y[,s]*expA*vA))   # 5 by Nsize
  Efrac <- sapply(1:Nsize, function(s) (Ecoefvv[,,s]*Ecoef[s]-Ecoefv[,s]%o%Ecoefv[,s])/(Ecoef[s]^2), simplify="array") 
  EE <- solve(apply(Efrac[,,as.logical(delta)], c(1,2), sum) / Nsize)
  
  IFbeta.temp <- array(NA, dim=c(Nsize, Nsize, dimv))
  for (i in 1:Nsize)
  {
    for (s in 1:Nsize)
    {
      IFbeta.temp[i,s,] <- (vA[i,]-Ecoefv[,s]/Ecoef[s])*dM[i,s]
    }
  }
  IFbeta <- t(EE %*% t(apply(IFbeta.temp, c(1,3), sum)))   # Nsize by 5
  
  phibeta <- colSums(Epbeta %*% t(IFbeta)/(Eden*Sgu))
  
  # part 8: Lambda0, partial derivative, length=Nevent
  pLambda0 <- matrix(NA, nrow=Nsize, ncol=Nevent)
  for (i in 1:Nsize)
  {
    for (u in 1:Nevent)
    {
      pLambda0[i,u] <- (1-w[i])*PT[i,u]*(-expg[i])*PC[u]*(expg[i]*dLambda0[u]*Sgu[u]+dSgu[u])
    }
  }
  EpLambda0 <- colMeans(pLambda0)
  
  # part 9: Lambda0, IF, Nsize by Nevent
  temp9.I <- t(t(dM)/Ecoef)
  IFLambda0.I <- t(apply(temp9.I[, 1:Nevent], 1, cumsum))   
  
  temp9.II <- t(t(Ecoefv)/(Ecoef^2)*delta)/Nsize
  temp9.II.fix <- apply(temp9.II[,1:Nevent], 1, cumsum)   # Nevent by 5
  
  IFLambda0.II <- -t(temp9.II.fix %*% t(IFbeta))
  IFLambda0 <- IFLambda0.I + IFLambda0.II   # Nsize by Nevent
  
  phiLambda0 <- colSums(t(IFLambda0)*EpLambda0/(Eden*Sgu))
  
  # part 10 & 11: dLambda0
  pdLambda0 <- IFdLambda0 <- matrix(NA, nrow=Nsize, ncol=Nevent)
  for (i in 1:Nsize)
  {
    for (u in 1:Nevent)
    {
      pdLambda0[i,u] <- (1-w[i])*PT[i,u]*PC[u]*expg[i]*Sgu[u]
    }
  }
  IFdLambda0 <- t(apply(IFLambda0, 1, function(x) diff(c(0, x)))) # surprisingly correct
  EpdLambda0 <- colMeans(pdLambda0)
  phidLambda0 <- colSums(t(IFdLambda0)*EpdLambda0/(Eden*Sgu))
  
  # part 12 & 13: LambdaC
  pLambdaC <- matrix(NA, nrow=Nsize, ncol=Nevent)
  for (i in 1:Nsize)
  {
    for (u in 1:Nevent)
    {
      pLambdaC[i,u] <- (1-w[i])*PT[i,u]*PC[u]*(-1)*(expg[i]*dLambda0[u]*Sgu[u]+dSgu[u])
    }
  }
  EpLambdaC <- colMeans(pLambdaC)
  
  temp13 <- t(t(dMC)/colMeans(Y))
  IFLambdaC <- (t(apply(temp13[,1:Nevent], 1, cumsum)))
  
  phiLambdaC <- colSums(t(IFLambdaC)*EpLambdaC/(Eden*Sgu))
  
  
  # part 14: combine
  phi <- phih + phitheta + phibeta + phiLambda0 + phidLambda0 + phiLambdaC
  
  Shat <- Sgu[Nevent] 
  surv.phi <- -Shat * phi
  
  
  # diff with trt 0: for subpopulation with A==0
  sub <- A == 0
  Nsubsize <- sum(sub)
  
  fit.sub <- survfit(Surv(Ttilde[sub], delta[sub]) ~ 1)
  Shat.sub <- fit.sub$surv[sum(fit.sub$time <= t0)]
  sub.Lambda <- -log(rep(fit.sub$surv, times = (fit.sub$n.event + fit.sub$n.censor)))
  sub.dLambda <- diff(c(0, sub.Lambda))
  expand.dLambda <- rep(0, Nsize)
  expand.dLambda[sub] <- sub.dLambda
  
  sub.dN <- sub.Y <- sub.num <- sub.den <- matrix(NA, nrow=Nsize, ncol=Nevent)
  for (i in 1:Nsize)
  {
    for (s in 1:Nevent)
    {
      sub.dN[i,s] <- ifelse(i==s, delta[i], 0)
      sub.Y[i,s] <- ifelse(i>=s, 1, 0)
      sub.num[i,s] <- sub[i] * (sub.dN[i,s]-sub.Y[i,s]*expand.dLambda[s])
      sub.den[i,s] <- sub[i] * sub.Y[i,s]
    }
  }
  E.sub.den <- colMeans(sub.den)
  psi <- rowSums(t(t(sub.num)/E.sub.den)) 
  
  surv.psi <- -Shat.sub * psi
  
  diff.se.surv.0 <- sqrt(mean((surv.phi - surv.psi)^2)/Nsize)
  diff.se.Lambda.0 <- sqrt(mean((phi-psi)^2)/Nsize)
  
  
  
  # diff with trt 1: for subpopulation with A==1
  sub <- A == 1
  Nsubsize <- sum(sub)
  
  fit.sub <- survfit(Surv(Ttilde[sub], delta[sub]) ~ 1)
  Shat.sub <- fit.sub$surv[sum(fit.sub$time <= t0)]
  sub.Lambda <- -log(rep(fit.sub$surv, times = (fit.sub$n.event + fit.sub$n.censor)))
  sub.dLambda <- diff(c(0, sub.Lambda))
  expand.dLambda <- rep(0, Nsize)
  expand.dLambda[sub] <- sub.dLambda
  
  sub.dN <- sub.Y <- sub.num <- sub.den <- matrix(NA, nrow=Nsize, ncol=Nevent)
  for (i in 1:Nsize)
  {
    for (s in 1:Nevent)
    {
      sub.dN[i,s] <- ifelse(i==s, delta[i], 0)
      sub.Y[i,s] <- ifelse(i>=s, 1, 0)
      sub.num[i,s] <- sub[i] * (sub.dN[i,s]-sub.Y[i,s]*expand.dLambda[s])
      sub.den[i,s] <- sub[i] * sub.Y[i,s]
    }
  }
  E.sub.den <- colMeans(sub.den)
  psi <- rowSums(t(t(sub.num)/E.sub.den)) 
  
  surv.psi <- -Shat.sub * psi
  
  diff.se.surv.1 <- sqrt(mean((surv.phi - surv.psi)^2)/Nsize)
  diff.se.Lambda.1 <- sqrt(mean((phi-psi)^2)/Nsize)
  
  
  return(c(diff.se.surv.0=diff.se.surv.0, diff.se.Lambda.0=diff.se.Lambda.0, 
           diff.se.surv.1=diff.se.surv.1, diff.se.Lambda.1=diff.se.Lambda.1))
}




###   compute t-year S at trt0 group and trt1 group 
Fsub <- function(DataList, t0)
{
  Ttilde <- DataList$Ttilde
  delta <- DataList$delta
  A <- DataList$A
  
  index0 <- A==0
  index1 <- A==1
  
  fit.A0 <- survfit(Surv(Ttilde[index0], delta[index0]) ~ 1)
  S0 <- fit.A0$surv[sum(fit.A0$time <= t0)]
  
  fit.A1 <- survfit(Surv(Ttilde[index1], delta[index1]) ~ 1)
  S1 <- fit.A1$surv[sum(fit.A1$time <= t0)]
  
  return(c(S1, S0))
}
