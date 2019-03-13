
#BP based survival analysis:
#Copyright(2017): Sujit K. Ghosh, NC State University

BPsurv=function(y,d,m=10){
n=length(y); m=min(m,floor(n/log(n)))
tau=max(y)

#hazard function:
h=function(t,gama){
m=length(gama)
return(sum(gama*dbeta(t/tau,1:m,m:1)/tau) + as.numeric(t>tau)*m*gama[m]/tau)
                   }
#Cumulative hazard function:
H=function(t,gama){
m=length(gama)
return(sum(gama*pbeta(t/tau,1:m,m:1)) + max(t-tau,0)*m*gama[m]/tau)
                  }
#Negative log-likelihood function:
nloglik=function(gamma){
-sum(d*log(sapply(y,h,gama=gamma))-sapply(y,H,gama=gamma))
                       }
fit=optim(par=rep(1,m),fn=nloglik, method="L-BFGS-B",lower=rep(1.0e-8,m))
gamma.est=fit$par
h.est=function(t){h(t,gama=gamma.est)}
hFun=function(t){sapply(t,h.est)}
H.est=function(t){H(t,gama=gamma.est)}
HFun=function(t){sapply(t,H.est)}
SFun=function(t){exp(-HFun(t))}
return(list(SFun=SFun,hFun=hFun,HFun=HFun,m=m))
}







