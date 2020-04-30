X_n<-function(n,theta){
  x_t<-rnorm(1,0,1)  
  for(t in 2:n) x_t[t]<-theta*x_t[t-1]+sqrt(1-theta^2)*rnorm(1,0,1) 
  return(x_t)
}


rho2rho_D<-function(rho,lambda){
  Sigma<-matrix(c(1,rho,rho,1),2,2)
  q<-qnorm(lambda)
  rho_D<-(pmvnorm(upper=c(q,q),mean=c(0,0),sigma=Sigma)[1]-lambda^2)/(lambda*(1-lambda))
  return(rho_D)  
}

nkt2rho_D<-function(n,kt){
  p<-mean(kt/n)
  ans<-(var(kt)-n*p*(1-p))/(n*(n-1)*p*(1-p))
  return(ans)
}

rho_dt2rho_Dt<-function(rho,p,n,dt){
ans<-(n*n*p*(1-p)*rho2rho_D(rho*dt,p))/(n*p*(1-p)+n*(n-1)*p*(1-p)*rho2rho_D(rho,p))
return(ans)
}


lambda_S<-function(lambda,rho,S) pnorm((qnorm(lambda)-sqrt(rho)*S)/sqrt(1-rho))

P_k_S<-function(k,n,lambda,rho,S) dbinom(k,n,lambda_S(lambda,rho,S))

P_k<-function(k,n,lambda,rho) integrate(f<-function(S)P_k_S(k,n,lambda,rho,S)*dnorm(S),-Inf,Inf,abs.tol=0.000000001)$value



theta2rho_D1<-function(lambda,rho,theta){
  Sigma<-matrix(c(1,theta,theta,1),2,2)
  q<-qnorm(lambda)
  f<-function(s) lambda_S(lambda,rho,s[1])*lambda_S(lambda,rho,s[2])*dmvnorm(s,c(0,0),Sigma)
  ans<-(adaptIntegrate(f,c(-Inf,-Inf),c(Inf,Inf))$integral-lambda^2)/(lambda*(1-lambda))
return(ans)
}

P_k1k2<-function(k,n,lambda,rho,theta){
Sigma<-matrix(c(1,theta,theta,1),2,2)
f<-function(s) P_k_S(k[1],n[1],lambda,rho,s[1])*P_k_S(k[2],n[2],lambda,rho,s[2])*dmvnorm(s,c(0,0),Sigma)
ans<-adaptIntegrate(f,c(-Inf,-Inf),c(Inf,Inf))$integral
return(ans)
}

P_k2_ck1<-function(k,n,lambda,rho,theta) P_k1k2(k,n,lambda,rho,theta)/P_k(k[1],n[1],lambda,rho)


lambda_S_Z<-function(lambda,nu,rho,S,Z) pnorm((sqrt(Z/nu)*qt(lambda,nu)-sqrt(rho)*S)/sqrt(1-rho))

P_k_S_Z<-function(k,n,lambda,rho,nu,S,Z) dbinom(k,n,lambda_S_Z(lambda,nu,rho,S,Z))
P_k2<-function(k,n,lambda,rho,nu=6) adaptIntegrate(f<-function(x)P_k_S_Z(k,n,lambda,rho,nu,x[1],x[2])*dnorm(x[1])*dchisq(x[2],nu),c(-Inf,0),c(Inf,Inf))$integral
P_k3<-function(k,n,lambda,rho,nu=6,N=1e4){
set.seed(1);S<-rnorm(N);Z<-rchisq(N,nu);ans<-0.0
for(i in 1:N) ans<-ans+P_k_S_Z(k,n,lambda,rho,nu,S[i],Z[i])/N
return(ans)
}


logL<-function(n,k_t){
  return(function(para){
    ans<-0.0
    for(i in 1:length(k_t)) ans<-ans+log(P_k(k_t[i],n,para[1],para[2]))
    return(ans)
    }
    )
  }

logL2<-function(n,k_t){
  return(function(lambda,rho){
    ans<-0.0
    for(i in 1:length(k_t)) ans<-ans+log(P_k(k_t[i],n,lambda,rho))
    return(ans)
    }
    )
  }

logLt<-function(n,k_t,nu){
  return(function(para){
    ans<-0.0
    for(i in 1:length(k_t)) ans<-ans+log(P_k3(k_t[i],n,para[1],para[2],nu,10000))
    return(ans)
    }
    )
  }


logLt2<-function(n,k_t,nu){
  return(function(lambda,rho){
    ans<-0.0
    for(i in 1:length(k_t)) ans<-ans+log(P_k3(k_t[i],n,lambda,rho,nu,10000))
    return(ans)
    }
    )
  }


k_t<-function(T=100,n=1000,lambda=0.1,rho=0,theta=0,seed=123){
set.seed(seed)
q<-qnorm(lambda)
S_t<-matrix(rep(X_n(T,theta),n),nrow=T,ncol=n)
Xi_it<-matrix(rnorm(T*n),nrow=T,ncol=n)
D_it<-ifelse(rho*S_t+sqrt(1-rho*rho)*Xi_it<q,1,0)
k_t<-apply(D_it,MARGIN=1,sum)
return(k_t)
}


P_k_t_S_t<-function(k_t,n,lambda,rho,S_t){
T<-length(k_t)
Ans<-1
for(t in 1:T){
Ans<-Ans*P_k_S(k_t[t],n,lambda,rho,S_t[t])
}
return(Ans)
}

log_P_k_t<-function(k_t,n,lambda,rho,theta){
T<-length(k_t)
Ans<-0.0
S<-100000
set.seed(123)
for(s in 1:S){
S_t<-X_n(T,theta)
Ans<-Ans+P_k_t_S_t(k_t,n,lambda,rho,S_t)
}
Ans<-log(Ans/S)
return(Ans)
}


logL3<-function(n,k_t){
  return(function(para) log_P_k_t(k_t,n,para[1],para[2],para[3]))
}



ms2pars<-function(ms){
m<-50
x<-ms$lambda;y<-ms$rho;z<-ms$theta
rd=kde3d(x=x,y=y,z=z,h=c(bandwidth.nrd(x),bandwidth.nrd(y),bandwidth.nrd(z)),n=m)
i<-which.max(rd$d)%%m
j<-(which.max(rd$d)%/%m)%%m+1
k<-(which.max(rd$d)%/%m)%/%m+1
lambda<-rd$x[i]
rho<-rd$y[j]
theta<-rd$z[k] 
return(c(lambda=lambda,rho=rho,theta=theta))  
}

ms2pars_EXP<-function(ms){
m<-50
x<-ms$p;y<-ms$rho_A;z<-ms$theta
rd=kde3d(x=x,y=y,z=z,h=c(bandwidth.nrd(x),bandwidth.nrd(y),bandwidth.nrd(z)),n=m)
i<-which.max(rd$d)%%m
j<-(which.max(rd$d)%/%m)%%m+1
k<-(which.max(rd$d)%/%m)%/%m+1
p<-rd$x[i]
rho_A<-rd$y[j]
theta<-rd$z[k] 
return(c(p=p,rho_A=rho_A,theta=theta))  
}


ms2pars_POW<-function(ms){
m<-50
x<-ms$p;y<-ms$rho_A;z<-ms$gamma
rd=kde3d(x=x,y=y,z=z,h=c(bandwidth.nrd(x),bandwidth.nrd(y),bandwidth.nrd(z)),n=m)
i<-which.max(rd$d)%%m
j<-(which.max(rd$d)%/%m)%%m+1
k<-(which.max(rd$d)%/%m)%/%m+1
p<-rd$x[i]
rho_A<-rd$y[j]
gamma<-rd$z[k] 
return(c(p=p,rho_A=rho_A,gamma=gamma))  
}



ms2pars2<-function(ms){
m<-50
x<-ms$lambda;y<-ms$rho
rd=kde2d(x=x,y=y,h=c(bandwidth.nrd(x),bandwidth.nrd(y)),n=m)
i<-which.max(rd$z)%%m
j<-(which.max(rd$z)%/%m)%%m+1
lambda<-rd$x[i]
rho<-rd$y[j]
return(c(lambda=lambda,rho=rho))  
}


ms2pars2_Merton<-function(ms){
m<-50
x<-ms$p;y<-ms$rho_A
rd=kde2d(x=x,y=y,h=c(bandwidth.nrd(x),bandwidth.nrd(y)),n=m)
i<-which.max(rd$z)%%m
j<-(which.max(rd$z)%/%m)%%m+1
p<-rd$x[i]
rho_A<-rd$y[j]
return(c(p=p,rho_A=rho_A))  
}


ms2pars3<-function(ms){
m<-50
x<-ms$rho;y<-ms$theta
rd=kde2d(x=x,y=y,h=c(bandwidth.nrd(x),bandwidth.nrd(y)),n=m)
i<-which.max(rd$z)%%m
j<-(which.max(rd$z)%/%m)%%m+1
rho<-rd$x[i]
theta<-rd$y[j]
return(c(rho=rho,theta=theta))  
}


pars2pars_B<-function(T,n,pars,seed){
k_t<-k_t(T,n,pars[1],pars[2],pars[3],seed)
data<-list(T=T,n=n,k_t=k_t)
fit<-sampling(model,data=data,seed=123,iter=1000,warmup=500)
ms<-rstan::extract(fit)
out<-ms2pars(ms)
return(out)
}




dLBBD<-function(x,size,prob,rho_D){
a<-(1-rho_D)*prob/rho_D  
b<-(1-rho_D)*(1-prob)/rho_D
ifelse(rho_D>=0.0,(lchoose(size,x)+lbeta(a+x,b+size-x)-lbeta(a,b)),dbinom(x,size,prob,log=TRUE))
}

dBBD<-function(x,size,prob,rho_D){
a<-(1-rho_D)*prob/rho_D  
b<-(1-rho_D)*(1-prob)/rho_D
ifelse(rho_D>=0.0,exp((lchoose(size,x)+lbeta(a+x,b+size-x)-lbeta(a,b))),dbinom(x,size,prob))
}

dbeta2<-function(x,theta,rhoD){
a<-(1-rhoD)*theta/rhoD
b<-(1-rhoD)*(1-theta)/rhoD
dbeta(x,a,b)  
}


dLBBD_d1<-function(x,size,prob,rho_D,d1,n,k){
if(rho_D>0 ){  
a<-((1-rho_D)*prob/rho_D)+d1*k  
b<-((1-rho_D)*(1-prob)/rho_D)+d1*(n-k)
}
ifelse(rho_D>0.0,(lchoose(size,x)+lbeta(a+x,b+size-x)-lbeta(a,b)),dbinom(x,size,prob,log=TRUE))
}


dBBD_d1<-function(x,size,prob,rho_D,d1,n,k){
a<-((1-rho_D)*prob/rho_D)+d1*k  
b<-((1-rho_D)*(1-prob)/rho_D)+d1*(n-k)
ifelse(rho_D>0.0,exp((lchoose(size,x)+lbeta(a+x,b+size-x)-lbeta(a,b))),dbinom(x,size,prob))
}


dLBBD_EXP<-function(t,kt,Nt,prob,rho_D,r){
a<-((1-rho_D)*prob/rho_D) 
b<-((1-rho_D)*(1-prob)/rho_D)
if(t>=2){
a<-a+sum(r**((t-1):1)*kt[1:(t-1)])
b<-b+sum(r**((t-1):1)*(Nt[1:(t-1)]-kt[1:(t-1)]))
}
ifelse(rho_D>0.0,(lchoose(Nt[t],kt[t])+lbeta(a+kt[t],b+Nt[t]-kt[t])-lbeta(a,b)),dbinom(kt[t],Nt[t],prob,log=TRUE))
}  


dLBBD_POW<-function(t,kt,Nt,prob,rho_D,r){
a<-((1-rho_D)*prob/rho_D) 
b<-((1-rho_D)*(1-prob)/rho_D)
if(t>=2){
a<-a+sum(kt[1:(t-1)]/(t:2)**r)
b<-b+sum((Nt[1:(t-1)]-kt[1:(t-1)])/(t:2)**r)
}
ifelse(rho_D>0.0,(lchoose(Nt[t],kt[t])+lbeta(a+kt[t],b+Nt[t]-kt[t])-lbeta(a,b)),dbinom(kt[t],Nt[t],prob,log=TRUE))
}  


dLBBD_d12<-function(x,size,prob,rho_D,d1,d2,n1,k1,n2,k2){
a<-((1-rho_D)*prob/rho_D)+d1*k1+d2*k2  
b<-((1-rho_D)*(1-prob)/rho_D)+d1*(n1-k1)+d2*(n2-k2)
ifelse(rho_D>0.0,(lchoose(size,x)+lbeta(a+x,b+size-x)-lbeta(a,b)),dbinom(x,size,prob))
}



