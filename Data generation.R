## Data generating function
#######################
##Joint model with interval censored data and two-phase longitudinal biomarker

##Using full likelihood and GH integration to calculate score function


#######################

## library(lme4)
#### Define variance process model here
Sigma_AR <- function(tseq, theta) {
  n <- length(tseq)
  
  # Check for correct parameter length
  if (length(theta) != n + 1) {
    stop("Length of theta must be n+1 (n variances + 1 rho)")
  }
  
  # Extract parameters
  variances <- theta[1:n]
  rho <- theta[n + 1]
  
  # D matrix (diagonal of standard deviations)
  sds <- sqrt(variances)
  D <- diag(sds)
  
  # R matrix (Correlation)
  R <- rho^abs(outer(tseq, tseq, "-"))
  
  # Sigma = D * R * D
  Sigma <- D %*% R %*% D
  
  return(as.matrix(Sigma))
}
Sigma_CS <- function(tseq, theta) {
  n <- length(tseq)
  
  # Check for correct parameter length
  if (length(theta) != n + 1) {
    stop("Length of theta must be n+1 (n variances + 1 rho)")
  }
  
  # Extract parameters
  variances <- theta[1:n]
  rho <- theta[n + 1]
  
  # D matrix (diagonal of standard deviations)
  sds <- sqrt(variances)
  D <- diag(sds)
  
  # R matrix (Correlation)
  # R_ii = rho^0 = 1
  # R_ij = rho^1 = rho
  R <- rho^(abs(outer(tseq, tseq, "-")) > 0)
  
  # Sigma = D * R * D
  Sigma <- D %*% R %*% D
  
  return(as.matrix(Sigma))
}

Lambda.0=function(t,par){
  if (t>t0) {return(((t-t0)/par[1])^par[2])} else {return(0)}
}


#Define inverse function of Lambda.0; 
## input cumu hazard value, output time t
inv.Lambda.0=function(lam, par){
  lam^(1/par[2])*par[1]+t0
}


######################################################
## change data.lb into wide format 
changeformat=function(data.lb)
{
  data.lb_base=data.lb[,c("id","visit","visits","befdiag","afterdiag",
                          "xval",paste("aval_",1:K,sep=""),"preetime","etime","ctime",
                          "followuptime","status",paste("quant",1:length(quant_d),sep=""),
                          "lambdaindex","b","group")]
  
  yname=paste(rep(paste("ynew",as.character(1:length(sigma.errlist)),sep=""), each=K),
              paste("_",1:K,sep=""),sep="")
  yerrname=paste(rep(paste("yerr",as.character(1:length(sigma.errlist)),sep=""), each=K),
                 paste("_",1:K,sep=""),sep="")
  
  valuename=c("visittime",yname,yerrname)
  for (i in 1:length(valuename))
  {
    data.lb_new=cbind(data.lb_base,data.lb[valuename[i]])
    data.lb_wide=spread(data.lb_new,"visit",valuename[i])
    newname=paste(valuename[i],"_vis",as.character(1:J),sep="")
    colnames(data.lb_wide)[which(colnames(data.lb_wide)%in%(1:J))]=newname
    if (i==1) {data.lb_all=data.lb_wide} else 
    {data.lb_all=cbind(data.lb_all,data.lb_wide[newname])}
  } ## end i
  return(data.lb_all)
}

####################################################
simexfun=function(data.intcensor,yterm,lambdavecterm,lambdaaddterm,groupterm)
{
  data.lb_wide_all=NULL
  ## get the columns containing the name yterm and "id"
  data.yterm=data.intcensor[grepl(paste(yterm,"|id",sep=""),
                                  colnames(data.intcensor))]
  for (lambda in lambdavecterm){
    for (bind in 1:B){
      
      dataall=NULL
      for (i in 1:n.id){
        datai=data.yterm[which(data.yterm$id==i),]
        for (sigma.err in sigma.errlist)
        {
          Sigma.err=Sigma_CS(1:K,c(sigma.err^2,0))
          simexerr=mvrnorm(length(which(data.yterm$id==i)),mu=rep(0,K),Sigma=Sigma.err)
          ## update yterm, column 1 is id
          datai[,1+(1:K)+(which(sigma.errlist==sigma.err)-1)*K]=
            datai[,1+(1:K)+(which(sigma.errlist==sigma.err)-1)*K]+
            sqrt(lambda+lambdaaddterm)*simexerr
        } ## end sigma.err
        dataall=rbind(dataall,datai)
      } ## end i
      
      data.lb=data.intcensor
      data.lb[colnames(dataall)[1+(1:(length(sigma.errlist)*K))]]=
        dataall[1+(1:(length(sigma.errlist)*K))]
      data.lb=data.frame(data.lb,lambdaindex=which(lambdavecterm==lambda),b=bind,
                         group=groupterm)
      
      data.lb_wide=changeformat(data.lb)
      
      data.lb_wide_all=rbind(data.lb_wide_all,data.lb_wide)
    } ## end bind
  }  ## end lambda
  return(data.lb_wide_all)
} ## end simexfun

################################################
## function for generating SIMEX data for 
## M*(no measurement error) and M (one measurement error)
## SIMEX  $M_i*(t, \lambda,b)=M_i*(t)+N(0, (\lambda+1)\Sigma)$
##        $M_i(t, \lambda,b)=M_i(t)+N(0, \lambda\Sigma)$
simex=function(data.intcensor)
{
  data.intcensor.lball=rbind(
    ## SIMEX for M_i^* (no measurement error)
    simexfun(data.intcensor,yterm="ynew",lambdavecterm=lambdavec1,
             lambdaaddterm=1,groupterm=0),
    
    ## SIMEX for M_i; 
    simexfun(data.intcensor,yterm="yerr",lambdavecterm=lambdavec2,
             lambdaaddterm=0,groupterm=1) )
  
  ## replace the missing value "NA" into 0
  data.intcensor.lball[is.na(data.intcensor.lball)]=c(0)
  
  return(data.intcensor.lball)
} ## end simex

## normal random effect a=N(0,1), Z_i(t)=1, Zi(t) does not change by time t
## M^*(t)=beta.0+beta.1*t+beta.2*x+a*1+Mvrnorm(0_K,Sigma_e)
## cumulative hazard Lambda(t)=0.2*t*exp(x+0.5*M(t))
## no theta.t*t for model identifiability issue
## we assume exp() is constant within (t_{ij-1},t_{ij}) 
## update M^*(t) based on T_i, M^*(t)=M^*(t)+gamma*I(t>T_i)*(t-T_i)
## when generate data, we assume we know T_i, when fit the data, we assume 
## we do not know T_i
###################################################
# Data generation function
# Lambda.0 is our assumed baseline hazard,status 
## par are parameters related to the baseline hazard
## n.id denotes sample size (# of individual)
## sigma_x denotes fixed effect x's variance in generating step
## sigma_a denotes random effect a's variance
## sigma_e denotes measurement error's variance
## beta.0 denotes coefficients of intercept
## beta.1 denotes coefficients of time t
## beta.2 denotes coefficients of fixed effect x
## gamma.1 denotes coefficient of time effect after event (t-T)
## theta.x denotes coefficient of fixed effect x in survival model
## theta.a denotes coefficient of random effect a in survival model
## theta.m denotes coefficient of biomarker effect M* in survival model
## lambda_0 denotes baseline hazard in survival model (constant)
## J denotes the number of observations for each patients
## p denotes dim of fixed effect x
## K denotes dim of biomarker M*

##generating data with 1-dim biomarker and 1-dim random effect
data.gen=function(lambdavec,B,beta.0,beta.1,beta.2,theta.x,theta.a,theta.m,
                  gamma.1,K,p,Sigma.x,Sigma.a,Sigma.e,sigma.errlist,n.id,J,par) 
{

  ## variable with n.id size:
  ## x,a,y.obsnew,y.obs,prey.obsnew
  ## x.time,trueevent.time,event.time,preevent.time,censor.time,status
  ## variable with n.id*J size:
  ## xval,aval,LambdasVal,y.err1-4,prey.err1-4,
  ## y,ynew,preynew,previsittime,visittime,trueendtime,preetime,etime,ctime,event,id
  
  ## Step 1: generate visitnum (include time 0)
  ## censored time is after time t0, and in (max(3,t0),J)
  ## censorevent=rmultinom(n.id, size=1, prob=seq(1,J-3,length.out=J-max(3,floor(t0))))
  ## censorevent=rmultinom(n.id, size=1, prob=c(1,2,5))
  censorevent=rmultinom(n.id, size=1, prob=c(0,0,1))
  ## censorevent=rmultinom(n.id, size=1, prob=c(rep(0,J-floor(t0)-2),1,5))
  
  ## visitnum=unlist(lapply(1:n.id,function(k) which(censorevent[,k]==1)+max(3,floor(t0)) ))
  visitnum=unlist(lapply(1:n.id,function(k) which(censorevent[,k]==1)+J-3))
  ## if all censor at Jth visit
  ## visitnum=rep(J,n.id)
  
  # Step 2: observational times for longitudinal markers: 
  ## xtime[i]={t_1=0,t_2,t_2,...,t_{J-1}} (first include all times even after censored)
  ## to avoid bias
  random.time=lapply(1:n.id,function(k) (1:(J-1))+runif(J-1,-0.2,0.2))
  x.time=lapply(1:n.id, function (k) c(0,random.time[[k]]))
  ## censor.time=unlist(lapply(1:n.id, function(k) tail(x.time[[k]],n=1)))
  censor.time=unlist(lapply(1:n.id, function(k) x.time[[k]][visitnum[k]]))
  
  ## Step 3: generate original lm values   
  # baseline covariate; Sigma.x is covariance
  x=t(mvrnorm(n.id,mu=rep(0,p),Sigma=Sigma.x))
  
  # random effect
  a=t(mvrnorm(n.id,mu=rep(0,K),Sigma=Sigma.a))
  
  ## lm value at time 0
  mean.1=beta.0%*%rep(1, n.id)+beta.2%*%x+a
  
  ## K rows, n.id*J columns; lm value
  ## y is original value (no simex) without measurement error
  ## ynew is no meansurement error with gamma update
  errnum=length(sigma.errlist) 
  y=matrix(rep(0, K*errnum*n.id*J),nrow=K*errnum)
  ynew=y
  rownames(y)=paste("y",rep(1:errnum,each=K),"_",1:K,sep="")
  rownames(ynew)=paste("ynew",rep(1:errnum,each=K),"_",1:K,sep="")
  
  ## y.err is original value (no simex) with measurement error 
  y.err=matrix(rep(0, K*errnum*n.id*J),nrow=K*errnum)
  rownames(y.err)=paste("yerr",rep(1:errnum,each=K),"_",1:K,sep="")
  
  LambdasVal=rep(0, n.id*J)
  ## covariate has p rows,n.id*J columns
  xval=matrix(rep(0, p*n.id*J),nrow=p)
  aval=matrix(rep(0, K*n.id*J),nrow=K)
  rownames(aval)=paste("aval_",1:K,sep="")
  id=rep(0, n.id*J)
  
  ## duplicate J times for one subject (x.time)
  visittime=rep(0,n.id*J)
  
  ## duplicate J times for one subject
  ctime=rep(0,n.id*J)
  
  ## event happening time
  trueevent.time=event.time=preevent.time=rep(0,n.id)
  
  ## duplicate J times for one subject
  ## followup.time=min(trueevent.time,censor.time)
  etime=preetime=followuptime=rep(0,n.id*J)
  
  ## status=I(trueetime<=C), event=0 (before diagnosis),1 (on),
  ## 2 (aftter diagnosis)
  status=event=rep(0, n.id*J)
  ## visit=0,1,...,J-1; visits=# of visits include 0; befdiag includes 0
  ## afterdiag=sum_j t_j\geq U; visits=beforediag+afterdiag;
  visit=visits=befdiag=afterdiag=rep(0, n.id*J)
  
  ## quant time
  quanttime=matrix(rep(quant_d,n.id*J),nrow=length(quant_d))
  rownames(quanttime)=paste("quant",1:length(quant_d),sep="")
  
  ## generate n.id samples
  for (i in 1:n.id)
  {
    # generate mean of longitudinal data (with time effect)
    ## dim rows=K domains, columns=visitnum[i]+1 times
    mean.u=beta.1%*%x.time[[i]]+matrix(rep(mean.1[,i],J),ncol=J)
    # generate longitudinal data with random variation
    y.obsnew=y.obs=mean.u+t(mvrnorm(J,mu=rep(0,K),Sigma=Sigma.e))
    
    ## Step 4: Compute the cumulative hazard function for each individual 
    ## given the longitudinal data
    # Compute the multiplicative factors; factors dim=visitnum[i]*1
    factors=exp(as.vector(theta.x%*%x[,i])+
                  as.vector(theta.a%*%a[,i])+as.vector(theta.m%*%y.obs))  
    
    # Compute cumulative hazard at each time points x.time[[i]]
    # assume the factors are piecewise constant (only change at visittime)
    # int_0^3 lambda_0(t)exp(theta*t)dt=int 0^2 lambda_0(t)k_1 dt+
    # int 2^3 lambda_0(t)k_2 dt=k_1*(Lambda(2)-Lambda(0))+k_2*(Lambda(3)-Lambda(2))
    Lambdas=rep(0,J)  
    for (j in 2:J){
      Lambdas[j]=Lambdas[j-1]+
        (Lambda.0(x.time[[i]][j],par)-Lambda.0(x.time[[i]][j-1],par))*factors[j-1]
    } ## end j
    
    # generate event time (any event time has survival probability in (0,1))
    ## (generate U~Unif(0,1), tmp=-log(1-U)~EXP(1) )
    tmp=-log(runif(1))
    if (tmp>=Lambdas[J]){
      trueevent.time[i]=inv.Lambda.0( (tmp-Lambdas[J])/factors[J]+Lambda.0(x.time[[i]][J], par) ,par)
    } else {
      j=min(which(Lambdas>tmp))-1
      trueevent.time[i]=inv.Lambda.0( (tmp-Lambdas[j])/factors[j]+Lambda.0(x.time[[i]][j], par) , par)
    }
    
    if (trueevent.time[i]>censor.time[i]){ 
      ## the subject is censored before trueevent.time
      preevent.time[i]=censor.time[i]
      event.time[i]=trueevent.time[i]
      befdiag[((i-1)*J+1):(i*J)]=visitnum[i]
    } else { 
      ## event happens before or on censored time
      ## Step 5: update lm value after event time
      j=min(which(Lambdas>=tmp))-1
      preevent.time[i]=x.time[[i]][j]
      event.time[i]=x.time[[i]][j+1]
      
      ## event=0: no event; event=1: event at this time; 
      ## event=2: time after event
      event[(i-1)*J+j+1]=1
      befdiag[((i-1)*J+1):(i*J)]=j
      afterdiag[((i-1)*J+1):(i*J)]=visitnum[i]-j
      
      # update longitudinal marker value on and after index j+1;
      
      y.obsnew[,(j+1):visitnum[i]]=y.obs[,(j+1):visitnum[i]]+
        gamma.1%*%(x.time[[i]][(j+1):visitnum[i]]-trueevent.time[i])
      ## +t(matrix(rep(gamma.2%*%x[,i],visitnum[i]-j+1),nrow=visitnum[i]-j+1,byrow=TRUE)*
      ## (x.time[[i]][j:visitnum[i]]-trueevent.time[i]))
      
      if (j+1<visitnum[i]) 
      {
        event[((i-1)*J+j+2):(i*J)]=2
      }
      
      ## the difference between prey.obsnew and y.obsnew is just column 0
      ## prey.obsnew=cbind(prey.obsnew[,1],y.obsnew[,1:(visitnum[i]-1)])
    } ## end if (trueevent.time[i]
    
    ## save y and ynew
    ## generate measurement error data based on sigma_errlist
    for (sigma.err in sigma.errlist)
    {
      y[(1:K)+(which(sigma.errlist==sigma.err)-1)*K,((i-1)*J+1):(i*J)]=y.obs
      ynew[(1:K)+(which(sigma.errlist==sigma.err)-1)*K,((i-1)*J+1):(i*J)]=y.obsnew
      Sigma.err=Sigma_CS(1:K,c(sigma.err^2,0))
      y.obserr=t(mvrnorm(J,mu=rep(0,K),Sigma=Sigma.err))
      ## y.err is original value (no simex) with measurement error added
      y.err[(1:K)+(which(sigma.errlist==sigma.err)-1)*K,((i-1)*J+1):(i*J)]=
        y.obsnew+y.obserr
    }
    
    xval[,((i-1)*J+1): (i*J)]=x[,i]
    aval[,((i-1)*J+1): (i*J)]=a[,i]
    id[((i-1)*J+1): (i*J)]=i
    visit[((i-1)*J+1): (i*J)]=seq(1,J,1)
    ## visits include time 0
    visits[((i-1)*J+1): (i*J)]=visitnum[i]
    visittime[((i-1)*J+1): (i*J)]=x.time[[i]]
    etime[((i-1)*J+1): (i*J)]=event.time[i]
    preetime[((i-1)*J+1): (i*J)]=preevent.time[i]
    ctime[((i-1)*J+1): (i*J) ]=censor.time[i]
    ## status: event happened=1, censored=0
    status[((i-1)*J+1):(i*J)]=ifelse(event.time[i]<=censor.time[i],1,0)
    followuptime[((i-1)*J+1):(i*J)]=min(trueevent.time[i],censor.time[i])
    LambdasVal[((i-1)*J+1): (i*J)]=Lambdas
  } ## end i
  
  ## data.intcensor=data.frame(id,xval,y,aval,LambdasVal,
  ## preynew,ynew,prey.err,y.err,previsittime,visittime,trueendtime,
  ## preetime,etime,ctime,event)
  
  ## do matrix transpose
  xval=t(xval)
  aval=t(aval)
  quanttime=t(quanttime)
  ynew=t(ynew)
  y.err=t(y.err)
  
  data.intcensor=data.frame(id,visit,visits,befdiag,afterdiag,xval,aval,
                            preetime,etime,ctime,followuptime,status,quanttime,
                            ynew,y.err,visittime,event)
  
  ## data.intcensor=data.intcensor[is.na(data.intcensor$ynew1_1)==FALSE,]
  ## data.surv=data.intcensor[data.intcensor$event<=1,]
  
  # Generate SIMEX file
  data.intcensor.all=simex(data.intcensor)
  long.data = data.intcensor.all[data.intcensor.all$b==1&
                                   data.intcensor.all$group==0&
                                   data.intcensor.all$lambdaindex==1,]
  long.data <- long.data %>%
    select(-c(b, group, lambdaindex))
  long.data <- long.data %>%
    pivot_longer(
      cols = -c(id, visits, befdiag, afterdiag, xval, aval_1, preetime, etime,
                ctime, followuptime, status, quant1, quant2),
      names_to = c('.value','occasion'),
      names_sep = '_vis'
    )
  long.data$occasion <- as.numeric(long.data$occasion)
  return(long.data)
  ## return(list(data.intcensor=data.intcensor,
  ## data.surv=data.surv,status=status))
  
  ## bootstrap
  ## for (boot in 1:bootn){
  ## data.boot=NULL
  ## bsample=sort(sample(1:n.id,replace=T))
  ## for (bindex in 1:length(bsample))
  ## { 
  ## data.now=data.intcensor[which(data.intcensor$id==bsample[bindex]),]
  ## data.now$id=bindex
  ## data.boot=rbind(data.boot,data.now)
  ## }
  
  ## data.intcensor.bootlb=simex(data.boot)
  ## path_out='C:\\Users\\ltong\\Documents\\interval data simex\\' 
  ## filename.1=paste("rep=",(simseed-1)*n.rep+ii," boot=",boot,".csv.gz",sep="")
  ## write_csv(data.intcensor.bootlb, filename.1)
  ## } ## end boot
  
} ## end data.gen

#####################################
## Generalized date generation, including multiple biomarkers and fixed effects
## 11/05/2025: x time-dependent
data.gen.vec=function(parameters,n.id,J,par,Sigma.x) 
{
  ## unpack parameters
  beta.0 = t(as.matrix(parameters %>% dplyr::select(contains("beta.0"))))
  # beta.0 = parameters[paste("beta.0",1:K,sep="_"),]
  # K=length(beta.0)
  beta.1 = t(as.matrix(parameters %>% dplyr::select(contains("beta.1"))))
  # beta.1 = parameters[paste("beta.1",1:K,sep="_"),]
  beta.2.vec = data.frame(parameters %>% dplyr::select(contains("beta.2")))
  # index_grid <- expand.grid(row = 1:K, col = 1:p)
  # beta.2.vec = parameters[paste("beta.2",index_grid$row, index_grid$col,sep='_'),]
  col_names <- colnames(beta.2.vec)
  if (length(col_names)==1){
    beta.2 = as.matrix(beta.2.vec)
    p = ncol(beta.2)
  } else{
  indices <- do.call(rbind, strsplit(col_names, "_"))[, 2:3]
  indices <- apply(indices, 2, as.integer)
  beta.2 = matrix(0, nrow=max(indices[,1]),ncol=max(indices[,2]))
  p=ncol(beta.2)
  for (i in seq_along(col_names)) {
    row_idx <- indices[i, 1]
    col_idx <- indices[i, 2]
    beta.2[row_idx, col_idx] <- beta.2.vec[[col_names[i]]]
  }
  }
  gamma = t(as.matrix(parameters %>% dplyr::select(contains("gamma"))))
  
  # lambda0 = as.numeric(parameters %>% dplyr::select(contains("lambda0")))
  theta.x = as.matrix(parameters %>% dplyr::select(contains("theta.x")))
  theta.a = as.matrix(parameters %>% dplyr::select(contains("theta.a")))
  theta.m = as.matrix(parameters %>% dplyr::select(contains("theta.m")))
  
  Sigma.a = matrix(0, nrow=K,ncol=K)
  upper_vec_a = as.matrix(parameters %>% dplyr::select(contains("Sigma.a")))
  Sigma.a[upper.tri(Sigma.a,diag = TRUE)] = upper_vec_a
  if (K>1){
    Sigma.a=Sigma.a+t(Sigma.a)-diag(diag(Sigma.a))
  }
  Sigma.e = matrix(0, nrow=K,ncol=K)
  upper_vec_e = as.matrix(parameters %>% dplyr::select(contains("Sigma.e")))
  Sigma.e[upper.tri(Sigma.e,diag = TRUE)] = upper_vec_e
  if (K>1){
  Sigma.e=Sigma.e+t(Sigma.e)-diag(diag(Sigma.e))
  }
  ## variable with n.id size:
  ## x,a,y.obsnew,y.obs,prey.obsnew
  ## x.time,trueevent.time,event.time,preevent.time,censor.time,status
  ## variable with n.id*J size:
  ## xval,aval,LambdasVal,y.err1-4,prey.err1-4,
  ## y,ynew,preynew,previsittime,visittime,trueendtime,preetime,etime,ctime,event,id
  
  ## Step 1: generate visitnum (include time 0)
  ## censored time is after time t0, and in (max(3,t0),J)
  ## censorevent=rmultinom(n.id, size=1, prob=seq(1,J-3,length.out=J-max(3,floor(t0))))
  ## censorevent=rmultinom(n.id, size=1, prob=c(1,2,5))
  censorevent=rmultinom(n.id, size=1, prob=c(0,0,1))
  ## censorevent=rmultinom(n.id, size=1, prob=c(rep(0,J-floor(t0)-2),1,5))
  
  ## visitnum=unlist(lapply(1:n.id,function(k) which(censorevent[,k]==1)+max(3,floor(t0)) ))
  visitnum=unlist(lapply(1:n.id,function(k) which(censorevent[,k]==1)+J-3))
  ## if all censor at Jth visit
  ## visitnum=rep(J,n.id)
  
  # Step 2: observational times for longitudinal markers: 
  ## xtime[i]={t_1=0,t_2,t_2,...,t_{J-1}} (first include all times even after censored)
  ## to avoid bias
  random.time=lapply(1:n.id,function(k) (1:(J-1))+runif(J-1,-0.2,0.2))
  x.time=lapply(1:n.id, function (k) c(0,random.time[[k]]))
  ## censor.time=unlist(lapply(1:n.id, function(k) tail(x.time[[k]],n=1)))
  censor.time=unlist(lapply(1:n.id, function(k) x.time[[k]][visitnum[k]]))
  
  ## Step 3: generate original lm values   
  # baseline covariate; Sigma.x is covariance
  # x time-dependent
  # x=t(mvrnorm(n.id*J,mu=rep(0,p),Sigma=Sigma.x))
  x=lapply(1:n.id, function(k) rep(mvrnorm(1, mu=rep(0,p), Sigma = Sigma.x),J))
  # random effect
  a=t(mvrnorm(n.id,mu=rep(0,K),Sigma=Sigma.a))
  
  ## lm value at time 0
  mean.1=beta.0%*%rep(1, n.id)+a #+beta.2%*%x
  
  ## K rows, n.id*J columns; lm value
  ## y is original value (no simex) without measurement error
  ## ynew is no meansurement error with gamma update
  # errnum=length(sigma.errlist) 
  y=matrix(rep(0, K*n.id*J),nrow=K)
  ynew=y
  rownames(y)=paste("y","_",1:K,sep="")
  rownames(ynew)=paste("ynew","_",1:K,sep="")
  
  ## y.err is original value (no simex) with measurement error 
  y.err=matrix(rep(0, K*n.id*J),nrow=K)
  rownames(y.err)=paste("yerr","_",1:K,sep="")
  
  LambdasVal=rep(0, n.id*J)
  ## covariate has p rows,n.id*J columns
  xval=matrix(rep(0, p*n.id*J),nrow=p)
  rownames(xval)=paste("xval_",1:p,sep='')
  aval=matrix(rep(0, K*n.id*J),nrow=K)
  rownames(aval)=paste("aval_",1:K,sep="")
  id=rep(0, n.id*J)
  
  ## duplicate J times for one subject (x.time)
  visittime=rep(0,n.id*J)
  
  ## duplicate J times for one subject
  ctime=rep(0,n.id*J)
  
  ## event happening time
  trueevent.time=event.time=preevent.time=rep(0,n.id)
  
  ## duplicate J times for one subject
  ## followup.time=min(trueevent.time,censor.time)
  etime=preetime=followuptime=rep(0,n.id*J)
  
  ## status=I(trueetime<=C), event=0 (before diagnosis),1 (on),
  ## 2 (aftter diagnosis)
  status=event=rep(0, n.id*J)
  ## visit=0,1,...,J-1; visits=# of visits include 0; befdiag includes 0
  ## afterdiag=sum_j t_j\geq U; visits=beforediag+afterdiag;
  visit=visits=befdiag=afterdiag=rep(0, n.id*J)
  
  ## quant time
  # quanttime=matrix(rep(quant_d,n.id*J),nrow=length(quant_d))
  # rownames(quanttime)=paste("quant",1:length(quant_d),sep="")
  
  ## generate n.id samples
  for (i in 1:n.id)
  {
    # generate mean of longitudinal data (with time effect)
    ## dim rows=K domains, columns=visitnum[i]+1 times
    mean.u=beta.1%*%x.time[[i]]+matrix(rep(mean.1[,i],J),ncol=J)+beta.2%*%t(x[[i]])
    # generate longitudinal data with random variation
    y.obsnew=y.obs=mean.u+t(mvrnorm(J,mu=rep(0,K),Sigma=Sigma.e))
    
    ## Step 4: Compute the cumulative hazard function for each individual 
    ## given the longitudinal data
    # Compute the multiplicative factors; factors dim=visitnum[i]*1
    factors=exp(as.vector(theta.x%*%t(x[[i]]))+
                  as.vector(theta.a%*%a[,i])+as.vector(theta.m%*%y.obs))  
    
    # Compute cumulative hazard at each time points x.time[[i]]
    # assume the factors are piecewise constant (only change at visittime)
    # int_0^3 lambda_0(t)exp(theta*t)dt=int 0^2 lambda_0(t)k_1 dt+
    # int 2^3 lambda_0(t)k_2 dt=k_1*(Lambda(2)-Lambda(0))+k_2*(Lambda(3)-Lambda(2))
    Lambdas=rep(0,J)  
    for (j in 2:J){
      Lambdas[j]=Lambdas[j-1]+
        (Lambda.0(x.time[[i]][j],par)-Lambda.0(x.time[[i]][j-1],par))*factors[j-1]
    } ## end j
    
    # generate event time (any event time has survival probability in (0,1))
    ## (generate U~Unif(0,1), tmp=-log(1-U)~EXP(1) )
    tmp=-log(runif(1))
    if (tmp>=Lambdas[J]){
      trueevent.time[i]=inv.Lambda.0( (tmp-Lambdas[J])/factors[J]+Lambda.0(x.time[[i]][J], par) , par)
    } else {
      j=min(which(Lambdas>tmp))-1
      trueevent.time[i]=inv.Lambda.0( (tmp-Lambdas[j])/factors[j]+Lambda.0(x.time[[i]][j], par) , par)
    }
    
    if (trueevent.time[i]>censor.time[i]){ 
      ## the subject is censored before trueevent.time
      preevent.time[i]=censor.time[i]
      event.time[i]=trueevent.time[i]
      befdiag[((i-1)*J+1):(i*J)]=visitnum[i]
    } else { 
      ## event happens before or on censored time
      ## Step 5: update lm value after event time
      j=min(which(Lambdas>=tmp))-1
      preevent.time[i]=x.time[[i]][j]
      event.time[i]=x.time[[i]][j+1]
      
      ## event=0: no event; event=1: event at this time; 
      ## event=2: time after event
      event[(i-1)*J+j+1]=1
      befdiag[((i-1)*J+1):(i*J)]=j
      afterdiag[((i-1)*J+1):(i*J)]=visitnum[i]-j
      
      # update longitudinal marker value on and after index j+1;
      
      y.obsnew[,(j+1):visitnum[i]]=y.obs[,(j+1):visitnum[i]]+
        gamma%*%t(x.time[[i]][(j+1):visitnum[i]]-trueevent.time[i])
      ## +t(matrix(rep(gamma.2%*%x[,i],visitnum[i]-j+1),nrow=visitnum[i]-j+1,byrow=TRUE)*
      ## (x.time[[i]][j:visitnum[i]]-trueevent.time[i]))
      
      if (j+1<visitnum[i]) 
      {
        event[((i-1)*J+j+2):(i*J)]=2
      }
      
      ## the difference between prey.obsnew and y.obsnew is just column 0
      ## prey.obsnew=cbind(prey.obsnew[,1],y.obsnew[,1:(visitnum[i]-1)])
    } ## end if (trueevent.time[i]
    
    ## save y and ynew
    ## generate measurement error data based on sigma_errlist
    

      y[(1:K),((i-1)*J+1):(i*J)]=y.obs
      ynew[(1:K),((i-1)*J+1):(i*J)]=y.obsnew
      # Sigma.err=Sigma_CS(1:K,c(sigma.err^2,0))
      # y.obserr=t(mvrnorm(J,mu=rep(0,K),Sigma=Sigma.e))
      ## y.err is original value (no simex) with measurement error added
      # y.err[(1:K),((i-1)*J+1):(i*J)]=y.obsnew+y.obserr

    
    xval[,((i-1)*J+1): (i*J)]=t(x[[i]])
    aval[,((i-1)*J+1): (i*J)]=a[,i]
    id[((i-1)*J+1): (i*J)]=i
    visit[((i-1)*J+1): (i*J)]=seq(1,J,1)
    ## visits include time 0
    visits[((i-1)*J+1): (i*J)]=visitnum[i]
    visittime[((i-1)*J+1): (i*J)]=x.time[[i]]
    etime[((i-1)*J+1): (i*J)]=event.time[i]
    preetime[((i-1)*J+1): (i*J)]=preevent.time[i]
    ctime[((i-1)*J+1): (i*J) ]=censor.time[i]
    ## status: event happened=1, censored=0
    status[((i-1)*J+1):(i*J)]=ifelse(event.time[i]<=censor.time[i],1,0)
    followuptime[((i-1)*J+1):(i*J)]=min(trueevent.time[i],censor.time[i])
    LambdasVal[((i-1)*J+1): (i*J)]=Lambdas
  } ## end i
  
  ## data.intcensor=data.frame(id,xval,y,aval,LambdasVal,
  ## preynew,ynew,prey.err,y.err,previsittime,visittime,trueendtime,
  ## preetime,etime,ctime,event)
  
  ## do matrix transpose
  xval=t(xval)
  aval=t(aval)
  # quanttime=t(quanttime)
  ynew=t(ynew)
  y.err=t(y.err)
  trueevent.time = kronecker(trueevent.time,rep(1,J))
  data.intercensor=data.frame(id,visit,visits,befdiag,afterdiag,xval,aval,
                            preetime,etime,ctime,followuptime,status,
                            ynew,visittime,event, trueevent.time)
  return(data.intercensor)
  ## data.intcensor=data.intcensor[is.na(data.intcensor$ynew1_1)==FALSE,]
  ## data.surv=data.intcensor[data.intcensor$event<=1,]
  
  # Generate SIMEX file
  # data.intcensor.all=simex(data.intcensor)
  # long.data = data.intcensor.all[data.intcensor.all$b==1&
  #                                  data.intcensor.all$group==0&
  #                                  data.intcensor.all$lambdaindex==1,]
  # long.data <- long.data %>%
  #   select(-c(b, group, lambdaindex))
  # long.data <- long.data %>%
  #   pivot_longer(
  #     cols = -c(id, visits, befdiag, afterdiag, xval, aval_1, preetime, etime,
  #               ctime, followuptime, status, quant1, quant2),
  #     names_to = c('.value','occasion'),
  #     names_sep = '_vis'
  #   )
  # long.data$occasion <- as.numeric(long.data$occasion)
  # return(long.data)
  ## return(list(data.intcensor=data.intcensor,
  ## data.surv=data.surv,status=status))
  
  ## bootstrap
  ## for (boot in 1:bootn){
  ## data.boot=NULL
  ## bsample=sort(sample(1:n.id,replace=T))
  ## for (bindex in 1:length(bsample))
  ## { 
  ## data.now=data.intcensor[which(data.intcensor$id==bsample[bindex]),]
  ## data.now$id=bindex
  ## data.boot=rbind(data.boot,data.now)
  ## }
  
  ## data.intcensor.bootlb=simex(data.boot)
  ## path_out='C:\\Users\\ltong\\Documents\\interval data simex\\' 
  ## filename.1=paste("rep=",(simseed-1)*n.rep+ii," boot=",boot,".csv.gz",sep="")
  ## write_csv(data.intcensor.bootlb, filename.1)
  ## } ## end boot
  
} ## end data.gen

#######################################
data.gen.vec.error=function(parameters,n.id,J,par,Sigma.x, Sigma.err, lambdavec) 
{
  ## unpack parameters
  beta.0 = t(as.matrix(parameters %>% dplyr::select(contains("beta.0"))))
  # beta.0 = parameters[paste("beta.0",1:K,sep="_"),]
  # K=length(beta.0)
  beta.1 = t(as.matrix(parameters %>% dplyr::select(contains("beta.1"))))
  # beta.1 = parameters[paste("beta.1",1:K,sep="_"),]
  beta.2.vec = data.frame(parameters %>% dplyr::select(contains("beta.2")))
  # index_grid <- expand.grid(row = 1:K, col = 1:p)
  # beta.2.vec = parameters[paste("beta.2",index_grid$row, index_grid$col,sep='_'),]
  col_names <- colnames(beta.2.vec)
  if (length(col_names)==1){
    beta.2 = as.matrix(beta.2.vec)
    p = ncol(beta.2)
  } else{
    indices <- do.call(rbind, strsplit(col_names, "_"))[, 2:3]
    indices <- apply(indices, 2, as.integer)
    beta.2 = matrix(0, nrow=max(indices[,1]),ncol=max(indices[,2]))
    p=ncol(beta.2)
    for (i in seq_along(col_names)) {
      row_idx <- indices[i, 1]
      col_idx <- indices[i, 2]
      beta.2[row_idx, col_idx] <- beta.2.vec[[col_names[i]]]
    }
  }
  gamma = t(as.matrix(parameters %>% dplyr::select(contains("gamma"))))
  
  # lambda0 = as.numeric(parameters %>% dplyr::select(contains("lambda0")))
  theta.x = as.matrix(parameters %>% dplyr::select(contains("theta.x")))
  theta.a = as.matrix(parameters %>% dplyr::select(contains("theta.a")))
  theta.m = as.matrix(parameters %>% dplyr::select(contains("theta.m")))
  
  Sigma.a = matrix(0, nrow=K,ncol=K)
  upper_vec_a = as.matrix(parameters %>% dplyr::select(contains("Sigma.a")))
  Sigma.a[upper.tri(Sigma.a,diag = TRUE)] = upper_vec_a
  if (K>1){
    Sigma.a=Sigma.a+t(Sigma.a)-diag(diag(Sigma.a))
  }
  Sigma.e = matrix(0, nrow=K,ncol=K)
  upper_vec_e = as.matrix(parameters %>% dplyr::select(contains("Sigma.e")))
  Sigma.e[upper.tri(Sigma.e,diag = TRUE)] = upper_vec_e
  if (K>1){
    Sigma.e=Sigma.e+t(Sigma.e)-diag(diag(Sigma.e))
  }
  ## variable with n.id size:
  ## x,a,y.obsnew,y.obs,prey.obsnew
  ## x.time,trueevent.time,event.time,preevent.time,censor.time,status
  ## variable with n.id*J size:
  ## xval,aval,LambdasVal,y.err1-4,prey.err1-4,
  ## y,ynew,preynew,previsittime,visittime,trueendtime,preetime,etime,ctime,event,id
  
  ## Step 1: generate visitnum (include time 0)
  ## censored time is after time t0, and in (max(3,t0),J)
  ## censorevent=rmultinom(n.id, size=1, prob=seq(1,J-3,length.out=J-max(3,floor(t0))))
  ## censorevent=rmultinom(n.id, size=1, prob=c(1,2,5))
  censorevent=rmultinom(n.id, size=1, prob=c(0,0,1))
  ## censorevent=rmultinom(n.id, size=1, prob=c(rep(0,J-floor(t0)-2),1,5))
  
  ## visitnum=unlist(lapply(1:n.id,function(k) which(censorevent[,k]==1)+max(3,floor(t0)) ))
  visitnum=unlist(lapply(1:n.id,function(k) which(censorevent[,k]==1)+J-3))
  ## if all censor at Jth visit
  ## visitnum=rep(J,n.id)
  
  # Step 2: observational times for longitudinal markers: 
  ## xtime[i]={t_1=0,t_2,t_2,...,t_{J-1}} (first include all times even after censored)
  ## to avoid bias
  random.time=lapply(1:n.id,function(k) (1:(J-1))+runif(J-1,-0.2,0.2))
  x.time=lapply(1:n.id, function (k) c(0,random.time[[k]]))
  ## censor.time=unlist(lapply(1:n.id, function(k) tail(x.time[[k]],n=1)))
  censor.time=unlist(lapply(1:n.id, function(k) x.time[[k]][visitnum[k]]))
  
  ## Step 3: generate original lm values   
  # baseline covariate; Sigma.x is covariance
  # x time-dependent
  # x=t(mvrnorm(n.id*J,mu=rep(0,p),Sigma=Sigma.x))
  x=lapply(1:n.id, function(k) rep(mvrnorm(1, mu=rep(0,p), Sigma = Sigma.x),J))
  # random effect
  a=t(mvrnorm(n.id,mu=rep(0,K),Sigma=Sigma.a))
  
  ## lm value at time 0
  mean.1=beta.0%*%rep(1, n.id)+a #+beta.2%*%x
  
  ## K rows, n.id*J columns; lm value
  ## y is original value (no simex) without measurement error
  ## ynew is no meansurement error with gamma update
  # errnum=length(sigma.errlist) 
  y=matrix(rep(0, K*n.id*J),nrow=K)
  ynew=y
  rownames(y)=paste("y","_",1:K,sep="")
  rownames(ynew)=paste("ynew","_",1:K,sep="")
  
  ## y.err is original value (no simex) with measurement error 
  y.err=matrix(rep(0, K*n.id*J),nrow=K)
  rownames(y.err)=paste("yerr","_",1:K,sep="")
  
  LambdasVal=rep(0, n.id*J)
  ## covariate has p rows,n.id*J columns
  xval=matrix(rep(0, p*n.id*J),nrow=p)
  rownames(xval)=paste("xval_",1:p,sep='')
  aval=matrix(rep(0, K*n.id*J),nrow=K)
  rownames(aval)=paste("aval_",1:K,sep="")
  id=rep(0, n.id*J)
  
  ## duplicate J times for one subject (x.time)
  visittime=rep(0,n.id*J)
  
  ## duplicate J times for one subject
  ctime=rep(0,n.id*J)
  
  ## event happening time
  trueevent.time=event.time=preevent.time=rep(0,n.id)
  
  ## duplicate J times for one subject
  ## followup.time=min(trueevent.time,censor.time)
  etime=preetime=followuptime=rep(0,n.id*J)
  
  ## status=I(trueetime<=C), event=0 (before diagnosis),1 (on),
  ## 2 (aftter diagnosis)
  status=event=rep(0, n.id*J)
  ## visit=0,1,...,J-1; visits=# of visits include 0; befdiag includes 0
  ## afterdiag=sum_j t_j\geq U; visits=beforediag+afterdiag;
  visit=visits=befdiag=afterdiag=rep(0, n.id*J)
  
  ## quant time
  # quanttime=matrix(rep(quant_d,n.id*J),nrow=length(quant_d))
  # rownames(quanttime)=paste("quant",1:length(quant_d),sep="")
  
  ## generate n.id samples
  for (i in 1:n.id)
  {
    # generate mean of longitudinal data (with time effect)
    ## dim rows=K domains, columns=visitnum[i]+1 times
    mean.u=beta.1%*%x.time[[i]]+matrix(rep(mean.1[,i],J),ncol=J)+beta.2%*%t(x[[i]])
    # generate longitudinal data with random variation
    y.obsnew=y.obs=mean.u+t(mvrnorm(J,mu=rep(0,K),Sigma=Sigma.e))
    
    ## Step 4: Compute the cumulative hazard function for each individual 
    ## given the longitudinal data
    # Compute the multiplicative factors; factors dim=visitnum[i]*1
    factors=exp(as.vector(theta.x%*%t(x[[i]]))+
                  as.vector(theta.a%*%a[,i])+as.vector(theta.m%*%y.obs))  
    
    # Compute cumulative hazard at each time points x.time[[i]]
    # assume the factors are piecewise constant (only change at visittime)
    # int_0^3 lambda_0(t)exp(theta*t)dt=int 0^2 lambda_0(t)k_1 dt+
    # int 2^3 lambda_0(t)k_2 dt=k_1*(Lambda(2)-Lambda(0))+k_2*(Lambda(3)-Lambda(2))
    Lambdas=rep(0,J)  
    for (j in 2:J){
      Lambdas[j]=Lambdas[j-1]+
        (Lambda.0(x.time[[i]][j],par)-Lambda.0(x.time[[i]][j-1],par))*factors[j-1]
    } ## end j
    
    # generate event time (any event time has survival probability in (0,1))
    ## (generate U~Unif(0,1), tmp=-log(1-U)~EXP(1) )
    tmp=-log(runif(1))
    if (tmp>=Lambdas[J]){
      trueevent.time[i]=inv.Lambda.0( (tmp-Lambdas[J])/factors[J]+Lambda.0(x.time[[i]][J], par) , par)
    } else {
      j=min(which(Lambdas>tmp))-1
      trueevent.time[i]=inv.Lambda.0( (tmp-Lambdas[j])/factors[j]+Lambda.0(x.time[[i]][j], par) , par)
    }
    
    if (trueevent.time[i]>censor.time[i]){ 
      ## the subject is censored before trueevent.time
      preevent.time[i]=censor.time[i]
      event.time[i]=trueevent.time[i]
      befdiag[((i-1)*J+1):(i*J)]=visitnum[i]
    } else { 
      ## event happens before or on censored time
      ## Step 5: update lm value after event time
      j=min(which(Lambdas>=tmp))-1
      preevent.time[i]=x.time[[i]][j]
      event.time[i]=x.time[[i]][j+1]
      
      ## event=0: no event; event=1: event at this time; 
      ## event=2: time after event
      event[(i-1)*J+j+1]=1
      befdiag[((i-1)*J+1):(i*J)]=j
      afterdiag[((i-1)*J+1):(i*J)]=visitnum[i]-j
      
      # update longitudinal marker value on and after index j+1;
      
      y.obsnew[,(j+1):visitnum[i]]=y.obs[,(j+1):visitnum[i]]+
        gamma%*%t(x.time[[i]][(j+1):visitnum[i]]-trueevent.time[i])
      ## +t(matrix(rep(gamma.2%*%x[,i],visitnum[i]-j+1),nrow=visitnum[i]-j+1,byrow=TRUE)*
      ## (x.time[[i]][j:visitnum[i]]-trueevent.time[i]))
      
      if (j+1<visitnum[i]) 
      {
        event[((i-1)*J+j+2):(i*J)]=2
      }
      
      ## the difference between prey.obsnew and y.obsnew is just column 0
      ## prey.obsnew=cbind(prey.obsnew[,1],y.obsnew[,1:(visitnum[i]-1)])
    } ## end if (trueevent.time[i]
    
    ## save y and ynew
    ## generate measurement error data based on sigma_errlist
    
    
    y[(1:K),((i-1)*J+1):(i*J)]=y.obs
    ynew[(1:K),((i-1)*J+1):(i*J)]=y.obsnew
    # Sigma.err=Sigma_CS(1:K,c(sigma.err^2,0))
    # y.obserr=t(mvrnorm(J,mu=rep(0,K),Sigma=Sigma.e))
    ## y.err is original value (no simex) with measurement error added
    # y.err[(1:K),((i-1)*J+1):(i*J)]=y.obsnew+y.obserr
    
    
    xval[,((i-1)*J+1): (i*J)]=t(x[[i]])
    aval[,((i-1)*J+1): (i*J)]=a[,i]
    id[((i-1)*J+1): (i*J)]=i
    visit[((i-1)*J+1): (i*J)]=seq(1,J,1)
    ## visits include time 0
    visits[((i-1)*J+1): (i*J)]=visitnum[i]
    visittime[((i-1)*J+1): (i*J)]=x.time[[i]]
    etime[((i-1)*J+1): (i*J)]=event.time[i]
    preetime[((i-1)*J+1): (i*J)]=preevent.time[i]
    ctime[((i-1)*J+1): (i*J) ]=censor.time[i]
    ## status: event happened=1, censored=0
    status[((i-1)*J+1):(i*J)]=ifelse(event.time[i]<=censor.time[i],1,0)
    followuptime[((i-1)*J+1):(i*J)]=min(trueevent.time[i],censor.time[i])
    LambdasVal[((i-1)*J+1): (i*J)]=Lambdas
  } ## end i
  
  ## data.intcensor=data.frame(id,xval,y,aval,LambdasVal,
  ## preynew,ynew,prey.err,y.err,previsittime,visittime,trueendtime,
  ## preetime,etime,ctime,event)
  
  ## do matrix transpose
  xval=t(xval)
  aval=t(aval)
  # quanttime=t(quanttime)
  ynew=t(ynew)
  y.err=t(y.err)
  trueevent.time = kronecker(trueevent.time,rep(1,J))
  data.intercensor=data.frame(id,visit,visits,befdiag,afterdiag,xval,aval,
                              preetime,etime,ctime,followuptime,status,
                              ynew,visittime,event, trueevent.time)
  return(data.intercensor)
  ## data.intcensor=data.intcensor[is.na(data.intcensor$ynew1_1)==FALSE,]
  ## data.surv=data.intcensor[data.intcensor$event<=1,]
  
  # Generate SIMEX file
  # data.intcensor.all=simex(data.intcensor)
  # long.data = data.intcensor.all[data.intcensor.all$b==1&
  #                                  data.intcensor.all$group==0&
  #                                  data.intcensor.all$lambdaindex==1,]
  # long.data <- long.data %>%
  #   select(-c(b, group, lambdaindex))
  # long.data <- long.data %>%
  #   pivot_longer(
  #     cols = -c(id, visits, befdiag, afterdiag, xval, aval_1, preetime, etime,
  #               ctime, followuptime, status, quant1, quant2),
  #     names_to = c('.value','occasion'),
  #     names_sep = '_vis'
  #   )
  # long.data$occasion <- as.numeric(long.data$occasion)
  # return(long.data)
  ## return(list(data.intcensor=data.intcensor,
  ## data.surv=data.surv,status=status))
  
  ## bootstrap
  ## for (boot in 1:bootn){
  ## data.boot=NULL
  ## bsample=sort(sample(1:n.id,replace=T))
  ## for (bindex in 1:length(bsample))
  ## { 
  ## data.now=data.intcensor[which(data.intcensor$id==bsample[bindex]),]
  ## data.now$id=bindex
  ## data.boot=rbind(data.boot,data.now)
  ## }
  
  ## data.intcensor.bootlb=simex(data.boot)
  ## path_out='C:\\Users\\ltong\\Documents\\interval data simex\\' 
  ## filename.1=paste("rep=",(simseed-1)*n.rep+ii," boot=",boot,".csv.gz",sep="")
  ## write_csv(data.intcensor.bootlb, filename.1)
  ## } ## end boot
  
} 