################################
## Try numerical differential method
diff.beta.0 = function(long.data, beta.0.hat=beta.0,
                       beta.1.hat=beta.1,
                       beta.2.hat=beta.2,
                       gamma.hat=gamma.1,
                       lambda0.hat=1/par[1],
                       theta.x.hat=theta.x,
                       theta.a.hat=theta.a,
                       theta.m.hat=theta.m,
                       sigma.e.hat=sigma.e,
                       sigma.a.hat=sigma.a,
                       delta = 1e-6){
  sigma.e.hat = as.matrix(sigma.e.hat)
  sigma.a.hat = as.matrix(sigma.a.hat)
  val1 = initial.likelihood(long.data,beta.0.hat+delta,beta.1.hat,beta.2.hat,gamma.hat,lambda0.hat,theta.x.hat,theta.a.hat,theta.m.hat,sigma.e.hat,sigma.a.hat)$ll.ind
  val2 = initial.likelihood(long.data,beta.0.hat-delta,beta.1.hat,beta.2.hat,gamma.hat,lambda0.hat,theta.x.hat,theta.a.hat,theta.m.hat,sigma.e.hat,sigma.a.hat)$ll.ind
  diff = (val1 - val2)/(2*delta)
  return(diff)
}

diff.beta.1 = function(long.data, beta.0.hat=beta.0,
                       beta.1.hat=beta.1,
                       beta.2.hat=beta.2,
                       gamma.hat=gamma.1,
                       lambda0.hat=1/par[1],
                       theta.x.hat=theta.x,
                       theta.a.hat=theta.a,
                       theta.m.hat=theta.m,
                       sigma.e.hat=sigma.e,
                       sigma.a.hat=sigma.a,
                       delta = 1e-6){
  sigma.e.hat = as.matrix(sigma.e.hat)
  sigma.a.hat = as.matrix(sigma.a.hat)
  val1 = initial.likelihood(long.data,beta.0.hat,beta.1.hat+delta,beta.2.hat,gamma.hat,lambda0.hat,theta.x.hat,theta.a.hat,theta.m.hat,sigma.e.hat,sigma.a.hat)$ll.ind
  val2 = initial.likelihood(long.data,beta.0.hat,beta.1.hat-delta,beta.2.hat,gamma.hat,lambda0.hat,theta.x.hat,theta.a.hat,theta.m.hat,sigma.e.hat,sigma.a.hat)$ll.ind
  diff = (val1 - val2)/(2*delta)
  return(diff)
}

diff.beta.2 = function(long.data, beta.0.hat=beta.0,
                       beta.1.hat=beta.1,
                       beta.2.hat=beta.2,
                       gamma.hat=gamma.1,
                       lambda0.hat=1/par[1],
                       theta.x.hat=theta.x,
                       theta.a.hat=theta.a,
                       theta.m.hat=theta.m,
                       sigma.e.hat=sigma.e,
                       sigma.a.hat=sigma.a,
                       delta = 1e-6){
  sigma.e.hat = as.matrix(sigma.e.hat)
  sigma.a.hat = as.matrix(sigma.a.hat)
  val1 = initial.likelihood(long.data,beta.0.hat,beta.1.hat,beta.2.hat+delta,gamma.hat,lambda0.hat,theta.x.hat,theta.a.hat,theta.m.hat,sigma.e.hat,sigma.a.hat)$ll.ind
  val2 = initial.likelihood(long.data,beta.0.hat,beta.1.hat,beta.2.hat-delta,gamma.hat,lambda0.hat,theta.x.hat,theta.a.hat,theta.m.hat,sigma.e.hat,sigma.a.hat)$ll.ind
  diff = (val1 - val2)/(2*delta)
  return(diff)
}

diff.gamma = function(long.data, beta.0.hat=beta.0,
                      beta.1.hat=beta.1,
                      beta.2.hat=beta.2,
                      gamma.hat=gamma.1,
                      lambda0.hat=1/par[1],
                      theta.x.hat=theta.x,
                      theta.a.hat=theta.a,
                      theta.m.hat=theta.m,
                      sigma.e.hat=sigma.e,
                      sigma.a.hat=sigma.a,
                      delta = 1e-6){
  sigma.e.hat = as.matrix(sigma.e.hat)
  sigma.a.hat = as.matrix(sigma.a.hat)
  val1 = initial.likelihood(long.data,beta.0.hat,beta.1.hat,beta.2.hat,gamma.hat+delta,lambda0.hat,theta.x.hat,theta.a.hat,theta.m.hat,sigma.e.hat,sigma.a.hat)$ll.ind
  val2 = initial.likelihood(long.data,beta.0.hat,beta.1.hat,beta.2.hat,gamma.hat-delta,lambda0.hat,theta.x.hat,theta.a.hat,theta.m.hat,sigma.e.hat,sigma.a.hat)$ll.ind
  diff = (val1 - val2)/(2*delta)
  return(diff)
}

diff.lambda0 = function(long.data, beta.0.hat=beta.0,
                        beta.1.hat=beta.1,
                        beta.2.hat=beta.2,
                        gamma.hat=gamma.1,
                        lambda0.hat=1/par[1],
                        theta.x.hat=theta.x,
                        theta.a.hat=theta.a,
                        theta.m.hat=theta.m,
                        sigma.e.hat=sigma.e,
                        sigma.a.hat=sigma.a,
                        delta = 1e-6){
  sigma.e.hat = as.matrix(sigma.e.hat)
  sigma.a.hat = as.matrix(sigma.a.hat)
  val1 = initial.likelihood(long.data,beta.0.hat,beta.1.hat,beta.2.hat,gamma.hat,lambda0.hat+delta,theta.x.hat,theta.a.hat,theta.m.hat,sigma.e.hat,sigma.a.hat)$ll.ind
  val2 = initial.likelihood(long.data,beta.0.hat,beta.1.hat,beta.2.hat,gamma.hat,lambda0.hat-delta,theta.x.hat,theta.a.hat,theta.m.hat,sigma.e.hat,sigma.a.hat)$ll.ind
  diff = (val1 - val2)/(2*delta)
  return(diff)
}

diff.theta.x = function(long.data, beta.0.hat=beta.0,
                        beta.1.hat=beta.1,
                        beta.2.hat=beta.2,
                        gamma.hat=gamma.1,
                        lambda0.hat=1/par[1],
                        theta.x.hat=theta.x,
                        theta.a.hat=theta.a,
                        theta.m.hat=theta.m,
                        sigma.e.hat=sigma.e,
                        sigma.a.hat=sigma.a,
                        delta = 1e-6){
  sigma.e.hat = as.matrix(sigma.e.hat)
  sigma.a.hat = as.matrix(sigma.a.hat)
  val1 = initial.likelihood(long.data,beta.0.hat,beta.1.hat,beta.2.hat,gamma.hat,lambda0.hat,theta.x.hat+delta,theta.a.hat,theta.m.hat,sigma.e.hat,sigma.a.hat)$ll.ind
  val2 = initial.likelihood(long.data,beta.0.hat,beta.1.hat,beta.2.hat,gamma.hat,lambda0.hat,theta.x.hat-delta,theta.a.hat,theta.m.hat,sigma.e.hat,sigma.a.hat)$ll.ind
  diff = (val1 - val2)/(2*delta)
  return(diff)
}

diff.theta.a = function(long.data, beta.0.hat=beta.0,
                        beta.1.hat=beta.1,
                        beta.2.hat=beta.2,
                        gamma.hat=gamma.1,
                        lambda0.hat=1/par[1],
                        theta.x.hat=theta.x,
                        theta.a.hat=theta.a,
                        theta.m.hat=theta.m,
                        sigma.e.hat=sigma.e,
                        sigma.a.hat=sigma.a,
                        delta = 1e-6){
  sigma.e.hat = as.matrix(sigma.e.hat)
  sigma.a.hat = as.matrix(sigma.a.hat)
  val1 = initial.likelihood(long.data,beta.0.hat,beta.1.hat,beta.2.hat,gamma.hat,lambda0.hat,theta.x.hat,theta.a.hat+delta,theta.m.hat,sigma.e.hat,sigma.a.hat)$ll.ind
  val2 = initial.likelihood(long.data,beta.0.hat,beta.1.hat,beta.2.hat,gamma.hat,lambda0.hat,theta.x.hat,theta.a.hat-delta,theta.m.hat,sigma.e.hat,sigma.a.hat)$ll.ind
  diff = (val1 - val2)/(2*delta)
  return(diff)
}

diff.theta.m = function(long.data, beta.0.hat=beta.0,
                        beta.1.hat=beta.1,
                        beta.2.hat=beta.2,
                        gamma.hat=gamma.1,
                        lambda0.hat=1/par[1],
                        theta.x.hat=theta.x,
                        theta.a.hat=theta.a,
                        theta.m.hat=theta.m,
                        sigma.e.hat=sigma.e,
                        sigma.a.hat=sigma.a,
                        delta = 1e-6){
  sigma.e.hat = as.matrix(sigma.e.hat)
  sigma.a.hat = as.matrix(sigma.a.hat)
  val1 = initial.likelihood(long.data,beta.0.hat,beta.1.hat,beta.2.hat,gamma.hat,lambda0.hat,theta.x.hat,theta.a.hat,theta.m.hat+delta,sigma.e.hat,sigma.a.hat)$ll.ind
  val2 = initial.likelihood(long.data,beta.0.hat,beta.1.hat,beta.2.hat,gamma.hat,lambda0.hat,theta.x.hat,theta.a.hat,theta.m.hat-delta,sigma.e.hat,sigma.a.hat)$ll.ind
  diff = (val1 - val2)/(2*delta)
  return(diff)
}

diff.sigma.e = function(long.data, beta.0.hat=beta.0,
                        beta.1.hat=beta.1,
                        beta.2.hat=beta.2,
                        gamma.hat=gamma.1,
                        lambda0.hat=1/par[1],
                        theta.x.hat=theta.x,
                        theta.a.hat=theta.a,
                        theta.m.hat=theta.m,
                        sigma.e.hat=sigma.e,
                        sigma.a.hat=sigma.a,
                        delta = 1e-6){
  sigma.e.hat = as.matrix(sigma.e.hat)
  sigma.a.hat = as.matrix(sigma.a.hat)
  val1 = initial.likelihood(long.data,beta.0.hat,beta.1.hat,beta.2.hat,gamma.hat,lambda0.hat,theta.x.hat,theta.a.hat,theta.m.hat,sigma.e.hat+delta,sigma.a.hat)$ll.ind
  val2 = initial.likelihood(long.data,beta.0.hat,beta.1.hat,beta.2.hat,gamma.hat,lambda0.hat,theta.x.hat,theta.a.hat,theta.m.hat,sigma.e.hat-delta,sigma.a.hat)$ll.ind
  diff = (val1 - val2)/(2*delta)
  return(diff)
}

diff.sigma.a = function(long.data, beta.0.hat=beta.0,
                        beta.1.hat=beta.1,
                        beta.2.hat=beta.2,
                        gamma.hat=gamma.1,
                        lambda0.hat=1/par[1],
                        theta.x.hat=theta.x,
                        theta.a.hat=theta.a,
                        theta.m.hat=theta.m,
                        sigma.e.hat=sigma.e,
                        sigma.a.hat=sigma.a,
                        delta = 1e-6){
  sigma.e.hat = as.matrix(sigma.e.hat)
  sigma.a.hat = as.matrix(sigma.a.hat)
  val1 = initial.likelihood(long.data,beta.0.hat,beta.1.hat,beta.2.hat,gamma.hat,lambda0.hat,theta.x.hat,theta.a.hat,theta.m.hat,sigma.e.hat,sigma.a.hat+delta)$ll.ind
  val2 = initial.likelihood(long.data,beta.0.hat,beta.1.hat,beta.2.hat,gamma.hat,lambda0.hat,theta.x.hat,theta.a.hat,theta.m.hat,sigma.e.hat,sigma.a.hat-delta)$ll.ind
  diff = (val1 - val2)/(2*delta)
  return(diff)
}

diff.total = function(long.data, beta.0.hat=beta.0,
                      beta.1.hat=beta.1,
                      beta.2.hat=beta.2,
                      gamma.hat=gamma.1,
                      lambda0.hat=1/par[1],
                      theta.x.hat=theta.x,
                      theta.a.hat=theta.a,
                      theta.m.hat=theta.m,
                      sigma.e.hat=sigma.e,
                      sigma.a.hat=sigma.a,
                      delta = 1e-6) {
  n.id = length(unique(long.data$id))
  Q = matrix(0, nrow = n.id, ncol = 10)
  diff.beta.0 = diff.beta.0(long.data, beta.0.hat,
                            beta.1.hat,
                            beta.2.hat,
                            gamma.hat,
                            lambda0.hat,
                            theta.x.hat,
                            theta.a.hat,
                            theta.m.hat,
                            sigma.e.hat,
                            sigma.a.hat)
  Q[, 1] = t(diff.beta.0)
  diff.beta.1 = diff.beta.1(long.data, beta.0.hat,
                            beta.1.hat,
                            beta.2.hat,
                            gamma.hat,
                            lambda0.hat,
                            theta.x.hat,
                            theta.a.hat,
                            theta.m.hat,
                            sigma.e.hat,
                            sigma.a.hat)
  Q[, 2] = t(diff.beta.1)
  diff.beta.2 = diff.beta.2(long.data, beta.0.hat,
                            beta.1.hat,
                            beta.2.hat,
                            gamma.hat,
                            lambda0.hat,
                            theta.x.hat,
                            theta.a.hat,
                            theta.m.hat,
                            sigma.e.hat,
                            sigma.a.hat)
  Q[, 3] = t(diff.beta.2)
  diff.gamma = diff.gamma(long.data, beta.0.hat,
                            beta.1.hat,
                            beta.2.hat,
                            gamma.hat,
                            lambda0.hat,
                            theta.x.hat,
                            theta.a.hat,
                            theta.m.hat,
                            sigma.e.hat,
                            sigma.a.hat)
  Q[, 4] = t(diff.gamma)
  diff.lambda0 = diff.lambda0(long.data, beta.0.hat,
                            beta.1.hat,
                            beta.2.hat,
                            gamma.hat,
                            lambda0.hat,
                            theta.x.hat,
                            theta.a.hat,
                            theta.m.hat,
                            sigma.e.hat,
                            sigma.a.hat)
  Q[, 5] = t(diff.lambda0)
  diff.theta.x = diff.theta.x(long.data, beta.0.hat,
                            beta.1.hat,
                            beta.2.hat,
                            gamma.hat,
                            lambda0.hat,
                            theta.x.hat,
                            theta.a.hat,
                            theta.m.hat,
                            sigma.e.hat,
                            sigma.a.hat)
  Q[, 6]= t(diff.theta.x)
  diff.theta.a = diff.theta.a(long.data, beta.0.hat,
                            beta.1.hat,
                            beta.2.hat,
                            gamma.hat,
                            lambda0.hat,
                            theta.x.hat,
                            theta.a.hat,
                            theta.m.hat,
                            sigma.e.hat,
                            sigma.a.hat)
  Q[, 7] = t(diff.theta.a)
  diff.theta.m = diff.theta.m(long.data, beta.0.hat,
                            beta.1.hat,
                            beta.2.hat,
                            gamma.hat,
                            lambda0.hat,
                            theta.x.hat,
                            theta.a.hat,
                            theta.m.hat,
                            sigma.e.hat,
                            sigma.a.hat)
  Q[, 8] = t(diff.theta.m)
  diff.sigma.e = diff.sigma.e(long.data, beta.0.hat,
                            beta.1.hat,
                            beta.2.hat,
                            gamma.hat,
                            lambda0.hat,
                            theta.x.hat,
                            theta.a.hat,
                            theta.m.hat,
                            sigma.e.hat,
                            sigma.a.hat)
  Q[, 9] = t(diff.sigma.e)
  diff.sigma.a = diff.sigma.a(long.data, beta.0.hat,
                            beta.1.hat,
                            beta.2.hat,
                            gamma.hat,
                            lambda0.hat,
                            theta.x.hat,
                            theta.a.hat,
                            theta.m.hat,
                            sigma.e.hat,
                            sigma.a.hat)
  Q[, 10] = t(diff.sigma.a)
  if (is.na(det(t(Q)%*%Q))||qr(t(Q)%*%Q)$rank<10){
    hessian = -inv(t(Q)%*%Q+1e-6*diag(10))
  }
  else {hessian = -inv(t(Q)%*%Q)}
  score = apply(Q, 2, sum)
  return(list(Score=score, Hessian=hessian))
}

initial.likelihood = function(long.data, beta.0.hat=beta.0,
                              beta.1.hat=beta.1,
                              beta.2.hat=beta.2,
                              gamma.hat=gamma.1,
                              lambda0.hat=lambda0,
                              theta.x.hat=theta.x,
                              theta.a.hat=theta.a,
                              theta.m.hat=theta.m,
                              sigma.e.hat=sigma.e,
                              sigma.a.hat=sigma.a)
{
  n.id = length(unique(long.data$id))
  Q=matrix(0,nrow = n.id, ncol = 1)
  # Qdbeta.0=matrix(0,n.id)
  # Qdbeta.1=matrix(0,n.id)
  # Qdbeta.2=matrix(0,n.id)
  # Qdgamma=matrix(0,n.id)
  # Qdlambda0=matrix(0,n.id)
  # Qdtheta.x=matrix(0,n.id)
  # Qdtheta.a=matrix(0,n.id)
  # Qdtheta.m=matrix(0,n.id)
  # Qdsigma.e=matrix(0,n.id)
  # Qdsigma.a=matrix(0,n.id)
  ## set GH nodes and weights(23 nodes)
  no.gauss=20
  out.gauss=gauss.quad(no.gauss, kind="hermit")
  nodes=as.matrix(out.gauss$nodes)
  weight=as.matrix(out.gauss$weights)
  sigma.e.hat = as.matrix(sigma.e.hat)
  sigma.a.hat = as.matrix(sigma.a.hat)
  a.re = nodes%*%sqrtm(2*sigma.a.hat)
  for (i in 1:n.id){
    ## for test
    # i = 1 
    data.i = long.data[long.data$id==i,]
    befdiag.i = data.i$befdiag[1]
    afterdiag.i = data.i$afterdiag[1]
    xval.i = data.i$xval[1]
    V.i = data.i$preetime[1]
    U.i = data.i$etime[1]
    status.i = data.i$status[1]
    visittime.i = data.i$visittime
    befdiagtime.i = visittime.i[data.i$occasion<=data.i$befdiag]
    ondiagtime.i = data.i$preetime[1]
    afterdiagtime.i = visittime.i[data.i$occasion>data.i$befdiag]
    visit_num.i = data.i$visits[1]
    y.i = data.i$ynew1_1
    befdiagy.i = y.i[data.i$occasion<data.i$befdiag]
    befondiagy.i = y.i[data.i$occasion<=data.i$befdiag]
    ondiagy.i = y.i[data.i$occasion==data.i$befdiag]
    afterdiagy.i = y.i[data.i$occasion>data.i$befdiag]
    A.i = rep(1, afterdiag.i)
    A.i.1 = rep(1, befdiag.i)
    if (status.i==1){
      a2 = apply(-(as.numeric(log(2*pi)+log(sigma.e.hat))*matrix(1,afterdiag.i,no.gauss)+
                     (matrix(1,afterdiag.i,no.gauss)*as.vector(afterdiagy.i-beta.0.hat-
                                                                 beta.1.hat*afterdiagtime.i-beta.2.hat*xval.i-gamma.hat*afterdiagtime.i)-A.i%*%t(a.re))^2/as.numeric(sigma.e.hat))/2,2,sum)+
        lambda0.hat*ondiagtime.i*exp(theta.x.hat*xval.i+theta.a.hat*t(a.re)+theta.m.hat*ondiagy.i)
      logp = log(lambda0.hat)+theta.x.hat*xval.i+theta.a.hat*t(a.re)+theta.m.hat*ondiagy.i-
                lambda0.hat*sum(diff(befdiagtime.i)*exp(theta.m.hat*befdiagy.i))*exp(theta.x.hat*xval.i+theta.a.hat*t(a.re))+a2
      boundlogp = rep(200, no.gauss)
      idx_na_p = which(logp[1,]==Inf)
      logp[, idx_na_p] = boundlogp[idx_na_p]
      b2 = apply(-(matrix(1,afterdiag.i,no.gauss)*(afterdiagy.i-beta.0.hat-beta.1.hat*afterdiagtime.i-beta.2.hat*xval.i-
                                                     gamma.hat*afterdiagtime.i)-A.i%*%t(a.re))*gamma.hat/as.numeric(sigma.e.hat),2,sum)-
        lambda0.hat*exp(theta.x.hat*xval.i+theta.a.hat*t(a.re)+theta.m.hat*ondiagy.i)
      c2 = apply(-matrix(1,afterdiag.i,no.gauss)*gamma.hat^2/2/as.numeric(sigma.e.hat),2,sum)
      optimal = -b2/(2*c2)
      idx_left = which(optimal[1,]<=V.i)
      idx_middle = which(optimal[1,]>V.i & optimal[1,]<U.i)
      idx_right = which(optimal[1,]>=U.i)
      bound1 = (b2*V.i+c2*V.i^2+logp)+log(U.i-V.i)
      bound2 = (logp-b2^2/(4*c2))+log(U.i-V.i)
      bound3 = (b2*U.i+c2*U.i^2+logp)+log(U.i-V.i)
      logpI_bound = rep(0, no.gauss)
      logpI_bound[idx_left] = bound1[,idx_left]
      logpI_bound[idx_middle] = bound2[,idx_middle]
      logpI_bound[idx_right] = bound3[, idx_right]
      logpI = log(sqrt(pi)*(pnorm(sqrt(-2*c2)*U.i-b2/sqrt(-2*c2),0,1)-pnorm(sqrt(-2*c2)*V.i-b2/sqrt(-2*c2),0,1)))+(logp-b2^2/(4*c2))-log(sqrt(-c2))
      logpI_apx = log(sqrt(pi))+(logp-b2^2/(4*c2))-log(sqrt(2*pi))+(-(sqrt(-2*c2)*(U.i+V.i)-2*b2/sqrt(-2*c2))^2/8)+log(sqrt(-2*c2)*(U.i-V.i))-log(sqrt(-c2))
      idx_zero <- which(!is.na(logpI[1, ]) & logpI[1, ] == -Inf)
      logpI[, idx_zero] <- logpI_apx[, idx_zero]
      idx_na <- which(is.na(logpI[1,])|is.infinite(logpI[1,]))
      logpI[, idx_na] <- logpI_bound[idx_na]
      bound_final = rep(-Inf, no.gauss)
      idx_na <- which(is.na(logpI[1,])|is.infinite(logpI[1,]))
      logpI[, idx_na] <- bound_final[idx_na]
      if (any(is.nan(exp(logpI)))){
        message('NaN detected')
      }
      
      logq = apply(-((log(2*pi)+as.numeric(log(sigma.e.hat)))*matrix(1,befdiag.i,no.gauss)+
                        (matrix(1,befdiag.i,no.gauss)*as.vector(befondiagy.i-beta.0.hat-beta.1.hat*befdiagtime.i-beta.2.hat*xval.i)-A.i.1%*%t(a.re))^2/as.numeric(sigma.e.hat))/2,2,sum)
      # f
      # f = p*I*q
      f = exp(logpI+logq)
      
      ## score of each component
      ## p
      # pdbeta.0 = apply((matrix(1,afterdiag.i,no.gauss)*(afterdiagy.i-beta.0.hat-beta.1.hat*afterdiagtime.i-beta.2.hat*xval.i-gamma.hat*afterdiagtime.i)-A.i%*%t(a.re))/as.numeric(sigma.e.hat),2,sum)*p
      # pdbeta.1 = apply((matrix(1,afterdiag.i,no.gauss)*(afterdiagy.i-beta.0.hat-beta.1.hat*afterdiagtime.i-beta.2.hat*xval.i-gamma.hat*afterdiagtime.i)*afterdiagtime.i-A.i%*%t(a.re)*afterdiagtime.i)/as.numeric(sigma.e.hat),2,sum)*p
      # pdbeta.2 = apply((matrix(1,afterdiag.i,no.gauss)*(afterdiagy.i-beta.0.hat-beta.1.hat*afterdiagtime.i-beta.2.hat*xval.i-gamma.hat*afterdiagtime.i)*xval.i-A.i%*%t(a.re)*xval.i)/as.numeric(sigma.e.hat),2,sum)*p
      # pdgamma = apply((matrix(1,afterdiag.i,no.gauss)*(afterdiagy.i-beta.0.hat-beta.1.hat*afterdiagtime.i-beta.2.hat*xval.i-gamma.hat*afterdiagtime.i)*afterdiagtime.i-A.i%*%t(a.re)*afterdiagtime.i)/as.numeric(sigma.e.hat),2,sum)*p
      # pdlambda0 = (1/lambda0.hat-sum(diff(befdiagtime.i)*exp(theta.m.hat*befdiagy.i))*exp(theta.x.hat*xval.i+theta.a.hat*t(a.re))+ondiagtime.i*exp(theta.x.hat*xval.i+theta.a.hat*t(a.re)+theta.m.hat*ondiagy.i))*p
      # pdtheta.x = (xval.i-lambda0.hat*xval.i*sum(diff(befdiagtime.i)*exp(theta.m.hat*befdiagy.i))*exp(theta.x.hat*xval.i+theta.a.hat*t(a.re))+lambda0.hat*ondiagtime.i*xval.i*exp(theta.x.hat*xval.i+theta.a.hat*t(a.re)+theta.m.hat*ondiagy.i))*p
      # pdtheta.a = (t(a.re)-lambda0.hat*sum(diff(befdiagtime.i)*befdiagy.i*exp(theta.m.hat*befdiagy.i))*exp(theta.x.hat*xval.i+theta.a.hat*t(a.re))+lambda0.hat*ondiagtime.i*t(a.re)*exp(theta.x.hat*xval.i+theta.a.hat*t(a.re)+theta.m.hat*ondiagy.i))*p
      # pdtheta.m = (ondiagy.i-lambda0.hat*xval.i*sum(diff(befdiagtime.i)*exp(theta.m.hat*befdiagy.i))*exp(theta.x.hat*xval.i+theta.a.hat*t(a.re))+lambda0.hat*ondiagtime.i*ondiagy.i*exp(theta.x.hat*xval.i+theta.a.hat*t(a.re)+theta.m.hat*ondiagy.i))*p
      # pdsigma.e = -apply(matrix(1,afterdiag.i,no.gauss)/as.numeric(sigma.e.hat)-(matrix(1,afterdiag.i,no.gauss)*(afterdiagy.i-beta.0.hat-beta.1.hat*afterdiagtime.i-beta.2.hat*xval.i-gamma.hat*afterdiagtime.i)-A.i%*%t(a.re))^2/(as.numeric(sigma.e.hat))^2,2,sum)*p/2
      # pdbeta.0I = apply((matrix(1,afterdiag.i,no.gauss)*(afterdiagy.i-beta.0.hat-beta.1.hat*afterdiagtime.i-beta.2.hat*xval.i-gamma.hat*afterdiagtime.i)-A.i%*%t(a.re))/as.numeric(sigma.e.hat),2,sum)*
      #   sqrt(pi)*exp(log(p)-b2^2/(4*c2))*(pnorm(sqrt(-2*c2)*U.i-b2/sqrt(-2*c2),0,1)-pnorm(sqrt(-2*c2)*V.i-b2/sqrt(-2*c2),0,1))/sqrt(-c2)
      # pdbeta.1I = apply((matrix(1,afterdiag.i,no.gauss)*(afterdiagy.i-beta.0.hat-beta.1.hat*afterdiagtime.i-beta.2.hat*xval.i-gamma.hat*afterdiagtime.i)*afterdiagtime.i-A.i%*%t(a.re)*afterdiagtime.i)/as.numeric(sigma.e.hat),2,sum)*
      #   sqrt(pi)*exp(log(p)-b2^2/(4*c2))*(pnorm(sqrt(-2*c2)*U.i-b2/sqrt(-2*c2),0,1)-pnorm(sqrt(-2*c2)*V.i-b2/sqrt(-2*c2),0,1))/sqrt(-c2)
      # pdbeta.2I = apply((matrix(1,afterdiag.i,no.gauss)*(afterdiagy.i-beta.0.hat-beta.1.hat*afterdiagtime.i-beta.2.hat*xval.i-gamma.hat*afterdiagtime.i)*xval.i-A.i%*%t(a.re)*xval.i)/as.numeric(sigma.e.hat),2,sum)*
      #   sqrt(pi)*exp(log(p)-b2^2/(4*c2))*(pnorm(sqrt(-2*c2)*U.i-b2/sqrt(-2*c2),0,1)-pnorm(sqrt(-2*c2)*V.i-b2/sqrt(-2*c2),0,1))/sqrt(-c2)
      # pdgammaI = apply((matrix(1,afterdiag.i,no.gauss)*(afterdiagy.i-beta.0.hat-beta.1.hat*afterdiagtime.i-beta.2.hat*xval.i-gamma.hat*afterdiagtime.i)*afterdiagtime.i-A.i%*%t(a.re)*afterdiagtime.i)/as.numeric(sigma.e.hat),2,sum)*
      #   sqrt(pi)*exp(log(p)-b2^2/(4*c2))*(pnorm(sqrt(-2*c2)*U.i-b2/sqrt(-2*c2),0,1)-pnorm(sqrt(-2*c2)*V.i-b2/sqrt(-2*c2),0,1))/sqrt(-c2)
      # pdlambda0I = (1/lambda0.hat-sum(diff(befdiagtime.i)*exp(theta.m.hat*befdiagy.i))*exp(theta.x.hat*xval.i+theta.a.hat*t(a.re))+ondiagtime.i*exp(theta.x.hat*xval.i+theta.a.hat*t(a.re)+theta.m.hat*ondiagy.i))*
      #   sqrt(pi)*exp(log(p)-b2^2/(4*c2))*(pnorm(sqrt(-2*c2)*U.i-b2/sqrt(-2*c2),0,1)-pnorm(sqrt(-2*c2)*V.i-b2/sqrt(-2*c2),0,1))/sqrt(-c2)
      # pdtheta.xI = (xval.i-lambda0.hat*xval.i*sum(diff(befdiagtime.i)*exp(theta.m.hat*befdiagy.i))*exp(theta.x.hat*xval.i+theta.a.hat*t(a.re))+lambda0.hat*ondiagtime.i*xval.i*exp(theta.x.hat*xval.i+theta.a.hat*t(a.re)+theta.m.hat*ondiagy.i))*
      #   sqrt(pi)*exp(log(p)-b2^2/(4*c2))*(pnorm(sqrt(-2*c2)*U.i-b2/sqrt(-2*c2),0,1)-pnorm(sqrt(-2*c2)*V.i-b2/sqrt(-2*c2),0,1))/sqrt(-c2)
      # pdtheta.aI = (t(a.re)-lambda0.hat*sum(diff(befdiagtime.i)*befdiagy.i*exp(theta.m.hat*befdiagy.i))*exp(theta.x.hat*xval.i+theta.a.hat*t(a.re))+lambda0.hat*ondiagtime.i*t(a.re)*exp(theta.x.hat*xval.i+theta.a.hat*t(a.re)+theta.m.hat*ondiagy.i))*
      #   sqrt(pi)*exp(log(p)-b2^2/(4*c2))*(pnorm(sqrt(-2*c2)*U.i-b2/sqrt(-2*c2),0,1)-pnorm(sqrt(-2*c2)*V.i-b2/sqrt(-2*c2),0,1))/sqrt(-c2)
      # pdtheta.mI = (ondiagy.i-lambda0.hat*xval.i*sum(diff(befdiagtime.i)*exp(theta.m.hat*befdiagy.i))*exp(theta.x.hat*xval.i+theta.a.hat*t(a.re))+lambda0.hat*ondiagtime.i*ondiagy.i*exp(theta.x.hat*xval.i+theta.a.hat*t(a.re)+theta.m.hat*ondiagy.i))*
      #   sqrt(pi)*exp(log(p)-b2^2/(4*c2))*(pnorm(sqrt(-2*c2)*U.i-b2/sqrt(-2*c2),0,1)-pnorm(sqrt(-2*c2)*V.i-b2/sqrt(-2*c2),0,1))/sqrt(-c2)
      # pdsigma.eI = -apply(matrix(1,afterdiag.i,no.gauss)/as.numeric(sigma.e.hat)-(matrix(1,afterdiag.i,no.gauss)*(afterdiagy.i-beta.0.hat-beta.1.hat*afterdiagtime.i-beta.2.hat*xval.i-gamma.hat*afterdiagtime.i)-A.i%*%t(a.re))^2/(as.numeric(sigma.e.hat))^2,2,sum)*
      #   sqrt(pi)*exp(log(p)-b2^2/(4*c2))*(pnorm(sqrt(-2*c2)*U.i-b2/sqrt(-2*c2),0,1)-pnorm(sqrt(-2*c2)*V.i-b2/sqrt(-2*c2),0,1))/sqrt(-c2)/2
      ## b2
      # b2dbeta.0 = apply(matrix(1,afterdiag.i,no.gauss)*gamma.hat/as.numeric(sigma.e.hat),2,sum)
      # b2dbeta.1 = apply(matrix(1,afterdiag.i,no.gauss)*gamma.hat*afterdiagtime.i/as.numeric(sigma.e.hat),2,sum)
      # b2dbeta.2 = apply(matrix(1,afterdiag.i,no.gauss)*gamma.hat*xval.i/as.numeric(sigma.e.hat),2,sum)
      # b2dgamma = apply((matrix(1,afterdiag.i,no.gauss)*(-(afterdiagy.i-beta.0.hat-beta.1.hat*afterdiagtime.i-beta.2.hat*xval.i-2*gamma.hat*afterdiagtime.i))+
      #                     A.i%*%t(a.re)*gamma.hat)/as.numeric(sigma.e.hat),2,sum)
      # b2dlambda0 = -exp(theta.x.hat*xval.i+theta.a.hat*t(a.re)+theta.m.hat*ondiagy.i)
      # b2dtheta.x = -lambda0.hat*xval.i*exp(theta.x.hat*xval.i+theta.a.hat*t(a.re)+theta.m.hat*ondiagy.i)
      # b2dtheta.a = -lambda0.hat*t(a.re)*exp(theta.x.hat*xval.i+theta.a.hat*t(a.re)+theta.m.hat*ondiagy.i)
      # b2dtheta.m = -lambda0.hat*ondiagy.i*exp(theta.x.hat*xval.i+theta.a.hat*t(a.re)+theta.m.hat*ondiagy.i)
      # b2dsigma.e = apply((matrix(1,afterdiag.i,no.gauss)*(afterdiagy.i-beta.0.hat-beta.1.hat*afterdiagtime.i-beta.2.hat*xval.i-gamma.hat*afterdiagtime.i)*gamma.hat
      #                     -A.i%*%t(a.re)*gamma.hat)/as.numeric(sigma.e.hat)^2,2,sum)
      ## c2
      # c2dgamma = apply(-(matrix(1,afterdiag.i,no.gauss)*gamma.hat/as.numeric(sigma.e.hat)),2,sum)
      # c2dsigma.e = apply((matrix(1,afterdiag.i,no.gauss)*gamma.hat^2/as.numeric(sigma.e.hat)^2),2,sum)
      
      ## I
      # Idbeta.0 = -2*b2/c2*b2dbeta.0*I+sqrt(pi)*exp(-b2^2/(4*c2))*(-b2dbeta.0/sqrt(-2*c2))*(dnorm(sqrt(-2*c2)*U.i-b2/sqrt(-2*c2))-dnorm(sqrt(-2*c2)*V.i-b2/sqrt(-2*c2)))/sqrt(-c2)
      # Idbeta.1 = -2*b2/c2*b2dbeta.1*I+sqrt(pi)*exp(-b2^2/(4*c2))*(-b2dbeta.1/sqrt(-2*c2))*(dnorm(sqrt(-2*c2)*U.i-b2/sqrt(-2*c2))-dnorm(sqrt(-2*c2)*V.i-b2/sqrt(-2*c2)))/sqrt(-c2)
      # Idbeta.2 = -2*b2/c2*b2dbeta.2*I+sqrt(pi)*exp(-b2^2/(4*c2))*(-b2dbeta.2/sqrt(-2*c2))*(dnorm(sqrt(-2*c2)*U.i-b2/sqrt(-2*c2))-dnorm(sqrt(-2*c2)*V.i-b2/sqrt(-2*c2)))/sqrt(-c2)
      #  Idgamma = ((sqrt(pi)/sqrt(-c2))*(-(2*b2*c2*b2dgamma-b2^2*c2dgamma)/(16*c2^2))*exp(-b2^2/(4*c2))+sqrt(pi)*(-1/(2*sqrt(-c2)^3))*(-c2dgamma)*exp(-b2^2/(4*c2)))*
      #.        (pnorm(sqrt(-2*c2)*U.i-b2/sqrt(-2*c2))-pnorm(sqrt(-2*c2)*V.i-b2/sqrt(-2*c2)))+
      #          sqrt(pi)*exp(-b2^2/(4*c2))*((-2*c2dgamma*U.i/(2*sqrt(-2*c2))-(b2dgamma/sqrt(-2*c2)-(-2*c2dgamma/(2*sqrt(-2*c2)^3))))*dnorm(sqrt(-2*c2)*U.i-b2/sqrt(-2*c2))-
      #                  (-2*c2dgamma*V.i/(2*sqrt(-2*c2))-(b2dgamma/sqrt(-2*c2)-(-2*c2dgamma/(2*sqrt(-2*c2)^3))))*dnorm(sqrt(-2*c2)*V.i-b2/sqrt(-2*c2)))/sqrt(-c2)
      # Idlambda0 = -2*b2/c2*b2dlambda0*I+sqrt(pi)*exp(-b2^2/(4*c2))*(-b2dlambda0/sqrt(-2*c2))*(dnorm(sqrt(-2*c2)*U.i-b2/sqrt(-2*c2))-dnorm(sqrt(-2*c2)*V.i-b2/sqrt(-2*c2)))/sqrt(-c2)
      # Idtheta.x = -2*b2/c2*b2dtheta.x*I+sqrt(pi)*exp(-b2^2/(4*c2))*(-b2dtheta.x/sqrt(-2*c2))*(dnorm(sqrt(-2*c2)*U.i-b2/sqrt(-2*c2))-dnorm(sqrt(-2*c2)*V.i-b2/sqrt(-2*c2)))/sqrt(-c2)
      # Idtheta.a = -2*b2/c2*b2dtheta.a*I+sqrt(pi)*exp(-b2^2/(4*c2))*(-b2dtheta.a/sqrt(-2*c2))*(dnorm(sqrt(-2*c2)*U.i-b2/sqrt(-2*c2))-dnorm(sqrt(-2*c2)*V.i-b2/sqrt(-2*c2)))/sqrt(-c2)
      # Idtheta.m = -2*b2/c2*b2dtheta.m*I+sqrt(pi)*exp(-b2^2/(4*c2))*(-b2dtheta.m/sqrt(-2*c2))*(dnorm(sqrt(-2*c2)*U.i-b2/sqrt(-2*c2))-dnorm(sqrt(-2*c2)*V.i-b2/sqrt(-2*c2)))/sqrt(-c2)
      #Idsigma.e = ((sqrt(pi)/sqrt(-c2))*(-(2*b2*c2*b2dsigma.e-b2^2*c2dsigma.e)/(16*c2^2))*exp(-b2^2/(4*c2))+sqrt(pi)*(-1/(2*sqrt(-c2)^3))*(-c2dsigma.e)*exp(-b2^2/(4*c2)))*
      # (pnorm(sqrt(-2*c2)*U.i-b2/sqrt(-2*c2))-pnorm(sqrt(-2*c2)*V.i-b2/sqrt(-2*c2)))+
      #  sqrt(pi)*exp(-b2^2/(4*c2))*((-2*c2dsigma.e*U.i/(2*sqrt(-2*c2))-(b2dsigma.e/sqrt(-2*c2)-(-2*c2dsigma.e/(2*sqrt(-2*c2)^3))))*dnorm(sqrt(-2*c2)*U.i-b2/sqrt(-2*c2))-
      #                               (-2*c2dsigma.e*V.i/(2*sqrt(-2*c2))-(b2dsigma.e/sqrt(-2*c2)-(-2*c2dsigma.e/(2*sqrt(-2*c2)^3))))*dnorm(sqrt(-2*c2)*V.i-b2/sqrt(-2*c2)))/sqrt(-c2)
      
      # pIdbeta.0 = -2*b2/c2*b2dbeta.0*pI+sqrt(pi)*exp(log(p)-b2^2/(4*c2))*(-b2dbeta.0/sqrt(-2*c2))*(dnorm(sqrt(-2*c2)*U.i-b2/sqrt(-2*c2))-dnorm(sqrt(-2*c2)*V.i-b2/sqrt(-2*c2)))/sqrt(-c2)
      # pIdbeta.1 = -2*b2/c2*b2dbeta.1*pI+sqrt(pi)*exp(log(p)-b2^2/(4*c2))*(-b2dbeta.1/sqrt(-2*c2))*(dnorm(sqrt(-2*c2)*U.i-b2/sqrt(-2*c2))-dnorm(sqrt(-2*c2)*V.i-b2/sqrt(-2*c2)))/sqrt(-c2)
      # pIdbeta.2 = -2*b2/c2*b2dbeta.2*pI+sqrt(pi)*exp(log(p)-b2^2/(4*c2))*(-b2dbeta.2/sqrt(-2*c2))*(dnorm(sqrt(-2*c2)*U.i-b2/sqrt(-2*c2))-dnorm(sqrt(-2*c2)*V.i-b2/sqrt(-2*c2)))/sqrt(-c2)
      # 
      # pIdgamma = ((sqrt(pi)/sqrt(-c2))*(-(2*b2*c2*b2dgamma-b2^2*c2dgamma)/(4*c2^2))*exp(log(p)-b2^2/(4*c2))+sqrt(pi)*(-1/(2*sqrt(-c2)^3))*(-c2dgamma)*exp(log(p)-b2^2/(4*c2)))*
      #   (pnorm(sqrt(-2*c2)*U.i-b2/sqrt(-2*c2))-pnorm(sqrt(-2*c2)*V.i-b2/sqrt(-2*c2)))+
      #   sqrt(pi)*exp(log(p)-b2^2/(4*c2))*((-2*c2dgamma*U.i/(2*sqrt(-2*c2))-(b2dgamma/sqrt(-2*c2)-(-2*c2dgamma/(2*sqrt(-2*c2)^3))))*dnorm(sqrt(-2*c2)*U.i-b2/sqrt(-2*c2))-
      #                                       (-2*c2dgamma*V.i/(2*sqrt(-2*c2))-(b2dgamma/sqrt(-2*c2)-(-2*c2dgamma/(2*sqrt(-2*c2)^3))))*dnorm(sqrt(-2*c2)*V.i-b2/sqrt(-2*c2)))/sqrt(-c2)
      # pIdlambda0 = -2*b2/c2*b2dlambda0*pI+sqrt(pi)*exp(log(p)-b2^2/(4*c2))*(-b2dlambda0/sqrt(-2*c2))*(dnorm(sqrt(-2*c2)*U.i-b2/sqrt(-2*c2))-dnorm(sqrt(-2*c2)*V.i-b2/sqrt(-2*c2)))/sqrt(-c2)
      # pIdtheta.x = -2*b2/c2*b2dtheta.x*pI+sqrt(pi)*exp(log(p)-b2^2/(4*c2))*(-b2dtheta.x/sqrt(-2*c2))*(dnorm(sqrt(-2*c2)*U.i-b2/sqrt(-2*c2))-dnorm(sqrt(-2*c2)*V.i-b2/sqrt(-2*c2)))/sqrt(-c2)
      # pIdtheta.a = -2*b2/c2*b2dtheta.a*pI+sqrt(pi)*exp(log(p)-b2^2/(4*c2))*(-b2dtheta.a/sqrt(-2*c2))*(dnorm(sqrt(-2*c2)*U.i-b2/sqrt(-2*c2))-dnorm(sqrt(-2*c2)*V.i-b2/sqrt(-2*c2)))/sqrt(-c2)
      # 
      # pIdtheta.m = -2*b2/c2*b2dtheta.m*pI+sqrt(pi)*exp(log(p)-b2^2/(4*c2))*(-b2dtheta.m/sqrt(-2*c2))*(dnorm(sqrt(-2*c2)*U.i-b2/sqrt(-2*c2))-dnorm(sqrt(-2*c2)*V.i-b2/sqrt(-2*c2)))/sqrt(-c2)
      # pIdsigma.e = ((sqrt(pi)/sqrt(-c2))*(-(2*b2*c2*b2dsigma.e-b2^2*c2dsigma.e)/(4*c2^2))*exp(log(p)-b2^2/(4*c2))+sqrt(pi)*(-1/(2*sqrt(-c2)^3))*(-c2dsigma.e)*exp(log(p)-b2^2/(4*c2)))*
      #   (pnorm(sqrt(-2*c2)*U.i-b2/sqrt(-2*c2))-pnorm(sqrt(-2*c2)*V.i-b2/sqrt(-2*c2)))+
      #   sqrt(pi)*exp(log(p)-b2^2/(4*c2))*((-2*c2dsigma.e*U.i/(2*sqrt(-2*c2))-(b2dsigma.e/sqrt(-2*c2)-(-2*c2dsigma.e/(2*sqrt(-2*c2)^3))))*dnorm(sqrt(-2*c2)*U.i-b2/sqrt(-2*c2))-
      #                                       (-2*c2dsigma.e*V.i/(2*sqrt(-2*c2))-(b2dsigma.e/sqrt(-2*c2)-(-2*c2dsigma.e/(2*sqrt(-2*c2)^3))))*dnorm(sqrt(-2*c2)*V.i-b2/sqrt(-2*c2)))/sqrt(-c2)
      ## q
      # qdbeta.0 = apply(matrix(1,befdiag.i,no.gauss)*(befondiagy.i-beta.0.hat-beta.1.hat*befdiagtime.i-beta.2.hat*xval.i)/as.numeric(sigma.e.hat),2,sum)*q
      # qdbeta.1 = apply(matrix(1,befdiag.i,no.gauss)*(befondiagy.i-beta.0.hat-beta.1.hat*befdiagtime.i-beta.2.hat*xval.i)*befdiagtime.i/as.numeric(sigma.e.hat),2,sum)*q
      # qdbeta.2 = apply(matrix(1,befdiag.i,no.gauss)*(befondiagy.i-beta.0.hat-beta.1.hat*befdiagtime.i-beta.2.hat*xval.i)*xval.i/as.numeric(sigma.e.hat),2,sum)*q
      # qdsigma.e = -apply(matrix(1,befdiag.i,no.gauss)/as.numeric(sigma.e.hat)-(matrix(1,befdiag.i,no.gauss)*(befondiagy.i-beta.0.hat-beta.1.hat*befdiagtime.i-beta.2.hat*xval.i))^2/(as.numeric(sigma.e.hat))^2,2,sum)*q/2
      
      
      ## score function
      # fdbeta.0 = pdbeta.0I*q*exp(d2*t(a.re)+e2*t(a.re)^2)+
      #   pIdbeta.0*q*exp(d2*t(a.re)+e2*t(a.re)^2)+
      #   pI*qdbeta.0*exp(d2*t(a.re)+e2*t(a.re)^2)+
      #   pI*q*(d2dbeta.0*t(a.re))*exp(d2*t(a.re)+e2*t(a.re)^2)
      # fdbeta.1 = pdbeta.1I*q*exp(d2*t(a.re)+e2*t(a.re)^2)+
      #   pIdbeta.1*q*exp(d2*t(a.re)+e2*t(a.re)^2)+
      #   pI*qdbeta.1*exp(d2*t(a.re)+e2*t(a.re)^2)+
      #   pI*q*(d2dbeta.1*t(a.re))*exp(d2*t(a.re)+e2*t(a.re)^2)
      # fdbeta.2 = pdbeta.2I*q*exp(d2*t(a.re)+e2*t(a.re)^2)+
      #   pIdbeta.2*q*exp(d2*t(a.re)+e2*t(a.re)^2)+
      #   pI*qdbeta.2*exp(d2*t(a.re)+e2*t(a.re)^2)+
      #   pI*q*(d2dbeta.2*t(a.re))*exp(d2*t(a.re)+e2*t(a.re)^2)
      # fdgamma = pdgammaI*q*exp(d2*t(a.re)+e2*t(a.re)^2)+
      #   pIdgamma*q*exp(d2*t(a.re)+e2*t(a.re)^2)
      # fdlambda0 = pdlambda0I*q*exp(d2*t(a.re)+e2*t(a.re)^2)+
      #   pIdlambda0*q*exp(d2*t(a.re)+e2*t(a.re)^2)
      # fdtheta.x = pdtheta.xI*q*exp(d2*t(a.re)+e2*t(a.re)^2)+
      #   pIdtheta.x*q*exp(d2*t(a.re)+e2*t(a.re)^2)
      # fdtheta.a = pdtheta.aI*q*exp(d2*t(a.re)+e2*t(a.re)^2)+
      #   pIdtheta.a*q*exp(d2*t(a.re)+e2*t(a.re)^2)
      # fdtheta.m = pdtheta.mI*q*exp(d2*t(a.re)+e2*t(a.re)^2)+
      #   pIdtheta.m*q*exp(d2*t(a.re)+e2*t(a.re)^2)
      # fdsigma.e = pdsigma.eI*q*exp(d2*t(a.re)+e2*t(a.re)^2)+
      #   pIdsigma.e*q*exp(d2*t(a.re)+e2*t(a.re)^2)+
      #   pI*qdsigma.e*exp(d2*t(a.re)+e2*t(a.re)^2)+
      #   pI*q*(d2dsigma.e*t(a.re)+e2dsigma.e*t(a.re)^2)*exp(d2*t(a.re)+e2*t(a.re)^2)
      # fdsigma.a = f*(-(1/as.numeric(sigma.a.hat)-t(a.re)^2/as.numeric(sigma.a.hat)^2)/2)
      ## GH integration
      Q[i,] = as.numeric(f%*%weight)
      # Qdbeta.0[i,] = as.numeric(fdbeta.0%*%weight)/Q[i,]
      # Qdbeta.1[i,] = as.numeric(fdbeta.1%*%weight)/Q[i,]
      # Qdbeta.2[i,] = as.numeric(fdbeta.2%*%weight)/Q[i,]
      # Qdgamma[i,] = as.numeric(fdgamma%*%weight)/Q[i,]
      # Qdlambda0[i,] = as.numeric(fdlambda0%*%weight)/Q[i,]
      # Qdtheta.x[i,] = as.numeric(fdtheta.x%*%weight)/Q[i,]
      # Qdtheta.a[i,] = as.numeric(fdtheta.a%*%weight)/Q[i,]
      # Qdtheta.m[i,] = as.numeric(fdtheta.m%*%weight)/Q[i,]
      # Qdsigma.e[i,] = as.numeric(fdsigma.e%*%weight)/Q[i,]
      # Qdsigma.a[i,] = as.numeric(fdsigma.a%*%weight)/Q[i,]
    }
    else if (status.i == 0){
      ## components
      q = exp(apply(-((log(2*pi)+as.numeric(log(sigma.e.hat)))*matrix(1,befdiag.i,no.gauss)+
                        (matrix(1,befdiag.i,no.gauss)*as.vector(befondiagy.i-beta.0.hat-beta.1.hat*befdiagtime.i-beta.2.hat*xval.i)-A.i.1%*%t(a.re))^2/as.numeric(sigma.e.hat))/2,2,sum)-
                lambda0.hat*sum(diff(befdiagtime.i)*exp(theta.m.hat*befdiagy.i))*exp(theta.x.hat*xval.i+theta.a.hat*t(a.re)))
      
      ## score
      # q
      # qdbeta.0 = apply(matrix(1,befdiag.i,no.gauss)*(befondiagy.i-beta.0.hat-beta.1.hat*befdiagtime.i-beta.2.hat*xval.i)/as.numeric(sigma.e.hat),2,sum)*q
      # qdbeta.1 = apply(matrix(1,befdiag.i,no.gauss)*(befondiagy.i-beta.0.hat-beta.1.hat*befdiagtime.i-beta.2.hat*xval.i)*befdiagtime.i/as.numeric(sigma.e.hat),2,sum)*q
      # qdbeta.2 = apply(matrix(1,befdiag.i,no.gauss)*(befondiagy.i-beta.0.hat-beta.1.hat*befdiagtime.i-beta.2.hat*xval.i)*xval.i/as.numeric(sigma.e.hat),2,sum)*q
      # qdsigma.e = -apply(matrix(1,befdiag.i,no.gauss)/as.numeric(sigma.e.hat)-(matrix(1,befdiag.i,no.gauss)*(befondiagy.i-beta.0.hat-beta.1.hat*befdiagtime.i-beta.2.hat*xval.i))^2/(as.numeric(sigma.e.hat))^2,2,sum)*q/2
      
      # f
      f = q
      ## score function
      # fdbeta.0 = qdbeta.0*exp(d2*t(a.re)+e2*t(a.re)^2+f2)+q*(d2dbeta.0*t(a.re))*exp(d2*t(a.re)+e2*t(a.re)^2+f2)
      # fdbeta.1 = qdbeta.1*exp(d2*t(a.re)+e2*t(a.re)^2+f2)+q*(d2dbeta.1*t(a.re))*exp(d2*t(a.re)+e2*t(a.re)^2+f2)
      # fdbeta.2 = qdbeta.2*exp(d2*t(a.re)+e2*t(a.re)^2+f2)+q*(d2dbeta.2*t(a.re))*exp(d2*t(a.re)+e2*t(a.re)^2+f2)
      # fdgamma = matrix(0,no.gauss,nrow=1)
      # fdlambda0 = q*(f2dlambda0)*exp(d2*t(a.re)+e2*t(a.re)^2+f2)
      # fdtheta.x = q*(f2dtheta.x)*exp(d2*t(a.re)+e2*t(a.re)^2+f2)
      # fdtheta.a = q*(f2dtheta.a)*exp(d2*t(a.re)+e2*t(a.re)^2+f2)
      # fdtheta.m = q*(f2dtheta.m)*exp(d2*t(a.re)+e2*t(a.re)^2+f2)
      # fdsigma.e = qdsigma.e*exp(d2*t(a.re)+e2*t(a.re)^2+f2)+q*(d2dsigma.e*t(a.re)+e2dsigma.e*t(a.re^2))*exp(d2*t(a.re)+e2*t(a.re)^2+f2)
      # fdsigma.a = f*(-(1/as.numeric(sigma.a.hat)-t(a.re)^2/as.numeric(sigma.a.hat)^2)/2)
      ## GH
      Q[i,] = as.numeric(f%*%weight)
      # Qdbeta.0[i,] = as.numeric(fdbeta.0%*%weight)/Q[i,]
      # Qdbeta.1[i,] = as.numeric(fdbeta.1%*%weight)/Q[i,]
      # Qdbeta.2[i,] = as.numeric(fdbeta.2%*%weight)/Q[i,]
      # Qdgamma[i,] = as.numeric(fdgamma%*%weight)/Q[i,]
      # Qdlambda0[i,] = as.numeric(fdlambda0%*%weight)/Q[i,]
      # Qdtheta.x[i,] = as.numeric(fdtheta.x%*%weight)/Q[i,]
      # Qdtheta.a[i,] = as.numeric(fdtheta.a%*%weight)/Q[i,]
      # Qdtheta.m[i,] = as.numeric(fdtheta.m%*%weight)/Q[i,]
      # Qdsigma.e[i,] = as.numeric(fdsigma.e%*%weight)/Q[i,]
      # Qdsigma.a[i,] = as.numeric(fdsigma.a%*%weight)/Q[i,]
    }
    
  }
  return(list(ll=sum(log(Q)), ll.ind=log(Q)))
}
########################


likelihood.vec = function(long.data, parameters)
{
  long.data$id = as.factor(long.data$id)
  id.matrix <- model.matrix(~ id - 1, data = long.data)
  sm = Matrix(t(id.matrix),sparse = T)
  n.id = length(unique(long.data$id))
  long.data.id = long.data[!duplicated(long.data$id),]
  beta.0.hat =  t(as.matrix(parameters %>% dplyr::select(contains("beta.0"))))
  K = length(beta.0.hat)
  
  beta.1.hat = t(as.matrix(parameters %>% dplyr::select(contains("beta.1"))))
  beta.2.vec = data.frame(parameters %>% dplyr::select(contains("beta.2")))
  col_names <- colnames(beta.2.vec)
  if (length(col_names)==1){
    beta.2.hat=as.matrix(beta.2.vec)
    p = ncol(beta.2.hat)
  } else{
  indices <- do.call(rbind, strsplit(col_names, "_"))[, 2:3]
  indices <- apply(indices, 2, as.integer)
  beta.2.hat = matrix(0, nrow=max(indices[,1]),ncol=max(indices[,2]))
  p=ncol(beta.2.hat)
  for (i in seq_along(col_names)){
    row_idx <- indices[i, 1]
    col_idx <- indices[i, 2]
    beta.2.hat[row_idx, col_idx] <- beta.2.vec[[col_names[i]]]
  }
  }
  gamma.hat = t(as.matrix(parameters %>% dplyr::select(contains("gamma"))))
  
  lambda0.hat = as.numeric(parameters %>% dplyr::select(contains("lambda0")))
  theta.x.hat = as.matrix(parameters %>% dplyr::select(contains("theta.x")))
  theta.a.hat = as.matrix(parameters %>% dplyr::select(contains("theta.a")))
  theta.m.hat = as.matrix(parameters %>% dplyr::select(contains("theta.m")))
  
  Sigma.a.hat = matrix(0, nrow=K,ncol=K)
  upper_vec_a = as.matrix(parameters %>% dplyr::select(contains("Sigma.a")))
  Sigma.a.hat[upper.tri(Sigma.a,diag = TRUE)] = upper_vec_a
  if (K>1){
    Sigma.a.hat=Sigma.a.hat+t(Sigma.a.hat)-diag(diag(Sigma.a.hat))
  }
  Sigma.e.hat = matrix(0, nrow=K,ncol=K)
  upper_vec_e = as.matrix(parameters %>% dplyr::select(contains("Sigma.e")))
  Sigma.e.hat[upper.tri(Sigma.e,diag = TRUE)] = upper_vec_e
  if (K>1){
    Sigma.e.hat=Sigma.e.hat+t(Sigma.e.hat)-diag(diag(Sigma.e.hat))
  }
  inv.Sigma.e.hat = solve(Sigma.e.hat)
  
  a.re = nodes%*%sqrtm(2*Sigma.a.hat)
  
  y = long.data %>% dplyr::select(contains("ynew"))
  xval = long.data %>% dplyr::select(contains("xval"))
  xval.id = long.data.id %>% dplyr::select(contains("xval"))
  time = long.data$visittime
  status = long.data.id$status
  V = long.data.id$preetime
  U = long.data.id$etime
  A.i.1 = as.numeric(long.data$visit>long.data$befdiag)
  A.i.2 = as.numeric(long.data$visit<long.data$befdiag)
  A.i.3 = as.numeric(long.data$visit<=long.data$befdiag)
  A.i.4 = as.numeric(long.data$visit==long.data$befdiag)
  theta.xx = theta.x.hat%*%t(xval.id) #400 for individual
  theta.xx.mat = data.frame(matrix(rep(t(theta.xx),no.gauss^K),ncol = no.gauss^K,byrow=FALSE))
  theta.aa = theta.a.hat%*%t(a.re) #400 for nodes
  theta.aa.mat = data.frame(matrix(rep(theta.aa,n.id),nrow = n.id,byrow = TRUE))
  theta.mm = A.i.4[A.i.4!=0]*theta.m.hat%*%t(y[A.i.4!=0,]) #400 for individual
  theta.mm.mat = data.frame(matrix(rep(t(theta.mm),no.gauss^K),ncol=no.gauss^K,byrow=FALSE))
  factor = as.matrix(exp(theta.xx.mat+theta.aa.mat+theta.mm.mat))
  if (K>1){
    part1 = as.matrix((y-rep(1,nrow(y))%*%t(beta.0.hat)-time%*%t(beta.1.hat)-t(beta.2.hat%*%t(xval))-time%*%t(gamma.hat))*A.i.1)
  } else {
    part1 = as.matrix((y-as.numeric(beta.0.hat)-time%*%t(beta.1.hat)-t(beta.2.hat%*%t(xval))-time%*%t(gamma.hat))*A.i.1)
  }
  part1_1 = rowSums(part1%*%inv.Sigma.e.hat*part1)
  # part1_2 = data.frame(id=long.data$id,part1_1)
  part2 = t(rowSums(a.re%*%inv.Sigma.e.hat*a.re)%*%t(A.i.1))
  # part2_1 = data.frame(id=long.data$id, part2)
  part3 = part1%*%inv.Sigma.e.hat%*%t(a.re)
  # part3_1 = data.frame(id=long.data$id,part3)
  # Part1 = part1_2 %>%
  #   group_by(id) %>%
  #   summarise(sum_part1=sum(part1_1,na.rm = TRUE))%>%
  #   select(-id)
  # Part2 = part2_1 %>%
  #   group_by(id) %>%
  #   summarise(across(starts_with("X"),sum,.names = "sum_{.col}"))%>%
  #   select(-id)
  # Part3 = part3_1 %>%
  #   group_by(id) %>%
  #   summarise(across(starts_with("X"),sum,.names = "sum_{.col}"))%>%
  #   select(-id)
  part4 = (K*log(2*pi)+log(det(Sigma.e.hat)))*A.i.1
  # part4_1 = data.frame(id=long.data$id,part4+part1_1+part2-2*part3)
  Part4 = sm%*%(part4+part1_1+part2-2*part3)
  # Part4 = part4_1 %>%
  #   group_by(id) %>%
  #   summarise(across(starts_with("X"),sum,.names = "sum_{.col}")) %>%
  #   dplyr::select(-id)
  a2 = -Part4/2+lambda0.hat*V*factor
  part5 = A.i.2*c(diff(time),0)*exp(theta.m.hat%*%t(y))
  Part5 = sm%*%t(part5)
  # rownames(part5)=paste("part5")
  # part5_1 = data.frame(id=long.data$id, t(part5))
  # Part5 = part5_1 %>%
  #   group_by(id) %>%
  #   summarise(sum_part5 = sum(part5,na.rm=TRUE)) %>%
  #   dplyr::select(-id)
  factor2 = as.matrix(exp(theta.xx.mat+theta.aa.mat))*Part5
  logp = log(lambda0.hat)+log(factor)-lambda0.hat*factor2+a2

  # logp[is.na(logp)] <- 0
  
  part6 = part1%*%inv.Sigma.e.hat%*%gamma.hat
  # part6_1 = data.frame(id=long.data$id, part6)
  # Part6 = part6_1 %>%
  #   group_by(id) %>%
  #   summarise(sum_part6 = sum(part6,na.rm=TRUE)) %>%
  #   select(-id)
  part7 = data.frame(A.i.1%*%t(a.re%*%inv.Sigma.e.hat%*%gamma.hat))
  Part7 = sm%*%as.matrix(part6-part7)
  # part7_1 = data.frame(id=long.data$id, part6-part7)
  # Part7 = part7_1 %>%
  #   group_by(id) %>%
  #   summarise(across(starts_with("X"),sum,.names = "sum_{.col}")) %>%
  #   dplyr::select(-id)
  b2 = -Part7-lambda0.hat*factor
  part8 = A.i.1%*%t(gamma.hat)%*%inv.Sigma.e.hat%*%gamma.hat
  Part8 = sm%*%part8
  # part8_1 = data.frame(id=long.data$id,part8)
  # Part8 = part8_1 %>% 
  #   group_by(id) %>%
  #   summarise(sum_part8 = sum(part8, na.rm = TRUE)) %>%
  #   dplyr::select(-id)
  # c2 = -data.frame(matrix(rep(Part8,no.gauss^K),ncol = no.gauss^K,byrow = FALSE))/2
  c2 = -Part8%*%rep(1,no.gauss^K)/2
  logpI = as.matrix(log(sqrt(pi))+log(pnorm(as.matrix(sqrt(-2*c2)*U-b2/sqrt(-2*c2)),0,1)-pnorm(as.matrix(sqrt(-2*c2)*V-b2/sqrt(-2*c2),0,1)))+(logp-b2^2/(4*c2))-log(sqrt(-c2)))
  logpI_apx = as.matrix(log(sqrt(pi))+(logp-b2^2/(4*c2))-log(sqrt(2*pi))+(-(sqrt(-2*c2)*(U+V)-2*b2/sqrt(-2*c2))^2/8)+log(sqrt(-2*c2)*(U-V))-log(sqrt(-c2)))
  logpI[is.infinite(logpI)] = logpI_apx[is.infinite(logpI)]
  logpI[is.na(logpI) | is.infinite(logpI)] = -700
  # logpI = data.frame(logpI)
  if (K>1){
    part9 = as.matrix((y-rep(1,nrow(y))%*%t(beta.0.hat)-time%*%t(beta.1.hat)-t(beta.2.hat%*%t(xval)))*A.i.3)
  } else {
    part9 = as.matrix((y-as.numeric(beta.0.hat)-time%*%t(beta.1.hat)-t(beta.2.hat%*%t(xval)))*A.i.3)
  }
  part9_1 = rowSums(part9%*%inv.Sigma.e.hat*part9)
  # part9_2 = data.frame(id=long.data$id,part9_1)
  # Part9 = part9_2 %>%
  #   group_by(id) %>%
  #   summarise(sum_part9=sum(part9_1, na.rm = TRUE)) %>%
  #   select(-id)
  part10 = t(rowSums(a.re%*%inv.Sigma.e.hat*a.re)%*%t(A.i.3))
  # part10_1 = data.frame(id=long.data$id, part10)
  # Part10 = part10_1 %>%
  #   group_by(id) %>%
  #   summarise(across(starts_with("X"),sum,.names = "sum_{.col}"))%>%
  #   select(-id)
  part11 = part9%*%inv.Sigma.e.hat%*%t(a.re)
  # part11_1 = data.frame(id=long.data$id,part11)
  # Part11 = part3_1 %>%
  #   group_by(id) %>%
  #   summarise(across(starts_with("X"),sum,.names = "sum_{.col}"))%>%
  #   select(-id)
  part12 = (K*log(2*pi)+log(det(Sigma.e.hat)))*A.i.3
  # part12_1 = data.frame(id=long.data$id,part12)
  # Part12 = part12_1 %>%
  #   group_by(id) %>%
  #   summarise(sum_part4=sum(part12,na.rm = TRUE)) %>%
  #   select(-id)
  Part13 = sm%*%(part12+part9_1+part10-2*part11)
  # part13_1 = data.frame(id=long.data$id, part12+part9_1+part10-2*part11)
  # Part13 = part13_1%>%
  #   group_by(id) %>%
  #   summarise(across(starts_with("X"),sum,.names = "sum_{.col}"))%>%
  #   dplyr::select(-id)
  logq1 = -Part13/2
  logq2 = logq1 - lambda0.hat*factor2
  f = exp(status*(logpI+logq1)+(1-status)*logq2)
  
  Q = as.matrix(t(f%*%weight))/(sqrt(pi)^K)
  
  return(list(ll=sum(log(Q)), ll.ind=log(Q)))
}
  
likelihood.vec2 = function(long.data, parameters)
{
  long.data$id = as.factor(long.data$id)
  id.matrix <- model.matrix(~ id - 1, data = long.data)
  sm = Matrix(t(id.matrix),sparse = T)
  n.id = length(unique(long.data$id))
  long.data.id = long.data[!duplicated(long.data$id),]
  parameters = data.frame(t(parameters))
  # colnames(parameters)[1:K] = paste("beta.0","_",1:K,sep="")
  # colnames(parameters)[(K+1):(2:K)] = paste("beta.1","_",1:K,sep="")
  # index_grid <- expand.grid(row = 1:K, col = 1:p)
  # colnames(parameters)[(2*K+1):(2*K+p*K)]=paste("beta.2",index_grid$row, index_grid$col,sep='_')
  # colnames(parameters)[(2*K+p*K+1):(3*K+p*K)]=paste("gamma",1:K,sep="_")
  # colnames(parameters)[(3*K+p*K+1)] = paste("lambda0")
  # colnames(parameters)[(3*K+p*K+2):(3*K+p*K+p+1)] = paste("theta.x",1:p,sep="_")
  # colnames(parameters)[(3*K+p*K+p+2):(4*K+p*K+p+1)] = paste("theta.a",1:K,sep="_")
  # colnames(parameters)[(4*K+p*K+p+2):(5*K+p*K+p+1)] = paste("theta.m",1:K,sep="_")
  # pos_labels_a <- sapply(1:K, function(i) {
  #   sapply(i:K, function(j) paste0("Sigma.a_", i, "_", j))
  # })
  # pos_labels_a <- unlist(pos_labels_a)
  # colnames(parameters)[(5*K+p*K+p+2):(5*K+p*K+p+2+K*(K+1)/2)] = pos_labels_a
  # pos_labels_e <- sapply(1:K, function(i) {
  #   sapply(i:K, function(j) paste0("Sigma.e_", i, "_", j))
  # })
  # pos_labels_e <- unlist(pos_labels_e)
  # colnames(parameters[(5*K+p*K+p+3+K*(K+1)/2):(5*K+p*K+p+2+K*(K+1))]) <- pos_labels_e
  
  beta.0.hat =  t(as.matrix(parameters %>% dplyr::select(contains("beta.0"))))
  # K = length(beta.0.hat)
  beta.1.hat = t(as.matrix(parameters %>% dplyr::select(contains("beta.1"))))
  beta.2.vec = data.frame(parameters %>% dplyr::select(contains("beta.2")))
  col_names <- colnames(beta.2.vec)
  if (length(col_names)==1){
    beta.2.hat=as.matrix(beta.2.vec)
    p = ncol(beta.2.hat)
  } else{
    indices <- do.call(rbind, strsplit(col_names, "_"))[, 2:3]
    indices <- apply(indices, 2, as.integer)
    beta.2.hat = matrix(0, nrow=max(indices[,1]),ncol=max(indices[,2]))
    p=ncol(beta.2.hat)
    for (i in seq_along(col_names)){
      row_idx <- indices[i, 1]
      col_idx <- indices[i, 2]
      beta.2.hat[row_idx, col_idx] <- beta.2.vec[[col_names[i]]]
    }
  }
  gamma.hat = t(as.matrix(parameters %>% dplyr::select(contains("gamma"))))
  
  lambda0.hat = as.numeric(parameters %>% dplyr::select(contains("lambda0")))
  theta.x.hat = as.matrix(parameters %>% dplyr::select(contains("theta.x")))
  theta.a.hat = as.matrix(parameters %>% dplyr::select(contains("theta.a")))
  theta.m.hat = as.matrix(parameters %>% dplyr::select(contains("theta.m")))
  
  Sigma.a.hat = matrix(0, nrow=K,ncol=K)
  upper_vec_a = as.matrix(parameters %>% dplyr::select(contains("Sigma.a")))
  Sigma.a.hat[upper.tri(Sigma.a,diag = TRUE)] = upper_vec_a
  if (K>1){
    Sigma.a.hat=Sigma.a.hat+t(Sigma.a.hat)-diag(diag(Sigma.a.hat))
  }
  Sigma.e.hat = matrix(0, nrow=K,ncol=K)
  upper_vec_e = as.matrix(parameters %>% dplyr::select(contains("Sigma.e")))
  Sigma.e.hat[upper.tri(Sigma.e,diag = TRUE)] = upper_vec_e
  if (K>1){
    Sigma.e.hat=Sigma.e.hat+t(Sigma.e.hat)-diag(diag(Sigma.e.hat))
  }
  inv.Sigma.e.hat = solve(Sigma.e.hat)
  
  a.re = nodes%*%sqrtm(2*Sigma.a.hat)
  
  y = long.data %>% dplyr::select(contains("ynew"))
  xval = long.data %>% dplyr::select(contains("xval"))
  xval.id = long.data.id %>% dplyr::select(contains("xval"))
  time = long.data$visittime
  status = long.data.id$status
  V = long.data.id$preetime
  U = long.data.id$etime
  A.i.1 = as.numeric(long.data$visit>long.data$befdiag)
  A.i.2 = as.numeric(long.data$visit<long.data$befdiag)
  A.i.3 = as.numeric(long.data$visit<=long.data$befdiag)
  A.i.4 = as.numeric(long.data$visit==long.data$befdiag)
  theta.xx = theta.x.hat%*%t(xval.id) #400 for individual
  theta.xx.mat = data.frame(matrix(rep(t(theta.xx),no.gauss^K),ncol = no.gauss^K,byrow=FALSE))
  theta.aa = theta.a.hat%*%t(a.re) #400 for nodes
  theta.aa.mat = data.frame(matrix(rep(theta.aa,n.id),nrow = n.id,byrow = TRUE))
  theta.mm = A.i.4[A.i.4!=0]*theta.m.hat%*%t(y[A.i.4!=0,]) #400 for individual
  theta.mm.mat = data.frame(matrix(rep(t(theta.mm),no.gauss^K),ncol=no.gauss^K,byrow=FALSE))
  factor = as.matrix(exp(theta.xx.mat+theta.aa.mat+theta.mm.mat))
  if (K>1){
    part1 = as.matrix((y-rep(1,nrow(y))%*%t(beta.0.hat)-time%*%t(beta.1.hat)-t(beta.2.hat%*%t(xval))-time%*%t(gamma.hat))*A.i.1)
  } else {
    part1 = as.matrix((y-as.numeric(beta.0.hat)-time%*%t(beta.1.hat)-t(beta.2.hat%*%t(xval))-time%*%t(gamma.hat))*A.i.1)
  }
  part1_1 = rowSums(part1%*%inv.Sigma.e.hat*part1)
  # part1_2 = data.frame(id=long.data$id,part1_1)
  part2 = t(rowSums(a.re%*%inv.Sigma.e.hat*a.re)%*%t(A.i.1))
  # part2_1 = data.frame(id=long.data$id, part2)
  part3 = part1%*%inv.Sigma.e.hat%*%t(a.re)
  # part3_1 = data.frame(id=long.data$id,part3)
  # Part1 = part1_2 %>%
  #   group_by(id) %>%
  #   summarise(sum_part1=sum(part1_1,na.rm = TRUE))%>%
  #   select(-id)
  # Part2 = part2_1 %>%
  #   group_by(id) %>%
  #   summarise(across(starts_with("X"),sum,.names = "sum_{.col}"))%>%
  #   select(-id)
  # Part3 = part3_1 %>%
  #   group_by(id) %>%
  #   summarise(across(starts_with("X"),sum,.names = "sum_{.col}"))%>%
  #   select(-id)
  part4 = (K*log(2*pi)+log(det(Sigma.e.hat)))*A.i.1
  # part4_1 = data.frame(id=long.data$id,part4+part1_1+part2-2*part3)
  Part4 = sm%*%(part4+part1_1+part2-2*part3)
  # Part4 = part4_1 %>%
  #   group_by(id) %>%
  #   summarise(across(starts_with("X"),sum,.names = "sum_{.col}")) %>%
  #   dplyr::select(-id)
  a2 = -Part4/2+lambda0.hat*V*factor
  part5 = A.i.2*c(diff(time),0)*exp(theta.m.hat%*%t(y))
  Part5 = sm%*%t(part5)
  # rownames(part5)=paste("part5")
  # part5_1 = data.frame(id=long.data$id, t(part5))
  # Part5 = part5_1 %>%
  #   group_by(id) %>%
  #   summarise(sum_part5 = sum(part5,na.rm=TRUE)) %>%
  #   dplyr::select(-id)
  factor2 = as.matrix(exp(theta.xx.mat+theta.aa.mat))*Part5
  logp = log(lambda0.hat)+log(factor)-lambda0.hat*factor2+a2
  
  logp[is.na(logp)] <- 0
  
  part6 = part1%*%inv.Sigma.e.hat%*%gamma.hat
  # part6_1 = data.frame(id=long.data$id, part6)
  # Part6 = part6_1 %>%
  #   group_by(id) %>%
  #   summarise(sum_part6 = sum(part6,na.rm=TRUE)) %>%
  #   select(-id)
  part7 = data.frame(A.i.1%*%t(a.re%*%inv.Sigma.e.hat%*%gamma.hat))
  Part7 = sm%*%as.matrix(part6-part7)
  # part7_1 = data.frame(id=long.data$id, part6-part7)
  # Part7 = part7_1 %>%
  #   group_by(id) %>%
  #   summarise(across(starts_with("X"),sum,.names = "sum_{.col}")) %>%
  #   dplyr::select(-id)
  b2 = -Part7-lambda0.hat*factor
  part8 = A.i.1%*%t(gamma.hat)%*%inv.Sigma.e.hat%*%gamma.hat
  Part8 = sm%*%part8
  # part8_1 = data.frame(id=long.data$id,part8)
  # Part8 = part8_1 %>% 
  #   group_by(id) %>%
  #   summarise(sum_part8 = sum(part8, na.rm = TRUE)) %>%
  #   dplyr::select(-id)
  # c2 = -data.frame(matrix(rep(Part8,no.gauss^K),ncol = no.gauss^K,byrow = FALSE))/2
  c2 = -Part8%*%rep(1,no.gauss^K)/2
  logpI = as.matrix(log(sqrt(pi))+log(pnorm(as.matrix(sqrt(-2*c2)*U-b2/sqrt(-2*c2)),0,1)-pnorm(as.matrix(sqrt(-2*c2)*V-b2/sqrt(-2*c2),0,1)))+(logp-b2^2/(4*c2))-log(sqrt(-c2)))
  logpI_apx = as.matrix(log(sqrt(pi))+(logp-b2^2/(4*c2))-log(sqrt(2*pi))+(-(sqrt(-2*c2)*(U+V)-2*b2/sqrt(-2*c2))^2/8)+log(sqrt(-2*c2)*(U-V))-log(sqrt(-c2)))
  logpI[is.infinite(logpI)] = logpI_apx[is.infinite(logpI)]
  logpI[is.na(logpI) | is.infinite(logpI)] = -700
  # logpI = data.frame(logpI)
  if (K>1){
    part9 = as.matrix((y-rep(1,nrow(y))%*%t(beta.0.hat)-time%*%t(beta.1.hat)-t(beta.2.hat%*%t(xval)))*A.i.3)
  } else {
    part9 = as.matrix((y-as.numeric(beta.0.hat)-time%*%t(beta.1.hat)-t(beta.2.hat%*%t(xval)))*A.i.3)
  }
  part9_1 = rowSums(part9%*%inv.Sigma.e.hat*part9)
  # part9_2 = data.frame(id=long.data$id,part9_1)
  # Part9 = part9_2 %>%
  #   group_by(id) %>%
  #   summarise(sum_part9=sum(part9_1, na.rm = TRUE)) %>%
  #   select(-id)
  part10 = t(rowSums(a.re%*%inv.Sigma.e.hat*a.re)%*%t(A.i.3))
  # part10_1 = data.frame(id=long.data$id, part10)
  # Part10 = part10_1 %>%
  #   group_by(id) %>%
  #   summarise(across(starts_with("X"),sum,.names = "sum_{.col}"))%>%
  #   select(-id)
  part11 = part9%*%inv.Sigma.e.hat%*%t(a.re)
  # part11_1 = data.frame(id=long.data$id,part11)
  # Part11 = part3_1 %>%
  #   group_by(id) %>%
  #   summarise(across(starts_with("X"),sum,.names = "sum_{.col}"))%>%
  #   select(-id)
  part12 = (K*log(2*pi)+log(det(Sigma.e.hat)))*A.i.3
  # part12_1 = data.frame(id=long.data$id,part12)
  # Part12 = part12_1 %>%
  #   group_by(id) %>%
  #   summarise(sum_part4=sum(part12,na.rm = TRUE)) %>%
  #   select(-id)
  Part13 = sm%*%(part12+part9_1+part10-2*part11)
  # part13_1 = data.frame(id=long.data$id, part12+part9_1+part10-2*part11)
  # Part13 = part13_1%>%
  #   group_by(id) %>%
  #   summarise(across(starts_with("X"),sum,.names = "sum_{.col}"))%>%
  #   dplyr::select(-id)
  logq1 = -Part13/2
  logq2 = logq1 - lambda0.hat*factor2
  f = exp(status*(logpI+logq1)+(1-status)*logq2)
  
  Q = as.matrix(t(f%*%weight))/(sqrt(pi)^K)
  
  return(sum(log(Q)))
}

likelihood.piecewise = function(long.data, parameters)
{
  long.data$id = as.factor(long.data$id)
  id.matrix <- model.matrix(~ id - 1, data = long.data)
  sm = Matrix(t(id.matrix),sparse = T)
  n.id = length(unique(long.data$id))
  long.data.id = long.data[!duplicated(long.data$id),]
  beta.0.hat =  t(as.matrix(parameters %>% dplyr::select(contains("beta.0"))))
  K = length(beta.0.hat)
  beta.1.hat = t(as.matrix(parameters %>% dplyr::select(contains("beta.1"))))
  beta.2.vec = data.frame(parameters %>% dplyr::select(contains("beta.2")))
  col_names <- colnames(beta.2.vec)
  if (length(col_names)==1){
    beta.2.hat=as.matrix(beta.2.vec)
    p = ncol(beta.2.hat)
  } else{
    indices <- do.call(rbind, strsplit(col_names, "_"))[, 2:3]
    indices <- apply(indices, 2, as.integer)
    beta.2.hat = matrix(0, nrow=max(indices[,1]),ncol=max(indices[,2]))
    p=ncol(beta.2.hat)
    for (i in seq_along(col_names)){
      row_idx <- indices[i, 1]
      col_idx <- indices[i, 2]
      beta.2.hat[row_idx, col_idx] <- beta.2.vec[[col_names[i]]]
    }
  }
  gamma.hat = t(as.matrix(parameters %>% dplyr::select(contains("gamma"))))
  
  # lambda0.hat = as.numeric(parameters %>% dplyr::select(contains("lambda0")))
  lambda.hat = as.matrix(parameters %>% dplyr::select(contains("hazard")))
  theta.x.hat = as.matrix(parameters %>% dplyr::select(contains("theta.x")))
  theta.a.hat = as.matrix(parameters %>% dplyr::select(contains("theta.a")))
  theta.m.hat = as.matrix(parameters %>% dplyr::select(contains("theta.m")))
  
  Sigma.a.hat = matrix(0, nrow=K,ncol=K)
  upper_vec_a = as.matrix(parameters %>% dplyr::select(contains("Sigma.a")))
  Sigma.a.hat[upper.tri(Sigma.a.hat,diag = TRUE)] = upper_vec_a
  if (K>1){
    Sigma.a.hat=Sigma.a.hat+t(Sigma.a.hat)-diag(diag(Sigma.a.hat))
  }
  Sigma.e.hat = matrix(0, nrow=K,ncol=K)
  upper_vec_e = as.matrix(parameters %>% dplyr::select(contains("Sigma.e")))
  Sigma.e.hat[upper.tri(Sigma.e.hat,diag = TRUE)] = upper_vec_e
  if (K>1){
    Sigma.e.hat=Sigma.e.hat+t(Sigma.e.hat)-diag(diag(Sigma.e.hat))
  }
  inv.Sigma.e.hat = solve(Sigma.e.hat)
  
  a.re = nodes%*%sqrtm(2*Sigma.a.hat)
  
  y = long.data %>% dplyr::select(contains("ynew"))
  xval = long.data %>% dplyr::select(contains("xval"))
  xval.id = long.data.id %>% dplyr::select(contains("xval"))
  time = long.data$visittime
  status = long.data.id$status
  V = long.data.id$preetime
  U = long.data.id$etime
  A.i.1 = as.numeric(long.data$visit>long.data$befdiag)
  A.i.2 = as.numeric(long.data$visit<long.data$befdiag)
  A.i.3 = as.numeric(long.data$visit<=long.data$befdiag)
  A.i.4 = as.numeric(long.data$visit==long.data$befdiag)
  theta.xx = theta.x.hat%*%t(xval.id) #400 for individual
  theta.xx.mat = data.frame(matrix(rep(t(theta.xx),no.gauss^K),ncol = no.gauss^K,byrow=FALSE))
  theta.aa = theta.a.hat%*%t(a.re) #400 for nodes
  theta.aa.mat = data.frame(matrix(rep(theta.aa,n.id),nrow = n.id,byrow = TRUE))
  theta.mm = A.i.4[A.i.4!=0]*theta.m.hat%*%t(y[A.i.4!=0,]) #400 for individual
  theta.mm.mat = data.frame(matrix(rep(t(theta.mm),no.gauss^K),ncol=no.gauss^K,byrow=FALSE))
  # factor = as.matrix(exp(theta.xx.mat+theta.aa.mat+theta.mm.mat))
  
  ## creating dataset including 
  interval = d[-length(d)]
  fixed_grid <- expand.grid(id = unique(long.data$id), visittime = interval)
  long.data_expanded <- full_join(long.data, fixed_grid, by = c("id", "visittime"))
  long.data_expanded <- long.data_expanded %>%
    arrange(id, visittime) %>%
    fill(everything(), .direction = "down")
  ## creating piecewise hazard in dataset
  
  hazard_intervals = data.frame(
    start = as.vector(interval),
    end = c(interval[-1], Inf),
    hazard = as.vector(lambda.hat)
  )
  long.data_expanded <- long.data_expanded %>%
    rowwise() %>%
    mutate(
      hazard = hazard_intervals$hazard[
        visittime >= hazard_intervals$start & visittime < hazard_intervals$end
      ]
    )
  ## defining quantum using new dataset
  A.i.5 = as.numeric(long.data_expanded$visittime<long.data_expanded$etime)
  A.i.6 = as.numeric(long.data_expanded$visittime<long.data_expanded$preetime)
  # A.i.7 = as.numeric(long.data_expanded$visit<=long.data_expanded$befdiag+1)
  
  long.data_expanded$id = as.factor(long.data_expanded$id)
  id.matrix_expanded <- model.matrix(~ id - 1, data = long.data_expanded)
  sm_expanded = Matrix(t(id.matrix_expanded),sparse = T)
  
  
  y_expanded = long.data_expanded %>% dplyr::select(contains("ynew"))
  xval_expanded = long.data_expanded %>% dplyr::select(contains("xval"))
  time_expanded = long.data_expanded$visittime
  hazard = long.data_expanded$hazard
  
  theta.xx_expanded = theta.x.hat%*%t(xval_expanded) #400 for individual
  theta.xx.mat_expanded = data.frame(matrix(rep(t(theta.xx_expanded),no.gauss^K),ncol = no.gauss^K,byrow=FALSE))
  theta.aa.mat_expanded = data.frame(matrix(rep(theta.aa,nrow(long.data_expanded)),nrow = nrow(long.data_expanded),byrow = TRUE))
  theta.mm_expanded = theta.m.hat%*%t(y_expanded) #400 for individual
  theta.mm.mat_expanded = data.frame(matrix(rep(t(theta.mm_expanded),no.gauss^K),ncol=no.gauss^K,byrow=FALSE))
  factor_expanded = exp(theta.xx.mat_expanded+theta.aa.mat_expanded+theta.mm.mat_expanded)
  factor_expanded = Matrix(as.matrix(factor_expanded), sparse = T)
  
  ##Calculating a2k, b2k, pk
  U_expanded <- long.data_expanded$etime
  time_re <- time_expanded
  time_re[time_re >= U_expanded] <- U_expanded[time_re >= U_expanded]
  time_diff <- c(diff(time_re),0)
  time_diff[time_diff < 0] <- 0 
  fe_full = -hazard*time_diff*factor_expanded
  fe_full = rbind(0,fe_full)
  fe_full = fe_full[-nrow(fe_full),]
  bk = -A.i.5*hazard*factor_expanded
  
  id2 = t(sm_expanded)%*%sm_expanded
  id2 = tril(id2)
  # time-consuming!
  fe <- id2%*%fe_full
  
  ak = fe - bk*time_expanded
  
  # ak = Matrix(as.matrix(ak),sparse = T)
  # bk = Matrix(as.matrix(bk),sparse = T)
  
  id_counts <- table(long.data_expanded$id)
  theta.mm.mat_expanded_1 = theta.mm.mat[rep(1:nrow(theta.mm.mat), times = id_counts), , drop = FALSE]
  factor_expanded_1 = exp(theta.xx.mat_expanded+theta.aa.mat_expanded+theta.mm.mat_expanded_1)
  factor_expanded_1 = Matrix(as.matrix(factor_expanded_1),sparse=T)
  pk = hazard*factor_expanded_1
  # pk[is.na(pk)] <- 1
  ##calculating a1, b1, c1 and expand them
  if (K>1){
    part1 = as.matrix((y-rep(1,nrow(y))%*%t(beta.0.hat)-time%*%t(beta.1.hat)-t(beta.2.hat%*%t(xval))-time%*%t(gamma.hat))*A.i.1)
  } else {
    part1 = as.matrix((y-as.numeric(beta.0.hat)-time%*%t(beta.1.hat)-t(beta.2.hat%*%t(xval))-time%*%t(gamma.hat))*A.i.1)
  }
  part1_1 = rowSums(part1%*%inv.Sigma.e.hat*part1)
  
  part2 = t(rowSums(a.re%*%inv.Sigma.e.hat*a.re)%*%t(A.i.1))
  
  part3 = part1%*%inv.Sigma.e.hat%*%t(a.re)
  
  part4 = (K*log(2*pi)+log(det(Sigma.e.hat)))*A.i.1
  
  Part4 = sm%*%(part4+part1_1+part2-2*part3)
  a1 = -Part4/2
  a1_expanded = a1[rep(1:nrow(a1), times = id_counts), , drop = FALSE]
  
  part6 = part1%*%inv.Sigma.e.hat%*%gamma.hat
  
  part7 = data.frame(A.i.1%*%t(a.re%*%inv.Sigma.e.hat%*%gamma.hat))
  Part7 = sm%*%as.matrix(part6-part7)
  
  b1 = -Part7
  b1_expanded = b1[rep(1:nrow(b1), times = id_counts), , drop = FALSE]
  
  part8 = A.i.1%*%t(gamma.hat)%*%inv.Sigma.e.hat%*%gamma.hat
  Part8 = sm%*%part8
  c1 = -Part8%*%rep(1,no.gauss^K)/2
  ck_1 = c1[rep(1:nrow(c1), times = id_counts), drop = FALSE]
  
  ak_1 = ak+a1_expanded
  bk_1 = bk+b1_expanded
  
  time_next <- ave(time_expanded, long.data_expanded$id, FUN = function(x) c(x[-1], x[length(x)]))
  U_1 = pmax(pmin(long.data_expanded$etime,time_next), long.data_expanded$preetime)
  V_1 = pmin(pmax(long.data_expanded$preetime,time_expanded), long.data_expanded$etime)
  
  pI = pk*sqrt(pi)*exp(ak_1-bk_1^2/(4*ck_1))/sqrt(-ck_1)*(pnorm((sqrt(-2*ck_1)*U_1-bk_1/sqrt(-2*ck_1))@x,0,1)-pnorm((sqrt(-2*ck_1)*V_1-bk_1/sqrt(-2*ck_1))@x,0,1))
  pI_apx = pk*sqrt(pi)*exp(ak_1-bk_1^2/(4*ck_1))/sqrt(-ck_1)*(1/sqrt(2*pi)*exp(-((sqrt(-2*ck_1)*(U_1+V_1)-2*bk_1/sqrt(-2*ck_1))/2)^2/2))*(sqrt(-2*ck_1)*(U_1-V_1))
  logpI = log(sm_expanded%*%pI)
  logpI_apx = log(sm_expanded%*%pI_apx)
  logpI[is.infinite(logpI)] = logpI_apx[is.infinite(logpI)]
  logpI[is.na(logpI) | is.infinite(logpI)] = -700
  #q
  if (K>1){
    part9 = as.matrix((y-rep(1,nrow(y))%*%t(beta.0.hat)-time%*%t(beta.1.hat)-t(beta.2.hat%*%t(xval)))*A.i.3)
  } else {
    part9 = as.matrix((y-as.numeric(beta.0.hat)-time%*%t(beta.1.hat)-t(beta.2.hat%*%t(xval)))*A.i.3)
  }
  part9_1 = rowSums(part9%*%inv.Sigma.e.hat*part9)
  
  part10 = t(rowSums(a.re%*%inv.Sigma.e.hat*a.re)%*%t(A.i.3))
  
  part11 = part9%*%inv.Sigma.e.hat%*%t(a.re)
  
  part12 = (K*log(2*pi)+log(det(Sigma.e.hat)))*A.i.3
  
  Part13 = sm%*%(part12+part9_1+part10-2*part11)
  
  logq1 = -Part13/2
  
  factor2_expanded = - A.i.6*hazard*c(diff(time_expanded),0)*factor_expanded
  factor2_1 = Matrix(as.matrix(factor2_expanded),sparse=T)
  logq2 = sm_expanded%*%factor2_1
  
  f = exp(status*(logpI+logq1)+(1-status)*(logq1+logq2))
  
  Q = as.matrix(t(f%*%weight))/(sqrt(pi)^K)
  
  return(list(ll=sum(log(Q)), ll.ind=log(Q)))
}

#############################
# 11/17/2025
# using M-spline and I-spline
evaluate_Mspline <- function(x, knots, alpha, boundary_knots, degree = 3) {
  # Step 1: Generate the B-spline basis matrix.
  beta = exp(alpha)
  basis_matrix <- mSpline(x, knots = knots, degree = degree, intercept = TRUE, Boundary.knots = boundary_knots)

  # Step 2: Validate the inputs.
  if (length(alpha) != ncol(basis_matrix)) {
    stop(paste(
      "Error: The number of alpha coefficients does not match the number of basis functions.",
      "\nProvided", length(alpha), "coefficients.",
      "\nExpected", ncol(basis_matrix), "basis functions for the given knots and degree."
    ))
  }

  # Step 3: Calculate the spline value.
  spline_value <- as.vector(basis_matrix %*% beta)

  return(spline_value)
}

evaluate_Mspline_matrix <- function(X, knots, alpha, boundary_knots, degree = 3) {
  # --- Step 1: Pre-computation and Validation ---
  beta = exp(alpha)
  # Store the original dimensions of the input matrix.
  original_dims <- dim(X)
  if (is.null(original_dims)) {
    stop("Error: Input 'X' must be a matrix or a data frame.")
  }

  # The bs() function requires a vector input. We unroll the matrix.
  x_vector <- as.vector(as.matrix(X))

  # Generate the basis matrix for the vector.
  basis_matrix <- mSpline(x_vector, knots = knots, degree = degree, intercept = TRUE, Boundary.knots = boundary_knots)

  # Validate that the number of alpha coefficients matches the basis functions.
  if (length(alpha) != ncol(basis_matrix)) {
    stop(paste(
      "Error: The number of alpha coefficients does not match the number of basis functions.",
      "\nProvided", length(alpha), "coefficients.",
      "\nExpected", ncol(basis_matrix), "basis functions for the given knots and degree."
    ))
  }

  # --- Step 2: Core Calculation ---

  # Calculate the spline values for the unrolled vector.
  spline_values_vector <- as.vector(basis_matrix %*% beta)

  # --- Step 3: Reshape the Output ---

  # Reshape the resulting vector back into a matrix with the original dimensions.
  output_matrix <- matrix(spline_values_vector, nrow = original_dims[1], ncol = original_dims[2])

  return(output_matrix)
}

evaluate_Ispline <- function(x, knots, alpha, boundary_knots, degree = 3) {
  # Step 1: Generate the B-spline basis matrix.
  beta = exp(alpha)
  basis_matrix <- iSpline(x, knots = knots, degree = degree, intercept = TRUE, Boundary.knots = boundary_knots)

  # Step 2: Validate the inputs.
  if (length(alpha) != ncol(basis_matrix)) {
    stop(paste(
      "Error: The number of alpha coefficients does not match the number of basis functions.",
      "\nProvided", length(alpha), "coefficients.",
      "\nExpected", ncol(basis_matrix), "basis functions for the given knots and degree."
    ))
  }

  # Step 3: Calculate the spline value.
  spline_value <- as.vector(basis_matrix %*% beta)
  # basis_at_zero <- ibs(0,
  #                     knots = knots,
  #                     degree = degree,
  #                     intercept = TRUE,
  #                     Boundary.knots = boundary_knots)
  # value_at_zero <- as.vector(basis_at_zero %*% beta)
  # anchored_hazard_value <- spline_value - value_at_zero
  return(spline_value)
}
###############################
# 09/25/2025
# using B-spline to approximate cumulative hazard, non-decreasing
spline_cumulative_hazard <- function(x, knots, alpha, boundary_knots, degree) {
  
  # --- 1. Transform Alphas to Monotonically Increasing Betas ---
  beta <- cumsum(exp(alpha))
  
  # --- 2. Define Boundary Knots ---
  # We explicitly set the lower boundary to 0 to calculate the offset correctly.
  # boundary_knots <- c(0, max(x, na.rm = TRUE))
  
  # --- 3. Generate the Standard B-spline Basis Matrix ---
  basis_matrix <- bs(x, 
                     knots = knots, 
                     degree = degree, 
                     intercept = TRUE, 
                     Boundary.knots = boundary_knots)
  
  # --- 4. Validation ---
  if (length(beta) != ncol(basis_matrix)) {
    stop(paste(
      "Error: The number of transformed beta coefficients does not match the number of basis functions.",
      "\nProvided", length(beta), "coefficients.",
      "\nExpected", ncol(basis_matrix), "basis functions for the given knots and degree."
    ))
  }
  
  # --- 5. Calculate the Un-anchored Cumulative Hazard ---
  cumulative_hazard_value <- as.vector(basis_matrix %*% beta)

  # --- 6. ANCHORING: Calculate the value at t=0 to find the offset ---
  # We create the basis at t=0 using the same boundary knots for consistency.
  basis_at_zero <- bs(0,
                      knots = knots,
                      degree = degree,
                      intercept = TRUE,
                      Boundary.knots = boundary_knots)
  value_at_zero <- as.vector(basis_at_zero %*% beta)

  # --- 7. ANCHORING: Subtract the offset to ensure H(0) = 0 ---
  # This preserves the non-decreasing shape of the function.
  anchored_hazard_value <- cumulative_hazard_value - value_at_zero
  
  # Return the anchored values, ensuring no small negative numbers due to floating point error.
  return(pmax(0, anchored_hazard_value))
}

spline_hazard <- function(x, knots, alpha, boundary_knots, degree) {
  # --- 1. Transform Alphas to Monotonically Increasing Betas ---
  # This must be identical to the transformation in the cumulative hazard function.
  beta <- cumsum(exp(alpha))
  
  # --- 2. Define Boundary Knots for Consistency ---
  # boundary_knots <- c(0, max(x, na.rm = TRUE))
  
  # --- 3. Generate the B-spline DERIVATIVE Basis Matrix ---
  # dbs() from the splines2 package gives the basis for the derivative.
  derivative_basis_matrix <- dbs(x, 
                                 knots = knots, 
                                 degree = degree, 
                                 intercept = TRUE, 
                                 Boundary.knots = boundary_knots)
  
  # --- 4. Calculate the Hazard Value ---
  # The derivative of the constant offset S(0) is zero, so we don't need anchoring.
  hazard_value <- as.vector(derivative_basis_matrix %*% beta)
  
  # The hazard must be non-negative.
  return(pmax(0, hazard_value))
}

spline_hazard_matrix <- function(X, knots, alpha, boundary_knots, degree) {
  
  original_dims <- dim(X)
  if (is.null(original_dims)) {
    stop("Error: Input 'X' must be a matrix or a data frame.")
  }
  
  # The bs() function requires a vector input. We unroll the matrix.
  x_vector <- as.vector(as.matrix(X))
  # --- 1. Transform Alphas to Monotonically Increasing Betas ---
  # This must be identical to the transformation in the cumulative hazard function.
  beta <- cumsum(exp(alpha))
  
  # --- 2. Define Boundary Knots for Consistency ---
  # boundary_knots <- c(0, max(x, na.rm = TRUE))
  
  # --- 3. Generate the B-spline DERIVATIVE Basis Matrix ---
  # dbs() from the splines2 package gives the basis for the derivative.
  derivative_basis_matrix <- dbs(x_vector, 
                                 knots = knots, 
                                 degree = degree, 
                                 intercept = TRUE, 
                                 Boundary.knots = boundary_knots)
  
  # --- 4. Calculate the Hazard Value ---
  # The derivative of the constant offset S(0) is zero, so we don't need anchoring.
  hazard_value <- as.vector(derivative_basis_matrix %*% beta)
  
  output_matrix <- matrix(pmax(0, hazard_value), nrow = original_dims[1], ncol = original_dims[2])
  # The hazard must be non-negative.
  return(output_matrix)
}
##############################
# use cumsum(exp(alpha)) to calculate cumulative hazard
likelihood.spline = function(long.data, parameters, knots){
  long.data$id = as.factor(long.data$id)
  id.matrix <- model.matrix(~ id - 1, data = long.data)
  sm = Matrix(t(id.matrix),sparse = T)
  n.id = length(unique(long.data$id))
  long.data.id = long.data[!duplicated(long.data$id),]
  beta.0.hat =  t(as.matrix(parameters %>% dplyr::select(contains("beta.0"))))
  K = length(beta.0.hat)
  beta.1.hat = t(as.matrix(parameters %>% dplyr::select(contains("beta.1"))))
  beta.2.vec = data.frame(parameters %>% dplyr::select(contains("beta.2")))
  col_names <- colnames(beta.2.vec)
  if (length(col_names)==1){
    beta.2.hat=as.matrix(beta.2.vec)
    p = ncol(beta.2.hat)
  } else{
    indices <- do.call(rbind, strsplit(col_names, "_"))[, 2:3]
    indices <- apply(indices, 2, as.integer)
    beta.2.hat = matrix(0, nrow=max(indices[,1]),ncol=max(indices[,2]))
    p=ncol(beta.2.hat)
    for (i in seq_along(col_names)){
      row_idx <- indices[i, 1]
      col_idx <- indices[i, 2]
      beta.2.hat[row_idx, col_idx] <- beta.2.vec[[col_names[i]]]
    }
  }
  gamma.hat = t(as.matrix(parameters %>% dplyr::select(contains("gamma"))))
  
  # lambda0.hat = as.numeric(parameters %>% dplyr::select(contains("lambda0")))
  # lambda.hat = as.matrix(parameters %>% dplyr::select(contains("lambda")))
  alpha.hat = t(as.matrix(parameters %>% dplyr::select(contains("alpha"))))
  theta.x.hat = as.matrix(parameters %>% dplyr::select(contains("theta.x")))
  theta.a.hat = as.matrix(parameters %>% dplyr::select(contains("theta.a")))
  theta.m.hat = as.matrix(parameters %>% dplyr::select(contains("theta.m")))
  
  Sigma.a.hat = matrix(0, nrow=K,ncol=K)
  upper_vec_a = as.matrix(parameters %>% dplyr::select(contains("Sigma.a")))
  Sigma.a.hat[upper.tri(Sigma.a.hat,diag = TRUE)] = upper_vec_a
  if (K>1){
    Sigma.a.hat=Sigma.a.hat+t(Sigma.a.hat)-diag(diag(Sigma.a.hat))
  }
  Sigma.e.hat = matrix(0, nrow=K,ncol=K)
  upper_vec_e = as.matrix(parameters %>% dplyr::select(contains("Sigma.e")))
  Sigma.e.hat[upper.tri(Sigma.e.hat,diag = TRUE)] = upper_vec_e
  if (K>1){
    Sigma.e.hat=Sigma.e.hat+t(Sigma.e.hat)-diag(diag(Sigma.e.hat))
  }
  inv.Sigma.e.hat = solve(Sigma.e.hat)
  
  a.re = nodes%*%sqrtm(2*Sigma.a.hat)
  
  y = long.data %>% dplyr::select(contains("ynew"))
  xval = long.data %>% dplyr::select(contains("xval"))
  xval.id = long.data.id %>% dplyr::select(contains("xval"))
  time = long.data$visittime
  status = long.data.id$status
  V = long.data.id$preetime
  U = long.data.id$etime
  
  A.i.1 = as.numeric(long.data$visit>long.data$befdiag)
  A.i.2 = as.numeric(long.data$visit<long.data$befdiag)
  A.i.3 = as.numeric(long.data$visit<=long.data$befdiag)
  A.i.4 = as.numeric(long.data$visit==long.data$befdiag)
  
  boundary_knots <- c(0, max(time))
  ##calculating a1, b1, c1 and expand them
  if (K>1){
    part1 = as.matrix((y-rep(1,nrow(y))%*%t(beta.0.hat)-time%*%t(beta.1.hat)-t(beta.2.hat%*%t(xval))-time%*%t(gamma.hat))*A.i.1)
  } else {
    part1 = as.matrix((y-as.numeric(beta.0.hat)-time%*%t(beta.1.hat)-t(beta.2.hat%*%t(xval))-time%*%t(gamma.hat))*A.i.1)
  }
  part1_1 = rowSums(part1%*%inv.Sigma.e.hat*part1)
  
  part2 = t(rowSums(a.re%*%inv.Sigma.e.hat*a.re)%*%t(A.i.1))
  
  part3 = part1%*%inv.Sigma.e.hat%*%t(a.re)
  
  part4 = (K*log(2*pi)+log(det(Sigma.e.hat)))*A.i.1
  
  Part4 = sm%*%(part4+part1_1+part2-2*part3)
  a1 = -Part4/2
  # a1_expanded = a1[rep(1:nrow(a1), times = id_counts), , drop = FALSE]
  
  part6 = part1%*%inv.Sigma.e.hat%*%gamma.hat
  
  part7 = data.frame(A.i.1%*%t(a.re%*%inv.Sigma.e.hat%*%gamma.hat))
  Part7 = sm%*%as.matrix(part6-part7)
  
  b1 = -Part7
  # b1_expanded = b1[rep(1:nrow(b1), times = id_counts), , drop = FALSE]
  
  part8 = A.i.1%*%t(gamma.hat)%*%inv.Sigma.e.hat%*%gamma.hat
  Part8 = sm%*%part8
  c1 = -Part8%*%rep(1,no.gauss^K)/2
  
  ## generate Monte Carlo time to integrate, batch size = 100, total sample = 1000
  intervals = cbind(V, U)
  MC_sample = 1000
  batch = 10
  batch_sample = MC_sample/batch
  
  ## compute \sum \int lambda(s)*exp(\theta_m*M)ds
  rows_to_insert <- long.data %>%
    filter(visit == befdiag)
  new_rows <- rows_to_insert %>%
    mutate(visit := visit + 0.1)
  theta.mm_expanded = theta.m.hat%*%t(y)
  
  #theta.mm at V
  theta.mm = t(A.i.4[A.i.4!=0]*theta.m.hat%*%t(y[A.i.4!=0,]))
  # \sum \int lambda(s)*exp(\theta_m*M)ds before V
  spline_time_int = spline_cumulative_hazard(time, knots, alpha.hat, boundary_knots, degree)
  lambda.mm1 = sm%*%(A.i.2*diff(c(spline_time_int,0))*t(exp(theta.mm_expanded)))
  
  # theta.xx and theta.aa
  theta.xx = theta.x.hat %*% t(xval.id)
  theta.aa = theta.a.hat %*% t(a.re)
  theta.xx_1 = kronecker(t(rep(1,no.gauss^K)),t(theta.xx))
  theta.aa_1 = kronecker(rep(1,n.id), theta.aa)
  # reformat
  theta.xx_full = theta.xx_1[, rep(1:ncol(theta.xx_1), each = batch_sample)]
  theta.aa_full = theta.aa_1[, rep(1:ncol(theta.aa_1), each = batch_sample)]
  c1_full = c1[, rep(1:ncol(c1), each = batch_sample)]
  b1_full = b1[, rep(1:ncol(b1), each = batch_sample)]
  int_with_t = 0
  for (k in 1:batch){
    u <- runif(n.id * batch_sample)   # uniform(0,1)
    a <- rep(intervals[, 1], each = batch_sample)
    b <- rep(intervals[, 2], each = batch_sample)
    sample.time <- t(matrix(a + (b - a) * u, nrow = batch_sample, ncol = n.id))
    
    #long.data.expanded <- bind_rows(long.data, new_rows) %>%
    # arrange(id, visit)
    #y_expanded <- long.data.expanded %>% dplyr::select(contains("ynew"))
    time_expanded = matrix(rows_to_insert$visittime)[,rep(1, each = batch_sample)]
    id_time = data.frame(id = rows_to_insert$id, visit = rows_to_insert$visit, time_expanded)
    id_MCtime = data.frame(id = long.data.id$id, visit = new_rows$visit, sample.time)
    id_time_expanded <- bind_rows(id_time, id_MCtime) %>%
      arrange(id, visit)
    
    
    # \sum \int lambda(s)*exp(\theta_m*M)ds from V to MC.time
    
    MC_time_expanded <- id_time_expanded %>% dplyr::select(contains("X"))
    MC_time.vector = as.vector(as.matrix(MC_time_expanded))
    MC_spline_int = matrix(spline_cumulative_hazard(MC_time.vector, knots, alpha.hat, boundary_knots, degree),nrow = n.id*2, ncol = batch_sample)
    int_matrix <- Matrix(kronecker(diag(n.id), t(c(-1,1))), sparse = T)
    MC_spline <- int_matrix%*%MC_spline_int
    lambda.mm2 = MC_spline * exp(theta.mm)
    
    #total lambda * exp(M)
    lambda.mm <- lambda.mm1 + status*lambda.mm2
    # lambda.mm[lambda.mm<0] = 1 # censor case
    
    # reformat
    # theta.xx_full = theta.xx_1[, rep(1:ncol(theta.xx_1), each = ncol(lambda.mm))]
    # theta.aa_full = theta.aa_1[, rep(1:ncol(theta.aa_1), each = ncol(lambda.mm))]
    MC_time_full = sample.time[, rep(1:ncol(sample.time), times = ncol(theta.xx_1))]
    lambda.mm_full = lambda.mm[, rep(1:ncol(lambda.mm), times = ncol(theta.xx_1))]
    factor1_full = exp(-lambda.mm_full*exp(theta.xx_full+theta.aa_full))
    lambda_t = spline_hazard_matrix(sample.time, knots, alpha.hat, boundary_knots, degree)
    lambda_t_full = lambda_t[, rep(1:ncol(lambda_t), times = ncol(theta.xx_1))]
    # c1_full = c1[, rep(1:ncol(c1), each = ncol(sample.time))]
    # b1_full = b1[, rep(1:ncol(b1), each = ncol(sample.time))]
    # calculate MC integrate
    MC_integrate = exp(b1_full*MC_time_full+c1_full*MC_time_full^2)*lambda_t_full*factor1_full
    averaging_matrix <- Matrix(kronecker(
      diag(no.gauss^K), rep(1/batch_sample, batch_sample)
    ), sparse = T)
    int_with_t = int_with_t + MC_integrate%*%averaging_matrix * (U-V)/batch
    gc()
  }
  
  # calculate pI
  theta.mm1 = kronecker(t(rep(1,no.gauss^K)),theta.mm)
  pI = int_with_t * exp(a1) * exp(theta.xx_1+theta.aa_1+theta.mm1)
  logpI = log(pI)
  logpI[is.infinite(logpI)] = -700
  # q1
  if (K>1){
    part9 = as.matrix((y-rep(1,nrow(y))%*%t(beta.0.hat)-time%*%t(beta.1.hat)-t(beta.2.hat%*%t(xval)))*A.i.3)
  } else {
    part9 = as.matrix((y-as.numeric(beta.0.hat)-time%*%t(beta.1.hat)-t(beta.2.hat%*%t(xval)))*A.i.3)
  }
  part9_1 = rowSums(part9%*%inv.Sigma.e.hat*part9)
  
  part10 = t(rowSums(a.re%*%inv.Sigma.e.hat*a.re)%*%t(A.i.3))
  
  part11 = part9%*%inv.Sigma.e.hat%*%t(a.re)
  
  part12 = (K*log(2*pi)+log(det(Sigma.e.hat)))*A.i.3
  
  Part13 = sm%*%(part12+part9_1+part10-2*part11)
  
  logq1 = -Part13/2
  # q2
  lambda.mm3 = sm%*%t(A.i.2*diff(c(spline_time_int,0))*exp(theta.mm_expanded))
  lambda.mm3_full = kronecker(t(rep(1,no.gauss^K)),lambda.mm3)
  logq2 = -lambda.mm3_full* exp(theta.xx_1+theta.aa_1)
  # f and Q
  f = exp(status*(logpI+logq1)+(1-status)*(logq1+logq2))
  
  Q = as.matrix(t(f%*%weight))/(sqrt(pi)^K)
  
  return(list(ll=sum(log(Q)), ll.ind=log(Q)))
}
# Gauss-Legendre integration
# 11/06/2025: x time-dependent
likelihood.spline2 = function(long.data, parameters, knots){
  long.data$id = as.factor(long.data$id)
  id.matrix <- model.matrix(~ id - 1, data = long.data)
  sm = Matrix(t(id.matrix),sparse = T)
  n.id = length(unique(long.data$id))
  long.data.id = long.data[!duplicated(long.data$id),]
  beta.0.hat =  t(as.matrix(parameters %>% dplyr::select(contains("beta.0"))))
  K = length(beta.0.hat)
  beta.1.hat = t(as.matrix(parameters %>% dplyr::select(contains("beta.1"))))
  beta.2.vec = data.frame(parameters %>% dplyr::select(contains("beta.2")))
  col_names <- colnames(beta.2.vec)
  if (length(col_names)==1){
    beta.2.hat=as.matrix(beta.2.vec)
    p = ncol(beta.2.hat)
  } else{
    indices <- do.call(rbind, strsplit(col_names, "_"))[, 2:3]
    indices <- apply(indices, 2, as.integer)
    beta.2.hat = matrix(0, nrow=max(indices[,1]),ncol=max(indices[,2]))
    p=ncol(beta.2.hat)
    for (i in seq_along(col_names)){
      row_idx <- indices[i, 1]
      col_idx <- indices[i, 2]
      beta.2.hat[row_idx, col_idx] <- beta.2.vec[[col_names[i]]]
    }
  }
  gamma.hat = t(as.matrix(parameters %>% dplyr::select(contains("gamma"))))
  
  # lambda0.hat = as.numeric(parameters %>% dplyr::select(contains("lambda0")))
  # lambda.hat = as.matrix(parameters %>% dplyr::select(contains("lambda")))
  alpha.hat = t(as.matrix(parameters %>% dplyr::select(contains("alpha"))))
  theta.x.hat = as.matrix(parameters %>% dplyr::select(contains("theta.x")))
  theta.a.hat = as.matrix(parameters %>% dplyr::select(contains("theta.a")))
  theta.m.hat = as.matrix(parameters %>% dplyr::select(contains("theta.m")))
  
  Sigma.a.hat = matrix(0, nrow=K,ncol=K)
  upper_vec_a = as.matrix(parameters %>% dplyr::select(contains("Sigma.a")))
  Sigma.a.hat[upper.tri(Sigma.a.hat,diag = TRUE)] = upper_vec_a
  if (K>1){
    Sigma.a.hat=Sigma.a.hat+t(Sigma.a.hat)-diag(diag(Sigma.a.hat))
  }
  Sigma.e.hat = matrix(0, nrow=K,ncol=K)
  upper_vec_e = as.matrix(parameters %>% dplyr::select(contains("Sigma.e")))
  Sigma.e.hat[upper.tri(Sigma.e.hat,diag = TRUE)] = upper_vec_e
  if (K>1){
    Sigma.e.hat=Sigma.e.hat+t(Sigma.e.hat)-diag(diag(Sigma.e.hat))
  }
  inv.Sigma.e.hat = solve(Sigma.e.hat)
  
  L = chol(Sigma.a.hat)
  a.re = sqrt(2)*nodes%*%t(L)
  
  y = long.data %>% dplyr::select(contains("ynew"))
  xval = long.data %>% dplyr::select(contains("xval"))
  # xval.id = long.data.id %>% dplyr::select(contains("xval"))
  time = long.data$visittime
  status = long.data.id$status
  V = long.data.id$preetime
  U = long.data.id$etime
  
  A.i.1 = as.numeric(long.data$visit>long.data$befdiag)
  A.i.2 = as.numeric(long.data$visit<long.data$befdiag)
  A.i.3 = as.numeric(long.data$visit<=long.data$befdiag)
  A.i.4 = as.numeric(long.data$visit==long.data$befdiag)
  
  # boundary_knots <- c(0, max(time))
  ##calculating a1, b1, c1 and expand them
  if (K>1){
    part1 = as.matrix((y-rep(1,nrow(y))%*%t(beta.0.hat)-time%*%t(beta.1.hat)-t(beta.2.hat%*%t(xval))-time%*%t(gamma.hat))*A.i.1)
  } else {
    part1 = as.matrix((y-as.numeric(beta.0.hat)-time%*%t(beta.1.hat)-t(beta.2.hat%*%t(xval))-time%*%t(gamma.hat))*A.i.1)
  }
  part1_1 = rowSums(part1%*%inv.Sigma.e.hat*part1)
  
  part2 = t(rowSums(a.re%*%inv.Sigma.e.hat*a.re)%*%t(A.i.1))
  
  part3 = part1%*%inv.Sigma.e.hat%*%t(a.re)
  
  part4 = (K*log(2*pi)+log(det(Sigma.e.hat)))*A.i.1
  
  Part4 = sm%*%(part4+part1_1+part2-2*part3)
  a1 = -Part4/2
  # a1_expanded = a1[rep(1:nrow(a1), times = id_counts), , drop = FALSE]
  
  part6 = part1%*%inv.Sigma.e.hat%*%gamma.hat
  
  part7 = data.frame(A.i.1%*%t(a.re%*%inv.Sigma.e.hat%*%gamma.hat))
  Part7 = sm%*%as.matrix(part6-part7)
  
  b1 = -Part7
  # b1_expanded = b1[rep(1:nrow(b1), times = id_counts), , drop = FALSE]
  
  part8 = A.i.1%*%t(gamma.hat)%*%inv.Sigma.e.hat%*%gamma.hat
  Part8 = sm%*%part8
  c1 = -Part8%*%rep(1,no.gauss^K)/2
  
  ## use Gauss-Legendre integration, 20 nodes
  
  t_nodes <- outer((U+V)/2, z, FUN = function(m, z) m + (U-V)/2 * z)      # (400 x 20) matrix
  t_weights <- outer((U-V)/2, w, `*`) 
  ## compute \sum \int lambda(s)*exp(\theta_m*M)ds
  rows_to_insert <- long.data %>%
    filter(visit == befdiag)
  new_rows <- rows_to_insert %>%
    mutate(visit := visit + 0.1)
  theta.mm_expanded = theta.m.hat%*%t(y)
  
  #theta.mm at V
  theta.mm = t(A.i.4[A.i.4!=0]*theta.m.hat%*%t(y[A.i.4!=0,]))
  #theta.xx at V
  theta.xx_V = t(A.i.4[A.i.4!=0]*theta.x.hat%*%t(xval[A.i.4!=0,]))
  # \sum \int lambda(s)*exp(\theta_m*M)ds before V
  spline_time_int = spline_cumulative_hazard(time, knots, alpha.hat, boundary_knots, degree)
  # lambda.mm1 = sm%*%(A.i.2*diff(c(spline_time_int,0))*t(exp(theta.mm_expanded)))
  
  # theta.xx and theta.aa
  theta.xx = theta.x.hat %*% t(xval)
  lambda.mx1 = sm%*%(A.i.2*diff(c(spline_time_int,0))*t(exp(theta.mm_expanded+theta.xx)))
  
  theta.aa = theta.a.hat %*% t(a.re)
  # theta.xx_1 = kronecker(t(rep(1,no.gauss^K)),t(theta.xx))
  theta.aa_1 = kronecker(rep(1,n.id), theta.aa)
  # reformat
  # theta.xx_full = theta.xx_1[, rep(1:ncol(theta.xx_1), each = no.gauss)]
  theta.aa_full = theta.aa_1[, rep(1:ncol(theta.aa_1), each = no.gauss)]
  c1_full = c1[, rep(1:ncol(c1), each = no.gauss)]
  b1_full = b1[, rep(1:ncol(b1), each = no.gauss)]
  # int_with_t = 0
  # for (k in 1:batch){
  #   u <- runif(n.id * batch_sample)   # uniform(0,1)
  #   a <- rep(intervals[, 1], each = batch_sample)
  #   b <- rep(intervals[, 2], each = batch_sample)
  #   sample.time <- t(matrix(a + (b - a) * u, nrow = batch_sample, ncol = n.id))
    
    #long.data.expanded <- bind_rows(long.data, new_rows) %>%
    # arrange(id, visit)
    #y_expanded <- long.data.expanded %>% dplyr::select(contains("ynew"))
    time_expanded = matrix(rows_to_insert$visittime)[,rep(1, each = no.gauss)]
    id_time = data.frame(id = rows_to_insert$id, visit = rows_to_insert$visit, time_expanded)
    id_GLtime = data.frame(id = long.data.id$id, visit = new_rows$visit, t_nodes)
    id_time_expanded <- bind_rows(id_time, id_GLtime) %>%
      arrange(id, visit)
    
    
    # \sum \int lambda(s)*exp(\theta_m*M)ds from V to MC.time
    
    GL_time_expanded <- id_time_expanded %>% dplyr::select(contains("X"))
    GL_time.vector = as.vector(as.matrix(GL_time_expanded))
    GL_spline_int = matrix(spline_cumulative_hazard(GL_time.vector, knots, alpha.hat, boundary_knots, degree),nrow = n.id*2, ncol = no.gauss)
    int_matrix <- Matrix(kronecker(diag(n.id), t(c(-1,1))), sparse = T)
    GL_spline <- int_matrix%*%GL_spline_int
    lambda.mx2 = GL_spline * exp(theta.mm+theta.xx_V)
    
    #total lambda * exp(M)
    lambda.mx <- lambda.mx1 + status*lambda.mx2
    # lambda.mm[lambda.mm<0] = 1 # censor case
    
    # reformat
    # theta.xx_full = theta.xx_1[, rep(1:ncol(theta.xx_1), each = ncol(lambda.mm))]
    theta.aa_full = theta.aa_1[, rep(1:ncol(theta.aa_1), each = ncol(lambda.mx))]
    GL_time_full = t_nodes[, rep(1:no.gauss, times = no.gauss^K)]
    lambda.mx_full = lambda.mx[, rep(1:ncol(lambda.mx), times = no.gauss^K)]
    factor1_full = exp(-lambda.mx_full*exp(theta.aa_full))
    lambda_t = spline_hazard_matrix(t_nodes, knots, alpha.hat, boundary_knots, degree)
    lambda_t_full = lambda_t[, rep(1:ncol(lambda_t), times = no.gauss^K)]
    # c1_full = c1[, rep(1:ncol(c1), each = ncol(sample.time))]
    # b1_full = b1[, rep(1:ncol(b1), each = ncol(sample.time))]
    # calculate MC integrate
    GL_integrate = exp(b1_full*GL_time_full+c1_full*GL_time_full^2)*lambda_t_full*factor1_full
    weighting_matrix <- kronecker(
      t(rep(1,no.gauss^K)),t_weights)
    sum_matrix = Matrix(kronecker(diag(no.gauss^K), rep(1, no.gauss)), sparse = T)
    int_GL = (GL_integrate*weighting_matrix) %*% sum_matrix 
  #   gc()
  # }
  
  # calculate pI
  theta.mm1 = kronecker(t(rep(1,no.gauss^K)),theta.mm)
  theta.xx1 = kronecker(t(rep(1,no.gauss^K)),theta.xx_V)
  # pI = int_GL * exp(a1) * exp(theta.xx1+theta.aa_1+theta.mm1)
  logpI = log(int_GL) + a1 + theta.xx1 + theta.aa_1 + theta.mm1
  # logpI[is.infinite(logpI)] = -700
  # q1
  if (K>1){
    part9 = as.matrix((y-rep(1,nrow(y))%*%t(beta.0.hat)-time%*%t(beta.1.hat)-t(beta.2.hat%*%t(xval)))*A.i.3)
  } else {
    part9 = as.matrix((y-as.numeric(beta.0.hat)-time%*%t(beta.1.hat)-t(beta.2.hat%*%t(xval)))*A.i.3)
  }
  part9_1 = rowSums(part9%*%inv.Sigma.e.hat*part9)
  
  part10 = t(rowSums(a.re%*%inv.Sigma.e.hat*a.re)%*%t(A.i.3))
  
  part11 = part9%*%inv.Sigma.e.hat%*%t(a.re)
  
  part12 = (K*log(2*pi)+log(det(Sigma.e.hat)))*A.i.3
  
  Part13 = sm%*%(part12+part9_1+part10-2*part11)
  
  logq1 = -Part13/2
  # q2
  lambda.mx3 = sm%*%t(A.i.2*diff(c(spline_time_int,0))*exp(theta.mm_expanded+theta.xx))
  lambda.mx3_full = kronecker(t(rep(1,no.gauss^K)),lambda.mx3)
  logq2 = -lambda.mx3_full* exp(theta.aa_1)
  # f and Q
  # f = exp(status*(logpI+logq1)+(1-status)*(logq1+logq2))
  # f[is.na(f)] = 0
  # 1. Initialize log_f as a standard matrix (as before)
  log_f <- matrix(NA, nrow = nrow(logpI), ncol = ncol(logpI))

  # 2. Ensure status is logical index
  is_observed <- as.logical(status)

  # 3. Calculate Observed Rows (status == 1)
  if(any(is_observed)) {
    # FIX: Use drop=FALSE to preserve dimensions and as.matrix() to convert type
    # We perform the Matrix addition, then convert to base matrix for assignment
    val_observed <- logpI[is_observed, , drop=FALSE] + logq1[is_observed, , drop=FALSE]
    log_f[is_observed, ] <- as.matrix(val_observed)
  }

  # 4. Calculate Censored Rows (status == 0)
  if(any(!is_observed)) {
    # FIX: Same here
    val_censored <- logq1[!is_observed, , drop=FALSE] + logq2[!is_observed, , drop=FALSE]
    log_f[!is_observed, ] <- as.matrix(val_censored)
  }

  # 5. Exponentiate
  f <- exp(log_f)
  Q = as.matrix(t(f%*%weight))/(sqrt(pi)^K)
  #logQ = log(Q)
  #logQ[is.na(logQ)] = -2000
  return(list(ll=sum(log(Q)), ll.ind=log(Q)))
}
# use alpha and M-spline to calculate hazard
likelihood.spline3 = function(long.data, parameters, knots){
  long.data$id = as.factor(long.data$id)
  id.matrix <- model.matrix(~ id - 1, data = long.data)
  sm = Matrix(t(id.matrix),sparse = T)
  n.id = length(unique(long.data$id))
  long.data.id = long.data[!duplicated(long.data$id),]
  beta.0.hat =  t(as.matrix(parameters %>% dplyr::select(contains("beta.0"))))
  K = length(beta.0.hat)
  beta.1.hat = t(as.matrix(parameters %>% dplyr::select(contains("beta.1"))))
  beta.2.vec = data.frame(parameters %>% dplyr::select(contains("beta.2")))
  col_names <- colnames(beta.2.vec)
  if (length(col_names)==1){
    beta.2.hat=as.matrix(beta.2.vec)
    p = ncol(beta.2.hat)
  } else{
    indices <- do.call(rbind, strsplit(col_names, "_"))[, 2:3]
    indices <- apply(indices, 2, as.integer)
    beta.2.hat = matrix(0, nrow=max(indices[,1]),ncol=max(indices[,2]))
    p=ncol(beta.2.hat)
    for (i in seq_along(col_names)){
      row_idx <- indices[i, 1]
      col_idx <- indices[i, 2]
      beta.2.hat[row_idx, col_idx] <- beta.2.vec[[col_names[i]]]
    }
  }
  gamma.hat = t(as.matrix(parameters %>% dplyr::select(contains("gamma"))))
  
  # lambda0.hat = as.numeric(parameters %>% dplyr::select(contains("lambda0")))
  # lambda.hat = as.matrix(parameters %>% dplyr::select(contains("lambda")))
  alpha.hat = t(as.matrix(parameters %>% dplyr::select(contains("alpha"))))
  theta.x.hat = as.matrix(parameters %>% dplyr::select(contains("theta.x")))
  theta.a.hat = as.matrix(parameters %>% dplyr::select(contains("theta.a")))
  theta.m.hat = as.matrix(parameters %>% dplyr::select(contains("theta.m")))
  
  Sigma.a.hat = matrix(0, nrow=K,ncol=K)
  upper_vec_a = as.matrix(parameters %>% dplyr::select(contains("Sigma.a")))
  Sigma.a.hat[upper.tri(Sigma.a.hat,diag = TRUE)] = upper_vec_a
  if (K>1){
    Sigma.a.hat=Sigma.a.hat+t(Sigma.a.hat)-diag(diag(Sigma.a.hat))
  }
  Sigma.e.hat = matrix(0, nrow=K,ncol=K)
  upper_vec_e = as.matrix(parameters %>% dplyr::select(contains("Sigma.e")))
  Sigma.e.hat[upper.tri(Sigma.e.hat,diag = TRUE)] = upper_vec_e
  if (K>1){
    Sigma.e.hat=Sigma.e.hat+t(Sigma.e.hat)-diag(diag(Sigma.e.hat))
  }
  inv.Sigma.e.hat = solve(Sigma.e.hat)
  
  #L = chol(Sigma.a.hat)
  #a.re = sqrt(2)*nodes%*%t(L)
  a.re = gh_aw$nodes
  weight_a = gh_aw$weights
  
  y = long.data %>% dplyr::select(contains("ynew"))
  xval = long.data %>% dplyr::select(contains("xval"))
  # xval.id = long.data.id %>% dplyr::select(contains("xval"))
  time = long.data$visittime
  status = long.data.id$status
  V = long.data.id$preetime
  U = long.data.id$etime
  
  A.i.1 = as.numeric(long.data$visit>long.data$befdiag)
  A.i.2 = as.numeric(long.data$visit<long.data$befdiag)
  A.i.3 = as.numeric(long.data$visit<=long.data$befdiag)
  A.i.4 = as.numeric(long.data$visit==long.data$befdiag)
  
  # boundary_knots <- c(0, max(time))
  ##calculating a1, b1, c1 and expand them
  if (K>1){
    part1 = as.matrix((y-rep(1,nrow(y))%*%t(beta.0.hat)-time%*%t(beta.1.hat)-t(beta.2.hat%*%t(xval))-time%*%t(gamma.hat))*A.i.1)
  } else {
    part1 = as.matrix((y-as.numeric(beta.0.hat)-time%*%t(beta.1.hat)-t(beta.2.hat%*%t(xval))-time%*%t(gamma.hat))*A.i.1)
  }
  part1_1 = rowSums(part1%*%inv.Sigma.e.hat*part1)
  
  part2 = t(rowSums(a.re%*%inv.Sigma.e.hat*a.re)%*%t(A.i.1))
  
  part3 = part1%*%inv.Sigma.e.hat%*%t(a.re)
  
  part4 = (K*log(2*pi)+log(det(Sigma.e.hat)))*A.i.1
  
  Part4 = sm%*%(part4+part1_1+part2-2*part3)
  a1 = -Part4/2
  # a1_expanded = a1[rep(1:nrow(a1), times = id_counts), , drop = FALSE]
  
  part6 = part1%*%inv.Sigma.e.hat%*%gamma.hat
  
  part7 = data.frame(A.i.1%*%t(a.re%*%inv.Sigma.e.hat%*%gamma.hat))
  Part7 = sm%*%as.matrix(part6-part7)
  
  b1 = -Part7
  # b1_expanded = b1[rep(1:nrow(b1), times = id_counts), , drop = FALSE]
  
  part8 = A.i.1%*%t(gamma.hat)%*%inv.Sigma.e.hat%*%gamma.hat
  Part8 = sm%*%part8
  c1 = -Part8%*%rep(1,no.gauss^K)/2
  
  ## use Gauss-Legendre integration, 20 nodes
  
  t_nodes <- outer((U+V)/2, z, FUN = function(m, z) m + (U-V)/2 * z)      # (400 x 20) matrix
  t_weights <- outer((U-V)/2, w, `*`) 
  ## compute \sum \int lambda(s)*exp(\theta_m*M)ds
  rows_to_insert <- long.data %>%
    filter(visit == befdiag)
  new_rows <- rows_to_insert %>%
    mutate(visit := visit + 0.1)
  theta.mm_expanded = theta.m.hat%*%t(y)
  
  #theta.mm at V
  theta.mm = t(A.i.4[A.i.4!=0]*theta.m.hat%*%t(y[A.i.4!=0,]))
  #theta.xx at V
  theta.xx_V = t(A.i.4[A.i.4!=0]*theta.x.hat%*%t(xval[A.i.4!=0,]))
  # \sum \int lambda(s)*exp(\theta_m*M)ds before V
  spline_time_int = evaluate_Ispline(time, knots, alpha.hat, boundary_knots, degree)
  # lambda.mm1 = sm%*%(A.i.2*diff(c(spline_time_int,0))*t(exp(theta.mm_expanded)))
  
  # theta.xx and theta.aa
  theta.xx = theta.x.hat %*% t(xval)
  lambda.mx1 = sm%*%(A.i.2*diff(c(spline_time_int,0))*t(exp(theta.mm_expanded+theta.xx)))
  
  theta.aa = theta.a.hat %*% t(a.re)
  # theta.xx_1 = kronecker(t(rep(1,no.gauss^K)),t(theta.xx))
  theta.aa_1 = kronecker(rep(1,n.id), theta.aa)
  # reformat
  # theta.xx_full = theta.xx_1[, rep(1:ncol(theta.xx_1), each = no.gauss)]
  theta.aa_full = theta.aa_1[, rep(1:ncol(theta.aa_1), each = no.gauss)]
  c1_full = c1[, rep(1:ncol(c1), each = no.gauss)]
  b1_full = b1[, rep(1:ncol(b1), each = no.gauss)]
  # int_with_t = 0
  # for (k in 1:batch){
  #   u <- runif(n.id * batch_sample)   # uniform(0,1)
  #   a <- rep(intervals[, 1], each = batch_sample)
  #   b <- rep(intervals[, 2], each = batch_sample)
  #   sample.time <- t(matrix(a + (b - a) * u, nrow = batch_sample, ncol = n.id))
  
  #long.data.expanded <- bind_rows(long.data, new_rows) %>%
  # arrange(id, visit)
  #y_expanded <- long.data.expanded %>% dplyr::select(contains("ynew"))
  time_expanded = matrix(rows_to_insert$visittime)[,rep(1, each = no.gauss)]
  id_time = data.frame(id = rows_to_insert$id, visit = rows_to_insert$visit, time_expanded)
  id_GLtime = data.frame(id = long.data.id$id, visit = new_rows$visit, t_nodes)
  id_time_expanded <- bind_rows(id_time, id_GLtime) %>%
    arrange(id, visit)
  
  
  # \sum \int lambda(s)*exp(\theta_m*M)ds from V to MC.time
  
  GL_time_expanded <- id_time_expanded %>% dplyr::select(contains("X"))
  GL_time.vector = as.vector(as.matrix(GL_time_expanded))
  GL_spline_int = matrix(evaluate_Ispline(GL_time.vector, knots, alpha.hat, boundary_knots, degree),nrow = n.id*2, ncol = no.gauss)
  int_matrix <- Matrix(kronecker(diag(n.id), t(c(-1,1))), sparse = T)
  GL_spline <- int_matrix%*%GL_spline_int
  lambda.mx2 = GL_spline * exp(theta.mm+theta.xx_V)
  
  #total lambda * exp(M)
  lambda.mx <- lambda.mx1 + status*lambda.mx2
  # lambda.mm[lambda.mm<0] = 1 # censor case
  
  # reformat
  # theta.xx_full = theta.xx_1[, rep(1:ncol(theta.xx_1), each = ncol(lambda.mm))]
  theta.aa_full = theta.aa_1[, rep(1:ncol(theta.aa_1), each = ncol(lambda.mx))]
  GL_time_full = t_nodes[, rep(1:no.gauss, times = no.gauss^K)]
  lambda.mx_full = lambda.mx[, rep(1:ncol(lambda.mx), times = no.gauss^K)]
  factor1_full = exp(-lambda.mx_full*exp(theta.aa_full))
  lambda_t = evaluate_Mspline_matrix(t_nodes, knots, alpha.hat, boundary_knots, degree)
  lambda_t_full = lambda_t[, rep(1:ncol(lambda_t), times = no.gauss^K)]
  # c1_full = c1[, rep(1:ncol(c1), each = ncol(sample.time))]
  # b1_full = b1[, rep(1:ncol(b1), each = ncol(sample.time))]
  # calculate MC integrate
  GL_integrate = exp(b1_full*GL_time_full+c1_full*GL_time_full^2)*lambda_t_full*factor1_full
  weighting_matrix <- kronecker(
    t(rep(1,no.gauss^K)),t_weights)
  sum_matrix = Matrix(kronecker(diag(no.gauss^K), rep(1, no.gauss)), sparse = T)
  int_GL = (GL_integrate*weighting_matrix) %*% sum_matrix 
  #   gc()
  # }
  
  # calculate pI
  theta.mm1 = kronecker(t(rep(1,no.gauss^K)),theta.mm)
  theta.xx1 = kronecker(t(rep(1,no.gauss^K)),theta.xx_V)
  # pI = int_GL * exp(a1) * exp(theta.xx1+theta.aa_1+theta.mm1)
  logpI = log(int_GL) + a1 + theta.xx1 + theta.aa_1 + theta.mm1
  # logpI[is.infinite(logpI)] = -700
  # q1
  if (K>1){
    part9 = as.matrix((y-rep(1,nrow(y))%*%t(beta.0.hat)-time%*%t(beta.1.hat)-t(beta.2.hat%*%t(xval)))*A.i.3)
  } else {
    part9 = as.matrix((y-as.numeric(beta.0.hat)-time%*%t(beta.1.hat)-t(beta.2.hat%*%t(xval)))*A.i.3)
  }
  part9_1 = rowSums(part9%*%inv.Sigma.e.hat*part9)
  
  part10 = t(rowSums(a.re%*%inv.Sigma.e.hat*a.re)%*%t(A.i.3))
  
  part11 = part9%*%inv.Sigma.e.hat%*%t(a.re)
  
  part12 = (K*log(2*pi)+log(det(Sigma.e.hat)))*A.i.3
  
  Part13 = sm%*%(part12+part9_1+part10-2*part11)
  
  logq1 = -Part13/2
  # q2
  lambda.mx3 = sm%*%t(A.i.2*diff(c(spline_time_int,0))*exp(theta.mm_expanded+theta.xx))
  lambda.mx3_full = kronecker(t(rep(1,no.gauss^K)),lambda.mx3)
  logq2 = -lambda.mx3_full* exp(theta.aa_1)
  # f and Q
  # f = exp(status*(logpI+logq1)+(1-status)*(logq1+logq2))
  # f[is.na(f)] = 0
  # 1. Initialize log_f as a standard matrix (as before)
  log_f <- matrix(NA, nrow = nrow(logpI), ncol = ncol(logpI))
  
  # 2. Ensure status is logical index
  is_observed <- as.logical(status)
  
  # 3. Calculate Observed Rows (status == 1)
  if(any(is_observed)) {
    # FIX: Use drop=FALSE to preserve dimensions and as.matrix() to convert type
    # We perform the Matrix addition, then convert to base matrix for assignment
    val_observed <- logpI[is_observed, , drop=FALSE] + logq1[is_observed, , drop=FALSE]
    log_f[is_observed, ] <- as.matrix(val_observed)
  }
  
  # 4. Calculate Censored Rows (status == 0)
  if(any(!is_observed)) {
    # FIX: Same here
    val_censored <- logq1[!is_observed, , drop=FALSE] + logq2[!is_observed, , drop=FALSE]
    log_f[!is_observed, ] <- as.matrix(val_censored)
  }
  
  # 5. Exponentiate
  f <- exp(log_f)
  f[is.na(f)|is.infinite(f)] = 0
  Q = as.matrix(t(f%*%weight_a))/(sqrt(pi)^K)
  #logQ = log(Q)
  #logQ[is.na(logQ)] = -2000
  return(list(ll=sum(log(Q)), ll.ind=log(Q)))
}
############################
perturb_likelihood <- function(index, parameters, long.data, delta=1e-6) {
  # Create perturbed versions of 'parameters'
  para.1 <- parameters
  para.2 <- parameters
  para.1[, index] <- parameters[, index] + delta
  para.2[, index] <- parameters[, index] - delta
  
  # Compute likelihoods
  if (mode == 1){
    likelihood.1 <- likelihood.vec(long.data, para.1)$ll.ind
    likelihood.2 <- likelihood.vec(long.data, para.2)$ll.ind
  }
  else if (mode == 2){
    likelihood.1 <- likelihood.piecewise(long.data, para.1)$ll.ind
    likelihood.2 <- likelihood.piecewise(long.data, para.2)$ll.ind
  }
  else if (mode == 3){
    likelihood.1 <- likelihood.spline3(long.data, para.1, knots)$ll.ind
    likelihood.2 <- likelihood.spline3(long.data, para.2, knots)$ll.ind
    gc()
  }
  # Compute the numerical difference
  return((likelihood.1 - likelihood.2) / (2 * delta))
}




diff.likelihood.vec = function(long.data, parameters){
  
  num_para = length(parameters)
  # parameters.all=rbind(parameters,rep(1,num_para)%*%as.matrix(parameters)+delta*diag(rep(1,num_para)),rep(1,num_para)%*%as.matrix(parameters)-delta*diag(rep(1,num_para)))
  n.id = length(unique(long.data$id))
  diff = t(sapply(1:num_para, perturb_likelihood, parameters=parameters,long.data=long.data))
  H <- diff %*% t(diff)
  H = (H+t(H))/2
  eig <- eigen(H)
  
  eig$values[eig$values < 1e-6] <- 1e-6   # enforce PD
  H_pd <- eig$vectors %*% diag(eig$values) %*% t(eig$vectors)
  
  hessian <- -solve(H_pd)
  score = apply(diff, 1, sum)
  return(list(Score=score, Hessian=hessian))
}

Est_sd = function(long.data, parameters){
  
  num_para = length(parameters)
  # parameters.all=rbind(parameters,rep(1,num_para)%*%as.matrix(parameters)+delta*diag(rep(1,num_para)),rep(1,num_para)%*%as.matrix(parameters)-delta*diag(rep(1,num_para)))
  n.id = length(unique(long.data$id))
  diff = t(sapply(1:num_para, perturb_likelihood, parameters=parameters,long.data=long.data))
  full_info = diff%*%t(diff)
  num_prim = length(parameters  %>% dplyr::select(!contains("hazard")))
  A11 = full_info[1:num_prim, 1:num_prim]
  A21 = full_info[-(1:num_prim), 1:num_prim]
  A22 = full_info[-(1:num_prim), -(1:num_prim)]
  prim_info = solve(A11 - t(A21)%*%ginv(A22)%*%A21)
  return(sqrt(diag(prim_info)))
}
############################
perturb_likelihood2 <- function(index,  parameters, long.data, delta=1e-6) {
  # Create perturbed versions of 'parameters'
  para.1 <- parameters
  para.2 <- parameters
  para.1[, index] <- parameters[, index] + delta
  para.2[, index] <- parameters[, index] - delta
  
  # Compute scores
  score.1 <- diff.likelihood.vec(long.data, para.1)$Score
  score.2 <- diff.likelihood.vec(long.data, para.2)$Score
  
  # Compute the numerical difference
  return((score.1 - score.2) / (2*delta))
}

#####
diff.likelihood.vec2 = function(long.data, parameters){
  
  num_para = length(parameters)
  n.id = length(unique(long.data$id))
  diff = t(sapply(1:num_para, perturb_likelihood2, parameters=parameters,long.data=long.data))
  
  
  if (is.na(det(diff))||qr(diff)$rank<num_para){
    hessian = -solve(diff+1e-6*diag(num_para))
  } else {
    hessian = -solve(diff)
  }
  return(hessian)
}

###use diff.likelihood.vec2(long.data, parameters) to replace diff$Hessian in NR algorithm (see bootstrap modification below).

####Also, although I agree with Dr. Zhang that the sandwich form is theoretically incorrect due to random knots, I don't think there shall be that large numerical difference. So just curious whether your previous sandwich code is accurate.
#####If compute sandwich variance, then the hessian1 from diff.likelihood.vec is -(nJ)^{-1}, the hessian2 from diff.likelihood.vec2 is -(nI)^{-1}. So the asymptotic variance n^{-1}I^{-1}JI^{-1}=(nI)^{-1}(nJ)(nI)^{-1} shall be 

###A=-diff.likelihood.vec2(long.data,parameters)
###B=-diff.likelihood.vec(long.data,parameters)$Hessian
####A%*%solve(B)%*%A





######################################
## check score if normally distributed
check.score = function(long.data, beta.0.hat = as.vector(beta.0),
                       beta.1.hat = as.vector(beta.1),
                       beta.2.hat = as.vector(beta.2),
                       gamma.hat = as.vector(gamma.1),
                       lambda0.hat = as.vector(1/par[1]),
                       theta.x.hat = as.vector(theta.x),
                       theta.a.hat = as.vector(theta.a),
                       theta.m.hat = as.vector(theta.m),
                       sigma.e.hat = as.matrix(Sigma.e),
                       sigma.a.hat = as.matrix(Sigma.a))
{
  error = FALSE
  Q=matrix(0,n.id)
  Qdbeta.0=matrix(0,n.id)
  Qdbeta.1=matrix(0,n.id)
  Qdbeta.2=matrix(0,n.id)
  Qdgamma=matrix(0,n.id)
  Qdlambda0=matrix(0,n.id)
  Qdtheta.x=matrix(0,n.id)
  Qdtheta.a=matrix(0,n.id)
  Qdtheta.m=matrix(0,n.id)
  Qdsigma.e=matrix(0,n.id)
  Qdsigma.a=matrix(0,n.id)
  Hessian = matrix(0, nrow = 10, ncol = 10)
  Score = rep(0, 10)
  ## transform GH nodes
  a.re = nodes%*%sqrtm(2*sigma.a.hat)$B
  for (i in 1:n.id){
    ## for test
    # i = 1 
    data.i = long.data[long.data$id==i,]
    befdiag.i = data.i$befdiag[1]
    afterdiag.i = data.i$afterdiag[1]
    xval.i = data.i$xval[1]
    V.i = data.i$preetime[1]
    U.i = data.i$etime[1]
    status.i = data.i$status[1]
    visittime.i = data.i$visittime
    befdiagtime.i = visittime.i[data.i$occasion<=data.i$befdiag]
    ondiagtime.i = data.i$preetime[1]
    afterdiagtime.i = visittime.i[data.i$occasion>data.i$befdiag]
    visit_num.i = data.i$visits[1]
    y.i = data.i$ynew1_1
    befdiagy.i = y.i[data.i$occasion<data.i$befdiag]
    befondiagy.i = y.i[data.i$occasion<=data.i$befdiag]
    ondiagy.i = y.i[data.i$occasion==data.i$befdiag]
    afterdiagy.i = y.i[data.i$occasion>data.i$befdiag]
    A.i = rep(1, afterdiag.i)
    A.i.1 = rep(1, befdiag.i)
    if (status.i==1){
      a2 = apply(-(as.numeric(log(2*pi)+log(sigma.e.hat))*matrix(1,afterdiag.i,no.gauss)+
                     (matrix(1,afterdiag.i,no.gauss)*as.vector(afterdiagy.i-beta.0.hat-
                                                                 beta.1.hat*afterdiagtime.i-beta.2.hat*xval.i-gamma.hat*afterdiagtime.i)-A.i%*%t(a.re))^2/as.numeric(sigma.e.hat))/2,2,sum)+
        lambda0.hat*ondiagtime.i*exp(theta.x.hat*xval.i+theta.a.hat*t(a.re)+theta.m.hat*ondiagy.i)
      p = exp(log(lambda0.hat)+theta.x.hat*xval.i+theta.a.hat*t(a.re)+theta.m.hat*ondiagy.i-
                lambda0.hat*sum(diff(befdiagtime.i)*exp(theta.m.hat*befdiagy.i))*exp(theta.x.hat*xval.i+theta.a.hat*t(a.re))+a2)
      b2 = apply(-(matrix(1,afterdiag.i,no.gauss)*(afterdiagy.i-beta.0.hat-beta.1.hat*afterdiagtime.i-beta.2.hat*xval.i-
                                                     gamma.hat*afterdiagtime.i)-A.i%*%t(a.re))*gamma.hat/as.numeric(sigma.e.hat),2,sum)-
        lambda0.hat*exp(theta.x.hat*xval.i+theta.a.hat*t(a.re)+theta.m.hat*ondiagy.i)
      c2 = apply(-matrix(1,afterdiag.i,no.gauss)*gamma.hat^2/2/as.numeric(sigma.e.hat),2,sum)
      # I = sqrt(pi)*exp(-b2^2/(4*c2))*(pnorm(sqrt(-2*c2)*U.i-b2/sqrt(-2*c2),0,1)-pnorm(sqrt(-2*c2)*V.i-b2/sqrt(-2*c2),0,1))/sqrt(-c2)
      apx = exp(-log(2*pi)/2-((sqrt(-2*c2)*U.i+sqrt(-2*c2)*V.i-2*b2/sqrt(-2*c2))^2/8)+log(sqrt(-2*c2)*(U.i-V.i)))
      pI = sqrt(pi)*exp(log(p)-b2^2/(4*c2))*(pnorm(sqrt(-2*c2)*U.i-b2/sqrt(-2*c2),0,1)-pnorm(sqrt(-2*c2)*V.i-b2/sqrt(-2*c2),0,1))/sqrt(-c2)
      pI_apx = sqrt(pi)*exp(log(p)-b2^2/(4*c2))*apx/sqrt(-c2)
      pI[,pI[1,]==0]=pI_apx[,pI[1,]==0]
      q = exp(apply(-((log(2*pi)+as.numeric(log(sigma.e.hat)))*matrix(1,befdiag.i,no.gauss)+
                        (matrix(1,befdiag.i,no.gauss)*as.vector(befondiagy.i-beta.0.hat-beta.1.hat*befdiagtime.i-beta.2.hat*xval.i)-A.i.1%*%t(a.re))^2/as.numeric(sigma.e.hat))/2,2,sum))
      # f
      # f = p*I*q*exp(d2*t(a.re)+e2*t(a.re^2))
      f = pI*q
      
      ## score of each component
      ## p
      # pdbeta.0 = apply((matrix(1,afterdiag.i,no.gauss)*(afterdiagy.i-beta.0.hat-beta.1.hat*afterdiagtime.i-beta.2.hat*xval.i-gamma.hat*afterdiagtime.i)-A.i%*%t(a.re))/as.numeric(sigma.e.hat),2,sum)*p
      # pdbeta.1 = apply((matrix(1,afterdiag.i,no.gauss)*(afterdiagy.i-beta.0.hat-beta.1.hat*afterdiagtime.i-beta.2.hat*xval.i-gamma.hat*afterdiagtime.i)*afterdiagtime.i-A.i%*%t(a.re)*afterdiagtime.i)/as.numeric(sigma.e.hat),2,sum)*p
      # pdbeta.2 = apply((matrix(1,afterdiag.i,no.gauss)*(afterdiagy.i-beta.0.hat-beta.1.hat*afterdiagtime.i-beta.2.hat*xval.i-gamma.hat*afterdiagtime.i)*xval.i-A.i%*%t(a.re)*xval.i)/as.numeric(sigma.e.hat),2,sum)*p
      # pdgamma = apply((matrix(1,afterdiag.i,no.gauss)*(afterdiagy.i-beta.0.hat-beta.1.hat*afterdiagtime.i-beta.2.hat*xval.i-gamma.hat*afterdiagtime.i)*afterdiagtime.i-A.i%*%t(a.re)*afterdiagtime.i)/as.numeric(sigma.e.hat),2,sum)*p
      # pdlambda0 = (1/lambda0.hat-sum(diff(befdiagtime.i)*exp(theta.m.hat*befdiagy.i))*exp(theta.x.hat*xval.i+theta.a.hat*t(a.re))+ondiagtime.i*exp(theta.x.hat*xval.i+theta.a.hat*t(a.re)+theta.m.hat*ondiagy.i))*p
      # pdtheta.x = (xval.i-lambda0.hat*xval.i*sum(diff(befdiagtime.i)*exp(theta.m.hat*befdiagy.i))*exp(theta.x.hat*xval.i+theta.a.hat*t(a.re))+lambda0.hat*ondiagtime.i*xval.i*exp(theta.x.hat*xval.i+theta.a.hat*t(a.re)+theta.m.hat*ondiagy.i))*p
      # pdtheta.a = (t(a.re)-lambda0.hat*sum(diff(befdiagtime.i)*befdiagy.i*exp(theta.m.hat*befdiagy.i))*exp(theta.x.hat*xval.i+theta.a.hat*t(a.re))+lambda0.hat*ondiagtime.i*t(a.re)*exp(theta.x.hat*xval.i+theta.a.hat*t(a.re)+theta.m.hat*ondiagy.i))*p
      # pdtheta.m = (ondiagy.i-lambda0.hat*xval.i*sum(diff(befdiagtime.i)*exp(theta.m.hat*befdiagy.i))*exp(theta.x.hat*xval.i+theta.a.hat*t(a.re))+lambda0.hat*ondiagtime.i*ondiagy.i*exp(theta.x.hat*xval.i+theta.a.hat*t(a.re)+theta.m.hat*ondiagy.i))*p
      # pdsigma.e = -apply(matrix(1,afterdiag.i,no.gauss)/as.numeric(sigma.e.hat)-(matrix(1,afterdiag.i,no.gauss)*(afterdiagy.i-beta.0.hat-beta.1.hat*afterdiagtime.i-beta.2.hat*xval.i-gamma.hat*afterdiagtime.i)-A.i%*%t(a.re))^2/(as.numeric(sigma.e.hat))^2,2,sum)*p/2
      pdbeta.0I = apply((matrix(1,afterdiag.i,no.gauss)*(afterdiagy.i-beta.0.hat-beta.1.hat*afterdiagtime.i-beta.2.hat*xval.i-gamma.hat*afterdiagtime.i)-A.i%*%t(a.re))/as.numeric(sigma.e.hat),2,sum)*pI
      
      pdbeta.1I = apply((matrix(1,afterdiag.i,no.gauss)*(afterdiagy.i-beta.0.hat-beta.1.hat*afterdiagtime.i-beta.2.hat*xval.i-gamma.hat*afterdiagtime.i)*afterdiagtime.i-A.i%*%t(a.re)*afterdiagtime.i)/as.numeric(sigma.e.hat),2,sum)*pI
      
      pdbeta.2I = apply((matrix(1,afterdiag.i,no.gauss)*(afterdiagy.i-beta.0.hat-beta.1.hat*afterdiagtime.i-beta.2.hat*xval.i-gamma.hat*afterdiagtime.i)*xval.i-A.i%*%t(a.re)*xval.i)/as.numeric(sigma.e.hat),2,sum)*pI
      
      pdgammaI = apply((matrix(1,afterdiag.i,no.gauss)*(afterdiagy.i-beta.0.hat-beta.1.hat*afterdiagtime.i-beta.2.hat*xval.i-gamma.hat*afterdiagtime.i)*afterdiagtime.i-A.i%*%t(a.re)*afterdiagtime.i)/as.numeric(sigma.e.hat),2,sum)*pI
      
      pdlambda0I = (1/lambda0.hat-sum(diff(befdiagtime.i)*exp(theta.m.hat*befdiagy.i))*exp(theta.x.hat*xval.i+theta.a.hat*t(a.re))+ondiagtime.i*exp(theta.x.hat*xval.i+theta.a.hat*t(a.re)+theta.m.hat*ondiagy.i))*pI
      
      pdtheta.xI = (xval.i-lambda0.hat*xval.i*sum(diff(befdiagtime.i)*exp(theta.m.hat*befdiagy.i))*exp(theta.x.hat*xval.i+theta.a.hat*t(a.re))+lambda0.hat*ondiagtime.i*xval.i*exp(theta.x.hat*xval.i+theta.a.hat*t(a.re)+theta.m.hat*ondiagy.i))*pI
      
      pdtheta.aI = (t(a.re)-lambda0.hat*sum(diff(befdiagtime.i)*befdiagy.i*exp(theta.m.hat*befdiagy.i))*exp(theta.x.hat*xval.i+theta.a.hat*t(a.re))+lambda0.hat*ondiagtime.i*t(a.re)*exp(theta.x.hat*xval.i+theta.a.hat*t(a.re)+theta.m.hat*ondiagy.i))*pI
      
      pdtheta.mI = (ondiagy.i-lambda0.hat*xval.i*sum(diff(befdiagtime.i)*exp(theta.m.hat*befdiagy.i))*exp(theta.x.hat*xval.i+theta.a.hat*t(a.re))+lambda0.hat*ondiagtime.i*ondiagy.i*exp(theta.x.hat*xval.i+theta.a.hat*t(a.re)+theta.m.hat*ondiagy.i))*pI
      
      pdsigma.eI = -apply(matrix(1,afterdiag.i,no.gauss)/as.numeric(sigma.e.hat)-(matrix(1,afterdiag.i,no.gauss)*(afterdiagy.i-beta.0.hat-beta.1.hat*afterdiagtime.i-beta.2.hat*xval.i-gamma.hat*afterdiagtime.i)-A.i%*%t(a.re))^2/(as.numeric(sigma.e.hat))^2,2,sum)*pI/2
      
      ## b2
      b2dbeta.0 = apply(matrix(1,afterdiag.i,no.gauss)*gamma.hat/as.numeric(sigma.e.hat),2,sum)
      b2dbeta.1 = apply(matrix(1,afterdiag.i,no.gauss)*gamma.hat*afterdiagtime.i/as.numeric(sigma.e.hat),2,sum)
      b2dbeta.2 = apply(matrix(1,afterdiag.i,no.gauss)*gamma.hat*xval.i/as.numeric(sigma.e.hat),2,sum)
      b2dgamma = apply((matrix(1,afterdiag.i,no.gauss)*(-(afterdiagy.i-beta.0.hat-beta.1.hat*afterdiagtime.i-beta.2.hat*xval.i-2*gamma.hat*afterdiagtime.i))+
                          A.i%*%t(a.re))/as.numeric(sigma.e.hat),2,sum)
      b2dlambda0 = -exp(theta.x.hat*xval.i+theta.a.hat*t(a.re)+theta.m.hat*ondiagy.i)
      b2dtheta.x = -lambda0.hat*xval.i*exp(theta.x.hat*xval.i+theta.a.hat*t(a.re)+theta.m.hat*ondiagy.i)
      b2dtheta.a = -lambda0.hat*t(a.re)*exp(theta.x.hat*xval.i+theta.a.hat*t(a.re)+theta.m.hat*ondiagy.i)
      b2dtheta.m = -lambda0.hat*ondiagy.i*exp(theta.x.hat*xval.i+theta.a.hat*t(a.re)+theta.m.hat*ondiagy.i)
      b2dsigma.e = apply((matrix(1,afterdiag.i,no.gauss)*(afterdiagy.i-beta.0.hat-beta.1.hat*afterdiagtime.i-beta.2.hat*xval.i-gamma.hat*afterdiagtime.i)*gamma.hat
                          -A.i%*%t(a.re)*gamma.hat)/as.numeric(sigma.e.hat)^2,2,sum)
      ## c2
      c2dgamma = apply(-(matrix(1,afterdiag.i,no.gauss)*gamma.hat/as.numeric(sigma.e.hat)),2,sum)
      c2dsigma.e = apply((matrix(1,afterdiag.i,no.gauss)*gamma.hat^2/as.numeric(sigma.e.hat)^2),2,sum)/2
      ## I
      # Idbeta.0 = -2*b2/c2*b2dbeta.0*I+sqrt(pi)*exp(-b2^2/(4*c2))*(-b2dbeta.0/sqrt(-2*c2))*(dnorm(sqrt(-2*c2)*U.i-b2/sqrt(-2*c2))-dnorm(sqrt(-2*c2)*V.i-b2/sqrt(-2*c2)))/sqrt(-c2)
      # Idbeta.1 = -2*b2/c2*b2dbeta.1*I+sqrt(pi)*exp(-b2^2/(4*c2))*(-b2dbeta.1/sqrt(-2*c2))*(dnorm(sqrt(-2*c2)*U.i-b2/sqrt(-2*c2))-dnorm(sqrt(-2*c2)*V.i-b2/sqrt(-2*c2)))/sqrt(-c2)
      # Idbeta.2 = -2*b2/c2*b2dbeta.2*I+sqrt(pi)*exp(-b2^2/(4*c2))*(-b2dbeta.2/sqrt(-2*c2))*(dnorm(sqrt(-2*c2)*U.i-b2/sqrt(-2*c2))-dnorm(sqrt(-2*c2)*V.i-b2/sqrt(-2*c2)))/sqrt(-c2)
      #  Idgamma = ((sqrt(pi)/sqrt(-c2))*(-(2*b2*c2*b2dgamma-b2^2*c2dgamma)/(16*c2^2))*exp(-b2^2/(4*c2))+sqrt(pi)*(-1/(2*sqrt(-c2)^3))*(-c2dgamma)*exp(-b2^2/(4*c2)))*
      #.        (pnorm(sqrt(-2*c2)*U.i-b2/sqrt(-2*c2))-pnorm(sqrt(-2*c2)*V.i-b2/sqrt(-2*c2)))+
      #          sqrt(pi)*exp(-b2^2/(4*c2))*((-2*c2dgamma*U.i/(2*sqrt(-2*c2))-(b2dgamma/sqrt(-2*c2)-(-2*c2dgamma/(2*sqrt(-2*c2)^3))))*dnorm(sqrt(-2*c2)*U.i-b2/sqrt(-2*c2))-
      #                  (-2*c2dgamma*V.i/(2*sqrt(-2*c2))-(b2dgamma/sqrt(-2*c2)-(-2*c2dgamma/(2*sqrt(-2*c2)^3))))*dnorm(sqrt(-2*c2)*V.i-b2/sqrt(-2*c2)))/sqrt(-c2)
      # Idlambda0 = -2*b2/c2*b2dlambda0*I+sqrt(pi)*exp(-b2^2/(4*c2))*(-b2dlambda0/sqrt(-2*c2))*(dnorm(sqrt(-2*c2)*U.i-b2/sqrt(-2*c2))-dnorm(sqrt(-2*c2)*V.i-b2/sqrt(-2*c2)))/sqrt(-c2)
      # Idtheta.x = -2*b2/c2*b2dtheta.x*I+sqrt(pi)*exp(-b2^2/(4*c2))*(-b2dtheta.x/sqrt(-2*c2))*(dnorm(sqrt(-2*c2)*U.i-b2/sqrt(-2*c2))-dnorm(sqrt(-2*c2)*V.i-b2/sqrt(-2*c2)))/sqrt(-c2)
      # Idtheta.a = -2*b2/c2*b2dtheta.a*I+sqrt(pi)*exp(-b2^2/(4*c2))*(-b2dtheta.a/sqrt(-2*c2))*(dnorm(sqrt(-2*c2)*U.i-b2/sqrt(-2*c2))-dnorm(sqrt(-2*c2)*V.i-b2/sqrt(-2*c2)))/sqrt(-c2)
      # Idtheta.m = -2*b2/c2*b2dtheta.m*I+sqrt(pi)*exp(-b2^2/(4*c2))*(-b2dtheta.m/sqrt(-2*c2))*(dnorm(sqrt(-2*c2)*U.i-b2/sqrt(-2*c2))-dnorm(sqrt(-2*c2)*V.i-b2/sqrt(-2*c2)))/sqrt(-c2)
      #Idsigma.e = ((sqrt(pi)/sqrt(-c2))*(-(2*b2*c2*b2dsigma.e-b2^2*c2dsigma.e)/(16*c2^2))*exp(-b2^2/(4*c2))+sqrt(pi)*(-1/(2*sqrt(-c2)^3))*(-c2dsigma.e)*exp(-b2^2/(4*c2)))*
      # (pnorm(sqrt(-2*c2)*U.i-b2/sqrt(-2*c2))-pnorm(sqrt(-2*c2)*V.i-b2/sqrt(-2*c2)))+
      #  sqrt(pi)*exp(-b2^2/(4*c2))*((-2*c2dsigma.e*U.i/(2*sqrt(-2*c2))-(b2dsigma.e/sqrt(-2*c2)-(-2*c2dsigma.e/(2*sqrt(-2*c2)^3))))*dnorm(sqrt(-2*c2)*U.i-b2/sqrt(-2*c2))-
      #                               (-2*c2dsigma.e*V.i/(2*sqrt(-2*c2))-(b2dsigma.e/sqrt(-2*c2)-(-2*c2dsigma.e/(2*sqrt(-2*c2)^3))))*dnorm(sqrt(-2*c2)*V.i-b2/sqrt(-2*c2)))/sqrt(-c2)
      pIdbeta.0 = -2*b2/c2*b2dbeta.0*pI+sqrt(pi)*exp(log(p)-b2^2/(4*c2))*(-b2dbeta.0/sqrt(-2*c2))*(dnorm(sqrt(-2*c2)*U.i-b2/sqrt(-2*c2))-dnorm(sqrt(-2*c2)*V.i-b2/sqrt(-2*c2)))/sqrt(-c2)
      pIdbeta.1 = -2*b2/c2*b2dbeta.1*pI+sqrt(pi)*exp(log(p)-b2^2/(4*c2))*(-b2dbeta.1/sqrt(-2*c2))*(dnorm(sqrt(-2*c2)*U.i-b2/sqrt(-2*c2))-dnorm(sqrt(-2*c2)*V.i-b2/sqrt(-2*c2)))/sqrt(-c2)
      pIdbeta.2 = -2*b2/c2*b2dbeta.2*pI+sqrt(pi)*exp(log(p)-b2^2/(4*c2))*(-b2dbeta.2/sqrt(-2*c2))*(dnorm(sqrt(-2*c2)*U.i-b2/sqrt(-2*c2))-dnorm(sqrt(-2*c2)*V.i-b2/sqrt(-2*c2)))/sqrt(-c2)
      
      pIdgamma = (-c2dgamma/c2/2-(2*b2*c2*b2dgamma-b2^2*c2dgamma)/(4*c2^2))*pI+
        ((-c2dgamma*U.i/sqrt(-2*c2)-(sqrt(-2*c2)*b2dgamma+b2*c2dgamma/sqrt(-2*c2))/(-2*c2))*dnorm(sqrt(-2*c2)*U.i-b2/sqrt(-2*c2))-
           (-c2dgamma*V.i/sqrt(-2*c2)-(sqrt(-2*c2)*b2dgamma+b2*c2dgamma/sqrt(-2*c2))/(-2*c2))*dnorm(sqrt(-2*c2)*V.i-b2/sqrt(-2*c2)))*sqrt(pi)*exp(log(p)-b2^2/(4*c2))/sqrt(-c2)
      pIdlambda0 = -2*b2/c2*b2dlambda0*pI+sqrt(pi)*exp(log(p)-b2^2/(4*c2))*(-b2dlambda0/sqrt(-2*c2))*(dnorm(sqrt(-2*c2)*U.i-b2/sqrt(-2*c2))-dnorm(sqrt(-2*c2)*V.i-b2/sqrt(-2*c2)))/sqrt(-c2)
      pIdtheta.x = -2*b2/c2*b2dtheta.x*pI+sqrt(pi)*exp(log(p)-b2^2/(4*c2))*(-b2dtheta.x/sqrt(-2*c2))*(dnorm(sqrt(-2*c2)*U.i-b2/sqrt(-2*c2))-dnorm(sqrt(-2*c2)*V.i-b2/sqrt(-2*c2)))/sqrt(-c2)
      pIdtheta.a = -2*b2/c2*b2dtheta.a*pI+sqrt(pi)*exp(log(p)-b2^2/(4*c2))*(-b2dtheta.a/sqrt(-2*c2))*(dnorm(sqrt(-2*c2)*U.i-b2/sqrt(-2*c2))-dnorm(sqrt(-2*c2)*V.i-b2/sqrt(-2*c2)))/sqrt(-c2)
      
      pIdtheta.m = -2*b2/c2*b2dtheta.m*pI+sqrt(pi)*exp(log(p)-b2^2/(4*c2))*(-b2dtheta.m/sqrt(-2*c2))*(dnorm(sqrt(-2*c2)*U.i-b2/sqrt(-2*c2))-dnorm(sqrt(-2*c2)*V.i-b2/sqrt(-2*c2)))/sqrt(-c2)
      pIdsigma.e = (-c2dsigma.e/c2/2-(2*b2*c2*b2dsigma.e-b2^2*c2dsigma.e)/(4*c2^2))*pI+
        ((-c2dsigma.e*U.i/sqrt(-2*c2)-(sqrt(-2*c2)*b2dsigma.e+b2*c2dsigma.e/sqrt(-2*c2))/(-2*c2))*dnorm(sqrt(-2*c2)*U.i-b2/sqrt(-2*c2))-
           (-c2dsigma.e*V.i/sqrt(-2*c2)-(sqrt(-2*c2)*b2dsigma.e+b2*c2dsigma.e/sqrt(-2*c2))/(-2*c2))*dnorm(sqrt(-2*c2)*V.i-b2/sqrt(-2*c2)))*sqrt(pi)*exp(log(p)-b2^2/(4*c2))/sqrt(-c2)
      ## q
      qdbeta.0 = apply((matrix(1,befdiag.i,no.gauss)*(befondiagy.i-beta.0.hat-beta.1.hat*befdiagtime.i-beta.2.hat*xval.i)-A.i.1%*%t(a.re))/as.numeric(sigma.e.hat),2,sum)*q
      qdbeta.1 = apply((matrix(1,befdiag.i,no.gauss)*(befondiagy.i-beta.0.hat-beta.1.hat*befdiagtime.i-beta.2.hat*xval.i)-A.i.1%*%t(a.re))*befdiagtime.i/as.numeric(sigma.e.hat),2,sum)*q
      qdbeta.2 = apply((matrix(1,befdiag.i,no.gauss)*(befondiagy.i-beta.0.hat-beta.1.hat*befdiagtime.i-beta.2.hat*xval.i)-A.i.1%*%t(a.re))*xval.i/as.numeric(sigma.e.hat),2,sum)*q
      qdsigma.e = -apply(matrix(1,befdiag.i,no.gauss)/as.numeric(sigma.e.hat)-(matrix(1,befdiag.i,no.gauss)*(befondiagy.i-beta.0.hat-beta.1.hat*befdiagtime.i-beta.2.hat*xval.i)-A.i.1%*%t(a.re))^2/(as.numeric(sigma.e.hat))^2,2,sum)*q/2
      
      ## score function
      fdbeta.0 = pdbeta.0I*q+
        pIdbeta.0*q+
        pI*qdbeta.0
      
      fdbeta.1 = pdbeta.1I*q+
        pIdbeta.1*q+
        pI*qdbeta.1
      
      fdbeta.2 = pdbeta.2I*q+
        pIdbeta.2*q+
        pI*qdbeta.2
      
      fdgamma = pdgammaI*q+
        pIdgamma*q
      
      fdlambda0 = pdlambda0I*q+
        pIdlambda0*q
      
      fdtheta.x = pdtheta.xI*q+
        pIdtheta.x*q
      fdtheta.a = pdtheta.aI*q+
        pIdtheta.a*q
      
      fdtheta.m = pdtheta.mI*q+
        pIdtheta.m*q
      
      fdsigma.e = pdsigma.eI*q+
        pIdsigma.e*q+
        pI*qdsigma.e
      
      fdsigma.a = f*(-(1/as.numeric(sigma.a.hat)-t(a.re)^2/as.numeric(sigma.a.hat)^2)/2)
      ## GH integration
      Q[i,] = as.numeric(f%*%weight)
      Qdbeta.0[i,] = as.numeric(fdbeta.0%*%weight)/Q[i,]
      Qdbeta.1[i,] = as.numeric(fdbeta.1%*%weight)/Q[i,]
      Qdbeta.2[i,] = as.numeric(fdbeta.2%*%weight)/Q[i,]
      Qdgamma[i,] = as.numeric(fdgamma%*%weight)/Q[i,]
      Qdlambda0[i,] = as.numeric(fdlambda0%*%weight)/Q[i,]
      Qdtheta.x[i,] = as.numeric(fdtheta.x%*%weight)/Q[i,]
      Qdtheta.a[i,] = as.numeric(fdtheta.a%*%weight)/Q[i,]
      Qdtheta.m[i,] = as.numeric(fdtheta.m%*%weight)/Q[i,]
      Qdsigma.e[i,] = as.numeric(fdsigma.e%*%weight)/Q[i,]
      Qdsigma.a[i,] = as.numeric(fdsigma.a%*%weight)/Q[i,]
    }
    else if (status.i == 0){
      ## components
      q = exp(apply(-((log(2*pi)+as.numeric(log(sigma.e.hat)))*matrix(1,befdiag.i,no.gauss)+
                        (matrix(1,befdiag.i,no.gauss)*as.vector(befondiagy.i-beta.0.hat-beta.1.hat*befdiagtime.i-beta.2.hat*xval.i)-A.i.1%*%t(a.re))^2/as.numeric(sigma.e.hat))/2,2,sum)-
                lambda0.hat*sum(diff(befdiagtime.i)*exp(theta.m.hat*befdiagy.i))*exp(theta.x.hat*xval.i+theta.a.hat*t(a.re)))
      
      ## score
      # q
      qdbeta.0 = apply((matrix(1,befdiag.i,no.gauss)*(befondiagy.i-beta.0.hat-beta.1.hat*befdiagtime.i-beta.2.hat*xval.i)-A.i.1%*%t(a.re))/as.numeric(sigma.e.hat),2,sum)*q
      qdbeta.1 = apply((matrix(1,befdiag.i,no.gauss)*(befondiagy.i-beta.0.hat-beta.1.hat*befdiagtime.i-beta.2.hat*xval.i)-A.i.1%*%t(a.re))*befdiagtime.i/as.numeric(sigma.e.hat),2,sum)*q
      qdbeta.2 = apply((matrix(1,befdiag.i,no.gauss)*(befondiagy.i-beta.0.hat-beta.1.hat*befdiagtime.i-beta.2.hat*xval.i)-A.i.1%*%t(a.re))*xval.i/as.numeric(sigma.e.hat),2,sum)*q
      qdsigma.e = -apply(matrix(1,befdiag.i,no.gauss)/as.numeric(sigma.e.hat)-(matrix(1,befdiag.i,no.gauss)*(befondiagy.i-beta.0.hat-beta.1.hat*befdiagtime.i-beta.2.hat*xval.i)-A.i.1%*%t(a.re))^2/(as.numeric(sigma.e.hat))^2,2,sum)*q/2
      
      
      qdlambda0 = -sum(diff(befdiagtime.i)*exp(theta.m.hat*befdiagy.i))*exp(theta.x.hat*xval.i+theta.a.hat*t(a.re))*q
      qdtheta.x = -lambda0.hat*xval.i*sum(diff(befdiagtime.i)*exp(theta.m.hat*befdiagy.i))*exp(theta.x.hat*xval.i+theta.a.hat*t(a.re))*q
      qdtheta.a = -lambda0.hat*t(a.re)*sum(diff(befdiagtime.i)*exp(theta.m.hat*befdiagy.i))*exp(theta.x.hat*xval.i+theta.a.hat*t(a.re))*q
      qdtheta.m = -lambda0.hat*sum(diff(befdiagtime.i)*exp(theta.m.hat*befdiagy.i)*befdiagy.i)*exp(theta.x.hat*xval.i+theta.a.hat*t(a.re))*q
      # f
      f = q
      ## score function
      fdbeta.0 = qdbeta.0
      fdbeta.1 = qdbeta.1
      fdbeta.2 = qdbeta.2
      fdgamma = matrix(0,no.gauss,nrow=1)
      fdlambda0 = qdlambda0
      fdtheta.x = qdtheta.x
      fdtheta.a = qdtheta.a
      fdtheta.m = qdtheta.m
      fdsigma.e = qdsigma.e
      fdsigma.a = f*(-(1/as.numeric(sigma.a.hat)-t(a.re)^2/as.numeric(sigma.a.hat)^2)/2)
      ## GH
      Q[i,] = as.numeric(f%*%weight)
      Qdbeta.0[i,] = as.numeric(fdbeta.0%*%weight)/Q[i,]
      Qdbeta.1[i,] = as.numeric(fdbeta.1%*%weight)/Q[i,]
      Qdbeta.2[i,] = as.numeric(fdbeta.2%*%weight)/Q[i,]
      Qdgamma[i,] = as.numeric(fdgamma%*%weight)/Q[i,]
      Qdlambda0[i,] = as.numeric(fdlambda0%*%weight)/Q[i,]
      Qdtheta.x[i,] = as.numeric(fdtheta.x%*%weight)/Q[i,]
      Qdtheta.a[i,] = as.numeric(fdtheta.a%*%weight)/Q[i,]
      Qdtheta.m[i,] = as.numeric(fdtheta.m%*%weight)/Q[i,]
      Qdsigma.e[i,] = as.numeric(fdsigma.e%*%weight)/Q[i,]
      Qdsigma.a[i,] = as.numeric(fdsigma.a%*%weight)/Q[i,]
    }
    # Use total hessian and score to do Newton Ralphson
    ## We can also use score and Hessian of each individual to do NR
    score.i = c(Qdbeta.0[i,],Qdbeta.1[i,],Qdbeta.2[i,],Qdgamma[i,],Qdlambda0[i,],Qdtheta.x[i,],Qdtheta.a[i,],Qdtheta.m[i,],Qdsigma.e[i,],Qdsigma.a[i,])
    Hessian.i = -score.i%*%t(score.i)
    if (is.nan(Q[i,])){
      error = TRUE
      break
    }
    
    
  }
  if (error == TRUE){
    break
  }
  hist(Qdbeta.0)
  hist(Qdbeta.1)
  hist(Qdbeta.2)
  hist(Qdgamma)
  hist(Qdlambda0)
  hist(Qdtheta.x)
  hist(Qdtheta.a)
  hist(Qdtheta.m)
  hist(Qdsigma.e)
  hist(Qdsigma.a)
}
# long.data.all = data.gen(lambdavec,B,beta.0,beta.1,beta.2,theta.x,theta.a,theta.m,
#                          gamma.1,K,p1,Sigma.x,Sigma.a,Sigma.e,sigma.errlist,n.id,J,par)
# check.score(long.data.all)
###############################
## check likelihood

# initial.likelihood(long.data.all)
# beta.0.test = seq(0,0.4,by=0.01)
# likelihood.test = rep(0,length(beta.0.test))
# for (i in 1:length(beta.0.test)){
#   likelihood.test[i]=initial.likelihood(long.data.all,ini.beta.0 = beta.0.test[i])
# }
# plot(beta.0.test,likelihood.test)
# 
# beta.1.test = seq(0,0.2,by=0.01)
# likelihood.test = rep(0,length(beta.1.test))
# for (i in 1:length(beta.1.test)){
#   likelihood.test[i]=initial.likelihood(long.data.all,ini.beta.1 = beta.1.test[i])
# }
# plot(beta.1.test,likelihood.test)
# 
# beta.2.test = seq(0,0.2,by=0.01)
# likelihood.test = rep(0,length(beta.2.test))
# for (i in 1:length(beta.2.test)){
#   likelihood.test[i]=initial.likelihood(long.data.all,ini.beta.2 = beta.2.test[i])
# }
# plot(beta.2.test,likelihood.test)
# 
# lambda0.test = seq(0.1,0.3,by=0.01)
# likelihood.test = rep(0,length(lambda0.test))
# for (i in 1:length(lambda0.test)){
#   likelihood.test[i]=initial.likelihood(long.data.all,ini.lambda0  = lambda0.test[i])
# }
# plot(lambda0.test,likelihood.test)
# 
# gamma.test = seq(0.01,0.2,by=0.01)
# likelihood.test = rep(0,length(gamma.test))
# for (i in 1:length(gamma.test)){
#   likelihood.test[i]=initial.likelihood(long.data.all,ini.gamma  = gamma.test[i])
# }
# plot(gamma.test,likelihood.test)
# 
# theta.x.test = seq(0,0.2,by=0.01)
# likelihood.test = rep(0,length(theta.x.test))
# for (i in 1:length(theta.x.test)){
#   likelihood.test[i]=initial.likelihood(long.data.all,ini.theta.x  = theta.x.test[i])
# }
# plot(theta.x.test,likelihood.test)
# 
# theta.a.test = seq(0,0.2,by=0.01)
# likelihood.test = rep(0,length(theta.a.test))
# for (i in 1:length(theta.a.test)){
#   likelihood.test[i]=initial.likelihood(long.data.all,ini.theta.a  = theta.a.test[i])
# }
# plot(theta.a.test,likelihood.test)
# 
# theta.m.test = seq(0,0.2,by=0.01)
# likelihood.test = rep(0,length(theta.m.test))
# for (i in 1:length(theta.m.test)){
#   likelihood.test[i]=initial.likelihood(long.data.all,ini.theta.m  = theta.m.test[i])
# }
# plot(theta.m.test,likelihood.test)
# 
# sigma.e.test = seq(0.8,1,by=0.01)
# likelihood.test = rep(0,length(sigma.e.test))
# for (i in 1:length(sigma.e.test)){
#   likelihood.test[i]=initial.likelihood(long.data.all,ini.sigma.e  = sigma.e.test[i])
# }
# plot(sigma.e.test,likelihood.test)
# 
# sigma.a.test = seq(0.7,0.9,by=0.01)
# likelihood.test = rep(0,length(sigma.a.test))
# for (i in 1:length(sigma.a.test)){
#   likelihood.test[i]=initial.likelihood(long.data.all,ini.sigma.a  = sigma.a.test[i])
# }
# plot(sigma.a.test,likelihood.test)





###############################
# Calculate score and Hessian
score_fun = function(long.data, beta.0.hat = as.vector(beta.0),
                     beta.1.hat = as.vector(beta.1),
                     beta.2.hat = as.vector(beta.2),
                     gamma.hat = as.vector(gamma.1),
                     lambda0.hat = as.vector(1/par[1]),
                     theta.x.hat = as.vector(theta.x),
                     theta.a.hat = as.vector(theta.a),
                     theta.m.hat = as.vector(theta.m),
                     sigma.e.hat = as.matrix(Sigma.e),
                     sigma.a.hat = as.matrix(Sigma.a)){
  error = FALSE
  Q=matrix(0,n.id)
  Qdbeta.0=matrix(0,n.id)
  Qdbeta.1=matrix(0,n.id)
  Qdbeta.2=matrix(0,n.id)
  Qdgamma=matrix(0,n.id)
  Qdlambda0=matrix(0,n.id)
  Qdtheta.x=matrix(0,n.id)
  Qdtheta.a=matrix(0,n.id)
  Qdtheta.m=matrix(0,n.id)
  Qdsigma.e=matrix(0,n.id)
  Qdsigma.a=matrix(0,n.id)
  Hessian = matrix(0, nrow = 10, ncol = 10)
  Score = rep(0, 10)
  ## transform GH nodes
  a.re = nodes%*%sqrtm(2*sigma.a.hat)$B
  for (i in 1:n.id){
    ## for test
    # i = 1 
    data.i = long.data[long.data$id==i,]
    befdiag.i = data.i$befdiag[1]
    afterdiag.i = data.i$afterdiag[1]
    xval.i = data.i$xval[1]
    V.i = data.i$preetime[1]
    U.i = data.i$etime[1]
    status.i = data.i$status[1]
    visittime.i = data.i$visittime
    befdiagtime.i = visittime.i[data.i$occasion<=data.i$befdiag]
    ondiagtime.i = data.i$preetime[1]
    afterdiagtime.i = visittime.i[data.i$occasion>data.i$befdiag]
    visit_num.i = data.i$visits[1]
    y.i = data.i$ynew1_1
    befdiagy.i = y.i[data.i$occasion<data.i$befdiag]
    befondiagy.i = y.i[data.i$occasion<=data.i$befdiag]
    ondiagy.i = y.i[data.i$occasion==data.i$befdiag]
    afterdiagy.i = y.i[data.i$occasion>data.i$befdiag]
    A.i = rep(1, afterdiag.i)
    A.i.1 = rep(1, befdiag.i)
    if (status.i==1){
      a2 = apply(-(as.numeric(log(2*pi)+log(sigma.e.hat))*matrix(1,afterdiag.i,no.gauss)+
                     (matrix(1,afterdiag.i,no.gauss)*as.vector(afterdiagy.i-beta.0.hat-
                                                                 beta.1.hat*afterdiagtime.i-beta.2.hat*xval.i-gamma.hat*afterdiagtime.i)-A.i%*%t(a.re))^2/as.numeric(sigma.e.hat))/2,2,sum)+
        lambda0.hat*ondiagtime.i*exp(theta.x.hat*xval.i+theta.a.hat*t(a.re)+theta.m.hat*ondiagy.i)
      p = exp(log(lambda0.hat)+theta.x.hat*xval.i+theta.a.hat*t(a.re)+theta.m.hat*ondiagy.i-
                lambda0.hat*sum(diff(befdiagtime.i)*exp(theta.m.hat*befdiagy.i))*exp(theta.x.hat*xval.i+theta.a.hat*t(a.re))+a2)
      b2 = apply(-(matrix(1,afterdiag.i,no.gauss)*(afterdiagy.i-beta.0.hat-beta.1.hat*afterdiagtime.i-beta.2.hat*xval.i-
                                                     gamma.hat*afterdiagtime.i)-A.i%*%t(a.re))*gamma.hat/as.numeric(sigma.e.hat),2,sum)-
        lambda0.hat*exp(theta.x.hat*xval.i+theta.a.hat*t(a.re)+theta.m.hat*ondiagy.i)
      c2 = apply(-matrix(1,afterdiag.i,no.gauss)*gamma.hat^2/2/as.numeric(sigma.e.hat),2,sum)
      # I = sqrt(pi)*exp(-b2^2/(4*c2))*(pnorm(sqrt(-2*c2)*U.i-b2/sqrt(-2*c2),0,1)-pnorm(sqrt(-2*c2)*V.i-b2/sqrt(-2*c2),0,1))/sqrt(-c2)
      apx = exp(-log(2*pi)/2-((sqrt(-2*c2)*U.i+sqrt(-2*c2)*V.i-2*b2/sqrt(-2*c2))^2/8)+log(sqrt(-2*c2)*(U.i-V.i)))
      pI = sqrt(pi)*exp(log(p)-b2^2/(4*c2))*(pnorm(sqrt(-2*c2)*U.i-b2/sqrt(-2*c2),0,1)-pnorm(sqrt(-2*c2)*V.i-b2/sqrt(-2*c2),0,1))/sqrt(-c2)
      pI_apx = sqrt(pi)*exp(log(p)-b2^2/(4*c2))*apx/sqrt(-c2)
      pI[,pI[1,]==0]=pI_apx[,pI[1,]==0]
      q = exp(apply(-((log(2*pi)+as.numeric(log(sigma.e.hat)))*matrix(1,befdiag.i,no.gauss)+
                        (matrix(1,befdiag.i,no.gauss)*as.vector(befondiagy.i-beta.0.hat-beta.1.hat*befdiagtime.i-beta.2.hat*xval.i)-A.i.1%*%t(a.re))^2/as.numeric(sigma.e.hat))/2,2,sum))
      # f
      # f = p*I*q*exp(d2*t(a.re)+e2*t(a.re^2))
      f = pI*q
      
      ## score of each component
      ## p
      # pdbeta.0 = apply((matrix(1,afterdiag.i,no.gauss)*(afterdiagy.i-beta.0.hat-beta.1.hat*afterdiagtime.i-beta.2.hat*xval.i-gamma.hat*afterdiagtime.i)-A.i%*%t(a.re))/as.numeric(sigma.e.hat),2,sum)*p
      # pdbeta.1 = apply((matrix(1,afterdiag.i,no.gauss)*(afterdiagy.i-beta.0.hat-beta.1.hat*afterdiagtime.i-beta.2.hat*xval.i-gamma.hat*afterdiagtime.i)*afterdiagtime.i-A.i%*%t(a.re)*afterdiagtime.i)/as.numeric(sigma.e.hat),2,sum)*p
      # pdbeta.2 = apply((matrix(1,afterdiag.i,no.gauss)*(afterdiagy.i-beta.0.hat-beta.1.hat*afterdiagtime.i-beta.2.hat*xval.i-gamma.hat*afterdiagtime.i)*xval.i-A.i%*%t(a.re)*xval.i)/as.numeric(sigma.e.hat),2,sum)*p
      # pdgamma = apply((matrix(1,afterdiag.i,no.gauss)*(afterdiagy.i-beta.0.hat-beta.1.hat*afterdiagtime.i-beta.2.hat*xval.i-gamma.hat*afterdiagtime.i)*afterdiagtime.i-A.i%*%t(a.re)*afterdiagtime.i)/as.numeric(sigma.e.hat),2,sum)*p
      # pdlambda0 = (1/lambda0.hat-sum(diff(befdiagtime.i)*exp(theta.m.hat*befdiagy.i))*exp(theta.x.hat*xval.i+theta.a.hat*t(a.re))+ondiagtime.i*exp(theta.x.hat*xval.i+theta.a.hat*t(a.re)+theta.m.hat*ondiagy.i))*p
      # pdtheta.x = (xval.i-lambda0.hat*xval.i*sum(diff(befdiagtime.i)*exp(theta.m.hat*befdiagy.i))*exp(theta.x.hat*xval.i+theta.a.hat*t(a.re))+lambda0.hat*ondiagtime.i*xval.i*exp(theta.x.hat*xval.i+theta.a.hat*t(a.re)+theta.m.hat*ondiagy.i))*p
      # pdtheta.a = (t(a.re)-lambda0.hat*sum(diff(befdiagtime.i)*befdiagy.i*exp(theta.m.hat*befdiagy.i))*exp(theta.x.hat*xval.i+theta.a.hat*t(a.re))+lambda0.hat*ondiagtime.i*t(a.re)*exp(theta.x.hat*xval.i+theta.a.hat*t(a.re)+theta.m.hat*ondiagy.i))*p
      # pdtheta.m = (ondiagy.i-lambda0.hat*xval.i*sum(diff(befdiagtime.i)*exp(theta.m.hat*befdiagy.i))*exp(theta.x.hat*xval.i+theta.a.hat*t(a.re))+lambda0.hat*ondiagtime.i*ondiagy.i*exp(theta.x.hat*xval.i+theta.a.hat*t(a.re)+theta.m.hat*ondiagy.i))*p
      # pdsigma.e = -apply(matrix(1,afterdiag.i,no.gauss)/as.numeric(sigma.e.hat)-(matrix(1,afterdiag.i,no.gauss)*(afterdiagy.i-beta.0.hat-beta.1.hat*afterdiagtime.i-beta.2.hat*xval.i-gamma.hat*afterdiagtime.i)-A.i%*%t(a.re))^2/(as.numeric(sigma.e.hat))^2,2,sum)*p/2
      pdbeta.0I = apply((matrix(1,afterdiag.i,no.gauss)*(afterdiagy.i-beta.0.hat-beta.1.hat*afterdiagtime.i-beta.2.hat*xval.i-gamma.hat*afterdiagtime.i)-A.i%*%t(a.re))/as.numeric(sigma.e.hat),2,sum)*pI
      
      pdbeta.1I = apply((matrix(1,afterdiag.i,no.gauss)*(afterdiagy.i-beta.0.hat-beta.1.hat*afterdiagtime.i-beta.2.hat*xval.i-gamma.hat*afterdiagtime.i)*afterdiagtime.i-A.i%*%t(a.re)*afterdiagtime.i)/as.numeric(sigma.e.hat),2,sum)*pI
      
      pdbeta.2I = apply((matrix(1,afterdiag.i,no.gauss)*(afterdiagy.i-beta.0.hat-beta.1.hat*afterdiagtime.i-beta.2.hat*xval.i-gamma.hat*afterdiagtime.i)*xval.i-A.i%*%t(a.re)*xval.i)/as.numeric(sigma.e.hat),2,sum)*pI
      
      pdgammaI = apply((matrix(1,afterdiag.i,no.gauss)*(afterdiagy.i-beta.0.hat-beta.1.hat*afterdiagtime.i-beta.2.hat*xval.i-gamma.hat*afterdiagtime.i)*afterdiagtime.i-A.i%*%t(a.re)*afterdiagtime.i)/as.numeric(sigma.e.hat),2,sum)*pI
      
      pdlambda0I = (1/lambda0.hat-sum(diff(befdiagtime.i)*exp(theta.m.hat*befdiagy.i))*exp(theta.x.hat*xval.i+theta.a.hat*t(a.re))+ondiagtime.i*exp(theta.x.hat*xval.i+theta.a.hat*t(a.re)+theta.m.hat*ondiagy.i))*pI
      
      pdtheta.xI = (xval.i-lambda0.hat*xval.i*sum(diff(befdiagtime.i)*exp(theta.m.hat*befdiagy.i))*exp(theta.x.hat*xval.i+theta.a.hat*t(a.re))+lambda0.hat*ondiagtime.i*xval.i*exp(theta.x.hat*xval.i+theta.a.hat*t(a.re)+theta.m.hat*ondiagy.i))*pI
      
      pdtheta.aI = (t(a.re)-lambda0.hat*sum(diff(befdiagtime.i)*befdiagy.i*exp(theta.m.hat*befdiagy.i))*exp(theta.x.hat*xval.i+theta.a.hat*t(a.re))+lambda0.hat*ondiagtime.i*t(a.re)*exp(theta.x.hat*xval.i+theta.a.hat*t(a.re)+theta.m.hat*ondiagy.i))*pI
      
      pdtheta.mI = (ondiagy.i-lambda0.hat*xval.i*sum(diff(befdiagtime.i)*exp(theta.m.hat*befdiagy.i))*exp(theta.x.hat*xval.i+theta.a.hat*t(a.re))+lambda0.hat*ondiagtime.i*ondiagy.i*exp(theta.x.hat*xval.i+theta.a.hat*t(a.re)+theta.m.hat*ondiagy.i))*pI
      
      pdsigma.eI = -apply(matrix(1,afterdiag.i,no.gauss)/as.numeric(sigma.e.hat)-(matrix(1,afterdiag.i,no.gauss)*(afterdiagy.i-beta.0.hat-beta.1.hat*afterdiagtime.i-beta.2.hat*xval.i-gamma.hat*afterdiagtime.i)-A.i%*%t(a.re))^2/(as.numeric(sigma.e.hat))^2,2,sum)*pI/2
      
      ## b2
      b2dbeta.0 = apply(matrix(1,afterdiag.i,no.gauss)*gamma.hat/as.numeric(sigma.e.hat),2,sum)
      b2dbeta.1 = apply(matrix(1,afterdiag.i,no.gauss)*gamma.hat*afterdiagtime.i/as.numeric(sigma.e.hat),2,sum)
      b2dbeta.2 = apply(matrix(1,afterdiag.i,no.gauss)*gamma.hat*xval.i/as.numeric(sigma.e.hat),2,sum)
      b2dgamma = apply((matrix(1,afterdiag.i,no.gauss)*(-(afterdiagy.i-beta.0.hat-beta.1.hat*afterdiagtime.i-beta.2.hat*xval.i-2*gamma.hat*afterdiagtime.i))+
                          A.i%*%t(a.re))/as.numeric(sigma.e.hat),2,sum)
      b2dlambda0 = -exp(theta.x.hat*xval.i+theta.a.hat*t(a.re)+theta.m.hat*ondiagy.i)
      b2dtheta.x = -lambda0.hat*xval.i*exp(theta.x.hat*xval.i+theta.a.hat*t(a.re)+theta.m.hat*ondiagy.i)
      b2dtheta.a = -lambda0.hat*t(a.re)*exp(theta.x.hat*xval.i+theta.a.hat*t(a.re)+theta.m.hat*ondiagy.i)
      b2dtheta.m = -lambda0.hat*ondiagy.i*exp(theta.x.hat*xval.i+theta.a.hat*t(a.re)+theta.m.hat*ondiagy.i)
      b2dsigma.e = apply((matrix(1,afterdiag.i,no.gauss)*(afterdiagy.i-beta.0.hat-beta.1.hat*afterdiagtime.i-beta.2.hat*xval.i-gamma.hat*afterdiagtime.i)*gamma.hat
                          -A.i%*%t(a.re)*gamma.hat)/as.numeric(sigma.e.hat)^2,2,sum)
      ## c2
      c2dgamma = apply(-(matrix(1,afterdiag.i,no.gauss)*gamma.hat/as.numeric(sigma.e.hat)),2,sum)
      c2dsigma.e = apply((matrix(1,afterdiag.i,no.gauss)*gamma.hat^2/as.numeric(sigma.e.hat)^2),2,sum)/2
      ## I
      # Idbeta.0 = -2*b2/c2*b2dbeta.0*I+sqrt(pi)*exp(-b2^2/(4*c2))*(-b2dbeta.0/sqrt(-2*c2))*(dnorm(sqrt(-2*c2)*U.i-b2/sqrt(-2*c2))-dnorm(sqrt(-2*c2)*V.i-b2/sqrt(-2*c2)))/sqrt(-c2)
      # Idbeta.1 = -2*b2/c2*b2dbeta.1*I+sqrt(pi)*exp(-b2^2/(4*c2))*(-b2dbeta.1/sqrt(-2*c2))*(dnorm(sqrt(-2*c2)*U.i-b2/sqrt(-2*c2))-dnorm(sqrt(-2*c2)*V.i-b2/sqrt(-2*c2)))/sqrt(-c2)
      # Idbeta.2 = -2*b2/c2*b2dbeta.2*I+sqrt(pi)*exp(-b2^2/(4*c2))*(-b2dbeta.2/sqrt(-2*c2))*(dnorm(sqrt(-2*c2)*U.i-b2/sqrt(-2*c2))-dnorm(sqrt(-2*c2)*V.i-b2/sqrt(-2*c2)))/sqrt(-c2)
      #  Idgamma = ((sqrt(pi)/sqrt(-c2))*(-(2*b2*c2*b2dgamma-b2^2*c2dgamma)/(16*c2^2))*exp(-b2^2/(4*c2))+sqrt(pi)*(-1/(2*sqrt(-c2)^3))*(-c2dgamma)*exp(-b2^2/(4*c2)))*
      #.        (pnorm(sqrt(-2*c2)*U.i-b2/sqrt(-2*c2))-pnorm(sqrt(-2*c2)*V.i-b2/sqrt(-2*c2)))+
      #          sqrt(pi)*exp(-b2^2/(4*c2))*((-2*c2dgamma*U.i/(2*sqrt(-2*c2))-(b2dgamma/sqrt(-2*c2)-(-2*c2dgamma/(2*sqrt(-2*c2)^3))))*dnorm(sqrt(-2*c2)*U.i-b2/sqrt(-2*c2))-
      #                  (-2*c2dgamma*V.i/(2*sqrt(-2*c2))-(b2dgamma/sqrt(-2*c2)-(-2*c2dgamma/(2*sqrt(-2*c2)^3))))*dnorm(sqrt(-2*c2)*V.i-b2/sqrt(-2*c2)))/sqrt(-c2)
      # Idlambda0 = -2*b2/c2*b2dlambda0*I+sqrt(pi)*exp(-b2^2/(4*c2))*(-b2dlambda0/sqrt(-2*c2))*(dnorm(sqrt(-2*c2)*U.i-b2/sqrt(-2*c2))-dnorm(sqrt(-2*c2)*V.i-b2/sqrt(-2*c2)))/sqrt(-c2)
      # Idtheta.x = -2*b2/c2*b2dtheta.x*I+sqrt(pi)*exp(-b2^2/(4*c2))*(-b2dtheta.x/sqrt(-2*c2))*(dnorm(sqrt(-2*c2)*U.i-b2/sqrt(-2*c2))-dnorm(sqrt(-2*c2)*V.i-b2/sqrt(-2*c2)))/sqrt(-c2)
      # Idtheta.a = -2*b2/c2*b2dtheta.a*I+sqrt(pi)*exp(-b2^2/(4*c2))*(-b2dtheta.a/sqrt(-2*c2))*(dnorm(sqrt(-2*c2)*U.i-b2/sqrt(-2*c2))-dnorm(sqrt(-2*c2)*V.i-b2/sqrt(-2*c2)))/sqrt(-c2)
      # Idtheta.m = -2*b2/c2*b2dtheta.m*I+sqrt(pi)*exp(-b2^2/(4*c2))*(-b2dtheta.m/sqrt(-2*c2))*(dnorm(sqrt(-2*c2)*U.i-b2/sqrt(-2*c2))-dnorm(sqrt(-2*c2)*V.i-b2/sqrt(-2*c2)))/sqrt(-c2)
      #Idsigma.e = ((sqrt(pi)/sqrt(-c2))*(-(2*b2*c2*b2dsigma.e-b2^2*c2dsigma.e)/(16*c2^2))*exp(-b2^2/(4*c2))+sqrt(pi)*(-1/(2*sqrt(-c2)^3))*(-c2dsigma.e)*exp(-b2^2/(4*c2)))*
      # (pnorm(sqrt(-2*c2)*U.i-b2/sqrt(-2*c2))-pnorm(sqrt(-2*c2)*V.i-b2/sqrt(-2*c2)))+
      #  sqrt(pi)*exp(-b2^2/(4*c2))*((-2*c2dsigma.e*U.i/(2*sqrt(-2*c2))-(b2dsigma.e/sqrt(-2*c2)-(-2*c2dsigma.e/(2*sqrt(-2*c2)^3))))*dnorm(sqrt(-2*c2)*U.i-b2/sqrt(-2*c2))-
      #                               (-2*c2dsigma.e*V.i/(2*sqrt(-2*c2))-(b2dsigma.e/sqrt(-2*c2)-(-2*c2dsigma.e/(2*sqrt(-2*c2)^3))))*dnorm(sqrt(-2*c2)*V.i-b2/sqrt(-2*c2)))/sqrt(-c2)
      pIdbeta.0 = -2*b2/c2*b2dbeta.0*pI+sqrt(pi)*exp(log(p)-b2^2/(4*c2))*(-b2dbeta.0/sqrt(-2*c2))*(dnorm(sqrt(-2*c2)*U.i-b2/sqrt(-2*c2))-dnorm(sqrt(-2*c2)*V.i-b2/sqrt(-2*c2)))/sqrt(-c2)
      pIdbeta.1 = -2*b2/c2*b2dbeta.1*pI+sqrt(pi)*exp(log(p)-b2^2/(4*c2))*(-b2dbeta.1/sqrt(-2*c2))*(dnorm(sqrt(-2*c2)*U.i-b2/sqrt(-2*c2))-dnorm(sqrt(-2*c2)*V.i-b2/sqrt(-2*c2)))/sqrt(-c2)
      pIdbeta.2 = -2*b2/c2*b2dbeta.2*pI+sqrt(pi)*exp(log(p)-b2^2/(4*c2))*(-b2dbeta.2/sqrt(-2*c2))*(dnorm(sqrt(-2*c2)*U.i-b2/sqrt(-2*c2))-dnorm(sqrt(-2*c2)*V.i-b2/sqrt(-2*c2)))/sqrt(-c2)
      
      pIdgamma = (-c2dgamma/c2/2-(2*b2*c2*b2dgamma-b2^2*c2dgamma)/(4*c2^2))*pI+
        ((-c2dgamma*U.i/sqrt(-2*c2)-(sqrt(-2*c2)*b2dgamma+b2*c2dgamma/sqrt(-2*c2))/(-2*c2))*dnorm(sqrt(-2*c2)*U.i-b2/sqrt(-2*c2))-
           (-c2dgamma*V.i/sqrt(-2*c2)-(sqrt(-2*c2)*b2dgamma+b2*c2dgamma/sqrt(-2*c2))/(-2*c2))*dnorm(sqrt(-2*c2)*V.i-b2/sqrt(-2*c2)))*sqrt(pi)*exp(log(p)-b2^2/(4*c2))/sqrt(-c2)
      pIdlambda0 = -2*b2/c2*b2dlambda0*pI+sqrt(pi)*exp(log(p)-b2^2/(4*c2))*(-b2dlambda0/sqrt(-2*c2))*(dnorm(sqrt(-2*c2)*U.i-b2/sqrt(-2*c2))-dnorm(sqrt(-2*c2)*V.i-b2/sqrt(-2*c2)))/sqrt(-c2)
      pIdtheta.x = -2*b2/c2*b2dtheta.x*pI+sqrt(pi)*exp(log(p)-b2^2/(4*c2))*(-b2dtheta.x/sqrt(-2*c2))*(dnorm(sqrt(-2*c2)*U.i-b2/sqrt(-2*c2))-dnorm(sqrt(-2*c2)*V.i-b2/sqrt(-2*c2)))/sqrt(-c2)
      pIdtheta.a = -2*b2/c2*b2dtheta.a*pI+sqrt(pi)*exp(log(p)-b2^2/(4*c2))*(-b2dtheta.a/sqrt(-2*c2))*(dnorm(sqrt(-2*c2)*U.i-b2/sqrt(-2*c2))-dnorm(sqrt(-2*c2)*V.i-b2/sqrt(-2*c2)))/sqrt(-c2)
      
      pIdtheta.m = -2*b2/c2*b2dtheta.m*pI+sqrt(pi)*exp(log(p)-b2^2/(4*c2))*(-b2dtheta.m/sqrt(-2*c2))*(dnorm(sqrt(-2*c2)*U.i-b2/sqrt(-2*c2))-dnorm(sqrt(-2*c2)*V.i-b2/sqrt(-2*c2)))/sqrt(-c2)
      pIdsigma.e = (-c2dsigma.e/c2/2-(2*b2*c2*b2dsigma.e-b2^2*c2dsigma.e)/(4*c2^2))*pI+
        ((-c2dsigma.e*U.i/sqrt(-2*c2)-(sqrt(-2*c2)*b2dsigma.e+b2*c2dsigma.e/sqrt(-2*c2))/(-2*c2))*dnorm(sqrt(-2*c2)*U.i-b2/sqrt(-2*c2))-
           (-c2dsigma.e*V.i/sqrt(-2*c2)-(sqrt(-2*c2)*b2dsigma.e+b2*c2dsigma.e/sqrt(-2*c2))/(-2*c2))*dnorm(sqrt(-2*c2)*V.i-b2/sqrt(-2*c2)))*sqrt(pi)*exp(log(p)-b2^2/(4*c2))/sqrt(-c2)
      ## q
      qdbeta.0 = apply((matrix(1,befdiag.i,no.gauss)*(befondiagy.i-beta.0.hat-beta.1.hat*befdiagtime.i-beta.2.hat*xval.i)-A.i.1%*%t(a.re))/as.numeric(sigma.e.hat),2,sum)*q
      qdbeta.1 = apply((matrix(1,befdiag.i,no.gauss)*(befondiagy.i-beta.0.hat-beta.1.hat*befdiagtime.i-beta.2.hat*xval.i)-A.i.1%*%t(a.re))*befdiagtime.i/as.numeric(sigma.e.hat),2,sum)*q
      qdbeta.2 = apply((matrix(1,befdiag.i,no.gauss)*(befondiagy.i-beta.0.hat-beta.1.hat*befdiagtime.i-beta.2.hat*xval.i)-A.i.1%*%t(a.re))*xval.i/as.numeric(sigma.e.hat),2,sum)*q
      qdsigma.e = -apply(matrix(1,befdiag.i,no.gauss)/as.numeric(sigma.e.hat)-(matrix(1,befdiag.i,no.gauss)*(befondiagy.i-beta.0.hat-beta.1.hat*befdiagtime.i-beta.2.hat*xval.i)-A.i.1%*%t(a.re))^2/(as.numeric(sigma.e.hat))^2,2,sum)*q/2
      
      ## score function
      fdbeta.0 = pdbeta.0I*q+
        pIdbeta.0*q+
        pI*qdbeta.0
      
      fdbeta.1 = pdbeta.1I*q+
        pIdbeta.1*q+
        pI*qdbeta.1
      
      fdbeta.2 = pdbeta.2I*q+
        pIdbeta.2*q+
        pI*qdbeta.2
      
      fdgamma = pdgammaI*q+
        pIdgamma*q
      
      fdlambda0 = pdlambda0I*q+
        pIdlambda0*q
      
      fdtheta.x = pdtheta.xI*q+
        pIdtheta.x*q
      
      fdtheta.a = pdtheta.aI*q+
        pIdtheta.a*q
      
      fdtheta.m = pdtheta.mI*q+
        pIdtheta.m*q
      
      fdsigma.e = pdsigma.eI*q+
        pIdsigma.e*q+
        pI*qdsigma.e
      
      fdsigma.a = f*(-(1/as.numeric(sigma.a.hat)-t(a.re)^2/as.numeric(sigma.a.hat)^2)/2)
      ## GH integration
      Q[i,] = as.numeric(f%*%weight)
      Qdbeta.0[i,] = as.numeric(fdbeta.0%*%weight)/Q[i,]
      Qdbeta.1[i,] = as.numeric(fdbeta.1%*%weight)/Q[i,]
      Qdbeta.2[i,] = as.numeric(fdbeta.2%*%weight)/Q[i,]
      Qdgamma[i,] = as.numeric(fdgamma%*%weight)/Q[i,]
      Qdlambda0[i,] = as.numeric(fdlambda0%*%weight)/Q[i,]
      Qdtheta.x[i,] = as.numeric(fdtheta.x%*%weight)/Q[i,]
      Qdtheta.a[i,] = as.numeric(fdtheta.a%*%weight)/Q[i,]
      Qdtheta.m[i,] = as.numeric(fdtheta.m%*%weight)/Q[i,]
      Qdsigma.e[i,] = as.numeric(fdsigma.e%*%weight)/Q[i,]
      Qdsigma.a[i,] = as.numeric(fdsigma.a%*%weight)/Q[i,]
    }
    else if (status.i == 0){
      ## components
      q = exp(apply(-((log(2*pi)+as.numeric(log(sigma.e.hat)))*matrix(1,befdiag.i,no.gauss)+
                        (matrix(1,befdiag.i,no.gauss)*as.vector(befondiagy.i-beta.0.hat-beta.1.hat*befdiagtime.i-beta.2.hat*xval.i)-A.i.1%*%t(a.re))^2/as.numeric(sigma.e.hat))/2,2,sum)-
                lambda0.hat*sum(diff(befdiagtime.i)*exp(theta.m.hat*befdiagy.i))*exp(theta.x.hat*xval.i+theta.a.hat*t(a.re)))
      
      ## score
      # q
      qdbeta.0 = apply((matrix(1,befdiag.i,no.gauss)*(befondiagy.i-beta.0.hat-beta.1.hat*befdiagtime.i-beta.2.hat*xval.i)-A.i.1%*%t(a.re))/as.numeric(sigma.e.hat),2,sum)*q
      qdbeta.1 = apply((matrix(1,befdiag.i,no.gauss)*(befondiagy.i-beta.0.hat-beta.1.hat*befdiagtime.i-beta.2.hat*xval.i)-A.i.1%*%t(a.re))*befdiagtime.i/as.numeric(sigma.e.hat),2,sum)*q
      qdbeta.2 = apply((matrix(1,befdiag.i,no.gauss)*(befondiagy.i-beta.0.hat-beta.1.hat*befdiagtime.i-beta.2.hat*xval.i)-A.i.1%*%t(a.re))*xval.i/as.numeric(sigma.e.hat),2,sum)*q
      qdsigma.e = -apply(matrix(1,befdiag.i,no.gauss)/as.numeric(sigma.e.hat)-(matrix(1,befdiag.i,no.gauss)*(befondiagy.i-beta.0.hat-beta.1.hat*befdiagtime.i-beta.2.hat*xval.i)-A.i.1%*%t(a.re))^2/(as.numeric(sigma.e.hat))^2,2,sum)*q/2
      
      
      qdlambda0 = -sum(diff(befdiagtime.i)*exp(theta.m.hat*befdiagy.i))*exp(theta.x.hat*xval.i+theta.a.hat*t(a.re))*q
      qdtheta.x = -lambda0.hat*xval.i*sum(diff(befdiagtime.i)*exp(theta.m.hat*befdiagy.i))*exp(theta.x.hat*xval.i+theta.a.hat*t(a.re))*q
      qdtheta.a = -lambda0.hat*t(a.re)*sum(diff(befdiagtime.i)*exp(theta.m.hat*befdiagy.i))*exp(theta.x.hat*xval.i+theta.a.hat*t(a.re))*q
      qdtheta.m = -lambda0.hat*sum(diff(befdiagtime.i)*exp(theta.m.hat*befdiagy.i)*befdiagy.i)*exp(theta.x.hat*xval.i+theta.a.hat*t(a.re))*q
      # f
      f = q
      ## score function
      fdbeta.0 = qdbeta.0
      fdbeta.1 = qdbeta.1
      fdbeta.2 = qdbeta.2
      fdgamma = matrix(0,no.gauss,nrow=1)
      fdlambda0 = qdlambda0
      fdtheta.x = qdtheta.x
      fdtheta.a = qdtheta.a
      fdtheta.m = qdtheta.m
      fdsigma.e = qdsigma.e
      fdsigma.a = f*(-(1/as.numeric(sigma.a.hat)-t(a.re)^2/as.numeric(sigma.a.hat)^2)/2)
      ## GH
      Q[i,] = as.numeric(f%*%weight)
      Qdbeta.0[i,] = as.numeric(fdbeta.0%*%weight)/Q[i,]
      Qdbeta.1[i,] = as.numeric(fdbeta.1%*%weight)/Q[i,]
      Qdbeta.2[i,] = as.numeric(fdbeta.2%*%weight)/Q[i,]
      Qdgamma[i,] = as.numeric(fdgamma%*%weight)/Q[i,]
      Qdlambda0[i,] = as.numeric(fdlambda0%*%weight)/Q[i,]
      Qdtheta.x[i,] = as.numeric(fdtheta.x%*%weight)/Q[i,]
      Qdtheta.a[i,] = as.numeric(fdtheta.a%*%weight)/Q[i,]
      Qdtheta.m[i,] = as.numeric(fdtheta.m%*%weight)/Q[i,]
      Qdsigma.e[i,] = as.numeric(fdsigma.e%*%weight)/Q[i,]
      Qdsigma.a[i,] = as.numeric(fdsigma.a%*%weight)/Q[i,]
    }
    # Use total hessian and score to do Newton Ralphson
    score.i = c(Qdbeta.0[i,],Qdbeta.1[i,],Qdbeta.2[i,],Qdgamma[i,],Qdlambda0[i,],Qdtheta.x[i,],Qdtheta.a[i,],Qdtheta.m[i,],Qdsigma.e[i,],Qdsigma.a[i,])
    Hessian.i = -score.i%*%t(score.i)
    Score = Score + score.i
    Hessian = Hessian + Hessian.i
  }
  return(list(Score=Score, Hessian=Hessian))
}

## check for score
# beta.0.test = seq(0,0.4,by=0.01)
# score.test = rep(0,length(beta.0.test))
# for (i in 1:length(beta.0.test)){
#   score.test[i]=score_fun(long.data.all, beta.0.hat = beta.0.test[i])$Score[1]
# }
# plot(beta.0.test,score.test)
# 
# beta.1.test = seq(0,0.2,by=0.01)
# score.test = rep(0,length(beta.1.test))
# for (i in 1:length(beta.1.test)){
#   score.test[i]=score_fun(long.data.all, beta.1.hat = beta.1.test[i])$Score[2]
# }
# plot(beta.1.test,score.test)
# 
# beta.2.test = seq(0,0.2,by=0.01)
# score.test = rep(0,length(beta.2.test))
# for (i in 1:length(beta.2.test)){
#   score.test[i]=score_fun(long.data.all, beta.2.hat = beta.2.test[i])$Score[3]
# }
# plot(beta.2.test,score.test)
# 
# gamma.test = seq(0.02,0.2,by=0.01)
# score.test = rep(0,length(gamma.test))
# for (i in 1:length(gamma.test)){
#   score.test[i]=score_fun(long.data.all, gamma.hat = gamma.test[i])$Score[4]
# }
# plot(gamma.test,score.test)
# 
# lambda0.test = seq(0.1,0.3,by=0.01)
# score.test = rep(0,length(lambda0.test))
# for (i in 1:length(lambda0.test)){
#   score.test[i]=score_fun(long.data.all, lambda0.hat = lambda0.test[i])$Score[5]
# }
# plot(lambda0.test,score.test)

