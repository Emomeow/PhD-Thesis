## Initialize parameter
## sample size
n.id=400
## assume before t0, the baseline hazard=0
t0=0
# visit time centered at=0,1,2,3,...,J-1 (if not censored)
## visit label: t1=0, t2, ...,t_J=J-1
J=11

par=c(5,1.5)
# True values of parameters
# K domains for biomarkers
K=2
## p dims for covarite X
p=1
Sigma.x=0.8*diag(p)
Sigma.a=0.8*diag(K)
Sigma.a.vec=matrix(Sigma.a[upper.tri(Sigma.a,diag = TRUE)])
pos_labels_a <- sapply(1:K, function(i) {
  sapply(i:K, function(j) paste0("Sigma.a_", i, "_", j))
})
pos_labels_a <- unlist(pos_labels_a)
rownames(Sigma.a.vec) <- pos_labels_a
Sigma.e=0.8*diag(K)
Sigma.e.vec=matrix(Sigma.e[upper.tri(Sigma.e,diag = TRUE)])
pos_labels_e <- sapply(1:K, function(i) {
  sapply(i:K, function(j) paste0("Sigma.e_", i, "_", j))
})
pos_labels_e <- unlist(pos_labels_e)
rownames(Sigma.e.vec) <- pos_labels_e

## fix effect parameter in lm
## beta.0 intercept, beta.1 for slope (t), beta.2 for X
beta.0=matrix(rep(0.1,K),nrow=K)
rownames(beta.0)=paste("beta.0","_",1:K,sep="")
beta.1=matrix(rep(0.1,K),nrow=K)
rownames(beta.1)=paste("beta.1","_",1:K,sep="")
beta.2=matrix(0.1,nrow = K,ncol= p) 
index_grid <- expand.grid(row = 1:K, col = 1:p)
beta.2.vec=matrix(as.vector(beta.2))
rownames(beta.2.vec)=paste("beta.2",index_grid$row, index_grid$col,sep='_')

## parameter in hazard exponential
theta.x=matrix(rep(0.1,p),nrow=p)
rownames(theta.x)=paste("theta.x",1:p,sep="_")
## dim p*1 
theta.a=matrix(rep(0.1,K),nrow=K)
rownames(theta.a)=paste("theta.a",1:K,sep="_")
## theta.a=rep(0.2,K)
theta.m=matrix(rep(0.1,K),nrow=K)
rownames(theta.m)=paste("theta.m",1:K,sep="_")

## parameter in the model for updating lm
## gamma.1 for slope (t), gamma.2 for X
gamma=matrix(0.1,nrow=K)
rownames(gamma)=paste("gamma",1:K,sep="_")

no.gauss = 20
out.gauss=gauss.quad(no.gauss, kind="hermit")
nodes=as.matrix(do.call(expand.grid,replicate(K,out.gauss$nodes,simplify = FALSE)))
if (K>1){
  weight=as.vector(do.call(outer,replicate(K,out.gauss$weights,simplify = FALSE)))
} else {
  weight=as.vector(out.gauss$weights)
}

## mode = 1 # constant hazard
# mode = 2 # piecewise constant
## mode = 3 # B-spline
if (mode == 1){

## baseline cumulative hazard parameter
lambda0 = matrix(1/par[1])
rownames(lambda0) = paste("lambda0")

## Combine all parameters into a dataframe
parameters = data.frame(t(beta.0),t(beta.1),t(beta.2.vec),t(gamma),t(lambda0),
                        t(theta.x),t(theta.a),t(theta.m),t(Sigma.a.vec),t(Sigma.e.vec))

}

##########################
if (mode == 2){
  knots = 3
  ## Combine all parameters into a dataframe
  parameters = data.frame(t(beta.0),t(beta.1),t(beta.2.vec),t(gamma),
                          t(theta.x),t(theta.a),t(theta.m),t(Sigma.a.vec),t(Sigma.e.vec))

}

##############################
if (mode == 3){
  num_knots = 6
  degree = 3
  parameters = data.frame(t(beta.0),t(beta.1),t(beta.2.vec),t(gamma),
                          t(theta.x),t(theta.a),t(theta.m),t(Sigma.a.vec),t(Sigma.e.vec))
  leg <- gauss.quad(no.gauss, kind = "legendre")
  z <- leg$nodes    # length = 20
  w <- leg$weights  # length = 20
}
