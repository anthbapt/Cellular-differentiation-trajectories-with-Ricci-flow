######################################################
#####Preliminaries: ##################################
#####x= n x T matrix of expression values, rows index genes, columns index time points
#####int.red= n x n adjacency matrix describing how the proteins corresponding to the genes interact:
#####
#####first compute the degree of each gene and the inverse

apply(int.red,1,sum)->deg
1/deg->D

#####generate a list to index the neighbours of each gene
lapply(1:nrow(int.red), function(i){
  print(i)
  return(which(int.red[i,]>0))
})->neighbs

#####define the Forman Ricci curvature using node weights W_i=1/deg(i); w_ij=1/xi*xj
F_Ricci2<-function(i,j,a){###input a is a weighted matrix where aij=xi*xj is (i,j) \in E, 0 else
  if(int.red[i,j]==0){return(NA)}else{
    (1/a[i,j])*(1/a[i,neighbs[[i]]])[-which(neighbs[[i]]==j)]->q
    (1/a[i,j])*(1/a[neighbs[[j]],j])[-which(neighbs[[j]]==i)]->q2
    D[i]+D[j]-(1/a[i,j])*(D[i]*sum((q)^(-1/2))+D[j]*sum((q2)^(-1/2)))->r
    return(r)
  }}

#####another function to output Ricci curvatures for all edges incident on a given node i
F_Ricci_row2<-function(i,a){
  return(mapply(function(j){return(F_Ricci2(i,j,a))},neighbs[[i]]))
}

#####computation of scaler Ricci curvature at each nodes, for each column of x
mcmapply(function(i){
  print(i)
  v<-x[,i]
  outer(v,v)->a##mass action princ
  a[which(int.red==0)]<-0##sparsify via network
  mapply(function(j){
    print(j)
    return(F_Ricci_row2(j,a))
  },1:ncol(int.red))->F_Ricci_mat
  mapply(function(k){
    return((sum(F_Ricci_mat[[k]])))
  },1:length(F_Ricci_mat))->F_n#
  return(F_n)},1:ncol(x),mc.cores = 8)->F_G_full_vec2

######compute stationary distribution for each column of x
mcmapply(function(i){
  print(i)
  v<-x[,i]
  outer(v,v)->a##mass action princ
  a[which(int.red==0)]<-0##sparsify via network
  apply(a,1,sum)->b###compute row sum
  #length(which(a!=t(a)))##good so symmetric random walk
  a/b->c##stoch mat
  b/sum(b)->mu
  return(mu)},1:ncol(x),mc.cores = 8)->mu_full_vec

#########compute total Ricci curvature for each column of x
mu_full_vec*D*F_G_full_vec2->F_G_pre_sum2
apply(F_G_pre_sum2,2,sum)->F_G2###this gives T values of total Ricci curvature, one for each column of x


#########compare with entropy rate
#########compute nodal entropies
mapply(function(i){
  print(i)
  v<-x[,i]
  outer(v,v)->a##mass action princ
  a[which(int.red==0)]<-0##sparsify via network
  apply(a,1,sum)->b###compute row sum
  #length(which(a!=t(a)))##good so symmetric random walk
  a/b->c##stoch mat
  apply(c,1,ent)->e
  return(e)},test)->ent_full_vec

##############compute signalling entropy
apply(mu_full_vec*ent_full_vec,2,sum)->SR

#######
plot(SR,F_G2)


#######iterate Ricci Flow
#######define starting and ending states of Ricci flow and corresponding edge curvatures:
########worked example for Chu et al., time course on average gene expression at 6 time points during differentiation
#######start state
v<-x[,1]
outer(v,v)->a
a[which(int.red==0)]<-0
a->a_1

########
mcmapply(function(j){
  print(j)
  return(F_Ricci_row2(j,a_1))
},1:ncol(int.red),mc.cores = 10)->F_Ricci_mat_1

#########end state
v<-x[,ncol(x)]
outer(v,v)->a
a[which(int.red==0)]<-0
a->a_2

##########
mcmapply(function(j){
  print(j)
  return(F_Ricci_row2(j,a_2))
},1:ncol(int.red),mc.cores = 10)->F_Ricci_mat_2

###########Pick largest Delta t (here called eps) s.t. keep positive weights post first iteration
seq(0.01,0.1,0.01)->eps_seq
mapply(function(k){
  R_flow_c<-function(a_t,F_Ricci_mat_t,F_Ricci_mat_2){return(
    mapply(function(i){
      mapply(function(j){
        return(1/(a_t[i,neighbs[[i]][j]])+eps_seq[k]*(F_Ricci_mat_t[[i]][j]-F_Ricci_mat_2[[i]][j])*1/(a_t[i,neighbs[[i]][j]]))},1:length(neighbs[[i]]))},1:ncol(a_t)))}
  ##############itterate weights one step
  R_flow_c(a_1,F_Ricci_mat_1,F_Ricci_mat_2)->a_t
  return(length(which(unlist(a_t)<0))/length(unlist(a_t)))},1:length(eps_seq))->prop_zero
cbind(eps_seq,prop_zero) ##this will output the optimal eps, here 0.06

####################
eps<-0.06
####################

############define ricci flow
R_flow_c<-function(a_t,F_Ricci_mat_t,F_Ricci_mat_2){return(
  mapply(function(i){
    mapply(function(j){
      return((1/a_t[i,neighbs[[i]][j]])+eps*(F_Ricci_mat_t[[i]][j]-F_Ricci_mat_2[[i]][j])*(1/a_t[i,neighbs[[i]][j]]))},1:length(neighbs[[i]]))},1:ncol(a_t)))}

############
#define the 4 intermediate states
v<-x[,2]
outer(v,v)->a
a[which(int.red==0)]<-0
a->a_1_1

###
v<-x[,3]
outer(v,v)->a
a[which(int.red==0)]<-0
a->a_1_2

###
v<-x[,4]
outer(v,v)->a
a[which(int.red==0)]<-0
a->a_1_3

###
v<-x[,4]
outer(v,v)->a
a[which(int.red==0)]<-0
a->a_1_4

#########set number of iterations to run and empty vectors to populate

K<-150
fts<-list()
ent_iters<-list()
mu_iters<-list()
cors_a1<-matrix(0,K,1)
cors_a1_1<-matrix(0,K,1)
cors_a1_2<-matrix(0,K,1)
cors_a1_3<-matrix(0,K,1)
cors_a1_4<-matrix(0,K,1)
cors_a2<-matrix(0,K,1)

##########initiate starting state
a_t<-a_1
F_Ricci_mat_t<-F_Ricci_mat_1

#############################run K iterations of the flow
for(l in 1:K){
  print(l)
  R_flow_c(a_t,F_Ricci_mat_t,F_Ricci_mat_2)->a_t
  #######
  a_t_2<-matrix(0,nrow = nrow(a_1),ncol=ncol(a_1))
  for(i in 1:length(a_t)){a_t_2[i,neighbs[[i]]]<-(1/a_t[[i]])}####invert a to get xixj weighted matrix back
  ###########output correlation with starting and ending states for a_t, can include further intermediate states if needed
  cors_a1[l]<-sum(sqrt(sum((a_1-a_t_2)^2)))
  cors_a1_1[l]<-sum(sqrt(sum((a_1_1-a_t_2)^2)))
  cors_a1_2[l]<-sum(sqrt(sum((a_1_2-a_t_2)^2)))
  cors_a1_3[l]<-sum(sqrt(sum((a_1_3-a_t_2)^2)))
  cors_a1_4[l]<-sum(sqrt(sum((a_1_4-a_t_2)^2)))
  cors_a2[l]<-sum(sqrt(sum((a_2-a_t_2)^2)))
  #####
  a_t_2->a_t
  ###output nodal entropies and stationary dist. for a_t
  apply(a_t,1,sum)->b
  a_t/b->c
  apply(c,1,ent)->e
  b/sum(b)->mu
  ###
  e->ent_iters[[l]]
  mu->mu_iters[[l]]
  #####output scalar Ricci curvature for each state
  mcmapply(function(j){
    print(j)
    return(F_Ricci_row2(j,a_t_2))
  },1:ncol(int.red),mc.cores = 8)->F_Ricci_mat_t
  mcmapply(function(k){return((sum(F_Ricci_mat_t[[k]])))},1:length(F_Ricci_mat_t),mc.cores = 8)->fts[[l]]
}

################################copmute scalar curvatures and signalling entropies for eack of K iterations
unlist(lapply(1:length(fts),function(x){return(sum(mu_iters[[x]]*D*fts[[x]]))}))->ft_as
c(F_G2[2],ft_as)->ft_as
# unlist(lapply(1:length(fts),function(x){return(sum(mu_iters[[x]]*ent_iters[[x]]))}))->SR_as
# c(SR[2],SR_as)->SR_as

####################################confirm convergence
plot(ft_as,type="l",lwd=2,xlab="k",
     col="blue",ylab="Total Forman-Ricci Curvature");abline(h=F_G2[c(1,ncol(x))],col="red")
plot(SR_as,type="l",lwd=2,ylim=c(6.707,6.717),xlab="k",
     col="blue",ylab="Signalling Entropy");abline(h=SR[c(1,ncol(x))],col="red")

#####################################
cbind(cors_a1,cors_a1_1,cors_a1_2,cors_a1_3,cors_a1_3,cors_a1_4,cors_a2)->cors_all
apply(cors_all_a,2,function(x){which(x==min(x))})###display first passing iteration for each intermediate data point

#############Euclidean trajectory null model
####convert to weighted graph
to_edge<-function(x){
  mapply(function(i){
    x[,i]->v
    outer(v,v)->a
    return(a[which(int.red!=0)])
  }, 1:ncol(x))->ed
  return(ed)}

####
to_edge(x)->x2
####

x2[,1]->x_1
x2[,ncol(x)]->x_n
####make straight line trajectory
y<-function(t){return(x_1+t*(x_n-x_1)/150)}

####compute distance between 6 intermediate time points
mapply(function(i){
  print(i)
  a<-y(i)
  return(c(sqrt(sum((a-x2[,1])^2)),
           sqrt(sum((a-x2[,2])^2)),
           sqrt(sum((a-x2[,3])^2)),
           sqrt(sum((a-x2[,4])^2)),
           sqrt(sum((a-x2[,5])^2)),
           sqrt(sum((a-x2[,6])^2))))
},1:150)->ds
t(ds)->ds
apply(ds,2,function(x){which(x==min(x))})
plot(apply(ds,2,function(x){which(x==min(x))}),c(0,12,24,36,72,96))
cor.test(apply(ds,2,function(x){which(x==min(x))})[-c(1,6)],c(0,12,24,36,72,96)[-c(1,6)])

################
times<-c(0,12,24,36,72,96)###true data point differentiation times

###RF, closest pass of trajectory plotted against true diff time
cor.test(times[-c(1,6)],apply(cors_all_a,2,function(x){which(x==min(x))})[-c(1,6)])->ct
plot(c(apply(cors_all_a,2,function(x){which(x==min(x))})[1:5],150),times,pch=16,col=
       c("darkred","red","pink","green","lightblue","blue"),
     main="RF",xlab="k",ylab="Sample time (hrs)",
     sub=paste(signif(ct$estimate,3),signif(ct$p.value,3)))#
abline(lm(times[-c(1,6)]~apply(cors_all_a,2,function(x){which(x==min(x))})[-c(1,6)]),col="black",lwd=1)

########Straight line, closest pass of trajectory plotted against true diff time
cor.test(times[-c(1,6)],apply(ds,2,function(x){which(x==min(x))})[-c(1,6)])->ct2
plot(apply(ds,2,function(x){which(x==min(x))}),times,pch=16,col=
       c("darkred","red","pink","green","lightblue","blue"),
     main="ED",xlab="k",ylab="Sample time (hrs)",
     sub=paste(signif(ct2$estimate,3),signif(ct2$p.value,3)))#
abline(lm(times[-c(1,6)]~apply(ds,2,function(x){which(x==min(x))})[-c(1,6)]),col="black",lwd=1)
dev.off()
