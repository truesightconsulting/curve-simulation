library(data.table)
setwd("d:\\BOA Optimizer\\2014 T3\\roll up final curves\\")
ex.curve=fread("opt_modelinput_curve_npv.csv")
ex.dim=fread("sim_input_dim.csv")
learn.rate.start=1e-8
x.factor=1 # the spend used to fit curve is actuall spend/x.factor
seed=1088 # for generate random number to tweak curves
##########################################################################################
# Code part
##########################################################################################
bdgt_dim=ex.dim$bdgt_dim[ex.dim$bdgt_dim!=0]
dim=ex.dim$dim[ex.dim$dim!=0]
if(sum(dim %in% bdgt_dim)==0) bdgt_dim1="all_count" else bdgt_dim1=dim[dim %in% bdgt_dim]

curve$b1=curve$b/curve$cpp

# agg current spend and decomp
setkey(curve,"bdgt_id")
for.sp=unique(curve)
sp_current=for.sp[,list(sp_current=sum(sp_current)),by=c(bdgt_dim1)]


# generate spend points
percent=c(seq(0,8,by=.1),9:20)
percent_mat=matrix(rep(percent,nrow(curve)),nr=nrow(curve),byrow=T)
sp_mat=curve$sp_current*percent_mat
colnames(sp_mat)=paste("sp_",1:length(percent),sep="")

# generate npv points
calc_npv=function(x){
  npv=curve$a*(1-exp(-curve$b1*x))
  npv
}
npv_mat=apply(sp_mat,2,calc_npv)
colnames(npv_mat)=paste("npv_",1:length(percent),sep="")

# combine spend and npv data
npv_mat=data.table(cbind(curve[,c(dim),with=F],npv_mat))
sp_mat=data.table(cbind(curve[,c(bdgt_dim,"all_count"),with=F],sp_mat))
sp_mat_shrink=sp_mat[!duplicated(sp_mat[,c(bdgt_dim),with=F]),]

sp_fit=sp_mat_shrink[,lapply(.SD,sum,na.rm=T),by=c(bdgt_dim1)]
npv_fit=npv_mat[,lapply(.SD,sum,na.rm=T),by=c(dim)]

curve_fit_final=merge(npv_fit,sp_fit,by=c(bdgt_dim1),all.x=T)
col_npv=which(names(curve_fit_final) %in% paste("npv_",1:length(percent),sep="")) 
col_sp=which(names(curve_fit_final) %in% paste("sp_",1:length(percent),sep="")) 

# fit curves
a=rep(0,nrow(curve_fit_final))
b=rep(0,nrow(curve_fit_final))

if(round(nrow(curve_fit_final)/20)==0) int=5 else int=round(nrow(curve_fit_final)/20)
for (i in 1:nrow(curve_fit_final)){
  # print(i)
  set.seed(seed+i)
  if (i%%int==0) print(paste("Curves: ",(round(i/nrow(curve_fit_final),digit=2) * 100),  "% Complete ", sep="",Sys.time()))
  x=as.vector(as.matrix(curve_fit_final[i,]))
  dataset=data.frame(d=x[col_npv],id=x[col_sp])
  dataset$d=dataset$d+rnorm(n=nrow(dataset),mean=0,sd=0.1)
  a.start <- max(dataset$d)
  b.start <- learn.rate.start
  control1 <- nls.control(maxiter= 1000000, minFactor= 1e-30, warnOnly= FALSE,tol=1e-06)
  nl.reg <- nls(d ~ a * (1-exp(-b * id/x.factor)),data=dataset,start= list(a=a.start,b=b.start),
                control= control1)
  b[i] <- coef(nl.reg)[2]
  a[i] <- coef(nl.reg)[1]
}

if (length(dim)==1&dim=="all") {

  match=ex.curve[,dim,with=F]
  match=match[unique(dim)] 

  }else {
  dim.name=paste(as.vector(do.call(cbind,strsplit(dim,"_count"))),"_name",sep="")
  match=ex.curve[,c(dim,dim.name),with=F] 
  match=match[!duplicated(match[,dim,with=F])]
}

final=data.frame(curve_fit_final[,dim,with=F],a=a,b=b)
final=merge(final,match,by=c(dim),all.x=T)
final=merge(final,sp_current,by=c(dim),all.x=T)
final$dim=do.call(paste, c(final[paste(as.vector(do.call(cbind,strsplit(dim,"_count"))),"_name",sep="")], sep = "_"))
final=data.table(final)
final=final[,':='(predicted=a*(1-exp(-b*sp_current)))]
write.csv(final,"sim_output_curve.csv",row.names=F)
