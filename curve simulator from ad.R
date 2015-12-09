##########################################################################################
# SETUP
##########################################################################################
# library(data.table);library(bit64)
# setwd("d:\\Archives\\Git\\curve-simulation\\")
# learn.rate.start=1e-5 # the start point of learn rate for fitting curves (vary from curve file to curve file)
ex.curve=fread("sim_modelinput_curve.csv")
ex.dim=fread("sim_input_setup.csv")
##########################################################################################
# Code part
##########################################################################################
# generate npv points
calc_ad=function(x){
  ad=function(x1,x2){(1 - ((1 - x1 * exp(log(0.5) / ex.curve$hl)) / (exp(1) ^ ((x2/ex.curve$cps) / 
                                                                                        (ex.curve$wks * ex.curve$max) * (-log(1 - ex.curve$hrf)) * 10))))}
  ad(x1=ad(x1=rep(0,length(ex.curve$wks)),x2=x),x2=x)
}
ad=calc_ad(ex.curve$spend)
ex.curve$beta=ex.curve$decomp/ad
# create id
dim.id=names(ex.curve)[grepl("_name",names(ex.curve))]
for (i in 1:length(dim.id)){
  tempid=dim.id[i]
  temp=data.table(a=unique(ex.curve[[tempid]]))
  temp$id=1:nrow(temp)
  setnames(temp,names(temp),c(dim.id[i],paste(unlist(strsplit(dim.id[i],"_name")),"_id",sep="")))
  ex.curve=merge(ex.curve,temp,by=dim.id[i],all.x=T)
}

# check sumup dim
bdgt_dim=ex.dim$bdgt_dim[ex.dim$bdgt_dim!=0]
dim=ex.dim$dim[ex.dim$dim!=0]
if(sum(dim %in% bdgt_dim)==0) bdgt_dim1="all_id" else bdgt_dim1=dim[dim %in% bdgt_dim]


# compute current spend and npv
ex.curve$npv=ex.curve$decomp
npv_current=ex.curve[,list(npv=sum(decomp)),by=c(dim)]

for.sp=ex.curve[!duplicated(ex.curve[,c(bdgt_dim),with=F]),]
sp_current=for.sp[,list(sp_current=sum(spend)),by=c(bdgt_dim1)]
curve_current=merge(npv_current,sp_current,by=c(bdgt_dim1),all.x=T)
curve_current=curve_current[curve_current$npv!=0,]

# generate spend points
curve_fit=ex.curve[ex.curve$npv!=0,]
percent=c(seq(0,8,by=.1),9:100)
percent_mat=matrix(rep(percent,nrow(curve_fit)),nr=nrow(curve_fit),byrow=T)
sp_mat=curve_fit$spend*percent_mat
colnames(sp_mat)=paste("sp_",1:length(percent),sep="")

# generate npv points
calc_npv=function(x){
  ad=function(x1,x2){(1 - ((1 - x1 * exp(1) ^ (log(0.5) / curve_fit$hl)) / (exp(1) ^ ((x2/curve_fit$cps) / 
                                                                                        (curve_fit$wks * curve_fit$max) * (-log(1 - curve_fit$hrf)) * 10))))}
  npv=curve_fit$beta*ad(x1=ad(x1=rep(0,length(curve_fit$wks)),x2=x),x2=x)
  npv
}


curve_fit$wks=as.numeric(curve_fit$wks)
curve_fit$max=as.numeric(curve_fit$max)
npv_mat=apply(sp_mat,2,calc_npv)
colnames(npv_mat)=paste("npv_",1:length(percent),sep="")

# fit curves
npv_mat=data.table(cbind(curve_fit[,c(dim),with=F],npv_mat))
sp_mat=data.table(cbind(curve_fit[,c(bdgt_dim,"all_id"),with=F],sp_mat))
sp_mat_shrink=sp_mat[!duplicated(sp_mat[,c(bdgt_dim),with=F]),]

sp_fit=sp_mat_shrink[,lapply(.SD,sum,na.rm=T),by=c(bdgt_dim1)]
npv_fit=npv_mat[,lapply(.SD,sum,na.rm=T),by=c(dim)]

curve_fit_final=merge(npv_fit,sp_fit,by=c(bdgt_dim1),all.x=T)
col_npv=which(names(curve_fit_final) %in% paste("npv_",1:length(percent),sep="")) 
col_sp=which(names(curve_fit_final) %in% paste("sp_",1:length(percent),sep="")) 


a=rep(0,nrow(curve_fit_final))
b=rep(0,nrow(curve_fit_final))

if(round(nrow(curve_fit_final)/20)==0) int=5 else int=round(nrow(curve_fit_final)/20)
iter=0
while(T & iter<10) {
  iter=iter+1
  if (iter %% 2==0) learn.rate.start=learn.rate.start/10 else learn.rate.start=learn.rate.start*10
  
  col_zero=which(a %in% 0)
  if (sum(col_zero)==0) {
    print("Curves Fitting is complete.")
    break
  } else {
    for(i in col_zero) {
      if (i%%int==0) print(paste("Curves: ",(round(i/nrow(curve_fit_final),digit=2) * 100),  "% Complete ", sep="",Sys.time()))
      tryCatch({
        x=as.vector(as.matrix(curve_fit_final[i,]))
        dataset=data.frame(id=x[col_sp],d=x[col_npv])
        dataset$d=dataset$d+rnorm(nrow(dataset),0,median(dataset$d)/100)
        a.start <- max(dataset$d)
        b.start <- learn.rate.start
        control1 <- nls.control(maxiter= 10000, minFactor= 1e-30, warnOnly= FALSE,tol=1e-05)
        nl.reg <- nls(d ~ a * (1-exp(-b * id)),data=dataset,start= list(a=a.start,b=b.start),
                      control= control1)
        a[i]=coef(nl.reg)[1]
        b[i]=coef(nl.reg)[2]
      },error=function(e){
        print(i)
      },finally=next
      )
    }
  }
}

final=data.frame(curve_fit_final[,dim,with=F],a=a,b=b)

# a adjustment
final=merge(final,curve_current,by=c(dim),all.x=T)
final$npv_p=final$a*(1-exp(-final$b*final$sp_current))
final$a=final$a*final$npv/final$npv_p

# output
final$bdgt_id=do.call(paste, c(final[bdgt_dim1], sep = "_"))
dim.name=paste(as.vector(do.call(cbind,strsplit(dim,"_id"))),"_name",sep="")
match=ex.curve[,c(dim,dim.name),with=F] 
match=match[!duplicated(match[,dim,with=F])]

final=merge(final,match,by=c(dim),all.x=T)
final$dim=do.call(paste, c(final[paste(as.vector(do.call(cbind,strsplit(dim,"_id"))),"_name",sep="")], sep = "_"))
setnames(final,"sp_current","sp")
if (for.opt){
  cps.table=data.frame(unique(ex.curve[,c(bdgt_dim,"cps"),with=F]))
  cps.table$bdgt_id=do.call(paste, c(cps.table[bdgt_dim1], sep = "_"))
  cps.table=cps.table[!duplicated(cps.table$bdgt_id),]
  final=merge(final,data.table(cps.table),by="bdgt_id",all.x=T)
  final$b=final$b*final$cps
  setnames(final,"cps","cpp")
}

# export files
write.csv(final,"sim_output_curve.csv",row.names=F)


  