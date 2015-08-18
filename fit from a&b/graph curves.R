#########################################################################################
# set up
#########################################################################################
library(ggplot2);library(scales);library(reshape2);library(data.table);library(RColorBrewer)
setwd("d:\\Archives\\R Code\\Curve simulator\\fit from a&b\\")
data=fread("output_curve.csv")
x.factor=2
y.factor=1.1

#########################################################################################
# code part
#########################################################################################
final=data
x=seq(0,x.factor*max(final$sp_current),length.out=500)
for.graph=matrix(0,nc=length(final$dim),nr=length(x))
colnames(for.graph)=final$dim
for (i in 1:nrow(final)){
  for.graph[,i]=final$a[i]*(1-exp(-final$b[i]*x))
}
for.graph=data.table(spend=x,for.graph)
#for.graph1=melt(for.graph,id="spend")

write.csv(for.graph,"sim_output_graph.csv",row.names=F)
