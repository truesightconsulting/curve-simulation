#########################################################################################
# set up
#########################################################################################
library(ggplot2);library(scales);library(reshape2);library(data.table);library(RColorBrewer)
#setwd("d:\\Archives\\R Code\\Curve simulator\\")
data=fread("sim_output_curve.csv")
x.factor=2
y.factor=1.1

#########################################################################################
# code part
#########################################################################################
final=data
x=seq(0,x.factor*max(final$spend),length.out=500)
for.graph=matrix(0,nc=length(final$dim),nr=length(x))
colnames(for.graph)=final$dim
for (i in 1:nrow(final)){
  for.graph[,i]=final$a[i]*(1-exp(-final$b[i]*x))
}
for.graph=data.table(spend=x,for.graph)

for.graph1=melt(for.graph,id="spend")

p=ggplot(data=for.graph1)+aes(x=spend,y=value)+geom_line(size=1,aes(color=variable))+
  geom_point(data=final,aes(x=spend,y=npv,color=dim),size=4)+
  labs(x="Spend",y="Response",title="Response Curves",color="Dimension")+
  theme(plot.title = element_text(size = 24,face="bold",vjust=1),
        axis.title=element_text(face="bold"))

print(p)

write.csv(for.graph,"sim_output_graph.csv",row.names=F)
