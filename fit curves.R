library(data.table);library(bit64)
setwd("c:\\Users\\XinZhou\\Documents\\GitHub\\kohls\\admin\\opt\\curve\\")
learn.rate.start=1e-6 # the start point of learn rate for fitting curves (vary from curve file to curve file)
for.opt=T # T if the fitting is for otptimization curves
is.refit=F # refit with a/b curve or use the adstock formula

source("c:\\Users\\XinZhou\\Documents\\GitHub\\curve-simulation\\curve simulator from ad.R")