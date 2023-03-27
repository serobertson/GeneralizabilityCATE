
# Code to run CASS example used in:
#``Regression-based estimation of heterogeneous treatment effects 
#when extending inferences from a randomized trial to a target population.''

#Note: some of the code in the source file
#is modified from code provided from ``Debiased machine learning of conditional
#average treatment effects and other causal functions.'' The Econometric journal. 2020.
#V Semenova and V Chernozhukov.

source('00_cate_fxns_source.R')

outpath=paste0("./FiguresCASS/")

#choose parameter options
lowy=-1
highy=1

max_knot=2
max_degree=3
#-------------

#SELECT first step (gam or glm)
#first step
mytype="glm"
#mytype="gam"


#SELECT second step (poly or splines)
#2nd step
#ss_method="poly"
ss_method="splines"

#-------------
library("ggplot2")
library("fda")
library("mgcv")


DFcass<-readRDS("DF_cass_for_ML_original_coding_for_pseudooutcome.rds")
head(DFcass)
str(DFcass)

#strategy: split from the beginning on MI, and analyze as 2 separate datasets 
DFcass_sub1=subset(DFcass, DFcass$prevmi==1) #MI=1
DFcass_sub0=subset(DFcass, DFcass$prevmi==0) #MI=0

#remove subgroup covariate from regression:
DFcass_sub1$prevmi=NULL
DFcass_sub0$prevmi=NULL

#-----------------------------------------------------------#
#Target population CATE

myestimator="phi" #original regression paper estimator

#First step: obtain pseudo-outcomes 
r2sub1fit<-DGM_Sim(DF=DFcass_sub1, regtype=mytype, estimator=myestimator)
head(r2sub1fit)
r2sub0fit<-DGM_Sim(DF=DFcass_sub0, regtype=mytype, estimator=myestimator)
head(r2sub0fit)

if(myestimator=="trialonly"){
r2sub1fit=subset(r2sub1fit, r2sub1fit$S==1)
r2sub0fit=subset(r2sub0fit, r2sub0fit$S==1)
}

cross_val_spline=F

#MI = 1
min=40
max=max(r2sub1fit$ejecfr)

#returns max t-statistic in each b run (here 200 runs)
str1_g1subMI1=wrappernew(tstat_provided=F,graph_xvarname="Ejection fraction",
                      titlename="History of MI", X=r2sub1fit$ejecfr,min=min, max=max,
                      Y=r2sub1fit$pseudo_est_all,outpath, 
                      lowy,highy,ss_method,max_degree,max_knot,cross_val_spline=cross_val_spline,eps=1)
#MI = 0
min=40
max=max(r2sub1fit$ejecfr)
str1_g1subMI0=wrappernew(tstat_provided=F,graph_xvarname="Ejection fraction",
                      titlename="No history of MI", X=r2sub0fit$ejecfr,min=min,max=max,
                      Y=r2sub0fit$pseudo_est_all,outpath, 
                      lowy,highy,ss_method,max_degree,max_knot,cross_val_spline=cross_val_spline,eps=1)

DF=data.frame(t1=str1_g1subMI1$max_tstatvec, t2=str1_g1subMI0$max_tstatvec)
head(DF)

#Obtain supremum t-statistic necessary for uniform inference

#max across the two t-values 
tstat_max=pmax(DF$t1, DF$t2)
alpha=0.05
tstat_order_statistic<-quantile(tstat_max,1-alpha)
print(tstat_order_statistic) 

#individual max t-statistic
tstat_order_statistic_MI1<-quantile(DF$t1,1-alpha) 
tstat_order_statistic_MI0<-quantile(DF$t2,1-alpha) 

##############
#Plotting the CATES in each subgroup

#to generate plots:
min=40
plot_str1_g1subMI1=wrappernew(tstat_provided=tstat_order_statistic_MI1,graph_xvarname="Ejection fraction",
                              titlename="History of MI", X=r2sub1fit$ejecfr,min=min, max=max,
                              Y=r2sub1fit$pseudo_est_all,outpath, 
                              lowy,highy,ss_method,max_degree,max_knot,cross_val_spline=cross_val_spline,eps=1,
                              plot=T)

#MI = 0
min=40
max=max(r2sub1fit$ejecfr)
plot_str1_g1subMI0=wrappernew(tstat_provided=tstat_order_statistic_MI0,graph_xvarname="Ejection fraction",
                              titlename="No history of MI", X=r2sub0fit$ejecfr,min=min,max=max,
                              Y=r2sub0fit$pseudo_est_all,outpath, 
                              lowy,highy,ss_method,max_degree,max_knot,cross_val_spline=cross_val_spline,eps=1,
                              plot=T)


library("cowplot")
cowplot::plot_grid(plot_str1_g1subMI1,plot_str1_g1subMI0, nrow=1, ncol=2)
setwd(outpath)
#gives 1:1 aspect ratio (square boxes each panel)
pdfFileName <- paste("Target_CATE_plot_step1_", mytype, "_step2_",ss_method, "_estimator_", myestimator,".pdf",sep="")
ggsave(pdfFileName, width = 20, height = 10)


#######
#to generate summary statistics of CI and CB widths
min=40
plot_str1_g1subMI1_df=wrappernew(tstat_provided=tstat_order_statistic_MI1,graph_xvarname="Ejection fraction",
                              titlename="History of MI", X=r2sub1fit$ejecfr,min=min, max=max,
                              Y=r2sub1fit$pseudo_est_all,outpath, 
                              lowy,highy,ss_method,max_degree,max_knot,cross_val_spline=cross_val_spline,eps=1,
                              plot=F)

table(plot_str1_g1subMI1_df$group)

#width, uniform:
UCI_MI1=subset(plot_str1_g1subMI1_df, plot_str1_g1subMI1_df$group=="UCI")
UCIL_MI1=subset(plot_str1_g1subMI1_df, plot_str1_g1subMI1_df$group=="UCIL")

width_UCI_MI1=UCIL_MI1$y-UCI_MI1$y

#width, pointwise:
PCI_MI1=subset(plot_str1_g1subMI1_df, plot_str1_g1subMI1_df$group=="PCI")
PCIL_MI1=subset(plot_str1_g1subMI1_df, plot_str1_g1subMI1_df$group=="PCIL")

width_PCI_MI1=PCIL_MI1$y-PCI_MI1$y

#MI = 0
min=40
max=max(r2sub1fit$ejecfr)
plot_str1_g1subMI0_df=wrappernew(tstat_provided=tstat_order_statistic_MI0,graph_xvarname="Ejection fraction",
                              titlename="No history of MI", X=r2sub0fit$ejecfr,min=min,max=max,
                              Y=r2sub0fit$pseudo_est_all,outpath, 
                              lowy,highy,ss_method,max_degree,max_knot,cross_val_spline=cross_val_spline,eps=1,
                              plot=F)

#width, uniform:
UCI_MI0=subset(plot_str1_g1subMI0_df, plot_str1_g1subMI0_df$group=="UCI")
UCIL_MI0=subset(plot_str1_g1subMI0_df, plot_str1_g1subMI0_df$group=="UCIL")

width_UCI_MI0=UCIL_MI0$y-UCI_MI0$y

#width, pointwise:
PCI_MI0=subset(plot_str1_g1subMI0_df, plot_str1_g1subMI0_df$group=="PCI")
PCIL_MI0=subset(plot_str1_g1subMI0_df, plot_str1_g1subMI0_df$group=="PCIL")

width_PCI_MI0=PCIL_MI0$y-PCI_MI0$y

DF_width_all=data.frame(estimator=myestimator, point=plot_str1_g1subMI0_df$x,CI_MI1=width_PCI_MI1,CB_MI1=width_UCI_MI1,
                        CI_MI0=width_PCI_MI0,CB_MI0=width_UCI_MI0)
csvFileName <- paste("DiffTarget_CATE_width_step1_", mytype, "_step2_",ss_method, "_estimator_", myestimator,".csv",sep="")
write.csv(DF_width_all,csvFileName)

summary(DF_width_all$CI_MI1)
#summary statistics
sum=summary(DF_width_all[,3:ncol(DF_width_all)])

sumCI_MI1=summary(DF_width_all$CI_MI1) 
sumCB_MI1=summary(DF_width_all$CB_MI1) 

sumCI_MI0=summary(DF_width_all$CI_MI0) 
sumCB_MI0=summary(DF_width_all$CB_MI0) 

summary_output=data.frame(estimator=myestimator,
                          statistic=c("Min", "1st Q",  "Median","Mean", "3rd Q","Max"), 
                          CI_MI1=as.vector(sumCI_MI1),CB_MI1=as.vector(sumCB_MI1),
                          CI_MI0=as.vector(sumCI_MI0),CB_MI0=as.vector(sumCB_MI0) )
csvFileName <- paste("DiffTarget_CATE_summarywidth_step1_", mytype, "_step2_",ss_method, "_estimator_", myestimator,".csv",sep="")
write.csv(summary_output,csvFileName)


