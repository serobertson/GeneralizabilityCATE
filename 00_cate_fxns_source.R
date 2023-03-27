wrapper<-function(graph_xvarname=graph_xvarname,titlename=titlename, X=X,Y=Y,outpath=outpath,
                  lowy=lowy,highy=highy,ss_method=ss_method,min=min, max=max,
                  max_degree=max_degree,
                  max_knot=max_knot,
                  cross_val_spline=cross_val_spline,eps){
  
  if (ss_method == "ortho_poly") {
    result<-least_squares_series(X=X,Y=Y,max_degree=max_degree,eps=eps,min=min, max=max)
    
  } 
  if (ss_method == "splines") {
    result<-least_squares_splines(X=X,Y=Y,max_knot=max_knot,eps=eps, min=min, max=max,
                                  norder=3,cross_val_spline=cross_val_spline) 
  }
  if (ss_method == "poly") {
    result<-least_squares_series_old(X=X,Y=Y,degree=max_degree,eps=eps,min=min, max=max)
  }
  len = length(result$X_grid)
  group <-c(rep("UCI",len), rep("PCI",len), rep("Estimate",len),rep("PCIL",len),rep("UCIL",len))
  group_type<- c(rep("CI",len), rep("CI",len), rep("Estimate",len),rep("CI",len),rep("CI",len))
  group_ci_type<-c(rep("Uniform",len), rep("Point",len), rep("Uniform",len),rep("Point",len),rep("Uniform",len))
  df<-data.frame(x=rep(result$X_grid,5), y = c(result$ghat.lower,result$ghat.lower.point,result$g.hat,result$ghat.upper.point,result$ghat.upper),group=group, group_col = group_type,group_line =group_ci_type )
  plot=make_plot(df,graph_xvarname=graph_xvarname,titlename=titlename,outpath=paste0(outpath,method,opt,ss_method,graph_xvarname,".png"),degree=result$degree,ss_method=ss_method,lowy=-1,highy=1)
  return(plot)
}

wrappernew<-function(tstat_provided=F,graph_xvarname=graph_xvarname,titlename=titlename, X=X,Y=Y,outpath=outpath,
                  lowy=lowy,highy=highy,ss_method=ss_method,min=min, max=max,
                  max_degree=max_degree,
                  max_knot=max_knot,
                  cross_val_spline=cross_val_spline,eps,plot=T){
  
  if (ss_method == "ortho_poly") {
    result<-least_squares_series(X=X,Y=Y,max_degree=max_degree,eps=eps,min=min, max=max)
    
  } 
  if (ss_method == "splines") {
    result<-least_squares_splines(tstat_provided,X=X,Y=Y,max_knot=max_knot,eps=eps, min=min, max=max,
                                  norder=3,cross_val_spline=cross_val_spline) 
  }
  if (ss_method == "poly") {
    result<-least_squares_series_old(tstat_provided,X=X,Y=Y,degree=max_degree,eps=eps,min=min, max=max)
  }
  
  if(tstat_provided != F ){
    #plot
   len = length(result$X_grid)
   group <-c(rep("UCI",len), rep("PCI",len), rep("Estimate",len),rep("PCIL",len),rep("UCIL",len))
   group_type<- c(rep("CI",len), rep("CI",len), rep("Estimate",len),rep("CI",len),rep("CI",len))
   group_ci_type<-c(rep("Uniform",len), rep("Point",len), rep("Uniform",len),rep("Point",len),rep("Uniform",len))
   df<-data.frame(x=rep(result$X_grid,5), y = c(result$ghat.lower,result$ghat.lower.point,result$g.hat,result$ghat.upper.point,result$ghat.upper),group=group, group_col = group_type,group_line =group_ci_type )
  plotfigure=make_plot(df,graph_xvarname=graph_xvarname,titlename=titlename,outpath=paste0(outpath,method,opt,ss_method,graph_xvarname,".png"),degree=result$degree,ss_method=ss_method,lowy=-1,highy=1)
     if(plot==T){
        object=plotfigure
     } else {object=df}
      
  return(object)} else {
    return(result)
  }
}

fit_models<-function(DF,regtype){
  #Model fitting
  DF_split1<-DF
  #---------------------------------- 
  #fit participation
  DF_split1$y<-NULL
  DF_split1$a<-NULL
  DF_split1$split<-NULL
  
  if(regtype=="glm"){
    knot=2
    norder=3
    psfit=glm(S ~ bsplineS(age, breaks = quantile(age, c(0:knot)/knot), norder = norder)[ ,-1] +
                bsplineS(ejecfr, breaks = quantile(ejecfr, c(0:knot)/knot), norder = norder)[ ,-1] +
                .-age -ejecfr
              , data=DF_split1, family="binomial")
    
    print("ps")
    print(coef(psfit))
    pshat<-predict(psfit, newdata=DF, type="response")
  }else if(regtype=="gam"){
    k=5
    w_reg <- gam(S~s(age,fx=TRUE, k=k) + s(ejecfr,fx=TRUE,k=k) +
               angina + beta + diuret + employ +job + lvscor +
                 nitro + sex +  smoke +actlim31 +recact31 +hyprtn31 + diabts31 +  lmca31 + prxves31 +  prxlad31  +  sysbp31 ,data=DF_split1,family=binomial)
    pshat<-predict(w_reg, type="response")
  }  
  #---------------------------------- 
  
  #fit pa
  DF_split1<-subset(DF, DF$S==1)
  pafit=glm(a ~ age + angina + ejecfr  + sysbp31 + prxlad31 + lvscor, data=DF_split1, family="binomial")
  print("pa")
  print(coef(pafit))
  pahat<-predict(pafit, newdata=DF, type="response")

  #---------------------------------- 
  #fit g1
  DFA1<-subset(DF,DF$a==1 & DF$S==1)
  DFA1_split2<-DFA1
  
  DFA1_split2$split<-NULL
  DFA1_split2$a<-NULL
  DFA1_split2$S<-NULL
  
  if(regtype=="glm"){
    mu1fit=glm(y ~ bsplineS(age, breaks = quantile(age, c(0:knot)/knot), norder = norder)[ ,-1] +
                 bsplineS(ejecfr, breaks = quantile(ejecfr, c(0:knot)/knot), norder = norder)[ ,-1] +
                 .-age-ejecfr, data=DFA1_split2,family="binomial")
    print("g1")
    print(coef(mu1fit))
    mu1hat<-predict(mu1fit, newdata=DF, type="response")
  } else if(regtype=="gam"){
    k=5
    mu1fit <- gam(y~s(age,fx=TRUE, k=k) + s(ejecfr,fx=TRUE,k=k) +
                    angina + beta + diuret + employ +job + lvscor +
                    nitro + sex +  smoke +actlim31 +recact31 +hyprtn31 + diabts31 +  lmca31 + prxves31 +  prxlad31  +  sysbp31,data=DFA1_split2,family=binomial)
    mu1hat<-predict(mu1fit, newdata=DF, type="response")
  }
  #---------------------------------- 
  #fit g0
  DFA0<-subset(DF,DF$a==0 & DF$S==1)
  DFA0_split2<-DFA0
  DFA0_split2$split<-NULL
  DFA0_split2$a<-NULL
  DFA0_split2$S<-NULL
  
  if(regtype=="glm"){ 
    mu0fit=glm(y ~ bsplineS(age, breaks = quantile(age, c(0:knot)/knot), norder = norder)[ ,-1] +
                 bsplineS(ejecfr, breaks = quantile(ejecfr, c(0:knot)/knot), norder = norder)[ ,-1] +
                 .-age-ejecfr, data=DFA0_split2,family="binomial")
    print("g0")
    print(coef(mu0fit))
    mu0hat<-predict(mu0fit, newdata=DF,type="response")
  } else if(regtype=="gam"){
    k=5
    mu0fit <- gam(y~s(age,fx=TRUE, k=k) + s(ejecfr,fx=TRUE,k=k) +
                    angina + beta + diuret + employ +job + lvscor +
                    nitro + sex +  smoke +actlim31 +recact31 +hyprtn31 + diabts31 +  lmca31 + prxves31 +  prxlad31  +  sysbp31,data=DFA0_split2,family=binomial)
    mu0hat<-predict(mu0fit, newdata=DF, type="response")
  }
  
  dataframe_predictions=data.frame(pahat=pahat,pshat=pshat,mu1hat=mu1hat, mu0hat=mu0hat)
  colnames(dataframe_predictions) <- c("pahat","pshat", "mu1hat","mu0hat")
  return(dataframe_predictions)
}



#rename to fit_nuisance_models
DGM_Sim<-function(DF, regtype, estimator="phi"){
  
  DF$y<-DF$Y
  DF$a<-DF$A
  
  DF$Y<-NULL
  DF$A<-NULL
  
  head(DF)
  
  
  if(estimator=="phi"){
  fitmod=fit_models(DF, regtype=regtype) 
  }else if(estimator=="trialonly"){
    DF=subset(DF, DF$S==1) #trial-only with covariates
    fitmod=fit_models(DF, regtype=regtype) 
  }
  #CALCULATE PSEUDO-OUTCOME
  a<-DF$a
  y<-DF$y
  S<-DF$S
  pa=fitmod$pahat
  ps=fitmod$pshat
  g1=fitmod$mu1hat
  g0=fitmod$mu0hat
  
  #trial-only pseudo
  pseudo_est <- ((a-pa)/(pa*(1-pa)))*(y-a*g1-(1-a)*g0) + g1-g0
  
  #generalizability pseudo
  if(estimator=="phi"){
    gA=a*g1+(1-a)*g0
    pseudo_est_all <- (S/ps)*((a-pa)/(pa*(1-pa)))*(y-a*g1-(1-a)*g0) + g1-g0
  }else if(estimator=="trialonly"){
    pseudo_est_all <- ((a-pa)/(pa*(1-pa)))*(y-a*g1-(1-a)*g0) + g1-g0
  }
  
  return(data.frame(DF, fitmod, pseudo_est, pseudo_est_all))
  
}


calculate_drl<-function(DF,pa, ps, g1,g0,type){
  a<-DF$A
  y<-DF$Y
  pseudo_est <- ((a-pa)/(pa*(1-pa)))*(y-a*g1-(1-a)*g0) + g1-g0
  DF$pseudo_est<-pseudo_est
  DF_split3<-subset(DF, DF$S==1)
  DF_split3$split<-NULL
  DF_split3$A<-NULL
  DF_split3$S<-NULL
  DF_split3$Y<-NULL
  if(type=="lasso"){
    drlfit=cv.glmnet(pseudo_est ~ ., data=DF_split3)
    coef(drlfit,s="lambda.min")
    drl<-predict(drlfit, newdata=DF,type="response",s="lambda.min")
  }else {
    drlfit=glm(pseudo_est ~ ., data=DF_split3)
    print("drl")
    print(coef(drlfit))
    drl<-predict(drlfit, newdata=DF)
  }
  return(drl)
}

calculate_drl_all<-function(DF,pa, ps, g1,g0,type){
  a<-DF$A
  y<-DF$Y
  S<-DF$S
  pseudo_est_all <- (S/ps)*((a-pa)/(pa*(1-pa)))*(y-a*g1-(1-a)*g0) + g1-g0
  DF$pseudo_est<-pseudo_est_all
  DF_split3<-DF
  DF_split3$split<-NULL
  DF_split3$A<-NULL
  DF_split3$S<-NULL
  DF_split3$Y<-NULL
  if(type=="lasso"){
    drlfit=cv.glmnet(pseudo_est ~ ., data=DF_split3)
    coef(drlfit,s="lambda.min")
    drl<-predict(drlfit, newdata=DF,type="response",s="lambda.min")
  }else {
    drlfit=glm(pseudo_est ~ ., data=DF_split3) # family="binomial"
    print("drl")
    print(coef(drlfit))
    drl<-predict(drlfit, newdata=DF)
  }
  return(drl)
}


least_squares_splinesnew<-function(X,Y,max_knot=9,eps=0.1,lower_q=0.025,upper_q=0.975,norder=4,nderiv=0,cross_val_spline=TRUE,...) {
  ## Choose a grid of points
  X_grid = seq(min(X),max(X),eps)
  ## Create technical regressors
 
  breaks<- quantile(X, c(0:knot)/knot)
  formula.bsp 	<- Y ~ smooth.spline(X)
  fit	<- lm(formula.bsp);
  cv_knot =1
  
  formula.bsp 	<- Y ~ smooth.spline(X)
  fit	<- lm(formula.bsp);
  regressors_grid<-cbind( rep(1,length(X_grid)), smooth.spline(X_grid)[ ,-1])
  
  ### Estimate fitted value of g(x) at the grid points
  g.hat<-regressors_grid%*%coef(fit)
  ###  Weighted Bootstrap for Coefficients
  regressors<-cbind( rep(1,length(X)), bsplineS(X, breaks =breaks, norder = norder, nderiv = nderiv)[ ,-1])
  bootstraped.coefs<-wboot(X=regressors,Y=Y,degree=dim(regressors)[2]-1)
  ### Estimated fitted value of g(x) at the grid points (bootstrap)
  bootstrapped.fitted.values<-regressors_grid%*%t(bootstraped.coefs)
  
  Omega.hat<-white_vcov(regressors,Y,b.hat=coef(fit))
  standard_error<-sqrt(diag(regressors_grid%*%Omega.hat%*%t(regressors_grid)))
  ### Lower Pointwise CI
  ghat.lower.point<-g.hat+qnorm(lower_q)*standard_error
  ### Upper Pointwise CI
  ghat.upper.point<-g.hat+qnorm(upper_q)*standard_error
  ## Max T-stat
  max_tstat<-tboot(bootstrapped.fitted.values,g.hat,standard_error,alpha = 1-(upper_q-lower_q))
  print("max_tstat:")
  print(max_tstat)
  ## Lower Uniform CI
  ghat.lower<-g.hat-max_tstat*standard_error
  ## Upper Uniform CI
  ghat.upper<-g.hat+max_tstat*standard_error
  return(list(ghat.lower=ghat.lower,g.hat=g.hat, ghat.upper=ghat.upper,X_grid=X_grid,fit=fit,degree=cv_knot,ghat.lower.point=ghat.lower.point,
              ghat.upper.point=ghat.upper.point))
}

least_squares_splines_rms<-function(X,Y,max_knot=9,eps=0.1,lower_q=0.025,upper_q=0.975,norder=4,nderiv=0,cross_val_spline=TRUE,...) {
  ## Choose a grid of points
  X_grid = seq(min(X),max(X),eps)
  ## Create technical regressors
  if (cross_val_spline) {
    cv.bsp<-rep(0,max_knot-1)
    for (knot in 2:max_knot) {
      breaks<- quantile(X, c(0:knot)/knot)
      formula.bsp 	<- Y ~ rcs(X, max_knot)[ ,-1]
      fit	<- lm(formula.bsp);
      cv.bsp[knot-1]		<- sum( (fit$res / (1 - hatvalues(fit)) )^2);
    }
    ## Number of knots chosen by cross-validation
    cv_knot<-which.min(cv.bsp)+1
    ## Breaks
    breaks<- quantile(X, c(0:cv_knot)/cv_knot)
  } else {
    cv_knot =1
    breaks = quantile(X, c(0:cv_knot)/cv_knot)
    
  }
  
  formula.bsp 	<-Y ~ rcs(X, cv_knot)[ ,-1]
  fit	<- lm(formula.bsp);
  regressors_grid<-cbind( rep(1,length(X_grid)),rcs(X_grid,cv_knot)[ ,-1])
  
  ### Estimate fitted value of g(x) at the grid points
  g.hat<-regressors_grid%*%coef(fit)
  ###  Weighted Bootstrap for Coefficients
  regressors<-cbind( rep(1,length(X)), rcs(X, cv_knot)[ ,-1])
  bootstraped.coefs<-wboot(X=regressors,Y=Y,degree=dim(regressors)[2]-1)
  ### Estimated fitted value of g(x) at the grid points (bootstrap)
  bootstrapped.fitted.values<-regressors_grid%*%t(bootstraped.coefs)
  
  Omega.hat<-white_vcov(regressors,Y,b.hat=coef(fit))
  standard_error<-sqrt(diag(regressors_grid%*%Omega.hat%*%t(regressors_grid)))
  ### Lower Pointwise CI
  ghat.lower.point<-g.hat+qnorm(lower_q)*standard_error
  ### Upper Pointwise CI
  ghat.upper.point<-g.hat+qnorm(upper_q)*standard_error
  ## Max T-stat
  max_tstat<-tboot(bootstrapped.fitted.values,g.hat,standard_error,alpha = 1-(upper_q-lower_q))
  ## Lower Uniform CI
  ghat.lower<-g.hat-max_tstat*standard_error
  ## Upper Uniform CI
  ghat.upper<-g.hat+max_tstat*standard_error
  return(list(ghat.lower=ghat.lower,g.hat=g.hat, ghat.upper=ghat.upper,X_grid=X_grid,fit=fit,degree=cv_knot,ghat.lower.point=ghat.lower.point,
              ghat.upper.point=ghat.upper.point))
}

least_squares_splines<-function(tstat_provided,X,Y,max_knot=9,eps=0.1,min,max,
                                lower_q=0.025,upper_q=0.975,
                                norder=4,nderiv=0,cross_val_spline=TRUE,...) {
  ## Choose a grid of points
  X_grid = seq(min,max,eps)
  
  ## Create technical regressors
  if (cross_val_spline) {
    cv.bsp<-rep(0,max_knot-1)
    for (knot in 2:max_knot) {
      breaks<- quantile(X, c(0:knot)/knot)
      formula.bsp 	<- Y ~ bsplineS(X, breaks =breaks, norder = norder, nderiv = nderiv)[ ,-1]
      fit	<- lm(formula.bsp);
      cv.bsp[knot-1]		<- sum( (fit$res / (1 - hatvalues(fit)) )^2);
    }
    ## Number of knots chosen by cross-validation
    cv_knot<-which.min(cv.bsp)+1
    ## Breaks
    breaks<- quantile(X, c(0:cv_knot)/cv_knot)
  } else {
    cv_knot=max_knot # supplied by user
    breaks = quantile(X, c(0:cv_knot)/cv_knot)
    
  }
  
  formula.bsp 	<- Y ~ bsplineS(X, breaks =breaks, norder = norder, nderiv = 0)[ ,-1]
  fit	<- lm(formula.bsp);
  regressors_grid<-cbind( rep(1,length(X_grid)), bsplineS(X_grid, breaks =breaks, norder = norder, nderiv = nderiv)[ ,-1])
  
  ### Estimate fitted value of g(x) at the grid points
  g.hat<-regressors_grid%*%coef(fit)
  ###  Weighted Bootstrap for Coefficients
  regressors<-cbind( rep(1,length(X)), bsplineS(X, breaks =breaks, norder = norder, nderiv = nderiv)[ ,-1])
  bootstraped.coefs<-wboot(X=regressors,Y=Y,degree=dim(regressors)[2]-1)
  ### Estimated fitted value of g(x) at the grid points (bootstrap)
  bootstrapped.fitted.values<-regressors_grid%*%t(bootstraped.coefs)
  
  Omega.hat<-white_vcov(regressors,Y,b.hat=coef(fit))
  standard_error<-sqrt(diag(regressors_grid%*%Omega.hat%*%t(regressors_grid)))
  ### Lower Pointwise CI
  ghat.lower.point<-g.hat+qnorm(lower_q)*standard_error
  ### Upper Uniform CI
  ghat.upper.point<-g.hat+qnorm(upper_q)*standard_error
  ## Max T-stat
  if(tstat_provided==F){
    max_tstatvec<-tboot(bootstrapped.fitted.values,g.hat,standard_error,alpha = 1-(upper_q-lower_q))
    alpha=0.05
    max_tstat<-quantile(max_tstatvec,1-alpha)
  } else{
    max_tstat<-tstat_provided #provided t-stat 
    max_tstatvec=NA
  }
  
  ghat.lower<-g.hat-max_tstat*standard_error
  # ## Upper Uniform CI
  ghat.upper<-g.hat+max_tstat*standard_error
  return(list(max_tstat=max_tstat,max_tstatvec=max_tstatvec,ghat.lower=ghat.lower,g.hat=g.hat, ghat.upper=ghat.upper,X_grid=X_grid,fit=fit,degree=cv_knot,ghat.lower.point=ghat.lower.point,
              ghat.upper.point=ghat.upper.point,standard_error=standard_error,bootstrapped.fitted.values=bootstrapped.fitted.values))
}


least_squares_series<-function(X, Y,max_degree,eps=0.1,min,max,
                               lower_q=0.025,upper_q=0.975,...) {
  cv.pol<-rep(0,max_degree)
  for (degree in 1:max_degree) {
    formula.pol 	<- Y ~ poly(X, degree)
    fit	<- lm(formula.pol );
    cv.pol[degree]		<- sum( (fit$res / (1 - hatvalues(fit)) )^2);
  }
  ## Number of knots chosen by cross-validation
  cv_degree<-which.min(cv.pol)
  ## Estimate coefficients
  formula.pol 	<- Y ~ poly(X, cv_degree)
  fit	<- lm(formula.pol);
  
  ## Choose a grid of points
  #X_grid = seq(min(X),max(X),eps)
  X_grid = seq(min,max,eps)
  ## Create technical regressors
  regressors_grid<-cbind( rep(1,length(X_grid)), poly(X_grid,cv_degree))
  ### Estimate fitted value of g(x) at the grid points
  g.hat<-regressors_grid%*%coef(fit)
  ###  Weighted Bootstrap for Coefficients
  regressors<-cbind( rep(1,length(X)), poly(X,cv_degree))
  bootstraped.coefs<-wboot(X=regressors,Y=Y,degree=dim(regressors)[2]-1)
  ### Estimated fitted value of g(x) at the grid points (bootstrap)
  bootstrapped.fitted.values<-regressors_grid%*%t(bootstraped.coefs)
  
  Omega.hat<-white_vcov(regressors,Y,b.hat=coef(fit))
  standard_error<-sqrt(diag(regressors_grid%*%Omega.hat%*%t(regressors_grid)))
  ### Lower Pointwise CI
  ghat.lower.point<-g.hat+qnorm(lower_q)*standard_error
  ### Upper Uniform CI
  ghat.upper.point<-g.hat+qnorm(upper_q)*standard_error
  ## Max T-stat
  max_tstat<-tboot(bootstrapped.fitted.values,g.hat,standard_error,alpha = 1-(upper_q-lower_q))
  ## Lower Uniform CI
  ghat.lower<-g.hat-max_tstat*standard_error
  ## Upper Uniform CI
  ghat.upper<-g.hat+max_tstat*standard_error
  
  
  return(list(ghat.lower=ghat.lower,g.hat=g.hat, ghat.upper=ghat.upper,X_grid=X_grid,fit=fit,degree=cv_degree,ghat.lower.point=ghat.lower.point,
              ghat.upper.point=ghat.upper.point))
}


least_squares_series_old<-function(tstat_provided,X, Y,degree,eps=0.1,min, max, lower_q=0.025,upper_q=0.975,...) {
  X_grid = seq(min,max,eps)
  
  regressors_grid = matrix(1,length(X_grid),degree+1)
  
  X_power = matrix(1,length(X),degree+1)
  for (p in 2:(degree+1)) {
    X_power[,p]<-X^(p-1)
    regressors_grid[,p]<-X_grid^(p-1)
  }
  ## Choose a grid of points
  fit<-lm(Y~X_power-1)
  ### Estimate fitted value of g(x) at the grid points
  g.hat<-regressors_grid%*%coef(fit)
  ###  Weighted Bootstrap for Coefficients
  bootstraped.coefs<-wboot(X=X_power,Y=Y,degree=dim(X_power)[2]-1)
  ### Estimated fitted value of g(x) at the grid points (bootstrap)
  bootstrapped.fitted.values<- regressors_grid%*%t(bootstraped.coefs)
  Omega.hat<-white_vcov(X_power,Y,b.hat=coef(fit))
  standard_error<-sqrt(diag(regressors_grid%*%Omega.hat%*%t(regressors_grid)))
  ### Lower Pointwise CI
  ghat.lower.point<-g.hat+qnorm(lower_q)*standard_error
  ### Upper Uniform CI
  ghat.upper.point<-g.hat+qnorm(upper_q)*standard_error
  ## Max T-stat
  if(tstat_provided==F){
  max_tstatvec<-tboot(bootstrapped.fitted.values,g.hat,standard_error,alpha = 1-(upper_q-lower_q))
  alpha=0.05
  max_tstat<-quantile(max_tstatvec,1-alpha)
  } else{
    max_tstat<-tstat_provided #provided t-stat 
    max_tstatvec=NA
  }
  
  # ## Lower Uniform CI
   ghat.lower<-g.hat-max_tstat*standard_error
  # ## Upper Uniform CI
   ghat.upper<-g.hat+max_tstat*standard_error
  # 
   return(list(max_tstat=max_tstat,max_tstatvec=max_tstatvec,ghat.lower=ghat.lower,g.hat=g.hat, ghat.upper=ghat.upper,X_grid=X_grid,fit=fit,degree=degree,ghat.lower.point=ghat.lower.point,
               ghat.upper.point=ghat.upper.point,standard_error=standard_error,bootstrapped.fitted.values=bootstrapped.fitted.values))
}


white_vcov<-function(X,y,b.hat,...) {
  # take subset of data
  # compute the estimated residuals
  N<-length(y)
  p.OLS<-dim(X)[2]
  error.sq<-(y-X%*%b.hat)^2
  # 
  Q<-t(X)%*%X/N
  Sigma<-matrix(0,p.OLS,p.OLS)
  for (i in 1:N) {
    Sigma<-Sigma+X[i,]%*%t(X[i,])*error.sq[i]
  }
  Sigma<-Sigma/N
  Omega.hat<-solve(Q)%*%Sigma%*%solve(Q)/N
  return(Omega.hat)
}

wboot<-function(X,Y,degree,B=200) {
  ## Sample size
  set.seed(88)
  N<-length(Y)
  b.hat<-matrix(0,B,degree+1)
  for (b in 1:B) {
    
    h<-rexp(N)
    h<-h/sum(h)
    b.hat[b,]<-coef(lm(Y~X-1,weights=h))
  }
  return(b.hat)
}

tboot<-function(bootstrapped.fitted.values,g.hat,standard_error,alpha) {
  B = dim(bootstrapped.fitted.values)[2]
  tstat<-abs((bootstrapped.fitted.values-matrix(rep(g.hat,B),ncol=B))/matrix(rep(standard_error,B),ncol=B))
  max_tstatvec<-apply(tstat,2,max) 
  return(max_tstatvec)
}



make_plot<-function(df,graph_xvarname, titlename, outpath,lowy,highy,degree,ss_method = "series",...) {
  show(degree)
  if (ss_method == "splines") {
    title<-paste0(titlename)
  } else {
    title<-paste0(titlename)
  }
  g<-ggplot(data=df)+
    aes(x=x,y=y,colour=group )+
    theme_classic()+
    xlab(graph_xvarname)+
    ylab("CATE")+  geom_hline(yintercept=0, linetype="solid", color = "black") +
    scale_colour_manual(values=c("black","grey","grey","grey","grey"))+
    theme(plot.title = element_text(size=40,hjust = 0.5,family="serif"),text=element_text(size=40,hjust = 0.5,family="serif"))+
    theme(legend.title=element_blank())+
    theme(legend.position="none")+
    ylim(low=lowy,high=highy)+
    geom_line(aes(linetype = group_line),size=1.5)+
    scale_linetype_manual(values=c("dashed","solid"))+ theme(axis.ticks.length = unit(0.5, "cm")) +
    ggtitle(title) 
  plot(g)
  return(g)
}
