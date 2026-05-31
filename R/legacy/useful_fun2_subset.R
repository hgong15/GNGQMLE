eta_lambda=function(x,lambda,w,params){
  #Cr0=moments0(r,w,params)
  #Cr=Cr0^(1/r)
  Cr=1/lambda
  #print(Cr0)
  if(w==0){
    return(log(dPE(x*Cr, nu =params[1])*Cr))
  }else if(w==1){
    t_df=params[1]
    return(log(dt(x*Cr, df = t_df)*Cr))
  }
  else if(w==2){
    m=params[1]
    nu=params[2]
    return(log(dpearsonIV(x,m= params[1], nu= params[2],location=0,scale=lambda)))
  }
  else if(w==3){
    return(dsgg(x,beta=params[1],nu= params[2],lambda=lambda,log=T ))
  }else if(w==4){
    return(log(dstable(x*Cr, alpha= params[1],beta= params[2],gamma= params[3],delta= params[4])*Cr))
  }
}

QMLEstimator_func=function(y,r,method=1)
{
  n=length(y)
  yt=y
  C11=1/(2*r^(1/r-1)*gamma(1/r))
  #print(C11)
  QMLEstimator=NULL
  log_likelihood=NULL
  if(method==1){
    L=function(lambda)
    { 
      theta=exp(lambda)
      result1=eta_theta(yt,theta)
      ht2=result1$ht
      eta2=result1$hat_eta
      return(mean(r/2*log(ht2)+abs(eta2)^r))
      # sum(log(exp(s)*sqrt(1+exp(a)*X^2))+((Y-p*X)/exp(s)/sqrt(1+exp(a)*X^2)/2))
    }
    opt=optim(log(c(0.01,0.2,0.9)),L)
    lambda=opt$par
    QMLEstimator=exp(lambda)
    log_likelihood=-opt$value
  }
  if(method==2){
    L=function(lambda)
    { 
      theta=c(20,5,5)*sigmoid(lambda)+c(0,1e-10,1e-10)
      #theta=sigmoid(lambda)
      result1=eta_theta(yt,theta)
      ht2=result1$ht
      eta2=result1$hat_eta
      return(mean(r/2*log(ht2)+abs(eta2)^r))
      # sum(log(exp(s)*sqrt(1+exp(a)*X^2))+((Y-p*X)/exp(s)/sqrt(1+exp(a)*X^2)/2))
    }
    #opt=optim(log(theta0),L)
    opt=optim(c(-2,-2,-2),L)
    #opt= nlminb(log(c(0.1,0.1,0.1)), L, lower =c(-18,-3,-3), upper = c(0,0.4,0.4))
    #opt= optim(c(0.1,0.1,0.1), L, NULL, method = "L-BFGS-B",lower =c(1e-8,0.05,0.05), upper = c(1,0.95,0.95))
    #opt= optim(log(c(0.1,0.1,0.1)), L, NULL, method = "BFGS")
    #opt=optim(c(0.25,0.15,0.6), L, NULL, method = "L-BFGS-B",lower =c(1e-8,0.05,0.1), upper = c(2,.95,0.95))
    #opt=optim(c(0.25,0.15,0.6), L, NULL,lower =c(1e-8,0.05,0.1), upper = c(2,.95,0.95))
    #QMLEstimator=exp(opt$par)
    lambda=opt$par
    #QMLEstimator=c(1,1,1)*sigmoid(lambda)+c(1e-5,1e-4,1e-4)
    QMLEstimator=c(20,5,5)*sigmoid(lambda)+c(0,1e-10,1e-10)
    log_likelihood=-opt$value
  }
  if(method==3){
    L=function(lambda)
    { 
      #theta=c(1,1,1)*sigmoid(lambda)+c(1e-5,1e-4,1e-4)
      #theta=sigmoid(lambda)
      theta=lambda
      result1=eta_theta(yt,theta)
      ht2=result1$ht
      eta2=result1$hat_eta
      return(mean(r/2*log(ht2)+abs(eta2)^r))
      # sum(log(exp(s)*sqrt(1+exp(a)*X^2))+((Y-p*X)/exp(s)/sqrt(1+exp(a)*X^2)/2))
    }
    #opt=optim(log(theta0),L)
    #opt=optim(c(-2,-2,-2),L)
    #opt= nlminb(log(c(0.1,0.1,0.1)), L, lower =c(-18,-3,-3), upper = c(0,0.4,0.4))
    opt= optim(theta0, L, NULL, method = "L-BFGS-B",lower =c(1e-10,1e-10,1e-10), upper = c(10,2,2))
    #opt= optim(log(c(0.1,0.1,0.1)), L, NULL, method = "BFGS")
    #opt=optim(c(0.25,0.15,0.6), L, NULL, method = "L-BFGS-B",lower =c(1e-8,0.05,0.1), upper = c(2,.95,0.95))
    #opt=optim(c(0.25,0.15,0.6), L, NULL,lower =c(1e-8,0.05,0.1), upper = c(2,.95,0.95))
    #QMLEstimator=exp(opt$par)
    lambda=opt$par
    QMLEstimator=lambda
    log_likelihood=-opt$value
    
  }
  log_likelihood1=n*log_likelihood/r+n*log(C11)
  return(list(QMLEstimator=QMLEstimator,log_likelihood=log_likelihood1,AIC=2*3-2*log_likelihood1,BIC=3*log(n)-2*log_likelihood1))
}
r=1
1/(2*r^(1/r-1)*gamma(1/r))

lambda_estimate_func=function(error_tilde,params)
{
  
  w=params[1]
  params1=params[-1]
  L2=function(lambda)
  {
    return(-mean(eta_lambda(error_tilde,exp(lambda),w,params1)))
  }
  lambda_hat=exp(optim(c(-1), L2, NULL, method = "L-BFGS-B",lower =c(-10), upper = c(10))$par)
  #opt0 = optim(c(0), L2)
  # opt0 =optimize(L2,c(-10,10), tol=1e-8)
  # #c(0, 1), tol = 0.0001,
  # lambda_hat=exp(opt0$objective)
  return(lambda_hat)
}



sigmoid=function(x){
  return(1/(1+exp(-x)))
}

lambda_estimate_func2 = function(error_tilde, method = 1,err=0.1)
{
  my_list <- list()
  obj_v = rep(1e10, 4)
  opt0 =NULL
  opt1=NULL
  opt2=NULL
  # s=sd(error_tilde)
  # #error_tilde=error_tilde+rnorm(length(error_tilde),0,err*s)
  # #error_tilde[error_tilde==0]=rnorm(sum(error_tilde==0),0,err*s)
  if (method == 1) {
    L0 = function(lambda)
    {
      theta1 = exp(lambda)
      #return(-mean(dPE(error_tilde / theta1[1], nu = theta1[2],log = T)) +log( theta1[1]))
      
      #error_tilde[error_tilde==0]=rnorm(sum(error_tilde==0),0,err*s)
      xx= dPE(error_tilde / theta1[1],  nu = theta1[2],log=T)
      #xx=xx[xx<log(10000)]
      return(-mean(xx)+log(theta1[1])
      )
    }
    opt0 = optim(c(0, 0), L0)
    
    L1 = function(lambda)
    {
      #theta1 = c(20, 29) * sigmoid(lambda) + c(1e-3, 1)
      theta1 = exp((lambda))
      s=1#sqrt(theta1[2]/(theta1[2]-2))
      return(-(mean(log(
        dt(error_tilde / theta1[1]*s, df = theta1[2]) / theta1[1]*s
      ))))
    }
    #opt1=nlminb(c(1,3), L1, lower = c(0.01,0.01), upper = c(20,15))
    opt1 = optim(c(0, 1), L1)
    
    
  }
  if (method == 2) {
    L0 = function(lambda)
    {
      #theta1 = c(20, 1.9) * sigmoid(lambda) + c(1e-2, 0.1)
      theta1 = c(20, 2) * sigmoid(lambda) + c(0.1, 0.5)
      
      xx= dPE(error_tilde / theta1[1],  nu = theta1[2],log=T)
      #xx=xx[xx<log(10000)]
      return(-(mean(xx)+log(theta1[1])
      ))
      
      # return(-mean(log(
      #   dPE(error_tilde / theta1[1], nu = theta1[2]) / theta1[1]
      # )))
    }
    opt0 = optim(c(1, 0), L0)
    #sigmoid(-1)*20
    #opt0= optim(c(0,0), L0)
    #opt0=nlminb(c(1,1), L0, lower = c(0.01,0.01), upper = c(20,3))
    #opt0= optim(c(1,1), L0, NULL, method = "L-BFGS-B",lower =c(0.01,0.01), upper = c(15,2.5))
    #opt0=nlminb(c(1,0,1,1), L0, lower = c(0.01,0,1,0.01), upper = c(20,0,1,2))
    L1 = function(lambda)
    {
      theta1 = c(20, 30) * sigmoid(lambda) + c(0.1, 2)
      return(-(mean(log(
        dt(error_tilde / theta1[1], df = theta1[2]) / theta1[1]
      ))))
    }
    opt1 = optim(c(1, -5), L1)
    #sigmoid((-1))
    # L2=function(lambda)
    # {
    #   theta1 = c(20, 200) * sigmoid(lambda[1:2]) + c(1e-3, 1)
    #   return(-(mean(log(dpearsonIV(error_tilde/theta1[1],m =theta1[2],nu=lambda[3],location=0, scale=1,log=T)-log(theta1[1])))))
    # }
    # 
    # opt2= optim(c(-3,-10,0), L2)
  }
  if (method== 3) {
    L0 = function(lambda)
    {
      theta1 = lambda
      return(-mean(log(
        dPE(error_tilde / theta1[1], nu = theta1[2]) / theta1[1]
      )))
    }
    #opt0 = optim(c(-3, -1), L0)
    #opt0= optim(c(0,0), L0)
    #opt0=nlminb(c(1,1), L0, lower = c(0.01,0.01), upper = c(20,3))
    opt0= optim(c(1,1), L0, NULL, method = "L-BFGS-B",lower =c(0.01,0.1), upper = c(15,2))
    #opt0=nlminb(c(1,0,1,1), L0, lower = c(0.01,0,1,0.01), upper = c(20,0,1,2))
    L1 = function(lambda)
    {
      theta1 = lambda
      return(-(mean(log(
        dt(error_tilde / theta1[1], df = theta1[2]) / theta1[1]
      ))))
    }
    opt1= optim(c(1,5), L1, NULL, method = "L-BFGS-B",lower =c(0.01,2.1), upper = c(15,20))
  }
  obj_v[1] = opt0$value
  obj_v[2] = opt1$value
  
  
  #print(opt0)
  #print(opt1)
  
  
  
  min_index <- which.min(obj_v)
  #print(obj_v)
  #print(min_index)
  #print(c(min_index-1,my_list[[min_index]]$par))
  theta1 = NULL
  # print(exp(opt0$par))
  # print(exp(opt1$par))
  
  #if (min_index - 1 == 0) {
  if (min_index - 1 == 0) {
    lambda = opt0$par
    if (method == 1) {
      theta1 =exp(lambda) 
    }
    if (method == 2) {
      theta1 = c(20, 2) * sigmoid(lambda) + c(0.1, 0.5)
    }
    if (method ==3) {
      theta1 = lambda
    }
  }
  if (min_index - 1 == 1) {
    lambda = opt1$par
    if (method == 1) {
      theta1 =exp(lambda) 
    }
    if (method == 2) {
      theta1 = c(20, 30) * sigmoid(lambda) + c(0.1, 2)
    }
    if (method==3) {
      theta1 = lambda
    }
  }
  #print(c(min_index - 1, theta1))
  
  return(c(min_index - 1, theta1))
  #print(c(min_index - 1, theta1))
}


lambda_estimate_func2 = function(error_tilde, method = 1,err=0.1)
{
  my_list <- list()
  obj_v = rep(1e10, 4)
  opt0 =NULL
  opt1=NULL
  opt2=NULL
  opt3=NULL
  # s=sd(error_tilde)
  # #error_tilde=error_tilde+rnorm(length(error_tilde),0,err*s)
  # #error_tilde[error_tilde==0]=rnorm(sum(error_tilde==0),0,err*s)
  if (method == 1) {
    L0 = function(lambda)
    {
      theta1 = exp(lambda)
      #return(-mean(dPE(error_tilde / theta1[1], nu = theta1[2],log = T)) +log( theta1[1]))
      
      #error_tilde[error_tilde==0]=rnorm(sum(error_tilde==0),0,err*s)
      xx= dPE(error_tilde / theta1[1],  nu = theta1[2],log=T)
      #xx=xx[xx<log(10000)]
      return(-mean(xx)+log(theta1[1])
      )
    }
    opt0 = optim(c(0, 0), L0)
    
    L1 = function(lambda)
    {
      #theta1 = c(20, 29) * sigmoid(lambda) + c(1e-3, 1)
      theta1 = exp((lambda))
      s=1#sqrt(theta1[2]/(theta1[2]-2))
      return(-(mean(log(
        dt(error_tilde / theta1[1]*s, df = theta1[2]) / theta1[1]*s
      ))))
    }
    #opt1=nlminb(c(1,3), L1, lower = c(0.01,0.01), upper = c(20,15))
    opt1 = optim(c(0, 1), L1)
    
    L2 = function(lambda)
    {
      #theta1 = c(20, 29) * sigmoid(lambda) + c(1e-3, 1)
      theta1 = c(exp(lambda[1]),exp(lambda[2]),lambda[3])
      s=1#sqrt(theta1[2]/(theta1[2]-2))
      return(-(mean(
        dpearsonIV(error_tilde ,m=theta1[2],nu=theta1[3], location=0,scale = theta1[1],log=T) 
      )))
    }
    #opt1=nlminb(c(1,3), L1, lower = c(0.01,0.01), upper = c(20,15))
    opt2 = optim(c(0, 1,0), L2)
    
    L3 = function(lambda)
    {
      #theta1 = c(20, 29) * sigmoid(lambda) + c(1e-3, 1)
      theta1 = c(exp(lambda[1]),exp(lambda[2]),lambda[3])
      s=1#sqrt(theta1[2]/(theta1[2]-2))
      return(-(mean(
        dsgg(error_tilde ,beta=theta1[2],nu=theta1[3],lambda = theta1[1],log=T) 
      )))
    }
    #opt1=nlminb(c(1,3), L1, lower = c(0.01,0.01), upper = c(20,15))
    opt3 = optim(c(0, 1,0), L3)
  }
  
  obj_v[1] = opt0$value
  obj_v[2] = opt1$value
  obj_v[3] = opt2$value
  obj_v[4] = opt3$value
  #print(opt0)
  #print(opt1)
  
  
  
  min_index <- which.min(obj_v)
 # print(obj_v)
  #print(min_index)
  #print(c(min_index-1,my_list[[min_index]]$par))
  theta1 = NULL
  # print(exp(opt0$par))
  # print(exp(opt1$par))
  
  if (min_index - 1 == 0) {
    lambda = opt0$par
    if (method == 1) {
      theta1 =exp(lambda) 
    }
    
  }
  if (min_index - 1 == 1) {
    lambda = opt1$par
    if (method == 1) {
      theta1 =exp(lambda) 
    }
  }
  if (min_index - 1 == 2) {
    lambda = opt2$par
    if (method == 1) {
      theta1 =  c(exp(lambda[1]),exp(lambda[2]),lambda[3])
    }
  }
  if (min_index - 1 == 3) {
    lambda = opt3$par
    if (method == 1) {
      theta1 =  c(exp(lambda[1]),exp(lambda[2]),lambda[3])
    }
  }
  print(c(min_index - 1, theta1))
  return(c(min_index - 1, theta1))
  #print(c(min_index - 1, theta1))
}



TS_NGQMLE_func = function(y, lambda, params, method = 1)
{
  n = length(y)
  yt = y
  w = params[1]
  params1 = params[-1]
  #cat("lambda:",lambda)
  #cat("\n parmas:",params)
  QMLEstimator = NULL
  log_likelihood = NULL
  #method=2
  if (method == 1) {
    L = function(theta1)
    {
      theta =exp(theta1)
      result1 = eta_theta(yt, theta)
      ht2 = result1$ht
      eta2 = result1$hat_eta
      return(- mean(eta_lambda(eta2, lambda, w, params1) - log(ht2) / 2))
    }
    opt = optim(log(theta0), L)
    theta1 = opt$par
    QMLEstimator =exp(theta1)
    log_likelihood = -opt$value
    
   
  }
  if (method == 2) {
    L = function(theta1)
    {
      theta = c(10,10,1)*sigmoid(theta1)+c(1e-5,1e-4,1e-4)
      result1 = eta_theta(yt, theta)
      ht2 = result1$ht
      eta2 = result1$hat_eta
      - mean(eta_lambda(eta2, lambda, w, params1) - log(ht2) / 2)
    }
    #cat("R#")
    opt = optim(c(-3.5, -3.5, -2), L)
    #opt= nlminb(log(c(0.1,0.1,0.1)), L, lower =c(-18,-3,-3), upper = c(0,0.4,0.4))
    #opt=optim(c(0.25,0.15,0.6), L, NULL, method = "L-BFGS-B",lower =c(1e-8,0.05,0.1), upper = c(2,.95,0.95))
    #QMLEstimator= sigmoid(-3.5)
    theta1 = opt$par
    #print(theta1)
    QMLEstimator =c(10,10,1)*sigmoid(theta1)+c(1e-5,1e-4,1e-4)
    log_likelihood = -opt$value
  }
  if (method == 3) {
    L = function(theta1)
    {
      theta = theta1
      result1 = eta_theta(yt, theta)
      ht2 = result1$ht
      eta2 = result1$hat_eta
      - mean(eta_lambda(eta2, lambda, w, params1) - log(ht2) / 2)
    }
    #cat("R#")
    #opt = optim(c(-2, -2, -2), L)
    #opt= nlminb(log(c(0.1,0.1,0.1)), L, lower =c(-18,-3,-3), upper = c(0,0.4,0.4))
    opt=optim(c(0.1,0.2,0.6), L, NULL, method = "L-BFGS-B",lower =c(0.0001,0.0001,0.005), upper = c(10,10,0.99))
    #QMLEstimator=
    theta1 = opt$par
    #print(theta1)
    QMLEstimator =theta1
    log_likelihood = -opt$value
  }
  return(list(QMLEstimator = QMLEstimator, log_likelihood = n*log_likelihood,AIC=2*3-2*log_likelihood*n,BIC=3*log(n)-2*log_likelihood*n))
}  




TS_NGQMLE_process=function(yt,r=2,NGQMLE_params=c(0,2),Cr,method=method,method_lambda=1,Adaptive=1){
  #The first step
  #cat("#The first step")
  result1=QMLEstimator_func(yt,r=r,method=method)
  tilde_theta=result1$QMLEstimator
  
  #print( result1)
  
  result2=eta_theta(yt,tilde_theta)
  
  #print( result2)
  
  hat_lambda=lambda_estimate_func(result2$hat_eta,NGQMLE_params)
  params2=NULL
  params3=NULL
  if(Adaptive == 1) {
    params2=lambda_estimate_func2(result2$hat_eta,method=method_lambda)
    params3 = params2[-2]
    
    #cat(",  params3:", params3, ", lambda_hat:", params2[2])
    
  }    #print("r2")
  #print(hat_lambda)
  ##The third step
  result3=TS_NGQMLE_func(yt,hat_lambda,NGQMLE_params,method=method)
  #result3=TS_NGQMLE_func(yt,1,NGQMLE_params)
  print(result3)
  result4= result3
  if(Adaptive==1){
    result4=TS_NGQMLE_func(yt,params2[2],params3,method=method)
  }
  #print(result4)
  #print("r4")
  hat_theta0=result3$QMLEstimator
  hat_theta1=result4$QMLEstimator
  return(list(GQMLE=list(tilde_theta=c(tilde_theta[1:2]/Cr^2,tilde_theta[3]),log_likelihood=result1$log_likelihood,AIC=result1$AIC,BIC=result1$BIC),
              TS_NGQMLE0=list(hat_theta=c(hat_theta0[1:2]/Cr^2,hat_theta0[3]),log_likelihood=result3$log_likelihood,AIC=result3$AIC,BIC=result3$BIC),
              TS_NGQMLE1=list(hat_theta=c(hat_theta1[1:2]/Cr^2,hat_theta1[3]),log_likelihood=result4$log_likelihood,AIC=result4$AIC,BIC=result4$BIC),
              lambda_hat0=hat_lambda, lambda_hat1=params2[2],qmle_params=params3
  ))
}



Simga2_hat=function(yt,theta,r=0.5,qmle_params){
  #theta=theta_hat1
  result=eta_theta(yt,theta)
  eta=result$hat_eta
  # cat("r-th moments:", mean(abs(eta)^r))
  # cat("2r-th moments:", mean(abs(eta)^(2*r)))
  ka_r=kappa_r(eta,r)
  # cat("qmle_params:", qmle_params,"\n")
  ka_f=kappa_f(eta,qmle_params[1],qmle_params[-1])
  
  v1=mean(theta[3]/(theta[2]*eta^2+theta[3]))
  v2=mean((theta[3]/(theta[2]*eta^2+theta[3]))^2)
  
  J=matrix(NA,2,2)
  J[1,1]=1/theta[2]^2
  J[1,2]=v1/theta[2]/theta[3]/(1-v1)
  J[2,1]=J[1,2]
  J[2,2]=(1+v1)*v2/theta[3]^2/(1-v1)/(1-v2)
  
  I_star=J
  I1=solve(I_star)
  
  e1=matrix(c(1,0),nrow=2)
  S20=ka_f*I1+(theta[2])^2*(ka_r-ka_f)*(e1%*%t(e1))
  S21=ka_r* I1
  
  
  #S_omega=ka_f*mean((psimga_ptheta$psimga_po)^2)+(theta[1])^2*(ka_r-ka_f)
  # cat("I_star^-1:",I1,"\n")
  return(list(S2G=S21,S2T=S20,I=I_star,I0=J))
}

Simga2_star=function(yt,theta,r=0.5,qmle_params){
  #theta=theta_hat1
  result=eta_theta_2(yt,theta)
  eta=result$hat_eta
  # cat("r-th moments:", mean(abs(eta)^r))
  # cat("2r-th moments:", mean(abs(eta)^(2*r)))
  ka_r=kappa_r(eta,r)
  # cat("qmle_params:", qmle_params,"\n")
  ka_f=kappa_f(eta,qmle_params[1],qmle_params[-1])
  # cat("kaf:",ka_f);
  #  cat(", kar:",ka_r,"; ")
  # # ##I_star
  psimga_ptheta=result$psimga_ptheta
  #print(psimga_ptheta$psimga_pa)
  s1=mean(psimga_ptheta$psimga_po^2)
  sigma_11=matrix(s1,1)
  sigma_23_1=matrix(c(mean(psimga_ptheta$psimga_po*psimga_ptheta$psimga_pa),mean(psimga_ptheta$psimga_po*psimga_ptheta$psimga_pb)),nrow=2)
  sigma_23_23=matrix(c(mean(psimga_ptheta$psimga_pa^2),mean(psimga_ptheta$psimga_pa*psimga_ptheta$psimga_pb),
                       0,mean(psimga_ptheta$psimga_pb^2)),nrow=2,byrow=T)
  sigma_23_23[2,1]= sigma_23_23[1,2]
  I_star=sigma_23_23-sigma_23_1%*%solve(sigma_11)%*%t(sigma_23_1)
  # print("sigma_23_23:")
  #  print(sigma_23_23)
  # # print(sigma_23_1%*%solve(sigma_11)%*%t(sigma_23_1))
  
  
  e1=matrix(c(1,0),nrow=2)
  I1=solve(I_star)
  S20=ka_f*I1+(theta[2])^2*(ka_r-ka_f)*(e1%*%t(e1))
  S21=ka_r* I1
  
  
  #S_omega=ka_f*mean((psimga_ptheta$psimga_po)^2)+(theta[1])^2*(ka_r-ka_f)
  # cat("I_star^-1:",I1,"\n")
  return(list(S2G=S21,S2T=S20,I=I_star,I0=sigma_23_23))
}

Simga1_hat=function(yt,theta,r=0.5,qmle_params){
  #theta=theta_hat1
  result=eta_theta_2(yt,theta)
  eta=result$hat_eta
  # cat("r-th moments:", mean(abs(eta)^r))
  # cat("2r-th moments:", mean(abs(eta)^(2*r)))
  ka_r=kappa_r(eta,r)
  # cat("qmle_params:", qmle_params,"\n")
  ka_f=kappa_f(eta,qmle_params[1],qmle_params[-1])
  # cat("kaf:",ka_f);
  #  cat(", kar:",ka_r,"; ")
  # # ##I_star
  psimga_ptheta=result$psimga_ptheta
  #print(psimga_ptheta$psimga_pa)
  s1=mean(psimga_ptheta$psimga_po^2)
  sigma_11=matrix(s1,1)
  sigma_23_1=matrix(c(mean(psimga_ptheta$psimga_po*psimga_ptheta$psimga_pa),mean(psimga_ptheta$psimga_po*psimga_ptheta$psimga_pb)),nrow=2)
  sigma_23_23=matrix(c(mean(psimga_ptheta$psimga_pa^2),mean(psimga_ptheta$psimga_pa*psimga_ptheta$psimga_pb),
                       0,mean(psimga_ptheta$psimga_pb^2)),nrow=2,byrow=T)
  sigma_23_23[2,1]= sigma_23_23[1,2]
  J=matrix(NA,3,3)
  J[1,1]=s1
  J[2:3,1]=sigma_23_1
  J[1,2:3]=t(sigma_23_1)
  J[2:3,2:3]=sigma_23_23
  b1=matrix(c(theta[1],theta[2],0),nrow=3)
  S1=ka_f*solve(J)+(ka_r-ka_f)*(b1%*%t(b1))
  SG=ka_r*solve(J)
  return(list(S1=S1,SG=SG,J=J))
}


kappa_r=function(eta,r){
  # cat("2r moments:",mean(abs(eta)^(2*r)))
  return((mean(abs(eta)^(2*r))-1)*4/r^2)
}


Test_2p=function(yt,theta_hat,n_samples=n_samples){
  result=eta_theta(yt,theta_hat)
  eta_hat=result$hat_eta
  x1=log(theta_hat[2]*eta_hat^2+theta_hat[3])
  gamma_hat=mean(x1)
  sigma_u_hat=mean(x1^2)-gamma_hat^2
  Tn=sqrt(n_samples)*gamma_hat/sqrt(sigma_u_hat)
  #cat("Tn: ",Tn,"\n")
  
  return(list(Tn=Tn,pvalue=pnorm(Tn)))
}

Test_1p=function(yt,theta_hat,n_samples=n_samples){
  result=eta_theta(yt,theta_hat)
  eta_hat=result$hat_eta
  x1=log(theta_hat[2]*eta_hat^2+theta_hat[3])
  gamma_hat=mean(x1)
  sigma_u_hat=mean(x1^2)-gamma_hat^2
  Tn=sqrt(n_samples)*gamma_hat/sqrt(sigma_u_hat)
  #cat("Tn: ",Tn,"\n")
  
  return(list(Tn=Tn,pvalue=1-pnorm(Tn)))
}


eta_theta=function(yt,theta){
  omega=theta[1]
  alpha=theta[2]
  beta=theta[3]
  n=length(yt)
  ht=rep(0,n)
  ht[1]=0.1
  eta=rep(0,n)
  psimga_po=rep(0,n)
  psimga_pa=rep(0,n)
  psimga_pb=rep(0,n)
  
  for(i in 2:length(yt)){
    ht[i]=omega+alpha*yt[i-1]^2+beta*ht[i-1]
    eta[i]=yt[i]/sqrt(ht[i])
    
    psimga_po[i]=psimga_po[i-1]+beta^(i-1)
    
    betav=exp(seq(i-2,0,by=-1)*log(beta))
    
    yt2=yt[1:(i-1)]^2
    psimga_pa[i]=sum(betav*yt2)
    #cat("i",i,"\n")
    #print(length(betav))
    #print(length(yt2))
    if(i>=3){
      betav2=exp(seq(i-3,0,by=-1)*log(beta))
      psimga_pb[i]=sum(seq(i-2,1,-1)* betav2*(omega+alpha*yt[1:(i-2)]^2))
    }
    
  }
  return(list(hat_eta=eta,ht=ht,psimga_ptheta=list(psimga_po=psimga_po/ht,psimga_pa=psimga_pa/ht,psimga_pb=psimga_pb/ht)))
}
