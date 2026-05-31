hgg_0 = function(x, beta, lambda) {
  C = (gamma(3 / beta) / gamma(1 / beta)) ^ (beta / 2)
  y = -(abs(x)/lambda)^beta*C-log(lambda)
  return(y)
}

hgg_1 = function(x, beta, lambda) {
  C = (gamma(3 / beta) / gamma(1 / beta)) ^ (beta / 2)
  dx0 = (abs(x) / lambda) ^ (beta - 1) * (beta * (abs(x) / lambda ^ 2)) * C - 1 /
    lambda
  #dx1 = -((-x / lambda) ^ (beta - 1) * (beta * (x / lambda ^ 2)) * C + 1 /
     #       lambda)
  #dx0[is.nan(dx0)] <- 0
  #dx1[is.nan(dx1)] <- 0
  #dx = dx0 + dx1
  return(dx0)
}

hgg_2=function(x,beta,lambda){
  C = (gamma(3 / beta) / gamma(1 / beta)) ^ (beta / 2)
  dx0=-(((abs(x) /lambda)^(beta - 1) * (beta * (abs(x)  * (2 * lambda)/(lambda^2)^2)) + 
           (abs(x) /lambda)^((beta - 1) - 1) * ((beta - 1) * (abs(x) /lambda^2)) * 
           (beta * (abs(x) /lambda^2))) * C - 1/lambda^2)
 # dx1=-(((-x/lambda)^((beta - 1) - 1) * ((beta - 1) * (x/lambda^2)) * 
  #         (beta * (x/lambda^2)) - (-x/lambda)^(beta - 1) * (beta * 
 #                                                              (x * (2 * lambda)/(lambda^2)^2))) * C - 1/lambda^2)
  #dx0[is.nan(dx0)]<-0
  #dx1[is.nan(dx1)]<-0
  #dx=dx0+dx1
  return(dx0)
}



ht_1 = function(x, df, lambda) {
  dx = -(-(df + 1)/2 * (2 * (x/lambda^2 * (x/lambda))/(df - 2)/(1 + (x/lambda)^2/(df - 2))) + 1/lambda)
  return(dx)
}

ht_2 = function(x, df, lambda) {
  dx =1/lambda^2 + -(df + 1)/2 * (2 * (x/lambda^2 * (x/lambda^2) + 
                                         x * (2 * lambda)/(lambda^2)^2 * (x/lambda))/(df - 2)/(1 + 
                                                                                                 (x/lambda)^2/(df - 2)) - 2 * (x/lambda^2 * (x/lambda))/(df - 
                                                                                                                                                           2) * (2 * (x/lambda^2 * (x/lambda))/(df - 2))/(1 + (x/lambda)^2/(df - 
                                                                                                                                                                                                                              2))^2)
  
  return(dx)
}


ht_1 = function(x, df, lambda) {
  dx =-(-(df + 1)/2 * (2 * (x/lambda^2 * (x/lambda))/(df)/(1 + (x/lambda)^2/(df))) + 
          1/lambda) 
    
  return(dx)
}

ht_2 = function(x, df, lambda) {
  dx =1/lambda^2 + -(df + 1)/2 * (2 * (x/lambda^2 * (x/lambda^2) + 
                                         x * (2 * lambda)/(lambda^2)^2 * (x/lambda))/(df)/(1 + (x/lambda)^2/(df)) - 
                                    2 * (x/lambda^2 * (x/lambda))/(df) * (2 * (x/lambda^2 * (x/lambda))/(df))/(1 + 
                                                                                                                 (x/lambda)^2/(df))^2)
  return(dx)
}
hPIV_1=function(x,m,nu,lambda){
  y=m * (2 * (x/lambda^2 * (x/lambda))/(1 + (x/lambda)^2)) + nu * 
    (x/lambda^2/(1 + (x/lambda)^2)) - 1/lambda
  return(y)
}
hPIV_2=function(x,m,nu,lambda){
  y=-(nu * (x * (2 * lambda)/(lambda^2)^2/(1 + (x/lambda)^2) - x/lambda^2 * 
              (2 * (x/lambda^2 * (x/lambda)))/(1 + (x/lambda)^2)^2) + m * 
        (2 * (x/lambda^2 * (x/lambda^2) + x * (2 * lambda)/(lambda^2)^2 * 
                (x/lambda))/(1 + (x/lambda)^2) - 2 * (x/lambda^2 * (x/lambda)) * 
           (2 * (x/lambda^2 * (x/lambda)))/(1 + (x/lambda)^2)^2) - 
        1/lambda^2)
  return(y)
}
hsgg_1=function(x,beta,nu,lambda){
  y=(abs(x)/lambda/(1 + nu*sign(x)))^(beta - 1) * (beta * (abs(x)/lambda^2/(1 + nu*sign(x))))/2 - 
    1/lambda
  return(y)
}
hsgg_2=function(x,beta,nu,lambda){
  y=-(((abs(x)/lambda/(1 + nu*sign(x)))^(beta - 1) * (beta * (abs(x) * (2 * lambda)/(lambda^2)^2/(1 + 
                                                                                  nu*sign(x)))) + (abs(x)/lambda/(1 + nu*sign(x)))^((beta - 1) - 1) * ((beta - 1) * 
                                                                                                                                    (abs(x)/lambda^2/(1 + nu*sign(x)))) * (beta * (abs(x)/lambda^2/(1 + nu*sign(x)))))/2 - 
        1/lambda^2)
  return(y)
}
