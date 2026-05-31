
library(doParallel)
stopImplicitCluster()
stopCluster(cl)
#showConnections(all = TRUE)


cl <- makeCluster(20)
registerDoParallel(cl)


replications=1000
#alpha0=0.06
alpha2_list=(0:7)*0.1

data_params_list=list(c(0,2),c(1,5),c(0,1),c(4,1.72,-0.01, 0.45, -0.01))
param_list=list(c(0.06,0.1096508,0.16),c(0.08,0.1201453,0.18),c(0.08,0.1206941,0.18),c(0.14,0.1800700,0.24))

result0<-NA
time1 <- Sys.time()
result0 <-  foreach(k0=4:4,.export=c('data_params_list',"alpha2_list"),.packages=c('gamlss.dist','stabledist'))%:% 
  foreach(l0=1:3)%:% 
  foreach(k1=1:8)%:% 
  foreach(l1= 1:2,.combine='cbind')%:% 
  foreach(i=1:replications,.combine='rbind')%dopar% {
    data_params=data_params_list[[k0]]
    alpha_list=param_list[k0]
    
    n_samples=500*2^(l1-1)
    alpha2=alpha2_list[k1]
    #print("=====================================")
    pvalues=matrix(0,1,1)
    eta = reta(n_samples + 1000, params = data_params)
    #theta0=c(0.1,0.2,0.9)
    #theta0=c(0.00001,0.1,0,0.9)
    #theta0=c(0.01,0.03,0.8,0.4)
    # theta0 = c(0.01, 0.03, 0.8, 0.2)
    # theta0 = c(0.01, 0.2, 0.9, 0.2)
    theta0=c(0.01, alpha_list[[1]][l0],0,alpha2)
    #theta0=c(0.01,0.03,0.8,0.2)
    #yt = rGarch(n_samples, theta0, eta, m = 1)
    yt=rGARCH2(n_samples, theta0, eta, m = 1)
    plot(yt,type="l")
    ngqmle_params = c(1,7)
    
    theta0=c(0.1,0.1,0.6,0.1)
    rlist=c(0.5)
    tryCatch({
      for (r0 in 1:1) {
        r=rlist[r0]
        results = TS_NGQMLE_process(yt,
                                    r,
                                    ngqmle_params,
                                    1,
                                    method = 1,
                                    Adaptive = 1)
        
        theta_tilde = results$GQMLE$tilde_theta
        theta_hat0 = results$TS_NGQMLE0$hat_theta
        theta_hat1 = results$TS_NGQMLE1$hat_theta
        lambda_hat0 = results$lambda_hat0
        lambda_hat1 = results$lambda_hat1
        
        qmle_params1 = results$qmle_params
        
        s = 1
        res2 = Sigma_h(m, s, r, yt, theta_hat1, qmle_params1, lambda_hat1)
        ss1 = (res2$pvalue < 0.05) + 0
        pvalues[1, 1*(r0-1)+1] = ss1[6]
        
        
      }
      
      
    }, error = function(e) {
      print(paste("Error:", e))
    })
    return(pvalues)
  }
time2 <- Sys.time()
print(time2-time1)
result0


res00=list()
#res0=matrix(NA,nrow=3*length(result0),ncol=2)

res0=matrix(NA,nrow=24*length(result0),ncol=2)
res0
for(i0 in 1:length(result0)) {
  result01 = result0[[i0]]
  for (i1 in 1:3) {
    result02 = result01[[i1]]
    for (i2 in 1:length(result02)) {
      res0[24 * (i0 - 1) + 8 * (i1 - 1) + i2, ] = colSums(result02[[i2]], na.rm = T) *
        100 / replications
    }
  }
}
res0
df2=data.frame(res0)
df2[is.na(df2)]=""
write.csv(df2, file = "test3_stb1_2.csv", row.names = FALSE) 

res_st5=res0
res_gg1=res0
res_stb1=res0

res0=res_st5

df2=read.csv(file = "C:\\Users\\gongh\\Documents\\test3_stb1.csv",header = TRUE) 
df2
res0=as.matrix(df2,ncol=16)



#res0
setwd("C:\\Users\\gongh\\Documents\\GNGQMLE -- JASA\\figures2")

df2=read.csv(file = "C:\\Users\\gongh\\Documents\\test3_gg2.csv",header = TRUE) 
res0=as.matrix(df2,ncol=16)
df2=read.csv(file = "C:\\Users\\gongh\\Documents\\test3_st5.csv",header = TRUE) 
res1=as.matrix(df2,ncol=16)
df2=read.csv(file = "C:\\Users\\gongh\\Documents\\test3_gg1.csv",header = TRUE) 
res2=as.matrix(df2,ncol=16)
df2=read.csv(file = "C:\\Users\\gongh\\Documents\\test3_stb1.csv",header = TRUE) 
res3=as.matrix(df2,ncol=16)
for(j in 1:3){
  pdf(paste0("test3_1_",j ,".pdf"),width = 3, height = 3)
  par(mfrow = c(1, 1),mar=c(2,2,0.5,0.5))
  pch_list=0:7
  X500=res0[(8*(j-1)+1):(8*(j-1)+8),1:8]/100
  X1000=res0[(8*(j-1)+1):(8*(j-1)+8),9:16]/100
  
  plot(alpha2_list,X500[1:8,1],type="l",ylim = c(0, 1),col="blue",lty = 1,xlab ="",ylab="" )
  #plot(alpha2_list,X500[1:8,1],type="p",col="blue",pch=7)
  lines(alpha2_list,X500[1:8,1],type="p",col="blue",pch=pch_list[1])
  
  for(i in 1:1){
    lines(alpha2_list,X1000[1:8,i],type="l",ylim = c(0, 1),col="red",lty = 2)
    lines(alpha2_list,X1000[1:8,i],type="p",col="red",pch=pch_list[i])
  }
  
  X500=res1[(8*(j-1)+1):(8*(j-1)+8),1:8]/100
  X1000=res1[(8*(j-1)+1):(8*(j-1)+8),9:16]/100
  lines(alpha2_list,X500[1:8,1],type="l",ylim = c(0, 1),col="blue",lty = 1 )
  lines(alpha2_list,X500[1:8,1],type="p",col="blue",pch=pch_list[2])
  lines(alpha2_list,X1000[1:8,1],type="l",ylim = c(0, 1),col="red",lty = 2)
  lines(alpha2_list,X1000[1:8,1],type="p",col="red",pch=pch_list[2])
  
  X500=res2[(8*(j-1)+1):(8*(j-1)+8),1:8]/100
  X1000=res2[(8*(j-1)+1):(8*(j-1)+8),9:16]/100
  lines(alpha2_list,X500[1:8,1],type="l",ylim = c(0, 1),col="blue",lty = 1 )
  lines(alpha2_list,X500[1:8,1],type="p",col="blue",pch=pch_list[3])
  lines(alpha2_list,X1000[1:8,1],type="l",ylim = c(0, 1),col="red",lty = 2)
  lines(alpha2_list,X1000[1:8,1],type="p",col="red",pch=pch_list[3])

  X500=res3[(8*(j-1)+1):(8*(j-1)+8),1:8]/100
  X1000=res3[(8*(j-1)+1):(8*(j-1)+8),9:16]/100
  lines(alpha2_list,X500[1:8,1],type="l",ylim = c(0, 1),col="blue",lty = 1 )
  lines(alpha2_list,X500[1:8,1],type="p",col="blue",pch=pch_list[4])
  lines(alpha2_list,X1000[1:8,1],type="l",ylim = c(0, 1),col="red",lty = 2)
  lines(alpha2_list,X1000[1:8,1],type="p",col="red",pch=pch_list[4])
  
  abline(h=0.05)
  dev.off()
}



df2=read.csv(file = "C:\\Users\\gongh\\Documents\\test3_gg2_2.csv",header = TRUE) 
res0=as.matrix(df2,ncol=16)
df2=read.csv(file = "C:\\Users\\gongh\\Documents\\test3_st5_2.csv",header = TRUE) 
res1=as.matrix(df2,ncol=16)
df2=read.csv(file = "C:\\Users\\gongh\\Documents\\test3_gg1_2.csv",header = TRUE) 
res2=as.matrix(df2,ncol=16)
df2=read.csv(file = "C:\\Users\\gongh\\Documents\\test3_stb1_2.csv",header = TRUE) 
res3=as.matrix(df2,ncol=16)
for(j in 1:3){
  pdf(paste0("test3_2",j ,".pdf"),width = 3, height = 3)
  par(mfrow = c(1, 1),mar=c(2,2,0.5,0.5))
  pch_list=0:7
  X500=res0[(8*(j-1)+1):(8*(j-1)+8),1:8]/100
  X1000=res0[(8*(j-1)+1):(8*(j-1)+8),9:16]/100
  
  plot(alpha2_list,X500[1:8,1],type="l",ylim = c(0, 1),col="blue",lty = 1,xlab ="",ylab="" )
  #plot(alpha2_list,X500[1:8,1],type="p",col="blue",pch=7)
  lines(alpha2_list,X500[1:8,1],type="p",col="blue",pch=pch_list[1])
  
  for(i in 1:1){
    lines(alpha2_list,X1000[1:8,i],type="l",ylim = c(0, 1),col="red",lty = 2)
    lines(alpha2_list,X1000[1:8,i],type="p",col="red",pch=pch_list[i])
  }
  
  X500=res1[(8*(j-1)+1):(8*(j-1)+8),1:8]/100
  X1000=res1[(8*(j-1)+1):(8*(j-1)+8),9:16]/100
  lines(alpha2_list,X500[1:8,1],type="l",ylim = c(0, 1),col="blue",lty = 1 )
  lines(alpha2_list,X500[1:8,1],type="p",col="blue",pch=pch_list[2])
  lines(alpha2_list,X1000[1:8,1],type="l",ylim = c(0, 1),col="red",lty = 2)
  lines(alpha2_list,X1000[1:8,1],type="p",col="red",pch=pch_list[2])
  
  X500=res2[(8*(j-1)+1):(8*(j-1)+8),1:8]/100
  X1000=res2[(8*(j-1)+1):(8*(j-1)+8),9:16]/100
  lines(alpha2_list,X500[1:8,1],type="l",ylim = c(0, 1),col="blue",lty = 1 )
  lines(alpha2_list,X500[1:8,1],type="p",col="blue",pch=pch_list[3])
  lines(alpha2_list,X1000[1:8,1],type="l",ylim = c(0, 1),col="red",lty = 2)
  lines(alpha2_list,X1000[1:8,1],type="p",col="red",pch=pch_list[3])
  
  X500=res3[(8*(j-1)+1):(8*(j-1)+8),1:8]/100
  X1000=res3[(8*(j-1)+1):(8*(j-1)+8),9:16]/100
  lines(alpha2_list,X500[1:8,1],type="l",ylim = c(0, 1),col="blue",lty = 1 )
  lines(alpha2_list,X500[1:8,1],type="p",col="blue",pch=pch_list[4])
  lines(alpha2_list,X1000[1:8,1],type="l",ylim = c(0, 1),col="red",lty = 2)
  lines(alpha2_list,X1000[1:8,1],type="p",col="red",pch=pch_list[4])
  
  abline(h=0.05)
  dev.off()
}

for(j in 1:3){
  #pdf(paste0("test3_2_",j ,".pdf"),width = 3, height = 3)
  #par(mfrow = c(1, 1),mar=c(2,2,0.5,0.5))
  pch_list=0:7
  X500=res0[(8*(j-1)+1):(8*(j-1)+8),1:8]/100
  X1000=res0[(8*(j-1)+1):(8*(j-1)+8),9:16]/100
  
  lines(alpha2_list,X500[1:8,1],type="l",col="blue",lty = 1 )
  #plot(alpha2_list,X500[1:8,1],type="p",col="blue",pch=7)
  lines(alpha2_list,X500[1:8,1],type="p",col="blue",pch=pch_list[2])
  
  for(i in 1:1){
    lines(alpha2_list,X1000[1:8,i],type="l",ylim = c(0, 1),col="red",lty = 2)
    lines(alpha2_list,X1000[1:8,i],type="p",col="red",pch=pch_list[2])
  }
 # abline(h=0.05)
  #dev.off()
}
df2=read.csv(file = "C:\\Users\\gongh\\Documents\\test3_gg1.csv",header = TRUE) 
#df2
res0=as.matrix(df2,ncol=16)
for(j in 1:3){
 # pdf(paste0("test3_2_",j ,".pdf"),width = 3, height = 3)
  #par(mfrow = c(1, 1),mar=c(2,2,0.5,0.5))
  pch_list=0:7
  X500=res0[(8*(j-1)+1):(8*(j-1)+8),1:8]/100
  X1000=res0[(8*(j-1)+1):(8*(j-1)+8),9:16]/100
  
  lines(alpha2_list,X500[1:8,1],type="l",col="blue",lty = 1 )
  #plot(alpha2_list,X500[1:8,1],type="p",col="blue",pch=7)
  lines(alpha2_list,X500[1:8,1],type="p",col="blue",pch=pch_list[3])
  
  for(i in 1:1){
    lines(alpha2_list,X1000[1:8,i],type="l",ylim = c(0, 1),col="red",lty = 2)
    lines(alpha2_list,X1000[1:8,i],type="p",col="red",pch=pch_list[3])
  }
  #abline(h=0.05)
  #dev.off()
}
df2=read.csv(file = "C:\\Users\\gongh\\Documents\\test3_stb1.csv",header = TRUE) 
#df2
res0=as.matrix(df2,ncol=16)
for(j in 1:3){
  # pdf(paste0("test3_2_",j ,".pdf"),width = 3, height = 3)
  #par(mfrow = c(1, 1),mar=c(2,2,0.5,0.5))
  pch_list=0:7
  X500=res0[(8*(j-1)+1):(8*(j-1)+8),1:8]/100
  X1000=res0[(8*(j-1)+1):(8*(j-1)+8),9:16]/100
  
  lines(alpha2_list,X500[1:8,1],type="l",col="blue",lty = 1 )
  #plot(alpha2_list,X500[1:8,1],type="p",col="blue",pch=7)
  lines(alpha2_list,X500[1:8,1],type="p",col="blue",pch=pch_list[4])
  
  for(i in 1:1){
    lines(alpha2_list,X1000[1:8,i],type="l",ylim = c(0, 1),col="red",lty = 2)
    lines(alpha2_list,X1000[1:8,i],type="p",col="red",pch=pch_list[4])
  }
  #abline(h=0.05)
  dev.off()
}



"test3_1_1.pdf"+str(1)
pdf("test3_1_2.pdf",width = 4, height = 4)
par(mfrow = c(1, 1),mar=c(2.5,2.5,1,1))
pch_list=0:7
X500=res0[,1:8]/100
X1000=res0[,9:16]/100
plot(alpha2_list,X500[1:8,1],type="l",ylim = c(0, 1),col="blue",lty = 1,xlab ="",ylab="" )
#plot(alpha2_list,X500[1:8,1],type="p",col="blue",pch=7)
lines(alpha2_list,X500[1:8,1],type="p",col="blue",pch=pch_list[1])
for(i in 2:4){
  lines(alpha2_list,X500[1:8,i],type="l",ylim = c(0, 1),col="blue",lty = 1)
  lines(alpha2_list,X500[1:8,i],type="p",col="blue",pch=pch_list[i])
}

for(i in 1:4){
  lines(alpha2_list,X1000[1:8,i],type="l",ylim = c(0, 1),col="red",lty = 2)
  lines(alpha2_list,X1000[1:8,i],type="p",col="red",pch=pch_list[i])
}
abline(h=0.05)
dev.off()




df2=data.frame(res0)
df2[is.na(df2)]=""
write.csv(df2, file = "test3_gg2.csv", row.names = FALSE) 



# 瀹氫箟鏁版嵁
x <- 1:10
y <- rnorm(10)

# 缁樺埗鍩烘湰鏁ｇ偣鍥?
plot(x, y, main = "Basic Scatter Plot", xlab = "X-axis", ylab = "Y-axis")

# 瀹氫箟棰滆壊鍜屽舰鐘?
colors <- 1:10
shapes <- c(1, 2, 3, 4, 5, 16, 17, 18, 19, 20)[1:10]

# 缁樺埗甯﹂鑹插拰褰㈢姸鐨勬暎鐐瑰浘
plot(x, y, col = colors, pch = 1, main = "Scatter Plot with Colors and Shapes")




# 瀹氫箟甯冨眬鐭╅樀
layout_matrix <- matrix(c(1, 3, 2, 4), nrow = 2, byrow = TRUE)
layout(layout_matrix)

# 缁樺埗鍥涗釜鍥撅紝鍏朵腑1鍜?涓轰竴涓粍鍚堬紝2鍜?涓哄彟涓€涓粍鍚?
plot(1:10, rnorm(10), main = "Plot 1")
plot(1:10, rnorm(10), main = "Plot 2")
plot(1:10, rnorm(10), main = "Plot 3")
plot(1:10, rnorm(10), main = "Plot 4")

