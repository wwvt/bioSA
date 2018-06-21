#25D-Dose Response Numerical Experiment, Global Sensitivity Analysis with data clustered into three time segments
# LHD Design
# Created on 4/24/2016, Last update on 4/24/2016
# 



rm(list=ls());
setwd("~/bioSA")

library(lhs)
library(mlegp)
library(R.matlab)


#-----------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------


 
load('M1_LP_norm_log_T2S10.Rdata');

N_SA = 10^4;

 
data_temp = readMat('X_SA_10k_norm_25.mat');
X_SA = data_temp$X;
X_SA = X_SA[,1:25];



 

# -----------------------------------------------------------------------

 
D_hat =array(0,dim=c(1,10)) # 1-dim: 10 time segments
 
 

for (i in 1:10)
{ 
  g0_hat = NULL;
  t_col = do.call(rbind, replicate(N_SA, i, simplify=FALSE));
  t_col = matrix(t_col, ncol = 1, byrow = FALSE);
  g0_hat  = predict(skriging_model_est,  cbind(X_SA,t_col));
  D_hat[i] = mean(g0_hat^2)-mean(g0_hat)^2;  # The estimated total variance
}
 
S1_hat=array(0,dim=c(10,25)) # 1-dim: 10 time segments; 2-dim: 25 different input parameters  

ST_hat=array(0,dim=c(10,25)) # 1-dim: 10 time segments; 2-dim: 25 different input parameters  

# 1 -----------------------------------------------------------------------
 

file_list = c("D1","D2", "D3","D4","D5","D6","D7","D8","D9", "D10", "D11","D12","D13","D14","D15","D16","D17","D18","D19","D20","D21","D22", "D23","D24","D25");

ptm <- proc.time();

for (i in 1:10)
{
  for (j in 1:25)
  {  
    D1_hat = NULL;
    Dt_hat = NULL;
    
    file_idx = j;
    file_name = paste("X_SA_", file_list[file_idx], "_10k_norm_25.mat", sep = ""); 
    data_temp = readMat(file_name);
    X_SA_D1 = data_temp$X.D1;
    X_SA_D2 = data_temp$X.D2;
    
    X_SA_T1 = data_temp$X.T1;
    X_SA_T2 = data_temp$X.T2;
    
    t_col = do.call(rbind, replicate(N_SA, i, simplify=FALSE));
    t_col = matrix(t_col, ncol = 1, byrow = FALSE);
    
    
   
    g0_hat_D1 = predict(skriging_model_est, cbind(X_SA_D1,t_col));
    g0_hat_D2 = predict(skriging_model_est, cbind(X_SA_D2,t_col));

    
    
    # D1_hat = mean(g0_hat_D1*g0_hat_D2)-mean(g0_hat)^2;  # The estimated total variance (method 1)
    D1_hat = D_hat[i] - 1/2*mean((g0_hat_D1-g0_hat_D2)^2);  #  method 2
    # D1_hat = ... ; # Method 3
    
    S1_hat[i,j] = D1_hat/D_hat[i];  # The estimated 1st-order Sobol's index for Xi
    
    
    
    g0_hat_T1 = predict(skriging_model_est, cbind(X_SA_T1,t_col));
    g0_hat_T2 = predict(skriging_model_est, cbind(X_SA_T2,t_col));
   
    
    # D_minus1_hat = mean(g0_hat_T1*g0_hat_T2)-mean(g0_hat)^2;  # The estimated total variance due to factors but the one of interest
    # ST_hat[j] = 1- D_minus1_hat/D_hat;  # The estimated total sensitivity index for Xi
    
    Dt_hat = 1/2*mean((g0_hat_T1-g0_hat_T2)^2);  # The estimated total variance  (guanranteed to be non-negative)
    
    ST_hat[i,j] =  Dt_hat/D_hat[i];  # The estimated total sensitivity index for Xi
    
  }
}

proc.time()-ptm

save.image("SA_M1_LP_norm_log_T2S10.Rdata")