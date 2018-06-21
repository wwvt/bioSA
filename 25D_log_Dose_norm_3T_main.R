#2D-Dose Response Numerical Experiment
# LHD Design
# Created on 4/12/2016, Last update on 2/8/2017
# 


# install.packages("lhs", lib="~/R/lib");
# install.packages("mlegp", lib="~/R/lib");
# install.packages("R.matlab", lib="~/R/lib");



library(lhs)
library(mlegp)
library(R.matlab)

#-----------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------
rm(list=ls());

 data = read.csv("M0_LP.csv",header=FALSE);
# Xpred = read.csv("Pred_points.csv",header=FALSE);
X = data[, 1:25];
#---------------Get the specific values for each design configuration---------------
factor_info = read.table("factor_info.txt", sep = "", header=FALSE, strip.white=TRUE)

DesignX = read.table(file = "OA_L128_4_25.csv", sep = ",", header = TRUE);
dim_info = dim(DesignX);
 
for (k in 1:25)
{
  factor_no = k;
  #  print( paste("Now is factor", as.character(factor_no), "--", factor_info[factor_no,1]) );
  
  A = DesignX[, factor_no+1];
  idx_L1 = which(A == 1);
  idx_L2 = which(A == 2);
  idx_L3 = which(A == 3);
  idx_L4 = which(A == 4);
  
  X[idx_L1,factor_no] = factor_info[factor_no,2];
  X[idx_L2,factor_no] = factor_info[factor_no,3];
  X[idx_L3,factor_no] = factor_info[factor_no,4];
  X[idx_L4,factor_no] = factor_info[factor_no,5];
}

ranges = matrix(, nrow = 2, ncol = 25 );

for (k in c(1:4, 17:25))
{  
  ranges[1,k] = 0;
  ranges[2,k] = 1;
}


for (k in c(5:6, 8:9, 11:12, 14:15))
{  
  ranges[1,k] = 0;
  ranges[2,k] = 10;
}

for (k in c(7, 10, 13, 16))
{  
  ranges[1,k] = 0;
  ranges[2,k] = 4;
}


for (k in 1:25)
{  
  X[,k] = (X[,k]-ranges[1,k])/(ranges[2,k]-ranges[1,k]);  # standadize the parameter to [0,1]
}



X_full = do.call(rbind, replicate(3, X, simplify=FALSE));   # replicate the design points for 3 times


# four observations made a day, 28 observations made for a week, 250 observations made for 9 weeks.
# Cluster the observations into three time segments, i.e., seg 1: week 1;  seg 2: weeks 2-6; seg 3: weeks 7-9
# {  
# w1=1:28   Week 1
# w2=29:56  Week 2
# w3=57:84  
# w4=85:112
# w5=113:140
# w6=141:168  Week 6
# w7=169:196
# w8=197:224
# w9=225:250 Week 9
# }

t_full = NULL;
y_mean = NULL;


for (k in c(1:3)) 
{
t = k; # period index, each period includes multiple days of interest
t_col = do.call(rbind, replicate(128, t, simplify=FALSE));
t_col = matrix(t_col, ncol = 1, byrow = FALSE);
t_full = cbind(t_full, t_col);


time_seg = switch(k, 1:28, 29:168, 169:250);

days = as.matrix(c(time_seg)); #specify the time index to look at, range of values of days; 
y_temp =  data[, ((26-1)+days)];
y_temp = as.matrix(y_temp, ncol = 1, byrow = FALSE);
y_mean_temp = as.matrix(rowMeans(y_temp)); 
y_mean = cbind(y_mean, y_mean_temp); 
}

dim(t_full) = c(128*3,1);
dim(y_mean) = c(128*3,1);

# y_mean_max = pmax(y_mean, matrix(1,128*3,1));  # transform the counts to reduce the output variation
y_mean_max =  y_mean + matrix(1,128*3,1);  # transform the counts to reduce the output variation
log_y_mean = log(y_mean_max); 

Xt_full = cbind(X_full, t_full);

 



 

# === >>> SK fit and predict
# stochastic kriging with constant trend (q=0), gammaP=2(use gauss correlation function)

skriging_model_est = mlegp(Xt_full,log_y_mean,nugget = 1, nugget.known = 0);


save.image('M0_LP_norm_log_3T_compare.Rdata')  # save the fitted SK metamodel for one output type, remember to rename it for different outputs differently
