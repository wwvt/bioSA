getwd()
setwd("C:/Users/V/Documents")
rm(list=ls())
origdata <- read.csv("M0_LP_First.csv",header = TRUE) # X1-X25, 25*3
temp <- t(origdata)
#colnames(temp)<- c('T11','T12','T2S1','T2S2','T2S3','T2S4','T2S5','T2S6','T2S7','T2S8','T2S9','T2S10','T31','T32','T33','T34','T35','T36')
colnames(temp)<- c('T1S1','T1S2','T2S1','T2S2','T2S3','T2S4','T2S5','T2S6','T2S7','T2S8','T2S9','T2S10','T3S1','T3S2','T3S3','T3S4','T3S5','T3S6')

sa <- temp # for all
#sa <- log2(temp) # for Th17_LP


#install.packages("heatmaply")
library(heatmaply)
heatmaply(sa, dendrogram="row",
          #showticklabels = c(TRUE, TRUE),
          colors = c( "blanchedalmond","brown1"),
          #grid_gap = 0.1, 
          showticklabels = c(FALSE, TRUE),
          row_text_angle = 0, 
          column_text_angle = 0,
          fontsize_row = 11, fontsize_col = 11,
      #    main = "M2_LP: estimates of the Sobol' total indices.",
          main = "M0_LP: estimates of the Sobol' first-order indices",
          k_col = NA, k_row = NA,
          margins = c(50,100,50,10)
         )





