# Code to select the r-k pairs that have the highest correlation between the F time-series between 3 species.
#1/19/2018
#############################################################################################################

# Clear workspace
rm(list = ls())

# Packages
library(plyr)
library(dplyr)
library(tidyr)
library(reshape2)
library(RColorBrewer)

setwd('C:/Users/kkleisner/Documents/ContractWork/Multi-Species_COMs/')
load('ms_cmsy.Rdata')

##Index the rk pairs:
a<-length(rv1)
b<-length(rv2)
c<-length(rv3)

# Calculate correlation matrix 
#utv <- rbind(utv1, utv2, utv3)
#cov_mat <- utv %*% t(utv)
corr_mat <- cov2cor(cov_mat)

# ##Check dim of supermatrix:
# dim(corr_mat)
# corr_mat2 <- corr_mat[1:a,(a+1):(a+b)]
# corr_mat3 <- corr_mat[1:a, (a+b+1):(a+b+c)]
# corr_mat4 <- corr_mat[(a+1):(a+b), 1:a]
# corr_mat6 <- corr_mat[(a+1):(a+b), (a+b+1):(a+b+c)]
# corr_mat7 <- corr_mat[(a+b+1):(a+b+c), 1:a]
# corr_mat8 <- corr_mat[(a+b+1):(a+b+c), (a+1):(a+b)]
# 
# ##Check dims of each matrix:
# dim(corr_mat2) #3536 5603
# dim(corr_mat3) #3536 5802
# dim(corr_mat4) #5603 3536
# dim(corr_mat6) #5603 5802
# dim(corr_mat7) #5802 3536
# dim(corr_mat8) #5802 5603


###This code will look at the normalized covariance matrix, 
###i.e., the 'supermatrix' and run each 'row' (1:3) pertaining to species 1-3.
###Each loop is for a 'row' in the supermatrix
###For row 1 we look at the correlation with species 1 and with species 2 (block 2) and with species 3 (block 3)
###For row 2 we look at the correlation with species 2 and with species 1 (block 4) and with species 3 (block 6)
###For row 3 we look at the correlation with species 3 and with species 1 (block 7) and with species 2 (block 8)
###We want the position of the maximum correlation in each row in each block
###Then, to calculate the correlation between species 1, 2 and 3, we sum the maximum correlations in each row and
###use some filter to select the maximum sums (here we use top 10 for illustration).
###The final outputs: Out1df_sort, Out2df_sort, Out3df_sort have the index positions for the r-k pairs that need to
###selected from the k and r vectors.

Max2<-NULL
Max3<-NULL
Idx2<-NULL
Idx3<-NULL
for (i in 1:a){
  max_cor<-0
  max_idx<-0
  for (j in (a+1):(a+b)){
    if(max_cor<corr_mat[i,j]){
      max_cor=corr_mat[i,j] 
      max_idx<-j
    }
  }
  Max2<-rbind(Max2,max_cor)
  Idx2<-rbind(Idx2, max_idx)
  
  max_cor<-0
  max_idx<-0
  for (j in (a+b+1):(a+b+c)){
    if(max_cor<corr_mat[i,j]){
      max_cor=corr_mat[i,j] 
      max_idx<-j
    }
  }
  Max3 <-rbind(Max3,max_cor)
  Idx3<-rbind(Idx3, max_idx)
}
Out1<-cbind(Max2+Max3, 1:a,Idx2-a,Idx3-(a+b))
Out1df<-data.frame(Out1, row.names = NULL)
##Sum of the maximum correlations, index for species 1, index for species 2, index for species 3: 
names(Out1df)<-c("sum", "index1","index2", "index3")
Out1df_sort<-Out1df[order(Out1df$sum, decreasing = T),]
###Decide how many max sums to select (this is top 10):
Out1df_sort<-Out1df_sort[1:10,]

Max4<-NULL
Max6<-NULL
Idx4<-NULL
Idx6<-NULL
for (i in 1:b){ 
  max_cor<-0
  max_idx<-0
  for (j in 1:a){
    if(max_cor<corr_mat[(a+i),j]){
      max_cor=corr_mat[(a+i),j] 
      max_idx<-j
    }
  }
  Max4<-rbind(Max4,max_cor)
  Idx4<-rbind(Idx4, max_idx)
  
  max_cor<-0
  max_idx<-0
  for (j in (a+b+1):(a+b+c)){
    if(max_cor<corr_mat[(a+i),j]){
      max_cor=corr_mat[(a+i),j] 
      max_idx<-j
    }
  }
  Max6<-rbind(Max6,max_cor)
  Idx6<-rbind(Idx6, max_idx)
}
Out2<-cbind(Max4+Max6, Idx4, 1:b,Idx6-(a+b))
Out2df<-data.frame(Out2, row.names = NULL)
##Sum of the maximum correlations, index for species 1, index for species 2, index for species 3: 
names(Out2df)<-c("sum", "index1", "index2", "index3")
Out2df_sort<-Out2df[order(Out2df$sum, decreasing = T),]
###Decide how many max sums to select (this is top 10):
Out2df_sort<-Out2df_sort[1:10,]

Max7<-NULL
Max8<-NULL
Idx7<-NULL
Idx8<-NULL
for (i in 1:c){ 
  max_cor<-0
  max_idx<-0
  for (j in 1:a){
    if(max_cor<corr_mat[(a+b+i),j]){
      max_cor=corr_mat[(a+b+i),j] 
      max_idx<-j
    }
  }
  Max7<-rbind(Max7,max_cor)
  Idx7<-rbind(Idx7, max_idx)
  
  max_cor<-0
  max_idx<-0
  for (j in (a+1):(a+b)){
    if(max_cor<corr_mat[(a+b+i),j]){
      max_cor=corr_mat[(a+b+i),j] 
      max_idx<-j
    }
  }
  Max8<-rbind(Max8,max_cor)
  Idx8<-rbind(Idx8, max_idx)
}
Out3<-cbind(Max7+Max8, Idx7,Idx8-a, (1):(c))
Out3df<-data.frame(Out3, row.names = NULL)
##Sum of the maximum correlations, index for species 1, index for species 2, index for species 3: 
names(Out3df)<-c("sum", "index1", "index2", "index3")
Out3df_sort<-Out3df[order(Out3df$sum, decreasing = T),]
###Decide how many max sums to select (this is top 10):
Out3df_sort<-Out3df_sort[1:10,]
