%  Code to obtain:
%  (1) A 95% confidence interval for the correlation between to unevenly 
%      spaced not coevally time series using the mGK method
%
%  (2) A 95 % confidence interval for the cross correlation at lag 1 between 
%     the two time series
%
%  The two time series correspond to the simulation 
%  scenario considered in the Monte Carlo experiment in the paper
%  where the true correlation is 0.2,
%  the inter-sampling-time intervals of each series have been generated from
%  a Gamma(400,1) distribution with skewness sk=0.1, the persistence of the two time 
%  series is tau=2 and the sample size  is n=200. T
%  
clear;clc;
x1=dlmread('x1.txt');
y1=dlmread('y1.txt');

tx1=dlmread('tx1.txt');
ty1=dlmread('ty1.txt');
%%

% 95% confidence interval for the correlation
[l1,u1]=bootstrap_ci_mGK(x1,y1,tx1,ty1);
int1=[l1,u1];
int1

% 95% confidence interval for the cross correlation at lag h=1
[l2,u2]= bootstrap_ci_mGK(x1,y1,tx1,ty1-1);
int2=[l2,u2];
int2
