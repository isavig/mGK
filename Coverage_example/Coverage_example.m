%  Code to obtain
%  the empirical coverage and average length of the confidence intervals 
%  for the correlation between two time series using the mGK method for one simulation 
%  scenario considered in the Monte Carlo experiment in the paper.
%  The scenario corresponds to the case where the true correlation is 0.2,
%  the inter-sampling-time intervals of each series have been generated from
%  a Gamma(400,1)  distribution with skewness sk=0.1, the persistence of the two time 
%  series is tau=2 and the sample size  is n=200.
% The empirical coverage  (EC)  and the average length (AL) should be close 
%  to the value reported in Tables 1  and 2 of the manuscript (0.903 for
%  coverage  and 0.352 for the average length).
%  
%
%  
%
%   
%  The code produces:
% (1) A matrix M of order 1000 x 2. Each row of the matrix corresponds to a
%                               confidence intervals with nominal level 0.95 
%                               for the correlation obtained applying the
%                               GK-modified method to one of the 1000 synthetic 
% (2) The empirical coverage (EC)
% (3) The average length (AL)
clear;clc;
x=dlmread('Xserie_sk1_tau3_200_Rho1.txt');
y=dlmread('Yserie_sk1_tau3_200_Rho1.txt');

tx=dlmread('tx_sk1_tau3_200_Rho1.txt');
ty=dlmread('ty_sk1_tau3_200_Rho1.txt');

x=x(1:1000,:)';y=y(1:1000,:)';tx=tx(1:1000,:)';ty=ty(1:1000,:)';
[a,b]=size(x);
for i=1:b;
    if i==1 || rem(i,10)==0;
        [i]
    end
    [lb,ub]=bootstrap_ci_mGK(x(:,i),y(:,i),tx(:,i),ty(:,i));
    M(i,:)=[lb,ub];
end
a1=M(:,1)<=0.2
a3=M(:,2)>=0.2
EC=sum(a1.*a3)/1000
AL=mean(M(:,2)-M(:,1))
