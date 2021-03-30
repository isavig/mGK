% BOOSTRAP_CI_mGK
%
% The function returns the 95% confidence interval for correlation between two 
% time series sampled (evenly or unevenly spaced) at different time points using 
% the mGK method.
%
% 
%   Inputs:
%      x     time series values (N x 1)
%      y     time series values (M x 1)
%      tx    time points associated with x
%      ty    time points associated with y
%            
%
%   Output:
%      bci1, bci2   upper and lower limits of the 95% confidence interval 
%
% The function can also be used to obtain  95% confidence interval for the cross correlation at lag h 
% between the two time series with the mGK method. To this end it is sufficient to replace the 
% input ty for (ty-h), that is
%  
% bootstrap_ci_mGK(x,y,tx,ty) produces 95% CI for Corr[X_{t}, Y_{t}] 
% bootstrap_ci_mGK(x,y,tx,ty-h) produces 95% CI for Corr[X_{t+h}, Y_{t}] 
%
% The function is an adaptation of the "boostrap_ci_m" function in Roberts et. al 2017.
% The main difference with respect to the original "boostrap_ci_m"  is the
% subroutine that generate the time points sequence for the two time
% series that we refer to as "MODIFIED SUBROUTINE 2.2" (lines 192-413).
% A Detailed description of the MODIFIED SUBROUTINE 2.2 and limitations of the original 
% subroutine can be found in Online Resource 1 
% of the manuscript. 


% Note that in the original boostrap_ci_m the  procedure to construct the
% confidence intervals is repeated 25 times producing 25 confidence intervals. 
% The lower (upper) bound of the  confidence interval output by the function
% is the median of the lower(upper) bounds of these 25 confidence intervals.
% In the bootstrap_ci_mGK.m function the procedure to construct the
% confidence intervals is run just once.


% The function evaluate a 95% confidence interval using B=2000 bootstrap
% replications. It is straightfowrard to modify the function to produce 
% confidence intervals with arbitrary confidence level (1-alpha) using 
% a arbitray  number of bootstrap
% replications B. 




function [lb,ub,taux,tauy]=bootstrap_ci_mGK(x,y,tx,ty)
max_iter = 2000;
  n = length(x);
  m = length(y);
  
  corr_bootstrap=nan(max_iter,1);
  
  xmean = nanmean(x);
  txmean = (tx(end)-tx(1)) / (n-1);
  ymean = nanmean(y);
  tymean = (ty(end)-ty(1)) / (m-1);
  delta_t = max([txmean,tymean]);

 
  temp = zeros(m+n,1);
  indx = zeros(m+n,2);
  seq = zeros(m+n+1,1);
  tstar = zeros(m+n,1);
  inds = zeros(m+n,1);
  nx = zeros(m+n,1);
  ny = zeros(m+n,1);
  orig = correlate_gaussian(x, y, tx, ty);
    % estimate persistence times
  txmean = (tx(n)-tx(1)) / (n-1);
  tymean = (ty(m)-ty(1)) / (m-1);
  tmean = max([txmean,tymean]);

  for i = 1:n/16
    tmp = tx(1:n) + i*txmean;
    taux = correlate_gaussian(x, x, tx, tmp(1:n));
    if abs(taux) < 0.368
      break
    end
  end
  taux = i*txmean;

  for j = 1:m/16
    tmp = ty(1:m) + j*tymean;
    tauy = correlate_gaussian(y, y, ty, tmp(1:m));
    if abs(tauy) < 0.368
      break
    end
  end
  tauy = j*tymean;

  p = 1.0 - tmean / (4*max([taux,tauy]));
  logp = log(p);
  % build integrated time-list
  temp = vertcat(tx(1:n),ty(1:m));

  % index elements for series x are positive integers, and negative for y
  tmp_indx = horzcat([1:n,-1:-1:-m]);

  [temp, tmp_indx] = joint_sort(temp, tmp_indx);
  j = 1;
  for i = 1:n+m
    tstar(j) = temp(i);
    if tmp_indx(i) > 0
      indx(j,1) = tmp_indx(i);
    else
      indx(j,2) = tmp_indx(i);
    end
    if i < n+m
      if temp(i) ~= temp(i+1)
	j = j + 1;
      end
    end
  end
  nm = j;
 % do the bootstraps
  n_low = 0;
  for iter = 1:max_iter
      iter;
      num_block=1;
      auxl=0;
      longitud_bloques=[0];
    
      j = 1;
      init = [0 0];
      rnd = rand;
      % random starting point
      i = floor((nm-1)*rnd) + 1;
      if indx(i,1) ~= 0
        seq(j) = indx(i,1);
	    j = j + 1;
        auxl=auxl+1;
        init(1) = 1;
      end
      if indx(i,2) ~= 0
	     seq(j) = indx(i,2);
	     j = j + 1;
         auxl=auxl+1;
         init(2) = 1;
      end
      
      while j <= n+m
        rnd = rand;
	    if i > 1
	       dt = tstar(i) - tstar(i-1);
        else
	       dt = tstar(2) - tstar(1);
        end  
	    if rnd > exp(dt/tmean*logp) 
           longitud_bloques(num_block)=auxl;
           num_block=num_block+1;
           auxl=0;
	       if ~all(init) 
	         j = 1;
	         init = [0 0];
             num_block=1;
             auxl=0;
             longitud_bloques=[0];
           end
           rnd = rand;
           i = floor((nm-1)*rnd) + 1;
           if indx(i,1) ~= 0
             seq(j) = indx(i,1);
             j = j + 1;
             auxl=auxl+1;
           end
           if indx(i,2) ~= 0
	         seq(j) = indx(i,2);
	         j = j + 1;
             auxl=auxl+1;
           end
        else 
	       i = i + 1;
	       if i > nm
	         i = 1;
           end
	       if indx(i,1) ~= 0
	         seq(j) = indx(i,1);
	         j = j + 1;
	         init(1) = 1;
             auxl=auxl+1;
           end
	       if indx(i,2) ~= 0
	         seq(j) = indx(i,2);
	         j = j + 1;
	         init(2) = 1;
             auxl=auxl+1;
           end
	    end
      end
       
    longitud_bloques(num_block)=auxl;
       
    % Modified subroutine 2.2 starts here
 
    q=length(longitud_bloques);
    block_seq=cell(1,q);
    block_seq{1}=seq(1:longitud_bloques(1));
    for i=2:q
       block_seq{i}=seq(sum(longitud_bloques(1:i-1))+1:sum(longitud_bloques(1:i)));
    end
   
    i=1;
    if tx(end)>ty(end);
      while i<=length(block_seq);
        if any(block_seq{i}==n) & any(block_seq{i}==-1);
           k=find(block_seq{i}==n);
           kk=find(block_seq{i}==-1);
           if k<kk;
             aux1=block_seq{i}(1:k);
             aux2=block_seq{i}(k+1:end);
             for j=length(block_seq):-1:i;
               block_seq{j+1}=block_seq{j};
             end
             block_seq{i}=aux1;
             block_seq{i+1}=aux2;
           end
         end
      i=i+1;
      end
    end

   i=1;
   if tx(end)<ty(end);
      while i<=length(block_seq);
        if any(block_seq{i}==-m) & any(block_seq{i}==1);
           k=find(block_seq{i}==-m);
           kk=find(block_seq{i}==1);
           if k<kk;
              aux1=block_seq{i}(1:k);
              aux2=block_seq{i}(k+1:end);

              for j=length(block_seq):-1:i;
                  block_seq{j+1}=block_seq{j};
              end
              block_seq{i}=aux1;
              block_seq{i+1}=aux2;
           end
        end
        i=i+1;
      end
   end

   i=1;
   if tx(end)==ty(end);
      while i<=length(block_seq);
         if (any(block_seq{i}==n) & any(block_seq{i}==-1)) | (any(block_seq{i}==-m) & any(block_seq{i}==1)) ;
            k1=find(block_seq{i}==n);
            k11=find(block_seq{i}==-1);    
            k2=find(block_seq{i}==-m);
            k22=find(block_seq{i}==1);
            if k1 > k2 & k1<k11;
               aux1=block_seq{i}(1:k1);
               aux2=block_seq{i}(k1+1:end);
               for j=length(block_seq):-1:i;
                    block_seq{j+1}=block_seq{j};
               end
               block_seq{i}=aux1;
               block_seq{i+1}=aux2;
               elseif k1 < k2 & k2<k22;
               aux1=block_seq{i}(1:k2);
               aux2=block_seq{i}(k2+1:end);

               for j=length(block_seq):-1:i;
                  block_seq{j+1}=block_seq{j};
               end
               block_seq{i}=aux1;
               block_seq{i+1}=aux2;
            end
         end
         i=i+1;
      end
   end

   q=length(block_seq);
   block_seq_x=cell(1,q);
   for i=1:q;
      block_seq_x{i}=block_seq{i}(block_seq{i}>0);
   end           
   block_seq_y=cell(1,q);
   for i=1:q;
      block_seq_y{i}=block_seq{i}(block_seq{i}<0);
   end          

   diff_x=cell(1,q);
   for i=1:q
     diff_x{i}=diff(tx(block_seq_x{i}));
     diff_x{i}(diff_x{i}<0)=delta_t;
   end
 
   diff_y=cell(1,q);
   aux=ty(2)-ty(1);
   for i=1:q
      diff_y{i}=diff(ty(-block_seq_y{i}));
      diff_y{i}(diff_y{i}<0)=delta_t;
   end           

   t_tilda=cell(q,2); 
   n_tilda=cell(q,2);

   if ~isempty(block_seq_x{1}) & ~isempty(block_seq_y{1}) 
     t_tilda{1,1}=tx(block_seq_x{1});
     t_tilda{1,2}=ty(-block_seq_y{1});
     elseif isempty(block_seq_x{1})
       t_tilda{1,1}=nan(length(block_seq_y{1}),1);
       t_tilda{1,2}=ty(-block_seq_y{1});
     elseif isempty(block_seq_y{1}) 
       t_tilda{1,1}=tx(block_seq_x{1});
       t_tilda{1,2}=nan(length(block_seq_x{1}),1);
   end

   if ~isempty(block_seq_x{1}) & ~isempty(block_seq_y{1}) 
     n_tilda{1,1}=x(block_seq_x{1});
     n_tilda{1,2}=y(-block_seq_y{1});
     elseif isempty(block_seq_x{1}) 
       n_tilda{1,1}=nan(length(block_seq_y{1}),1);
       n_tilda{1,2}=y(-block_seq_y{1});
     elseif isempty(block_seq_y{1}) 
       n_tilda{1,1}=x(block_seq_x{1});
       n_tilda{1,2}=nan(length(block_seq_x{1}),1);
   end

   ntx=t_tilda{1,1};
   nty=t_tilda{1,2};
   nx=n_tilda{1,1};
   ny=n_tilda{1,2};

   for i=2:q;
      if max(ntx)>=max(nty) || all(isnan(nty))  
         if ~isempty(block_seq_x{i}) & ~isempty(block_seq_y{i}) 
            if block_seq_x{i}(1)==1;
               dt=delta_t;
            else
               dt=tx(block_seq_x{i}(1))-tx(block_seq_x{i}(1)-1);
            end
            a=dt+max(ntx);
            t_tilda{i,1}=[a;a+cumsum(diff_x{i})];
            b=t_tilda{i,1}(1)+ty(-block_seq_y{i}(1))-tx(block_seq_x{i}(1));
            t_tilda{i,2}=[b;b+cumsum(diff_y{i})];
            n_tilda{i,1}=x(block_seq_x{i});
            n_tilda{i,2}=y(-block_seq_y{i});
         elseif isempty(block_seq_x{i}) 
           t_tilda{i,1}=nan(length(block_seq_y{i}),1);  
           if -block_seq_y{i}(1)==1;
              dt=delta_t;
           else
               dt=ty(-block_seq_y{i}(1))-ty(-block_seq_y{i}(1)-1);
           end
           
           b=dt+max(ntx);
           t_tilda{i,2}=[b;b+cumsum(diff_y{i})];
           
           n_tilda{i,2}=y(-block_seq_y{i});
           n_tilda{i,1}=nan(length(block_seq_y{i}),1);
         elseif isempty(block_seq_y{i}) 
           t_tilda{i,2}=nan(length(block_seq_x{i}),1);
           if block_seq_x{i}(1)==1;
              dt=delta_t;
           else
              dt=tx(block_seq_x{i}(1))-tx(block_seq_x{i}(1)-1);
           end
           a=dt+max(ntx);
           t_tilda{i,1}=[a;a+cumsum(diff_x{i})];
           n_tilda{i,1}=x(block_seq_x{i});
           n_tilda{i,2}=nan(length(block_seq_x{i}),1);
         end
         
      else 
        if ~isempty(block_seq_x{i}) & ~isempty(block_seq_y{i}) 
           if -block_seq_y{i}(1)==1;
              dt=delta_t;
           else
              dt=ty(-block_seq_y{i}(1))-ty(-block_seq_y{i}(1)-1);
           end
           a=dt+max(nty);
           t_tilda{i,2}=[a;a+cumsum(diff_y{i})];
           b=t_tilda{i,2}(1)+tx(block_seq_x{i}(1))-ty(-block_seq_y{i}(1));
           t_tilda{i,1}=[b;b+cumsum(diff_x{i})];
           n_tilda{i,1}=x(block_seq_x{i});
           n_tilda{i,2}=y(-block_seq_y{i});
        elseif isempty(block_seq_x{i})
           t_tilda{i,1}=nan(length(block_seq_y{i}),1);    
           if -block_seq_y{i}(1)==1;
              dt=delta_t;
           else
              dt=ty(-block_seq_y{i}(1))-ty(-block_seq_y{i}(1)-1);
           end
           a=dt+max(nty);
           t_tilda{i,2}=[a;a+cumsum(diff_y{i})];   
           n_tilda{i,2}=y(-block_seq_y{i});
           n_tilda{i,1}=nan(length(block_seq_y{i}),1);
       elseif isempty(block_seq_y{i})
           t_tilda{i,2}=nan(length(block_seq_x{i}),1);           
           if block_seq_x{i}(1)==1;
              dt=delta_t;
           else
              dt=tx(block_seq_x{i}(1))-tx(block_seq_x{i}(1)-1);
           end
           b=dt+max(nty);
           t_tilda{i,1}=[b;b+cumsum(diff_x{i})]; 
           n_tilda{i,1}=x(block_seq_x{i});
           n_tilda{i,2}=nan(length(block_seq_x{i}),1);
       end
      end
    ntx=[ntx;t_tilda{i,1}];
    nty=[nty;t_tilda{i,2}];
    nx=[nx;n_tilda{i,1}];
    ny=[ny;n_tilda{i,2}];
  end
  ntx=ntx(~isnan(ntx));
  nty=nty(~isnan(nty));
  nx=nx(~isnan(nx));
  ny=ny(~isnan(ny));
  
% Modified subroutine 2.2 ends here
 
  bse(iter) = correlate_gaussian(nx, ny, ntx, nty);

  if bse(iter) < orig
        n_low = n_low + 1;
  end
      
 end

 % call shell(bse(1:max_iter),max_iter)
 % calculate bias parameter
 % call PNI(real(n_low)/max_iter,1.0-real(n_low)/max_iter,real(n_low)/max_iter-0.5,z0,ierr)
 z0 = icpdf(real(n_low)/max_iter);

 % calculate acceleration parameter
 asum = sum(bse);
 astar = (sum(asum-bse)/(max_iter-1)) / max_iter;
 anum = sum((astar-(asum-bse)/(max_iter-1)).^3);
 aden = sum((astar-(asum-bse)/(max_iter-1)).^2);
 a = anum / (6.0*exp(1.5*log(aden)));

 % calculate interval points
 b0 = z0 + (z0-1.960) / (1-a*(z0-1.960));
 b1 = z0 + (z0+1.960) / (1-a*(z0+1.960));

 ci1 = nanmin([nanmax([floor(max_iter*cpdf(b0)),1]),max_iter]);
 ci2 = nanmax([nanmin([ceil(max_iter*cpdf(b1)),max_iter]),1]);

 bse = sort(bse);
 lb = bse(ci1);
 ub = bse(ci2);

end