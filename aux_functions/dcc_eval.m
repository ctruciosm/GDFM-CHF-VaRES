function [Qt, Ht, stdresid, a, b, Qbar]=dcc_eval(params, data, H)
% PURPOSE:
%   returns the estimated conditional covariances
% 
% USAGE:
%   [Qtfull, Htfull, stdresidfull,a,b,Qbar] = dcc_eval(dccparameters, data, htfull)
% 
% INPUTS:
%   parameters  - dcc parameters
%   data        - returns
%   H           - univarite conditional variances (stage 1)
% 
% OUTPUTS:
%   Qt          - Full Qt matrices
%   Ht          - Full Ht matrices (conditional covariances)
%   stdresid
%   a,b           - DCC parameters
%   Qbar
%
% COMMENTS:
%   Function adapted from MFE Toolbox of Kevin Sheppard
% 
% Author: Carlos Trucios
% Date: 2020-10-09


[t,k] = size(data);

stdresid = data./H.^.5;
Qbar = cov(stdresid);
Hstd = H.^0.5;
stdresid = [ones(1,k);stdresid];

a = params(1);
b = params(2);

Qt = zeros(k,k,t);
Rt = zeros(k,k,t);
Ht = zeros(k,k,t);
stdresid = zeros(t,k);

Qt(:,:,1) = Qbar;
Rt(:,:,1) = Qt(:,:,1)./(sqrt(diag(Qt(:,:,1)))*sqrt(diag(Qt(:,:,1)))');
Ht(:,:,1) = diag(Hstd(1,:))*Rt(:,:,1)*diag(Hstd(1,:));
stdresid(1,:) = data(1,:)*Ht(:,:,1)^(-0.5);  

for j=2:t
   Qt(:,:,j) = Qbar*(1-a-b) + a*(stdresid(j-1,:)'*stdresid(j-1,:)) + b*Qt(:,:,j-1);
   Rt(:,:,j) = Qt(:,:,j)./(sqrt(diag(Qt(:,:,j)))*sqrt(diag(Qt(:,:,j)))');
   Ht(:,:,j) = diag(Hstd(j,:))*Rt(:,:,j)*diag(Hstd(j,:));
   stdresid(j,:) = data(j,:)*Ht(:,:,j)^(-0.5);  
end;


