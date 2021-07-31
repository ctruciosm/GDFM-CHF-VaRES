% ML_ols - OLS Estimation
% 
% [beta,u,v,esu,r2,espqr]=ML_ols(y,x,det);
%   Inputs:
%       y   = Endogenous Variable (vector)
%       x   = Exogenous Variables (matrix)
%       det = 0 noconstant
%       det = 1 constant
%       det = 2 time trend
%       det = 3 constant + time trend
%   Outputs:
%       beta  = estimated coefficient
%       u     = residuals
%       v     = varcov matrix of the estimates
%       esu   = Residual Variance
%       r2    = R-Squared
%       espqr = estimates standard errors
%

% Written by Matteo Luciani
% This is a modified version of the codes available on Fabio Canova webpage

function [beta,u,v,esu,r2,espar]=ML_ols(y,x,det)
[T k] = size(x);

cons=ones(T,1); trend=(1:1:T)';
if      det==1; x=[cons x];
elseif  det==2; x=[trend x];
elseif  det==3; x=[cons trend x];
end;

xx=inv(x'*x);
beta=xx*x'*y;                   % ols coeff
u=y-x*beta;                     % residuals
uu=u'*u;                        % SSR
esu=uu/(T-k);                   % Residual Variance
r2=1-(uu)/(y'*y);               % R2
% r2c=1-(1-r2)*(T-1)/(T-k);       % adjusted R2
v=esu*xx;                       % varcov matrix of the estimates
espar=sqrt(diag(v));            % standard errors


