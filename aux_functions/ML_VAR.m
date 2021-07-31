% ML_VAR - Estimates a VAR(k) for a vector of variables y
% if y is a single variable it estimates an AR model with OLS
%
% [A,u]=ML_VAR(y,k,jj);
%   y  = vector of endogenous variables
%   k  = number of lags
%   jj = 0 noconstant
%   jj = 1 constant
%   jj = 2 time trend
%   jj = 3 constant + time trend
%   A  = Matrix of coefficients for the reduced form 
%   u  = Vector of Residuals
% ----------------------------------------------------------

% written by Matteo Luciani (matteo.luciani@ulb.ac.be)

function [A,u]=ML_VAR(y,k,jj)
[T N] = size(y);

%%% Building Up the vector for OLS regression %%%
for i=1:N,
 yy(:,i)=y(k+1:T,i);
 for j=1:k,
  xx(:,k*(i-1)+j)=y(k+1-j:T-j,i);
 end; 
end;


%%% OLS Equation-By-Equation %%%
for i=1:N;
    [A(:,i),u(:,i)]=ML_ols(yy(:,i),xx,jj);
end;
