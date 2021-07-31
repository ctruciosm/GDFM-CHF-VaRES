% ML_edynfactors2 - Estimation of Dynamic Factors Innovation 
% 
% Principal Components of the Residuals of a VAR(p) estimated on r Static Factors
%
% [etahat, G]=ML_edynfactors2(y,q)
% The Model:
%       X(t) = lambda*F(t) + xsi(t)
%       F(t) = A(L)*F(t-1) + epsilon
%       epsilon = G*eta
%       where G is r by q
% Outputs:
%       etahat = estimates of dynamic Factors Innovations
%       G = estimates of G
% Inputs:
%       y = epsilon
%       q = number of dynamic factors
%

% Written by Matteo Luciani (matteo.luciani@ulb.ac.be)

function [eta, G]=ML_edynfactors2(y,q)
opt.disp=0;
N=size(y,2);
sigma=cov(y);                                                               % Variance Covariance Matrix of VAR Residuals
if q<N; [P,M]=eigs(sigma,q,'LM',opt); else [P,M]=eig(sigma); end            % Eigenvalue eigenvectors decomposition
eta=y*P*(M^-.5);                                                            % Dynamic Factor Innovations            
G=P(:,1:q)*(M^.5);


