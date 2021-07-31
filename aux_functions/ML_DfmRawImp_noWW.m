% ML_DfmRawImp - Not Identified IRF for Factor Model
% ML_DfmRawImp(x, q, r, p, s,tb)
%   x - Data
%   q - n dynfactors
%   r - n static factors
%   p - n lags VAR Static Factors
%   s - n lags in the MA representation
%  tb - break point(optonal)
% 
% tb is 2x1 first indicate where is the break, second is which subsample
% 

% Written by Mario Forni
% Modified and commented by Matteo Luciani (matteo.luciani@ulb.ac.be)

function [imp, Chi, eta, lambda, G]  = ML_DfmRawImp_noWW(x, q, r, p, s)
%p=k; s=nlagsimp; 
[T N] = size(x);
y = ML_center(x);
[F,lambda,chi]=ML_efactors2(y,r,2);                                         % Estimating Static Factors with the normalization lambda'*lamnda/N=I

CL=zeros(N,q,s);                                                            % preallocate
Chi = chi + ones(T,1)*mean(x);                                           % Common Component

%p=BIC_function(F,T,r); 
[A u]=ML_VAR(F,p,0);                                                        % Estimating the Law of Motion fo the Static Factors   
[eta G]=ML_edynfactors2(u,q);                                               % Estimating Dynamic Factors Innovation
B = ML_MA(s,A,0);                                                           % MA Representation of the static factors
for ii=1:s
CL(:,:,ii)=lambda*B(:,:,ii)*G; % MA Representation of common components 
end     

[imp, v] = ML_DfmCholIdent(CL, 1:q, eta); % Choleski decomposition
