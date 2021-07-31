% ML_DfmCholIdent - Given an MA representation compute the Structural-Choleski form
%
% [imp, u] = ML_DfmCholIdent(rawimp, idvar, rawu)
%   imp    - structural impulse responses
%   u      - structural residuals
%   rawimp - MA coefficients
%   idvar  - order of variables for choleski identification
%   rawu   - residuals
%
% By Mario Forni modified by Matteo Luciani (matteo.luciani@ulb.ac.be), 
% do not use it for VAR
%

function [imp, u, imp1, H] = ML_DfmCholIdent(rawimp, idvar, rawu,rawimp1)

B0 = rawimp(idvar,:,1);
C = chol(B0*B0')';
H = inv(B0)*C;
k = size(rawimp,3);

for j =  1:k; imp(:,:,j) = rawimp(:,:,j)*H; end;

if nargin>2; u = rawu*H; end;
if nargin>3; imp1 = rawimp1*H; end;