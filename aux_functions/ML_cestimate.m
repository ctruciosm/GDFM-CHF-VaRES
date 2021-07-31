% ML_cestimate - computes the cross and auto covariances of the common component
%
% S = ML_cestimate(x,nfactors,w)
%   x        - T by N data matrix
%   nfactors - number of dynamic Factor
%   w        - number of lagged covariances used in the Bartlett lag-window estimation of the spectral matrix
%
% Written by Mario Forni modified by Matteo Luciani (mluciani@ulb.ac.be)
% See: Forni, Hallin, Lippi and Reichlin 2000 "The Generalized Factor Model: Identification and Estimation", 
%             The Review of Economics and Statistics 82, 540-54

function S = ML_cestimate(x,nfactors,w)

%%% define some useful quantities %%%
[T,N] = size(x);
W = 2*w+1;
B = triang(W);

%%% compute covariances %%%
S = zeros(N,N,W);
for k = 1:w+1,
     S(:,:,w+k) = B(w+k)*ML_center(x(k:T,:))'*ML_center(x(1:T+1-k,:))/(T-k);
     S(:,:,w-k+2) = S(:,:,w+k)';
end

%%% compute the spectral matrix in W points (S) %%%
Factor = exp(-sqrt(-1)*(-w:w)'*(0:2*pi/W:4*pi*w/W));
for j = 1:N
   S(j,:,:) = squeeze(S(j,:,:))*Factor;
end

%%% compute the egenvectors  for all points (E) %%%
opt.disp = 0;
[A,D] = eigs(S(:,:,1),nfactors,'LM',opt);
S(:,:,1) = A*D*A'; 
for j = 2:w+1,
   [A,D] = eigs(S(:,:,j),nfactors,'LM',opt);
   S(:,:,j) = A*D*A'; 
   S(:,:,W+2-j) = conj(S(:,:,j));
 end
 for j = 1:N
    S(:,j,:) = real(squeeze(S(:,j,:))*conj(Factor).'/W);
end
