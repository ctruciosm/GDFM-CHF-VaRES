% ML_MA - Compute the Wold Representation of a VAR(p)
% Y(t) = A(L)Y(t-1) + u(t)                      VAR(p)
% Y(t) = (I + B(1) + B(2) + ... )u(t) = B1*u(t) MA(inf)
%
% [B, B1] = ML_MA(s,beta,det)
%  s = number of periods for which computing the Wold Representatio
%  beta = matrix containing the coefficients of the VAR
%  deterministic part in the VAR

% Written by Matteo Luciani - matteo.luciani@.ulb.ac.be

function  [phi, B1,A] = ML_MA(s,beta,det)
if nargin>3; error('elimina i residui'); end;
[g N]=size(beta); n=N;

%%% Retrieving the Number of Lags in the VAR %%%
if      det==0; p=(g/N);            
elseif  det==3; p=((g-2)/N); 
else            p=((g-1)/N); end

%%% These lines are necessary to organize VAR coefficients in a way %%%
%%% such that computation of IRFs is feasible                       %%%
if      det==0; temp1 = beta; 
elseif  det==3; temp1 = beta(3:g,:);
else            temp1 = beta(2:g,:); end
temp2=temp1'; A = zeros(N, N, p);
for i=1:p; k=i; for j=1:N; A(:,j,i)= temp2(:,k); k=k+p;  end; end;


%%%     Computing The Wold Representation  %%%
phi = zeros(N,n,s);
phi(:,:, 1) = eye(N,N);
for jj=2:p;
    phi(:,:,jj) = 0;
    for ii = 1:jj-1;
        temp3=A(:,:,ii)*phi(:,:,jj-ii);        
        phi(:,:,jj)=phi(:,:,jj)+temp3;
    end;
end;
for jj=p+1:s;
    phi(:,:,jj) = 0;
    for ii = 1:p;
        temp3=A(:,:,ii)*phi(:,:,jj-ii);        
        phi(:,:,jj)=phi(:,:,jj)+temp3;
    end;
end; 

%%% Computing the Long-Run Multipliers Matrix %%%
C1=eye(N); for i=1:p; C1=C1-A(:,:,i); end; B1=inv(C1); 