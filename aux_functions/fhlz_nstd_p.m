% FHLZ - **** **** *****
% This code is a fork of the one in the Matteo Barigozzi website
% Some minor modifications were made
% X  - data matrix
% qq - should correspond to (q+1) in the paper
% k  - number of minimum lag for the polynomial A(L)
% w  - number of lag for the covariance matrix estimation
% nlagsimp
%
%   H(L)*x(t)=R*v(t)+H(L)*xsi(t)
%   chi(t)=C(L)v(t) --> H(L)*chi(t)=R*v(t) 
%                                   if C(L) is zeroless
%           Aj(L)*chi(t)=Rj*vj(t) for j=1:m
%
%X= r(1:end-1,:);
%qq = q+1;
%w = m;
%nlagsimp = K;
%idvar=1:q;
%nrepli=30;
function [cstar,CL,v,C1,eta1,Rfinal,ufinal] = fhlz_nstd_p(X,qq,k,w,nlagsimp,idvar,nrepli)

T = size(X,1);
q = length(idvar);                                                         % number of shocks

sigma = std(X);
mu = mean(X);
z = (X - ones(T,1)*mu);
% ./(ones(T,1)*sigma);

covmat = ML_cestimate(z,q,w);                                              % cross and auto covariances of the common component
[n , o , h] = size(covmat);
H =(h+1)/2;
m = floor(n/qq);                                                            % n=(q+1)m

for h = 1:nrepli
    imp = nan*ones(n,q,nlagsimp);
    imp1 = nan*ones(n,q);
    impR = nan*ones(n,q);
    riord = randperm(n);                                                    % genero una permutazione delle variabili
    while not(isempty(intersect(idvar ,riord(qq*m+1:end))))
        riord = randperm(n);                                                % Assicura che nessuna delle variabili in id var sia esclusa nella permutazione
    end    
    x = z(:,riord);                                                         % permuto le variabili
    covmatrix = covmat(riord,riord,:);                                      % Permutation of the covariance matrix;    
    [HL GL GL1 y]=StaticRepresentation(x,covmatrix,qq,m,nlagsimp,T,k,H);        % Estimation of H(L) - G(L)=inv(H(L)) - w(t)
    S = diag(std(y)); yy = ML_center(y)*(S^-1);  %Para garantir condicoes do Theorema B in Forni (2015)
    yy1 = ML_center(y);
    eta(:,:,h)=yy1;
    opt.disp = 0; [V, MM] = eigs(cov(yy), q,'LM',opt);                      % MM = Lambda - V = P
    M = diag(sqrt(diag(MM)));                                               % M = sqrt(Lambda)                                                        % v(t)=sqrt(inv(Lambda))P'w(t)
    uu = yy*V*inv(M);
    uf(:,:,h) = uu;  
    %RR(:,:,h) = V*inv(M);  %versao anterio ao 14 de junio mudou para S*V*M
    Rf(:,:,h) = S*V*M;   
     
               
    for lag = 1 : nlagsimp; 
        BBB(:, :, lag) = GL(:, :, lag)*S*V*M;         % C(L)     
    end;

    imp(riord(1:qq*m),:,:) = BBB;                                           % Riordina le variabili
    imp1(riord(1:qq*m),:,:) = GL1*S*V*M; 
    [beta(:,:,:,h), u(:,:,h) beta1(:,:,h)] = ML_DfmCholIdent(imp,idvar,uu,imp1);
    impR(riord(1:qq*m),:,h)=Rf(:,:,h);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %P1=S*V*M; P2=nan*ones(n,q); P2(riord(1:qq*m),:)=P1; P0=P2(idvar,:);
    %Pc = chol(P0*P0')';
    %Ph = inv(P0)*Pc;
    %for lag = 1 : nlagsimp ;Pbbb(:, :, lag) = BB(:, :, lag)*S*V*M*Ph;  end
    %Pbeta = nan*ones(n,q,nlagsimp); Pbeta(riord(1:qq*s),:,:) = Pbbb;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
end

CL = nanmean(beta,4);
CL2 = nanmean(BBB,4);
v = nanmean(u,3);
ufinal = nanmean(uf,3);
C1 = nanmean(beta1,3);
eta1 = nanmean(eta,3);
Rfinal = nanmean(impR,3);


vv=[zeros(nlagsimp-1,q); v];
cstar=zeros(n,T-k);
for ii=1:T-k;
    for jj=1:nlagsimp;
        cstar(:,ii)=cstar(:,ii)+CL(:,:,jj)*vv(ii+nlagsimp-jj,:)';
    end;
end;
cstar=cstar';
end

% bbeta = permute(beta,[4 2 3 1]);
% stdbeta = permute(nanstd(bbeta),[4 2 3 1]);
% cumbeta = CumImp(meanbeta,CodeData);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% =================================================================== %%%
%%% =================================================================== %%%

% VARcov estimate a VAR from autocovariance of the variables
% C - autoregressive coefficients 
% y - residuals

function [HL GL GL1 w]=StaticRepresentation(x,covmatrix,qq,m,nlagsimp,T,k,H);

GL = zeros(qq*m,qq*m,nlagsimp); HL = zeros(qq*m,qq*m,k+1); w = zeros(T-k,qq*m);

for mm = 1:m;                
    [C CC B u]=VARcov(x,covmatrix,T,k,qq,mm,H,nlagsimp);                    % Estimation of Aj(L) - Gj(L) - vj(t)    
    HL((mm-1)*qq+1:qq*mm,(mm-1)*qq+1:qq*mm,:)=CC;                           % Building the H(L) Matrix
    w(:,(mm-1)*qq+1:qq*mm) = u;                                             % Residual w(t)   
    GL((mm-1)*qq+1:qq*mm,(mm-1)*qq+1:qq*mm,:)=B;                            % Building G(L)=inv(H(L))
    GL1((mm-1)*qq+1:qq*mm,(mm-1)*qq+1:qq*mm,:)=inv(sum(CC,3));              % long run impact matrix   
end
end
%%% =================================================================== %%%

%%% =================================================================== %%%
%%% =================================================================== %%%
 
% VARcov estimate a VAR from autocovariance of the variables
% C - autoregressive coefficients 
% B - Moving average Represetation
% u - residuals
function [C CC B u]=VARcov(x,covmatrix,T,k,qq,mm,H,nlagsimp);

A = zeros(qq*k,qq*k); B = zeros(qq*k,qq); xx = zeros(T-k,qq*k);
for j = 1:k
    for i = 1:k;    
        A((j-1)*qq+1:qq*j,(i-1)*qq+1:qq*i) = covmatrix((mm-1)*qq+1:mm*qq, (mm-1)*qq+1:mm*qq,H + i - j);   
    end    
    B((j-1)*qq+1:qq*j,:) = covmatrix((mm-1)*qq+1:mm*qq, (mm-1)*qq+1:mm*qq,H - j);     
    xx(:,(j-1)*qq+1:qq*j) =  x(k+1-j:T-j ,(mm-1)*qq+1:mm*qq);
end


C = inv(A)*B;                                                               % Aj(L)
u = x(k+1:T , (mm-1)*qq+1:mm*qq) - xx*C;                                    % Residual wj(t)
CC(:,:,1) = eye(qq); CC(:,:,2:k + 1) = - reshape(C',qq,qq,k);               % Reshaping Aj(L) s.t. it is invertible
B = InvPolMatrix(CC,nlagsimp);                                              % Building Gj(L)= inv(Aj(L))
end
%%% =================================================================== %%%

%%% =================================================================== %%%
%%% =================================================================== %%%

% inversion of a matrix of polynomials in the lag operator
function inverse = InvPolMatrix(poly,nlags); 
n = size(poly,1); k = size(poly,3) - 1;
for s = 1:k+1; newpoly(:,:,s)=inv(poly(:,:,1))*poly(:,:,s); end;            % Guarantees that poly(:,:,1)=eye(n)
polynomialmatrix = - newpoly(:,:,2:k+1);
A = zeros(n*k,n*k);
A(n+1:n*k,1:n*(k-1)) = eye(n*(k-1));
for j = 1:k
   A(1:n,(j-1)*n+1:j*n) = polynomialmatrix(:,:,j);
end
inverse = zeros(n,n,nlags);
D = eye(n*k);
for j = 1:nlags 
   inverse(:,:,j) = D(1:n,1:n)*inv(poly(:,:,1));
   D = A*D;
end
%%% =================================================================== %%%
end


