%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%% GDFM-CHF: DCC with GJR and idiosyncratic GJR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
clc
warning off

%addpath(genpath('/home/ctruciosm/Downloads/gdfm')) % GDFM codes from Barigozzi website
addpath(genpath('/home/ctruciosm/Downloads/mfe-toolbox-master')) % MFE Toolbox
addpath(genpath('/home/ctruciosm/GDFM-CHF-VaRES/aux_functions')) % Auxiliary functions provided in this repo
addpath(genpath('/home/ctruciosm/GDFM-CHF-VaRES/Data')) % Dataset provided in this repo
data = importdata('retornos_GDFM_CHF_VaR_APP2.txt');

[T N] = size(data);
WR = 1386;  %Out-of-Sample Period
W  = T-WR;  %In-Sample

H_day_ahead = zeros(WR,N*N);
weights(1:N,1) = 1/N;

% Number of common shocks in the first panel
AUX_data = data(1:1+W-1,:);
mu = mean(AUX_data);
datatemp = bsxfun(@minus,AUX_data,mu);

q_max = 10;
nmin = round(3*N/4);
q_aux = HL2(datatemp,q_max,2,nmin,'p1')
q = q_aux(2);
% end numbers of common shocks


for l = 1:WR
display(['Estimation GDFM-CHF: ',num2str(l),' out of ',num2str(WR)])

AUX_data = data(l:l+W-1,:);
mu = mean(AUX_data);
datatemp = bsxfun(@minus,AUX_data,mu);

% in each replication, we need to determine k, q and nfactors
K = 10;                 % Lag B(L)
m = floor(sqrt(T));     % Lag Spectral density matrix
k = 1;                  % Lags for the VAR
nrepli = 30;            % number os permutations

[chi, CL, v] = fhlz_nstd_p(datatemp(1:end,:),q+1,k,m,K,1:q,nrepli);
idioest = datatemp(k+1:end,:)-chi;
[Htfull, H_one] = DCC_full(v);

for i = 1:N
    [pari, ~, hidio(:,i)] = tarch(idioest(:,i),1,1,1,'STUDENTST');
    h_one(i) = pari(1) + pari(2)*idioest(end,i)^2 + pari(3)*idioest(end,i)^2*(idioest(end,i)<0) + pari(4)*hidio(end,i);
end
% New Bootstrap Scheme (second approach in the paper):
% The common shock ~ MGARCH: 
% u_t = Hu_t^{1/2} \eta_t  
% ==> Hu_{t}^{-1/2} u_t = \eta_t
% ==> u^{\ast}_{T+1} = Hu_{T+1}^{1/2} \eta^{\ast}_t
% The idiosincratic components ~ MGARCH 
% \epsilon_t = Hi_t^{1/2} \eta_t
% ==> Hi_t^{-1/2} \epsilon_t = \eta_t
% ==> \epsilon^{\ast}_{T+1} = Hi_{T+1}^{1/2} \eta^{\ast}_t
% \chi_{T+1} = C(L) u^{\ast}_{T+1}

for j = 2:W
    idio_one(j,:) = diag(sqrt(h_one))*diag(1./sqrt(hidio(j-1,:)))*idioest(j-1,:)';
    u_one(j,:) = chol(H_one, 'lower')*inv(chol(Htfull(:,:,j-1), 'lower'))* v(j-1,:)';
    vv=[zeros(K-1,q); v; u_one(j,:)];
    chi_one=zeros(N,1);
    ii=W-k+1;
    for jj=1:K
        chi_one(:,1)=chi_one(:,1)+CL(:,:,jj)*vv(ii+K-jj,:)';
    end;
    common_one(j,:) = chi_one';
    ret_one(j,:) = common_one(j,:) + idio_one(j,:);
    rp(j,l) = weights'*ret_one(j,:)';
end
end
save('rp_GDFM-CHF_Boot.txt', 'rp', '-ASCII');