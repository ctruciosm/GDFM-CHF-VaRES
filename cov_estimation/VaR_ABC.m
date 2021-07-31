%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%% ABC: DCC with GJR and idiosyncratic GJR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
clc
warning off

addpath(genpath('/Volumes/CTRUCIOS_SD/Research/VaR_GDFM_CHF/CovarianceEstimationCodes/GDFM_VaR')) % GDFM codes from Barigozzi website
addpath(genpath('/Volumes/CTRUCIOS_SD/Research/VaR_GDFM_CHF/CovarianceEstimationCodes/GDFM_VaR/MFE')) % MFE Toolbox
addpath(genpath('/Users/ctruciosm/Desktop/GDFM-CHF-VaRES/aux_functions')) % Auxiliary functions provided in this repo
addpath(genpath('/Users/ctruciosm/Desktop/GDFM-CHF-VaRES/data')) % Dataset provided in this repo
data = importdata('retornos_GDFM_CHF_VaR_APP2.txt');

[T N] = size(data);
WR = 1386;  %Out-of-Sample Period
W  = T-WR;  %In-Sample

H_day_ahead = zeros(WR,N*N);
weights(1:N,1) = 1/N;

% Number of common shocks and common factors in the first panel
AUX_data = data(1:1+W-1,:);
mu = mean(AUX_data);
datatemp = bsxfun(@minus,AUX_data,mu);

q_max = 10;
nmin = round(3*N/4);
q_aux = HL2(datatemp,q_max,2,nmin,'p1')
q = q_aux(2);
nfactors = BN_crit2(datatemp,q_max);
% end numbers of common shocks and common factors


for l = 1:WR
display(['Estimation ABC: ',num2str(l),' out of ',num2str(WR)])

AUX_data = data(l:l+W-1,:);
mu = mean(AUX_data);
datatemp = bsxfun(@minus,AUX_data,mu);

% in each replication, we need to determine k, q and nfactors
K = 10;                 % Lag B(L)
m = floor(sqrt(T));     % Lag Spectral density matrix
k = 1;                  % Lags for the VAR

[CL, chi, uest,lambda,G]  = ML_DfmRawImp_noWW(datatemp, q, nfactors, k, K);
idioest = datatemp-chi;

% uest has one less element (by the VAR(1) inside the GDFM), 
% so the first Htfull is actually the second one when estimate the whole H
[Htfull, H_one] = DCC_full(uest);
for i = 1:N
    [pari, ~, hidio(:,i)] = tarch(idioest(:,i),1,1,1,'STUDENTST');
    h_one(i) = pari(1) + pari(2)*idioest(end,i)^2 + pari(3)*idioest(end,i)^2*(idioest(end,i)<0) + pari(4)*hidio(end,i);
end

% As we do not have Hu for the common components at time 1, we start in 2.
Haux = lambda*G*H_one*G'*lambda'+diag(h_one(1,:));
Hone = 0.5*(Haux+Haux');
for j = 2:W
    Haux = lambda*G*Htfull(:,:,j-1)*G'*lambda'+diag(hidio(j,:));
    H = 0.5*(Haux+Haux');
    rp(j,l) = weights'*chol(Hone, 'lower') * inv(chol(H, 'lower'))* datatemp(j,:)'; 
end
end
save('rp_ABC.txt', 'rp', '-ASCII');