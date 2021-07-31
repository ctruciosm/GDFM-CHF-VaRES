%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%% GDFM-CHF: DCC with GJR and idiosyncratic GJR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
clc
warning off

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
nrepli = 30;           % number os permutations

[chi, CL, v] = fhlz_nstd_p(datatemp(1:end,:),q+1,k,m,K,1:q,nrepli);
idioest = datatemp(k+1:end,:)-chi;
[Htfull, H_one] = DCC_full(v);

for i = 1:N
    [pari, ~, hidio(:,i)] = tarch(idioest(:,i),1,1,1,'STUDENTST');
    h_one(i) = pari(1) + pari(2)*idioest(end,i)^2 + pari(3)*idioest(end,i)^2*(idioest(end,i)<0) + pari(4)*hidio(end,i);
end
Haux = CL(:,:,1)*H_one*CL(:,:,1)'+diag(h_one(1,:));
Hone = 0.5*(Haux+Haux');
for j = 2:W
    Haux = CL(:,:,1)*Htfull(:,:,j-1)*CL(:,:,1)'+diag(hidio(j-1,:)); % Htfull(:,:,j-1) because the first Htfull corresponds to the secondtime
    H = 0.5*(Haux+Haux');
    rp(j,l) = weights'*chol(Hone, 'lower') * inv(chol(H, 'lower'))* datatemp(j,:)';    
end
end
save('rp_GDFM-CHF.txt', 'rp', '-ASCII');