%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%% GDFM-CHF: DCC with GJR and idiosyncratic GJR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
clc
warning off

addpath(genpath('/home/alunos/10/ra109078/GDFM_VaR'))
%addpath(genpath('/Volumes/CTRUCIOS_SD/VaR_GDFM_CHF/CovarianceEstimationCodes'))


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

% uest has one less element (by the VAR(1) inside the GDFM), 
% so the first Htfull is actually the second one when estimate the whole H
[Htfull, H_one] = DCC_full(v);

for i = 1:N
    [pari, ~, hidio(:,i)] = tarch(idioest(:,i),1,1,1,'STUDENTST');
    h_one(i) = pari(1) + pari(2)*idioest(end,i)^2 + pari(3)*idioest(end,i)^2*(idioest(end,i)<0) + pari(4)*hidio(end,i);
end
% As we do not have H for the common components at time 1, we start in 2.
% Htfull(:,:,j-1) because the first Htfull corresponds to the secondtime
% u_t

% New Two lines: 23-02-2021
Haux = CL(:,:,1)*H_one*CL(:,:,1)'+diag(h_one(1,:));
Hone = 0.5*(Haux+Haux');

for j = 2:W
    Haux = CL(:,:,1)*Htfull(:,:,j-1)*CL(:,:,1)'+diag(hidio(j-1,:));
    H = 0.5*(Haux+Haux');
    sigma_p = sqrt(weights'*H*weights);
    e(j,l) = datatemp(j,:)*weights/sigma_p;
    % new line 23-02-2021
    rp(j,l) = weights'*chol(Hone, 'lower') * inv(chol(H, 'lower'))* datatemp(j,:)';    
end


%Haux = CL(:,:,1)*H_one*CL(:,:,1)'+diag(h_one(1,:));
%H = 0.5*(Haux+Haux');
%sigma_p_day_ahead(l) = sqrt(weights'*H*weights);
%H_day_ahead(l,:) = H(:);
sigma_p_day_ahead(l) = sqrt(weights'*Hone*weights);
end



%save('H_GDFM-CHF.txt', 'H_day_ahead', '-ASCII');
save('rp_GDFM-CHF.txt', 'rp', '-ASCII');
save('s_GDFM-CHF.txt', 'sigma_p_day_ahead', '-ASCII');
save('e_GDFM-CHF.txt', 'e', '-ASCII');






