%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%% DCC Composite likelihood: GJR Student-T marginals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
clc
warning off

%addpath(genpath('/home/alunos/10/ra109078/GDFM_VaR'))
addpath(genpath('/Users/ctruciosm/Dropbox/Academico/VaR-GDFM-CHF/Codes'))


data = importdata('retornos_GDFM_CHF_VaR_APP2.txt');

[T N] = size(data);

WR = 1386;  %Out-of-Sample Period
W  = T-WR;  %In-Sample

H_day_ahead = zeros(WR,N*N);
weights(1:N,1) = 1/N;

for l = 1:WR
display(['Estimation DCCc: ',num2str(l),' out of ',num2str(WR)])

AUX_data = data(l:l+W-1,:);
mu = mean(AUX_data);
datatemp = bsxfun(@minus,AUX_data,mu);

[Htfull, H_one] = DCC_composite(datatemp);

for j = 1:W
    Haux = Htfull(:,:,j);
    H = 0.5*(Haux+Haux');
    sigma_p = sqrt(weights'*H*weights);
    e(j,l) = datatemp(j,:)*weights/sigma_p;
end

Haux = H_one;
H = 0.5*(Haux+Haux');
sigma_p_day_ahead(l) = sqrt(weights'*H*weights);
H_day_ahead(l,:) = H(:);
end

save('H_DCCc.txt', 'H_day_ahead', '-ASCII');
save('s_DCCc.txt', 'sigma_p_day_ahead', '-ASCII');
save('e_DCCc.txt', 'e', '-ASCII');







