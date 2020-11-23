%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%% RM2006
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
clc
warning off

addpath(genpath('/home/alunos/10/ra109078/GDFM_VaR'))
%addpath(genpath('/Users/ctruciosm/Dropbox/Academico/VaR-GDFM-CHF/Codes'))


data = importdata('retornos_GDFM_CHF_VaR_APP2.txt');

[T N] = size(data);



WR = 1386;  %Out-of-Sample Period
W  = T-WR;  %In-Sample

H_day_ahead = zeros(WR,N*N);
weights(1:N,1) = 1/N;

for l = 1:WR
display(['Estimation RM2006: ',num2str(l),' out of ',num2str(WR)])

AUX_data = data(l:l+W-1,:);
OoSdata(l,:) = data(l+W-1+1,:);
mu = mean(AUX_data);
OoSmu(l,:) = mu;
datatemp = bsxfun(@minus,AUX_data,mu);
[Htfull, ~]=riskmetrics2006(datatemp);

for j = 1:W
    Haux = Htfull(:,:,j);
    H = 0.5*(Haux+Haux');
    sigma_p = sqrt(weights'*H*weights);
    e(j,l) = datatemp(j,:)*weights/sigma_p;
end

Haux = Htfull(:,:,end);
H = 0.5*(Haux+Haux');
sigma_p_day_ahead(l) = sqrt(weights'*H*weights);
H_day_ahead(l,:) = H(:);
end

save('OoSdata.txt', 'OoSdata', '-ASCII');
save('OoSmu.txt', 'OoSmu', '-ASCII');
save('H_RM2006.txt', 'H_day_ahead', '-ASCII');
save('s_RM2006.txt', 'sigma_p_day_ahead', '-ASCII');
save('e_RM2006.txt', 'e', '-ASCII');






