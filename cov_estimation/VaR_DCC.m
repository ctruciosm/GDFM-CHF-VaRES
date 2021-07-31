%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%% DCC Composite likelihood: GJR Student-T marginals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
clc
warning off

%addpath(genpath('/home/ctruciosm/Downloads/gdfm')) % GDFM codes from Barigozzi website
addpath(genpath('/home/ctruciosm/Downloads/mfe-toolbox-master')) % MFE Toolbox
addpath(genpath('/home/ctruciosm/GDFM-CHF-VaRES/aux_functions')) % Auxiliary functions provided in this repo
addpath(genpath('/home/ctruciosm/GDFM-CHF-VaRES/Data')) % Dataset provided in this repo
data = importdata('retornos_GDFM_CHF_VaR_APP2.txt');

[T N] = size(data)
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
Haux = H_one;
Hone = 0.5*(Haux+Haux');

for j = 1:W
    Haux = Htfull(:,:,j);
    H = 0.5*(Haux+Haux');
    rp(j,l) = weights'*chol(Hone, 'lower') * inv(chol(H, 'lower'))* datatemp(j,:)'; 
end
end
save('rp_DCCc.txt', 'rp', '-ASCII');