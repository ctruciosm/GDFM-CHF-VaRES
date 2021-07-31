%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%% RM2006
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

for l = 1:WR
display(['Estimation RM2006: ',num2str(l),' out of ',num2str(WR)])

AUX_data = data(l:l+W-1,:);
OoSdata(l,:) = data(l+W-1+1,:);
mu = mean(AUX_data);
OoSmu(l,:) = mu;
datatemp = bsxfun(@minus,AUX_data,mu);
[Htfull, ~]=riskmetrics2006_fore(datatemp);
Haux = Htfull(:,:,end);
Hone = 0.5*(Haux+Haux');
CholHone = chol(Hone, 'lower');
for j = 1:W
    Haux = Htfull(:,:,j);
    H = 0.5*(Haux+Haux');
    rp(j,l) = weights'*CholHone* inv(chol(H, 'lower'))* datatemp(j,:)'; 
end
end
save('OoSdata.txt', 'OoSdata', '-ASCII');
save('OoSmu.txt', 'OoSmu', '-ASCII');
save('rp_RM2006.txt', 'rp', '-ASCII');