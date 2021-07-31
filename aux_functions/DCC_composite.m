function [Htfull, H_one]=DCC_composite(data)

[t,k]  = size(data);
htfull = zeros(t,k);

for i=1:k
    [pari, ~, httemp] = tarch(data(:,i),1,1,1,'STUDENTST');
    h_one(i) = pari(1) + pari(2)*data(end,i)^2 + pari(3)*data(end,i)^2*(data(end,i)<0) + pari(4)*httemp(end,1);
    htfull(:,i)=httemp;
end
D_one = diag(sqrt(h_one));    
    
ht = htfull;
Hstd = ht.^.5; % Transform univariate variances into standard deviations
Hstdfull = htfull.^0.5;
stdData = data./Hstd; % Obtain standardised residuals


options  =  optimset('fmincon');
options  =  optimset(options , 'Display'     , 'off');
options  =  optimset(options , 'Diagnostics' , 'off');
options  =  optimset(options , 'LargeScale'  , 'on', 'Algorithm','active-set');

dccstarting = [0.01 0.95];

% Estimate de DCC parameters
dccparameters = fmincon('likelihood_dcc_composite',dccstarting,ones(size(dccstarting)),[1-2*options.TolCon],[],[],zeros(size(dccstarting))+2*options.TolCon,[],[],options,stdData);


% Obtain the conditional variances
[Qtfull, Htfull, stdresidfull,a,b,Qbar] = dcc_eval(dccparameters, data, htfull);

Q_one = Qbar*(1-a-b) + a*stdresidfull(end,:)'*stdresidfull(end,:) + b*Qtfull(:,:,end);
q_one = sqrt(diag(Q_one));
R_one = Q_one./ (q_one*q_one');
H_one = D_one*R_one*D_one;