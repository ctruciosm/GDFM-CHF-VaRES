function [CL, Rt, likelihoods, Qt]= likelihood_dcc_composite(params, stdresid)

[t,K] = size(stdresid);
a = params(1);
b = params(2);

k = 2;                                                                      % based on pairs
CL = 0;
Qt = zeros(k,k,t);
Rt = zeros(k,k,t);


for index=1:K-1                                                             % total of K-1 contiguos pairs
    Qbar = cov(stdresid(:,[index index+1]));                                % select the pair
    Qt(:,:,1) = Qbar;
    Rt(1,1,1) = 1; Rt(2,2,1) = 1; 
    Rt(1,2,1) = Qt(1,2,1)/sqrt(Qt(1,1,1)*Qt(2,2,1)); 
    Rt(2,1,1) = Rt(1,2,1); 
    for j=2:t
       Qt(:,:,j) = Qbar*(1-a-b) + a*(stdresid(j-1,[index index+1])'*stdresid(j-1,[index index+1])) +  b*Qt(:,:,j-1);
       Rt(1,1,j) = 1; Rt(2,2,j) = 1; 
       Rt(1,2,j) = Qt(1,2,j)/sqrt(Qt(1,1,j)*Qt(2,2,j)); 
       Rt(2,1,j) = Rt(1,2,j);
       s11 = Rt(1,1,j); s22 = Rt(2,2,j); s12 = Rt(1,2,j);
       d = s11*s22-s12^2;
       temp = log(d)+ (s22*stdresid(j,index)^2+s11*stdresid(j,index+1)^2-2*s12*stdresid(j,index)*stdresid(j,index+1))/d;
       CL = CL+temp/(t*K); % Before: CL = CL+temp/2;
    end;
end

if isreal(CL) 
else
   disp('Imag')
   params
   CL=10E+8;
end

if isinf(CL)
   disp('Inf')
   params
   CL=10E+8;
end