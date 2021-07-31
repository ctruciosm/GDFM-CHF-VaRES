% Compute the number of static factors in a static factor model
% x=carichi*f+csi following the criterium in Bai Ng (2003)
% the static factors are estimated using simple static principal components

function [n4,n5,n6]=BN_crit2(w,rmax)
% INPUT         w           is a H*J data matrix
%               rmax        is the max number of static factors for which the
%                           criterium is computed (should be high e.g. 30)
% OUTPUT        three possibilities according to 3 different penalty functions
%               PC_1 PC_2 and PC_3 are the functions that we want to minimize
%               n1 ,n2, n3 are the number of static factors suggested by the 3 criteria

[H,J] = size(w);

G=cov(w');
[autovet,autoval]=eigs(G,rmax);
[autoval,IXX] = sort((diag(autoval)));
autoval = flipud(autoval);
IXX = flipud(IXX);
autovet = autovet(:,IXX);

for k=1:rmax
    R = autovet(:,1:k);
    F = R*sqrt(H);
    lambda = (F'*w)/H;                           % see lines 7-12 of page 9 (198) of Bai-Ng 2002
%     newF = F*((F'*F)/H)^(1/2);
%     newlambda = ((F'*F)\F'*w)'; 
    SQ=0;
%     newF=newF';
    for h=1:H
        for j=1:J
            SQ=SQ+(w(h,j)-F(h,:)*lambda(:,j))^2;
        end
    end
    V(k,1)=(SQ)/(H*J);
    IC_1(k,1)=log(V(k,1))+k*((J+H)/(J*H))*log((J*H)/(J+H));
    IC_2(k,1)=log(V(k,1))+k*((J+H)/(J*H))*log(min(sqrt(H),sqrt(J))^2);
    IC_3(k,1)=log(V(k,1))+k*((log(min(sqrt(H),sqrt(J))^2))/(min(sqrt(H),sqrt(J))^2));
    clear V SQ
end

[n4,j]=find(IC_1==min(IC_1));
[n5,j]=find(IC_2==min(IC_2));
[n6,j]=find(IC_3==min(IC_3));
