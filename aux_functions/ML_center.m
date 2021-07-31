% ML_center - Demean variables
% CENTER XC = center(X)
%	Centers each column of X.

%	J. Rodrigues 26/IV/97, jrodrig@ulb.ac.be
function XC = ML_center(X)
[T n] = size(X);
XC = X - ones(T,1)*(sum(X)/T); % Much faster than MEAN with a FOR loop