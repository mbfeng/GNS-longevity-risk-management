function [A,B,K] = lc_fit_svd(E,D)
% LC_FIT_SVD  Lee-Carter estimates by singular value decomposition (start values for lc_fit_mle).
%   log m(t,x) = A(x) + B(x)*K(t), with m = deaths/exposures, sum(B) = 1 and
%   A(x) = average of log m over the years.
%
% Inputs:
%   E : years-by-ages matrix of exposures
%   D : years-by-ages matrix of deaths (same size as E)
% Outputs:
%   A : ages-by-1 vector of alpha_x
%   B : ages-by-1 vector of beta_x
%   K : years-by-1 vector of kappa_t

[m,~] = size(E);               % number of years

CDR = D ./ E;                  % central death rates
A = mean(log(CDR))';

mat = log(CDR) - repmat(A',m,1);
[K,s,B] = svds(mat,1);         % first singular vectors

temp = sum(B);                 % scale so that sum(B) = 1
B = B ./ temp;
K = K .* s .* temp;
end
