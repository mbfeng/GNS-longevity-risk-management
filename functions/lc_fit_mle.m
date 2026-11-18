function [a,b,k] = lc_fit_mle(E,D,a,b,k)
% LC_FIT_MLE  Poisson maximum likelihood fit of the Lee-Carter model (Brouhns et al., 2002).
%   Deaths D(t,x) ~ Poisson(E(t,x)*exp(a_x + b_x*k_t)). The parameters are updated
%   one at a time by Newton steps until the log-likelihood changes by less than
%   1e-9, with the constraints sum(b) = 1 and mean(k) = 0.
%
% Inputs:
%   E, D    : years-by-ages matrices of exposures and deaths
%   a, b, k : start values (e.g. from lc_fit_svd): ages-by-1, ages-by-1, years-by-1
% Outputs:
%   a, b, k : fitted parameters

E = E';
D = D';
[numage,numyear] = size(D);

loglike = sum(sum((D.*(a*ones(1,numyear)+b*k')...
                    - E.*exp(a*ones(1,numyear)+b*k'))));
dif = 1;
while dif > 0.000000001

    for i = 1 : numage
       a(i) = a(i) + sum((D(i,:)-E(i,:).*exp(a(i)+b(i)*k'))) ...
                   ./sum((E(i,:).*exp(a(i)+b(i)*k')));
    end

    for i = 1 : numage
       b(i) = b(i) + sum((k'.*(D(i,:)-E(i,:).*exp(a(i)+b(i)*k')))) ...
                   ./sum((k'.^2.*E(i,:).*exp(a(i)+b(i)*k')));
    end

    for i = 1 : numyear
       k(i) = k(i) + sum((b.*(D(:,i)-E(:,i).*exp(a+b*k(i))))) ...
                   ./sum((b.^2.*E(:,i).*exp(a+b*k(i))));
    end

    % impose sum(b) = 1 and mean(k) = 0
    sumb = sum(b);
    meank = mean(k);
    k = sumb*(k - meank);
    a = a + b*meank;
    b = b/sumb;

    new_loglike = sum(sum((D.*(a*ones(1,numyear)+b*k')...
                            - E.*exp(a*ones(1,numyear)+b*k'))));
    dif = abs(new_loglike-loglike);
    loglike = new_loglike;
end
end
