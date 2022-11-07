function [X,varargout] = adelogm(A,m,tol,trunc)
%DELOGM Computation of log(A) based on the DE formula.
%   INPUT: A matrix
%          m number of terms in the formula (equals number of inversions)
%          tol tolerance for the interval truncation error
%          trunc a tolerance for the trapezoidal truncation error.
%   OUTPUT: X = log(A) approximation of the matrix logarithm
% This is an implementation of Algorithm 2 from:
% Tatsuoka, Fuminori; Sogabe, Tomohiro; Miyatake, Yuto; Zhang, Shao-Liang.
% Algorithms for the computation of the matrix logarithm based on the
% double exponential formula. J. Comput. Appl. Math. 373 (2020), 112396,
% 11 pp.

if issparse(A)
    A = full(A);
end
n = size(A,1);
I = eye(n);
kmax = floor(log2(n));

FDE = @(x) cosh(x)*sech(sinh(x))^2*inv((1 + tanh(sinh(x)))*(A - I) + 2*I);

nrma    = norm(A - I,2);                % Compute ||A - I||_2
nrmainv = 1./svds(A,1,"smallest");      % Compute ||A^{-1}||_2
rhoa    = abs(eigs(A,1,"largestabs"));  % Compute rho(A)
rhoainv = abs(eigs(A,1,"smallestabs")); % Compute rho(A^{-1})

th = max(abs(log(rhoa)),abs(log(rhoainv)));
epsmax = (3*nrma*nrmainv)/(th*(1+nrmainv));
if tol >= epsmax
    tol = epsmax/2;
end
a = min(th*tol/(3*nrma),1/(2*nrma));
b = max(1 - th*tol/(3*nrma*nrmainv),2*nrmainv/(2*nrmainv+1));
l = asinh(atanh(2*a-1));
r = asinh(atanh(2*b-1));
h = (r-l)/(m-1);
Told = h/2*(FDE(l) + FDE(r)) ...
    + h*sumcell(arrayfun(FDE,l+(1:m-2)*h,'UniformOutput',false));
for k = 0:kmax
    h = h/2;
    Tnew = 0.5*Told + ...
        h*sumcell(arrayfun(FDE,l+(2*(1:m-1)-1)*h,'UniformOutput',false));
    m = 2*m - 1;
    if norm(Tnew-Told)/(3*th) <= trunc
        break
    else
        Told = Tnew;
    end
end
if nargout == 2
    varargout{1} = m;
end
X = (A-I)*Tnew;
end

function X = sumcell(C)

if isempty(C)
    X = 0;
    return
end

X = zeros(size(C{1}));
for i=1:length(C)
    X = X + C{i};
end

end