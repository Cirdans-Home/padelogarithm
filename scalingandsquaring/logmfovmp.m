function X = logmfovmp(A,varargin)
%LOGMFOV uses the theoretical bound to compute log(A) via inverse scaling
%and squaring and the Gauss-Legendre quadrature poles.
%   INPUT: A matrix with eigenvalues with positive real parts
%          'fov' f(theta) boundary of the FOV of A
%          'k' k pre-selected number of Gauss-Legendre poles
%          's' s pre-selected number of inverse-scaling step
% Description: If A is diagonal, or has a diagonal Schur decomposition
% we just compute the logarithm on the diagonal. It does not make sense
% doing any fancy thing. If the Schur form is actually triangular, we
% employ the bound to determine how-many triangular square root we have to
% compute, and how many poles our Gauss-Legendre formula have to use.

%% A is not square
if size(A,1) ~= size(A,2)
    error(message('MATLAB:logmfov:inputMustBeSquare'));
end
%% Parsing the inputs
checkfov = @(x) any([ishandle(x),isequal(class(x),'chebfun'),isnumeric(x)]);
p = inputParser;
p.KeepUnmatched = true;
addRequired(p,'A',@ismatrix);
addParameter(p,'fov',[],@(x) checkfov(x));
addParameter(p,'s',-1,@isnumeric);
addParameter(p,'k',-1,@isnumeric);

parse(p,A,varargin{:});

if ~isempty(fieldnames(p.Unmatched))
    disp('Extra inputs:')
    disp(p.Unmatched)
end
if ~isempty(p.UsingDefaults)
    disp('Using defaults: ')
    disp(p.UsingDefaults)
end

%% Do we have to compute a FOV?
if isempty(p.Results.fov)
    f = fov(A);
    theta = linspace(0,1,100);
    ft = f(theta);
else
    f = p.Results.fov;
    if ishandle(f)
        theta = linspace(0,1,1000);
        ft = f(theta);
    else
        ft = f;
    end
end
s = p.Results.s;
k = p.Results.k;
%% Are you giving us a sparse matrix?
if issparse(A)
    A = full(A);
end
%% Is it already triangular?
schurInput = matlab.internal.math.isschur(A);
if schurInput
    S = A;
else
    % Assume A has finite elements: we should probably check somewhere that
    % this is indeed the case
    [Q, S] = matlab.internal.math.nofinitecheck.schur(A);
end

if isdiag(S)      % Check if S is diagonal, in this case we don't have
    d = diag(S);  % much to do: just compute log on the diagonal entries
    % and bring it home.
    if any(real(d) <= 0 & imag(d) == 0)
        warning(message('MATLAB:logmfov:nonPosRealEig'))
    end
    if schurInput
        X = diag(log(d));
    else
        logd = log(d);
        X = (Q.*logd.')*Q';
        if isreal(logd)
            X = (X+X')/2;
        end
    end
else
    %% We use the information given by the bound:
    bound = @(x,s,k) 2*pi*(1+sqrt(2))*max(abs( (1 - x.^(1./(2.^(s+1))))...
        ./(1 + x.^(1./(2.^(s+1))))).^(2*k+1));
    [sval,kval] = meshgrid(0:10,1:16);
    % to compute the values of s and k to use:
    if s > -1 && k == -1
        % We know how to scale but not how many poles
        bval = arrayfun(@(k) bound(ft,s,k),1:16);
        k = find(bval < eps,1,"first");
    elseif s == -1 && k > -1
        % We know how many poles we want to use, but not how to scale
        bval = arrayfun(@(s) bound(ft,s,k),0:10);
        s = find(bval < eps,1,"first"); s = s-1;
    elseif s == -1 && k == -1
        % We know nothing
        bval = arrayfun(@(s,k) bound(ft,s,k),sval,kval);
        [k,s] = find(bval < eps); 
        s = s-1;
        [~,costind] = min(2*k/3 + 28*s/3);
        k = k(costind);
        s = s(costind);
    end
    fprintf("s = %d k = %d bound = %e cost = %d x n^3\n",...
        s,k,bound(ft,s,k),round(2*k/3 +28*s/3));
    %% Inverse scaling
    for i=1:s
        S = sqrtm_tri(S);
    end
    %% Padè based on Gauss-Legendre:
    mp.OverrideDoubleBasicArrays(true);
    N = size(S,1);
    I = eye(N,N);
    T = S - I;
    [x,omega] = mp.GaussLegendre(k);
    X = zeros(N,N);
    for i=1:k
        X = X + omega(i)*((mp(1+x(i))*T + mp(2)*I)\T);
    end
    mp.OverrideDoubleBasicArrays(false);
    %% and squaring
    X = 2^s*X;
    %% Multiply back by Q:
    if ~schurInput
        X = Q*X*Q';
    end
end

end