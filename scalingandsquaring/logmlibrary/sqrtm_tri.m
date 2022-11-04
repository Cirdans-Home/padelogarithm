function R = sqrtm_tri(T)
%SQRTM_TRI Square root of quasi-upper triangular matrix.
%   Works entirely in real arithmetic when possible.

%   Nicholas J. Higham and Samuel D. Relton
%   Copyright 1984-2018 The MathWorks, Inc.
n = size(T,1);
switch n
    case 1 % 1x1 block
        R = sqrt(T);
    case 2 % 2x2 block
        R = sqrtm_tbt(T);
    otherwise % Larger block - divide and conquer.
        n2 = floor(n/2);
        if T(n2+1,n2) ~= 0
            n2 = n2+1;
        end
        T11 = T(1:n2,1:n2);
        T22 = T(n2+1:end, n2+1:end);
        T12 = T(1:n2, n2+1:end);
        
        R11 = sqrtm_tri(T11);
        R22 = sqrtm_tri(T22);
        R12 = matlab.internal.math.sylvester_tri(R11, R22, T12, 'I', 'I', 'notransp');

        R(n2+1:n, n2+1:n) = R22;
        R(1:n2,1:n2) = R11;
        R(1:n2, n2+1:n) = R12;
end
