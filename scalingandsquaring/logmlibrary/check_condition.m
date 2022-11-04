function c = check_condition(z)
%CHECK_CONDITION Check whether z is in safe range.
tol = eps^(1/8);
rtol = 1/tol;
c = abs(z) > rtol || abs(z - 1) < tol || abs(z + 1) < tol;
