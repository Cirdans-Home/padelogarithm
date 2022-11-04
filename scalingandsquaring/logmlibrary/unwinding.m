function u = unwinding(z)
%UNWINDING Unwinding number of z.
u = ceil( (imag(z) - pi)/(2.*pi) );
