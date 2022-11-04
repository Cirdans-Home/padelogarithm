function X = make_real(X)
%MAKE_REAL Remove imaginary part comprising roundoff in real case.
if norm(imag(X),1) <= 10*size(X,1)*eps*norm(X,1)
   X = real(X);
end