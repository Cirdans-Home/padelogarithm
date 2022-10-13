%% Toeplitz Matrix Example

clear; clc; close all;

a = @(t) 2.5 -t + t.^(-5);
th = linspace(-pi,pi,1000);


%% Visualize:
N = 100;
e = ones(N,1);
A = spdiags([-e,2.5*e,e],[-1,0,5],N,N);
if N <= 100
    vals = a(exp(1i*th));
    ev = eig(full(A));
    z = fov(full(A));
    figure(1)
    clf;
    hold on
    plot(real(vals),imag(vals),'r--',real(ev),imag(ev),'x');
    plot(z,'k','LineWidth',2)
    axis square
    hold off
end

%% Logarithm computation
N = 500;
e = ones(N,1);
A = spdiags([-e,2.5*e,e],[-1,0,5],N,N);
I = eye(N,N);
T = A-I;
ev = eig(T);
lmax = max(abs(a(exp(1i*th))-1));
lmin = min(abs(a(exp(1i*th))-1));
nA = norm(A,2);

bound = @(z,k) 2*pi*abs(((1 - sqrt(1+z))./(1+sqrt(1+z))).^(2*k+1));

truelog = logm(full(I+T));
for k=1:15
    [x,omega] = legpts(k);
    leglog = zeros(N,N);
    for i=1:k
        leglog = leglog + omega(i)*(((1+x(i))*T + 2*I)\I);
    end
    leglog = T*leglog;
    err(k) = norm(truelog-leglog,2);
    fprintf('k = %d Err %e\n',k,err(k));
    boundmin = bound(lmin,k);
    boundmax = bound(lmax,k);
    boundval(k) = (1+sqrt(2))*max(boundmin,boundmax);
    boundval2(k) = bound();
end

figure(2)
semilogy(1:length(err),err,'o-',1:length(boundval),boundval,'r--','LineWidth',2);
legend('Error w.r.t logm','Bound')
xlabel('k')