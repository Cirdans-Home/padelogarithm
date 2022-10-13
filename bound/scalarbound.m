%% Scalar bound for the Legendre based approximant

clear; clc; close all;

k = 3;
[x,omega] = legpts(k);
R = @(z,x,omega) arrayfun(@(t) sum(omega./(t.*(1+x) + 2)),z);
bound = @(z,k) arrayfun(@(t) max(2*pi*abs(((1 - sqrt(1+t))...
    ./(1+sqrt(1+t))).^(2*k+1)),eps),z);
err = @(z) abs(log(1+z) - z.*R(z,x,omega.'));
z = linspace(-1,3,100);
semilogy(z,err(z),'-',z,bound(z,k),'--','LineWidth',2)
grid on, axis tight, title(sprintf('k = %d',k))
legend({'Error','Estimate'})

for k=3:3:18
    figure(2)
    [x,omega] = legpts(k);
    R = @(z,x,omega) arrayfun(@(t) sum(omega./(t.*(1+x) + 2)),z);
    bound = @(z,k) arrayfun(@(t) max(2*pi*abs(((1 - sqrt(1+t))...
        ./(1+sqrt(1+t))).^(2*k+1)),eps),z);
    err = @(z) abs(log(1+z) - z.*R(z,x,omega.'));
    z = linspace(-1,3,100);
    subplot(2,3,k/3)
    semilogy(z,err(z),'-',z,bound(z,k),'--','LineWidth',2)
    grid on, axis tight, title(sprintf('k = %d',k))
    legend({'Error','Estimate'})
end