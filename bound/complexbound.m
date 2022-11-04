%% Figures for the complex bound
% Level sets in log10(.) scale of the absolute value of the bound for
% different values of k.

bound = @(z,k) arrayfun(@(t) 2*pi*abs(((1 - sqrt(1+t))...
    ./(1+sqrt(1+t))).^(2*k+1)),z);
x = linspace(-1,10);
y = linspace(-5,6);
[X,Y] = meshgrid(x,y);

for k=3:3:18
    figure(1)
    subplot(2,3,k/3)
    contourf(X,Y,log10(bound(X+1i*Y,k)))
    colorbar
    colormap("gray")
    title(sprintf('k = %d',k))
    axis square
end