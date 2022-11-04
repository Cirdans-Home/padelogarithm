%% TEST FOR SCALING AND SQUARING
% Numerical examples for the selection of the values s and k

clear; clc; close all;
addpath("logmlibrary/");
addpath("../chebfun/");

for testmat = [1,5,11]
    
    if testmat == 1 || testmat == 4
        nsize = 10;
    else
        nsize = 100;
    end

    [A, set, id, nmats] = logm_testmats(testmat,nsize);
    L = logminfo(full(A));

    f = fov(A);
    ev = eig(A);
    figure(1)
    ntheta = 100;
    theta = linspace(0,2*pi,ntheta);
    plot(real(f(theta)),imag(f(theta)),'k-',...
        real(ev),imag(ev),'.',...
        "LineWidth",2);

    bound = @(x,s,k) 2*pi*(1+sqrt(2))*abs( (1 - x.^(1/(2.^(s+1))))...
        ./(1 + x.^(1/(2.^(s+1))))).^(2*k+1);

    figure(2)
    s = 0;
    k = 6;
    subplot(1,3,1)
    plot3(real(f(theta)),imag(f(theta)),bound(f(theta),s,k),'r-',...
        real(f(theta)),imag(f(theta)),zeros(ntheta,1),'k-',...
        'LineWidth',2);
    hline(0,'k--');
    vline(0,'k--');
    xlabel('Re')
    ylabel('Im')
    axis square
    title('s = 0, k = 6')
    s = 1;
    k = 6;
    subplot(1,3,2)
    plot3(real(f(theta)),imag(f(theta)),bound(f(theta),s,k),'r-',...
        real(f(theta)),imag(f(theta)),zeros(ntheta,1),'k-',...
        'LineWidth',2);
    hline(0,'k--');
    vline(0,'k--')
    axis square
    title('s = 1, k = 6')
    xlabel('Re')
    ylabel('Im')
    subplot(1,3,3)
    s = 2;
    k = 6;
    plot3(real(f(theta)),imag(f(theta)),bound(f(theta),s,k),'r-',...
        real(f(theta)),imag(f(theta)),zeros(ntheta,1),'k-',...
        'LineWidth',2);
    hline(0,'k--');
    vline(0,'k--')
    xlabel('Re')
    ylabel('Im')
    axis square
    title('s = 2, k = 6')

    %% Test of the routine
    X = logmfov(A);

    fprintf("Absolute error is: %1.2e\n",norm(X - L,2));
    fprintf("Relative error is: %1.2e\n",norm(X - L,2)./norm(L,2));

    disp('press enter to continue')
    pause()
end