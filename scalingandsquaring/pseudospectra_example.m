%% Example of the bound with the Pseudospectra
% This codes needs both the EIGTOOL and the ADVANTPIX packages to be run.
% The EIGTOOL is used to obtain an estimate of the pseudospectra, while
% ADVANTPIX is used to overcome the limitation caused by the condition
% number of the test matrices.

clear; clc; close all;
mp.OverrideDoubleBasicArrays(false);
mp.ExtendConstAccuracy(false);
%% Test matrices
% These are test matrices from the paper:
% Cardoso, João R.; Silva Leite, F.
% Padé and Gregory error estimates for the logarithm of block triangular matrices.
% Appl. Numer. Math.56(2006), no.2, 253–267.
T = @(a,b,c) exp(a)*[1 b b^2/2+c; 0 1 b; 0 0 1];
logT = @(a,b,c) [a b c; 0 a b; 0 0 a];

a = 0.5; % Change it to 0.05, 0.1, 0.3, 0.5
b = 1e3;
c = 1e-3;
%% Pseudospectra
% We use here the `eigtool` software to evaluate the epsilon pseudospectra
% of the matrix we are interested in.
opts.npts=20;
opts.levels=-10:1:-3;
opts.ax=[-2 3 -1.5 1.5];
opts.proj_lev=Inf;
opts.colour=1;
opts.thick_lines=1;
opts.scale_equal=1;
opts.grid=1;
opts.dim=1;
opts.no_ews=0;
opts.no_psa=0;
opts.fov=0;
opts.colourbar=1;
opts.imag_axis=1;
opts.unit_circle=0;
opts.direct=1;
opts.print_plot=0;
[X,Y,SIGS] = eigtool(T(a,b,c),opts);
M = contour(X,Y,log10(SIGS));
epsilon = 1e-6;
ind = find(M(1,:) == log10(epsilon));
numvert = M(2,ind);
xval = M(1,ind+1:ind+numvert);
yval = M(2,ind+1:ind+numvert);
le = chebfun(xval.'+1i*yval.');
%% Evaluate bound
lengeps = sum(arcLength(le));
bound = @(x,s,k) (lengeps/epsilon)*max(abs( (1 - x.^(1./(2.^(s+1))))...
    ./(1 + x.^(1./(2.^(s+1))))).^(2*k+1),[],"all");
[sval,kval] = meshgrid(0:10,1:16);
ft = xval.'+1i*yval.';
for i = 0:10
    for j = 1:16
        bval(i+1,j) = bound(ft,i,j);
    end
end
[k,s] = find(bval < eps);
s = s-1;
[~,costind] = min(2*k/3 + 28*s/3);
k = k(costind);
s = s(costind);
%% Visualize the contour
handlefig = figure(2);
hold on
contour(X,Y,log10(SIGS),'ShowText','on','LineWidth',2)
plot(le,'k-','Linewidth',2)
plot(ft,'or','MarkerSize',5)
hold off
axis square
filename = char(sprintf("pseudoa%03d.tex",round(100*a)));
try
    matlab2tikz('filename',filename, ...
        'figurehandle',handlefig,'externalData',true);
catch
    warning("Plots for the paper are made using matlab2tikz")
end
%% Compute the logarithm
mp.Digits(34) % Set quadruple precision
a = mp(0.5);
b = mp(1e3);
c = mp(1e-3);
X = logmfovmp(mp(T(a,b,c)),'k',k,'s',s);
%% Evaluate the error
fprintf("Error is %1.2e bound is %1.2e\n",norm(X-logT(a,b,c),2),...
    bval(s+1,k));

%% Why the bound doesn't work?
f = @(x,s,k) (lengeps/epsilon)*abs( (1 - x.^(1./(2.^(s+1))))...
    ./(1 + x.^(1./(2.^(s+1))))).^(2*k+1);
figure(3)
subplot(1,2,1)
x = linspace(0,3,100);
y = linspace(-1.5,1.5,100);
xcheb = linspace(-1,1,100);
[Xg,Yg] = meshgrid(x,y);
Z = Xg + 1i*Yg;
mesh(Xg,Yg,log10(f(Z,s,k)));
subplot(1,2,2)
hold on
plot(le,'k-','Linewidth',2)
plot3(real(le(xcheb)),imag(le(xcheb)),log10(f(le(xcheb),s,k)));
view(20,30)
set(gca, 'Zdir', 'reverse')
hold off

