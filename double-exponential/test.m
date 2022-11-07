%% TEST OF THE DOUBLE EXPONENTIAL FORMULA
% Example of Section 5.2 comparing the Adaptive Doube Exponential formula
% and the approach based on Padé and the new error analysis.

clear; clc; close all;
addpath('../scalingandsquaring/');

%% PARTER

A = gallery('parter',10);
lA = logm(A);

[X,made] = adelogm(A,10,1e-15,1e-9);
errabs = norm(lA - X,2); 
errrel = errabs/norm(lA);

Xfov = logmfov(A);
errabsfov = norm(lA - Xfov,2); 
errrelfov = errabsfov/norm(lA);

fprintf('%d & %1.2e & %1.2e & & & %1.2e & %1.2e \\\\\n',made,errabs,errrel,...
    errabsfov,errrelfov);



%% DORR

A = gallery('dorr',10,0.05); A = full(A);
lA = logm(A);

[X,made] = adelogm(A,10,1e-14,1e-15);
errabs = norm(lA - X,2); 
errrel = errabs/norm(lA);

Xfov = logmfov(A);
errabsfov = norm(lA - Xfov,2); 
errrelfov = errabsfov/norm(lA);

fprintf('%d & %1.2e & %1.2e & & & %1.2e & %1.2e \\\\\n',made,errabs,errrel,...
    errabsfov,errrelfov);

%% HANOWA

A = -gallery('hanowa',10);
lA = logm(A);

[X,made] = adelogm(A,10,1e-15,1e-16);
errabs = norm(lA - X,2); 
errrel = errabs/norm(lA);

Xfov = logmfov(A);
errabsfov = norm(lA - Xfov,2); 
errrelfov = errabsfov/norm(lA);

fprintf('%d & %1.2e & %1.2e & & & %1.2e & %1.2e \\\\\n',made,errabs,errrel,...
    errabsfov,errrelfov);