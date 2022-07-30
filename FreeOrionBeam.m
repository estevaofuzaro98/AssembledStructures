function dy = FreeOrionBeam(t,y,par)
% FREE ORION BEAM
% Estevao Fuzaro de Almeida - July, 2022

% ALIAS
m = par.m;
c = par.c;
k = par.k;
beta = par.beta;
alpha = par.alpha;

% DEFINITION
dy(1,1) = y(2);
dy(2,1) = 1/m*(- c*y(2) - k*y(1) - alpha*y(1)^3 - beta*y(2)*y(1)^2);