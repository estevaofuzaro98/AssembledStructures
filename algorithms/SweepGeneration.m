function [time,signal] = SweepGeneration(fmin,fmax,SweepTime)
%% Signal processing settings
fs = 10*fmax; % Hz
dt = 1/fs; % s
T = SweepTime; % s

time = 0:dt:T; % s

%% Sweep sine settings
S = 1;
rate = (fmax-fmin)/T;
phi = 2 * pi * (fmin*time + rate*time.^2/2);

%% signal generation
signal = sin(phi);
