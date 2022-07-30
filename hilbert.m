function yh=hilbert(y)
%
% function hilbert.m allows to compute an analytical signal from a real one
%
% INPUTS
% y : signal
%
% Gaël CHEVALLIER 2021.1
%

Y=fft(y);
N = length(y);
if isrow(Y)
    Y = Y.';
end
if N/2-floor(N/2)==0
    YH = [Y(1:N/2+1,1);zeros(N/2-1,1)];
else
    YH = [Y(1:(N+1)/2,1);zeros((N-1)/2,1)];
end
yh=2*ifft(YH);