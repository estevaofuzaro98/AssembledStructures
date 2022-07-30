function [f]=hanning(N)
t=linspace(0,1,N);
f=0.5-0.5*cos(2*pi*t);
end