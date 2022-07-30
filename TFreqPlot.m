function TFreqPlot(t,y,par)
% Frequency analysis
% Estevao Fuzaro 2022

dt = t(2)-t(1);     % Time Increment [s]
Fs = 1/dt;          % Sample Frequency [Hz]
duration = t(end);  % Signal Duration [s]
N = Fs*duration;    % Number of Samples [-]

y = y(1:end-1);
t = t(1:end-1);

y = y.*hanning(N)';

y = [y; zeros(2000,1)];
N2 = length(y);

Y = fft(y);
Y_1Side = Y(1:N2/2);
f = Fs*(0:N2/2-1)/N2;
Y_meg = abs(Y_1Side)/(N/4);

fig = figure;
fig.WindowState = 'maximized';
plot(f,1E+3/par.f*Y_meg,'k-','LineWidth',3)
xlabel('Frequency [Hz]','FontSize',16)
ylabel('Amplitude [mm/N]','FontSize',16)
set(gca,'YScale','linear','FontSize',16)

end