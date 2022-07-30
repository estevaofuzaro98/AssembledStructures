function TimeFreqPlot(t,v)
%
% Time Frequency analysis
%
% Gaël CHEVALLIER 2022

Ns = length(t);   % Size the sample to analyse
Nk = 1*256;         % Size of the samples for the time frequency analysis
dt = t(2)-t(1);    % Sample time
Fs = 1/dt;         % Sampling frequency

per = linspace(0,1,Nk);
f = Fs/2*linspace(0,1,Nk/2+1);

dNk = Nk/4;
ix = Nk/2;
ii = 1;
while ix < Ns-Nk/2
    index = [ix-Nk/2+1 : ix+Nk/2] ;
    vk = v(index);
    vk = vk.*hanning(Nk);
    Y(:,ii) = fft(vk,Nk)/Nk;
    ix = ix + dNk;
    ii = ii+1;
end

tk = t(Nk/2:dNk:Ns-Nk/2);
Y = (abs(Y(1:end/2+1,:)));

fig = figure;
fig.WindowState = 'maximized';
su = surf(tk,f,Y);
view(0,90)
set(su,'Linestyle','None','FaceAlpha',0.6)
xlabel('Time [s]','Fontsize',14)
ylabel('Frequency [Hz]','Fontsize',14)
map = colormap(jet(64));
colormap(map);
set(gca,'FontSize',14)
grid on;grid minor;

end