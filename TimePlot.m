function TimePlot(t,y,label)
%
% function TimePlot allows to save figure into .pdf or .png and .fig
%
% INPUT
% t : absisse
% y : ordonnée
% label : ylabel
%
% Gaël C. 2021.1

fig = figure;
fig.WindowState = 'maximized';
plot(t,y,'-','LineWidth',2)
xlabel('Time [s]')
ylabel(label)
set(gca,'FontSize',14,'XGrid','on','YGrid','on')
end