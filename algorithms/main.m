% PROJECT - VIBRATIONS OF ASSEMBLED STRUCTURES
% SIMULATIONS OF ORION'S DATA - Prof. Gael Chevallier
% Author: Estevao Fuzaro de Almeida
% Ilha Solteira, Brazil – July, 2022

%% INITIALIZATION
clc; clear all; 
close all;

% == This script changes all interpreters from tex to latex and sets all grids to on and minor.
list_factory = fieldnames(get(groot,'factory'));
index_interpreter = find(contains(list_factory,'Interpreter'));
for i = 1:length(index_interpreter)
    default_name = strrep(list_factory{index_interpreter(i)},'factory','default');
    set(groot, default_name,'latex');
end
set(groot,'defaultAxesXMinorGrid','on','defaultAxesXMinorGridMode','manual');
set(groot,'defaultAxesYMinorGrid','on','defaultAxesYMinorGridMode','manual');
set(groot,'defaultAxesZMinorGrid','on','defaultAxesZMinorGridMode','manual');
set(groot,{'DefaultAxesXColor','DefaultAxesYColor','DefaultAxesZColor'},{'k','k','k'});


%% PARAMETERS OF THE STUDY

par.f = 200E-3;     % Excitation Amplitude Force [N]
    
% == Glued Parameters
mass = 0.0917;      % Mass [kg]
damp = 2.7946;      % Damping [Ns/m]
stiff = 1.1204E7;   % Stiffness [N/m]
beta = 1.4862E11;   % Non-linear Damping [Ns/m^3]
alpha = -2.0090E15; % Non-linear Stiffness [N/m^3]

% == Aprroximation from Normalized Chart
par.m = mass*[1 1.03 0.99 0.95 0.91];           % Mass [kg]
par.c = damp*[1 1.09 1.22 1.70 1.92];           % Damping [Ns/m]
par.k = stiff*[1 1.02 0.97 0.92 0.86];          % Stiffness [N/m]
par.beta = beta*[1 2.33 4.77 11.10 31.04];      % Non-linear Damping [Ns/m^3]
par.alpha = alpha*[1 3.40 6.32 13.34 24.79];    % Non-linear Stiffness [N/m^3]

conditions = {'Excitation','Glued','80cNm','30cNm','20cNm','10cNm'};

fmin = 1680;        % Freq. Min. [Hz]
fmax = 1820;        % Freq. Max. [Hz]


%% PART 1 : COMPUTATION UNDER SWEEPSINE EXCITATION
SweepTime = 10;     % Time Duration of Sweep [s]
[par.time,par.signal] = SweepGeneration(fmin,fmax,SweepTime);
par.force = par.f * par.signal;

% == Plot Excitation - It takes some time

% TimePlot(par.time,par.force,'Excitation Force')
% title([conditions(1)],'FontWeight','normal')
% TimeFreqPlot(par.time,par.force)
% title([conditions(1)],'FontWeight','normal')
% 
% params.force = par.force;
% params.time = par.time;
% 
% % Computation for each condition
% for i=1:1:5
%     params.m = par.m(i);
%     params.c = par.c(i);
%     params.k = par.k(i);
%     params.beta = par.beta(i);
%     params.alpha = par.alpha(i);
% 
%     ops = odeset('OutputFcn',@odetpbar);
%     y0 = [0 0];
%     [t,y] = ode45(@(t,y) ForcedOrionBeam(t,y,params),par.time,y0,ops);
%     CompDisp = y(:,1)';
% 
%     TimePlot(t,CompDisp,'Computed Displacement')
%     title([conditions(i+1)],'FontWeight','normal')
%     TimeFreqPlot(t,CompDisp)
%     title([conditions(i+1)],'FontWeight','normal')
% end

% == Provide a list for user choose which condition he/she wants
spaceList = {'Glued','80cNm','30cNm','20cNm','10cNm'}; 
[idx, tf] = listdlg('ListString',spaceList,'SelectionMode','Single','PromptString','Which condition do you want?','Initialvalue',1,'ListSize',[250,150],'Name','Make a choice');
par.cond = idx;     % Evaluated Condition

% == Retrieve data from generated figures and plot Time, Freq-Time and Freq

% figNameT = strcat('figures/Part1_',string(spaceList(idx)),'_Time.fig');
% open(figNameT);
% a = get(gca,'Children');
% t = get(a, 'XData');
% CompDisp = get(a, 'YData');
% 
% TimeFreqPlot(t,CompDisp)
% title([conditions(par.cond+1)],'FontWeight','normal')
% 
% TFreqPlot(t,CompDisp',par)
% title([conditions(par.cond+1)],'FontWeight','normal')

figNameTF = strcat('figures/Part1_',string(spaceList(idx)),'_Time_Freq.fig');
open(figNameTF);
a = get(gca,'Children');
xdata = get(a, 'XData');
ydata = get(a, 'YData');
zdata = get(a, 'ZData');
[valMax,idxMax] = max(max(zdata));

fig = figure;
fig.WindowState = 'maximized';
plot(ydata,1E+3/par.f*zdata(:,idxMax),'k-','LineWidth',3)
xlabel('Frequency [Hz]','FontSize',16)
ylabel('Amplitude [mm/N]','FontSize',16)
set(gca,'YScale','linear','FontSize',16)
% xlim([fmin fmax])
title([conditions(par.cond+1)],'FontWeight','normal','FontSize',16)


%% PART 2 : COMPUTATION USING HARMONIC BALANCE
par.f = 200E-3;     % Excitation Amplitude Force [N]
par.phi = 0;        % Excitation Phase [rad]

clear Iter Step

% == Provide a list for user choose which condition he/she wants
spaceList = {'Glued','80cNm','30cNm','20cNm','10cNm'}; 
[idx, tf] = listdlg('ListString',spaceList,'SelectionMode','Single','PromptString','Which condition do you want?','Initialvalue',1,'ListSize',[250,150],'Name','Make a choice');
par.cond = idx;     % Evaluated Condition

% == Computing Parameters
Comp.FreqMin = fmin;
Comp.FreqMax = fmax;
Comp.direction = sign(Comp.FreqMax - Comp.FreqMin);

prMes = sprintf('How much harmonic fields do you want?');
button = questdlg(prMes, 'Question', '1', '2', '1');
switch button
    case '1'
        Comp.H = 1;     % One Harmonic Field
        U0 = [0;0;0];
        TanRef = [0,0,0,1];
        hComp = 4;
    case '2'
        Comp.H = 2;     % Two Harmonic Field
        U0 = [0;0;0;0;0];
        TanRef = [0,0,0,0,0,1];
        hComp = 6;
end

Comp.options = optimoptions('fsolve','Display','iter');

Comp.ds = 0.1;
Comp.dsmin = 1E-5;      % Initial Continuation Step
Comp.dsmax = 4;

% == Initial Computation
W0 = Comp.FreqMin*2*pi;
X0 = [U0;W0];

fun = @(X) HBM_Model(X,par,Comp,'R');
options = optimset(optimset(@fsolve),'Display','on','Jacobian','off','MaxIter',5000,'TolFun', 1.0000e-30,'TolX', 1.0000e-30);
[Xc,fval,exitflag,output] = fsolve(fun,X0,options);

% == Continuation Cycle
kkk = 1;
Fc = Xc(hComp,1)/2/pi;
Qc = vecnorm(Xc(2:hComp-1,1));
Xold = Xc;
X0 = Xc;
NItFail = 0;
NItStab = 0;
ds = Comp.ds;
dsmax = Comp.dsmax;
dsmin = Comp.dsmin;
dir = Comp.direction;

figure
while (Fc <= max([Comp.FreqMin,Comp.FreqMax]))&&(Fc >= min([Comp.FreqMin,Comp.FreqMax]))
    if (NItFail > 4)||(ds < 1E-9)
        break
    end
    
    % == Prediction
    [Xpre,t] = Prediction(X0,Xold,TanRef,par,Comp,ds,dir);
    Xp(:,kkk) = Xpre;
    Qp = vecnorm(Xpre(2:hComp-1,1));
    Fp = Xpre(hComp,1)/2/pi;
    subplot(2,1,1)
    plot(Fp,Qp,'o')
    plot([Fc,Fp],[Qc,Qp],':')
    
    % == Calculation Phase - Correction
    [Xcor,IterK,exitflag,fval] = Correction(Xpre,X0,par,Comp,ds);
    
    if kkk > 2
        dIter = IterK-Iter(1,kkk-1);
    else
        dIter = 0;
    end
    
    if exitflag==0
        ds = max([dsmin,ds/2]);
        NItFail = NItFail+1;
    else
        if (dIter > 3)
            ds = max([dsmin,ds/2]);
            NItStab = 0;
        elseif dIter == 0
               if NItStab > 5
                    ds = min([dsmax,ds*2]);
               end
            NItStab = NItStab+1;
        end
        NItFail = 0;
        Xold = X0;
        X0 = Xcor;
        Xc(:,kkk) = Xcor;
        Step(:,kkk) = ds;
        Qc = vecnorm(Xcor(2:hComp-1,1));
        Fc = Xcor(hComp,1)/2/pi;
        Iter(1,kkk)=IterK;
        TanRef = t';
        kkk = kkk + 1;
    end
    subplot(3,1,1)
    if Comp.H == 1
        Amp = 1E+3/par.f*sqrt(Xc(2,:).^2 + Xc(3,:).^2);
        plot(Xc(hComp,:),Amp,'x')
    elseif Comp.H == 2
        Amp = 1E+3/par.f*sqrt(Xc(2,:).^2 + Xc(3,:).^2 + Xc(4,:).^2 + Xc(5,:).^2);
        plot(Xc(hComp,:),Amp,'x')
    end
    hold on
    if Comp.H == 1
        Amp = 1E+3/par.f*sqrt(Xc(2,end).^2 + Xc(3,end).^2);
        plot(Xc(hComp,end),Amp,'xr')
    elseif Comp.H == 2
        Amp = 1E+3/par.f*sqrt(Xc(2,end).^2 + Xc(3,end).^2 + Xc(4,:).^2 + Xc(5,:).^2);
        plot(Xc(hComp,end),Amp,'xr')
    end
    hold off
    drawnow
    subplot(3,1,2)
    plot(Xc(hComp,:),Iter,'x')
    drawnow
    subplot(3,1,3)
    plot(Xc(hComp,:),Step,'x')
    drawnow
end

% == After
fig = figure;
fig.WindowState = 'maximized';
if Comp.H == 1
    Amp = 1E+3/par.f*(sqrt(Xc(2,:).^2 + Xc(3,:).^2));
    plot(Xc(hComp,:)/2/pi,Amp,'-','LineWidth',3)
elseif Comp.H == 2
    Amp = 1E+3/par.f*(sqrt(Xc(2,:).^2 + Xc(3,:).^2 + Xc(4,:).^2 + Xc(5,:).^2));
    plot(Xc(hComp,:)/2/pi,Amp,'-','LineWidth',3)
end
xlabel('Frequency [Hz]','FontSize',16)
ylabel('Amplitude [mm/N]','FontSize',16)
set(gca,'YScale','log','FontSize',16)
title([conditions(par.cond+1)],'FontWeight','normal','FontSize',16)


%% PART 3 : COMPUTATION OF FREE VIBRATIONS
% == Computation
fs = 10*fmax;   % Sample Frequency [Hz]
dt = 1/fs;      % Time Increment [s]
T = 2;          % Time Range [s]

par.time = 0:dt:T;  % Time Vector [s]

for i=1:1:5
    params.m = par.m(i);
    params.c = par.c(i);
    params.k = par.k(i);
    params.beta = par.beta(i);
    params.alpha = par.alpha(i);

    ops = odeset('OutputFcn',@odetpbar);
    y0 = [3E-6 0];

    [t,y] = ode23tb(@(t,y) FreeOrionBeam(t,y,params),par.time,y0,ops);
    CompDisp = y(:,1)';
%     TimePlot(t,CompDisp,'Computed Displacement')
%     title([strcat(conditions(i+1), ' - ode23tb')],'FontWeight','normal')
%     xlim([0,0.4])

%     [t,y] = ode45(@(t,y) FreeOrionBeam(t,y,params),par.time,y0,ops);
%     CompDisp = y(:,1)';
%     TimePlot(t,CompDisp,'Computed Displacement')
%     title([strcat(conditions(i+1), ' - ode45')],'FontWeight','normal')
%     xlim([0,0.4])

    %     TimeFreqPlot(t,CompDisp)
    %     title([conditions(i+1)],'FontWeight','normal')

%     TFreqPlot(t,CompDisp',par)
%     title([conditions(i+1)],'FontWeight','normal')
%     a = get(gca,'Children');
%     xdata = get(a, 'XData');
%     ydata = get(a, 'YData');
%     [valMax,idxMax] = max(ydata);
%     wn = xdata(idxMax);


    % PART 4 : IDENTIFICATION OF PROPERTIES FROM THE FREE VIBRATIONS
    CompDispH = hilbert(CompDisp);

%     TimePlot(t,CompDisp,'Computed Displacement')
%     hold on
%     plot(t,abs(CompDispH),'LineWidth',3)
%     hold off
%     title([strcat(conditions(i+1))],'FontWeight','normal')
%     xlim([0,0.4])

    %     firstPoint = [1200 1200 1000 800 800];
    %     lastPoint = [4100 4100 3200 2000 1500];
    %
    %     fig = figure;
    %     fig.WindowState = 'maximized';
    %     tp = t(firstPoint(i):lastPoint(i),1);
    %     CompDispHp = CompDispH(firstPoint(i):lastPoint(i),1);
    %     plot(t,log(abs(CompDispH)),'k','linewidth', 2), hold on
    %     plot(tp,log(abs(CompDispHp)),'r','linewidth', 2.2), hold on
    %     title([conditions(i+1)],'FontWeight','normal')
    %     xlabel('Time [s]','FontSize',16)
    %     ylabel('$\ln |\mathcal{H} [x(t)]|$','FontSize',16)
    %     set(gca,'YScale','linear','FontSize',16)
    %
    %     fig = figure;
    %     fig.WindowState = 'maximized';
    %     plot(t,log(abs(CompDispH)),'k','linewidth', 2), hold on
    %     plot(tp,log(abs(CompDispHp)),'r','linewidth', 2.2), hold on
    %     title([conditions(i+1)],'FontWeight','normal')
    %     xlim([0 0.4])
    %     xlabel('Time [s]','FontSize',16)
    %     ylabel('$\ln |\mathcal{H} [x(t)]|$','FontSize',16)
    %     set(gca,'YScale','linear','FontSize',16)
    %
    %     grad = polyfit(tp,log(CompDispHp),1);
    %     Cc(i) = -real(grad(1))/wn;

    % Instantaneous amplitude
    instAmp = abs(CompDispH);

    % Instantaneous frequency
    [instFreq,tIF] = instfreq(CompDisp,fs);
    InterpIF = interp1(tIF,instFreq,t);

    % Instantaneous damping
    CsiInst = (par.c(i) + (par.beta(i)*instAmp.^2)/4)./(pi*par.m(i)*InterpIF);
    
%     fig = figure;
%     fig.WindowState = 'maximized';
%     plot(t,InterpIF,'k','linewidth', 2)
%     xlim([0 0.4])
%     title([strcat(conditions(i+1),' - Inst. Freq.')],'FontWeight','normal')
%     xlabel('Time [s]','FontSize',16)
%     ylabel('Frequency [Hz]','FontSize',16)
%     set(gca,'YScale','linear','FontSize',16)

    if i==1
        color = 'k';
    elseif i==2
        color = 'm';
    elseif i==3
        color = 'g';
    elseif i==4
        color = 'b';
    elseif i==5
        color = 'r';
    end
    
    fig1 = figure(1);
    fig1.WindowState = 'maximized';
    subplot(1,2,1)
    plot(t,instAmp,color,'linewidth', 3), hold on
    xlim([0 0.4])
    title('Hilbert Envelope','FontWeight','normal')
    xlabel('Time [s]','FontSize',24)
    ylabel('Amplitude [mm/N]','FontSize',24)
    set(gca,'YScale','linear','FontSize',24)

    subplot(1,2,2)
    plot(t,InterpIF,color,'linewidth', 3), hold on
    xlim([0 0.6])
    title('Instantaneous Frequency','FontWeight','normal')
    xlabel('Time [s]','FontSize',24)
    ylabel('Frequency [Hz]','FontSize',24)
    set(gca,'YScale','linear','FontSize',24)
    
    fig2 = figure(2);
    fig2.WindowState = 'maximized';
    if i==1
        subplot(5,2,1)
        plot(instAmp,InterpIF,color,'linewidth', 3), hold on
        title('Backbone Curves','FontWeight','normal')
        set(gca,'xticklabel',[])
        xlim([0 2.7E-6])
        ylim([1753 1756])
        set(gca,'YScale','linear','FontSize',20)
    elseif i==2
        subplot(5,2,3)
        plot(instAmp,InterpIF,color,'linewidth', 3), hold on
        set(gca,'xticklabel',[])
        xlim([0 2.7E-6])
        ylim([1742 1748])
        set(gca,'YScale','linear','FontSize',20)
    elseif i==3
        subplot(5,2,5)
        plot(instAmp,InterpIF,color,'linewidth', 3), hold on
        set(gca,'xticklabel',[])
        xlim([0 2.7E-6])
        ylim([1732 1738])
        ylabel('Frequency [Hz]','FontSize',20)
        set(gca,'YScale','linear','FontSize',20)
    elseif i==4
        subplot(5,2,7)
        plot(instAmp,InterpIF,color,'linewidth', 3), hold on
        set(gca,'xticklabel',[])
        xlim([0 2.7E-6])
        ylim([1719 1728])
        set(gca,'YScale','linear','FontSize',20)
    elseif i==5
        subplot(5,2,9)
        plot(instAmp,InterpIF,color,'linewidth', 3), hold on
        xlim([0 2.7E-6])
        ylim([1695 1707])
        xlabel('Amplitude [mm/N]','FontSize',20)
        set(gca,'YScale','linear','FontSize',20)
    end
    
    if i==1
        subplot(5,2,2)
        plot(instAmp,abs(CsiInst),color,'linewidth', 3), hold on
        title('Damping Curves','FontWeight','normal')
        set(gca,'xticklabel',[])
        xlim([0 2.7E-6])
        ylim([5.4E-3 6.1E-3])
        set(gca,'YScale','linear','FontSize',20)
    elseif i==2
        subplot(5,2,4)
        plot(instAmp,abs(CsiInst),color,'linewidth', 3), hold on
        set(gca,'xticklabel',[])
        xlim([0 2.7E-6])
        ylim([5.8E-3 7.1E-3])
        set(gca,'YScale','linear','FontSize',20)
    elseif i==3
        subplot(5,2,6)
        plot(instAmp,abs(CsiInst),color,'linewidth', 3), hold on
        set(gca,'xticklabel',[])
        xlim([0 2.7E-6])
        ylim([6.8E-3 9.1E-3])
        ylabel('Damping $[\%]$','FontSize',20)
        set(gca,'YScale','linear','FontSize',20)
    elseif i==4
        subplot(5,2,8)
        plot(instAmp,abs(CsiInst),color,'linewidth', 3), hold on
        set(gca,'xticklabel',[])
        xlim([0 2.7E-6])
        ylim([0.9E-2 1.5E-2])
        set(gca,'YScale','linear','FontSize',20)
    elseif i==5
        subplot(5,2,10)
        plot(instAmp,abs(CsiInst),color,'linewidth', 3), hold on
        xlim([0 2.7E-6])
        ylim([1E-2 2.1E-2])
        xlabel('Amplitude [mm/N]','FontSize',20)
        set(gca,'YScale','linear','FontSize',20)
    end

end
hold off
% fprintf('Csi: %.4f\n',Cc);
