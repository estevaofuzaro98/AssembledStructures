%% PART 2 : COMPUTATION USING HARMONIC BALANCE

par.force = linspace(2E-30,2E-1,10);

% == Provide a list for user choose which condition he/she wants
spaceList = {'Glued','80cNm','30cNm','20cNm','10cNm'};
[idx, tf] = listdlg('ListString',spaceList,'SelectionMode','Single','PromptString','Which condition do you want?','Initialvalue',1,'ListSize',[250,150],'Name','Make a choice');
par.cond = idx;     % Evaluated Condition

fig = figure;
fig.WindowState = 'maximized';

for ii = 1:size(par.force,2)

    fprintf('Stage %i of %i',ii,size(par.force,2));

    par.f = par.force(ii);     % Excitation Amplitude Force [N]
    par.phi = 0;        % Excitation Phase [rad]

    clear Iter Step

    % == Computing Parameters
    Comp.FreqMin = fmin;
    Comp.FreqMax = fmax;
    Comp.direction = sign(Comp.FreqMax - Comp.FreqMin);

    Comp.H = 1;     % One Harmonic Field
    U0 = [0;0;0];
    TanRef = [0,0,0,1];
    hComp = 4;

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

%     figure
    while (Fc <= max([Comp.FreqMin,Comp.FreqMax]))&&(Fc >= min([Comp.FreqMin,Comp.FreqMax]))
        if (NItFail > 4)||(ds < 1E-9)
            break
        end

        % == Prediction
        [Xpre,t] = Prediction(X0,Xold,TanRef,par,Comp,ds,dir);
        Xp(:,kkk) = Xpre;
        Qp = vecnorm(Xpre(2:hComp-1,1));
        Fp = Xpre(hComp,1)/2/pi;
%         subplot(2,1,1)
%         plot(Fp,Qp,'o')
%         plot([Fc,Fp],[Qc,Qp],':')

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
%         subplot(3,1,1)
        if Comp.H == 1
            Amp = 1E+3/par.f*sqrt(Xc(2,:).^2 + Xc(3,:).^2);
%             plot(Xc(hComp,:),Amp,'x')
        end
%         hold on
        if Comp.H == 1
            Amp = 1E+3/par.f*sqrt(Xc(2,end).^2 + Xc(3,end).^2);
%             plot(Xc(hComp,end),Amp,'xr')
        end
%         hold off
%         drawnow
%         subplot(3,1,2)
%         plot(Xc(hComp,:),Iter,'x')
%         drawnow
%         subplot(3,1,3)
%         plot(Xc(hComp,:),Step,'x')
%         drawnow
    end

    % == After
%     fig = figure;
%     fig.WindowState = 'maximized';
    if Comp.H == 1
        Amp = 1E+3/par.f*(sqrt(Xc(2,:).^2 + Xc(3,:).^2));
        plot(Xc(hComp,:)/2/pi,Amp,'-','LineWidth',3)
        hold on
    end
    xlabel('Frequency [Hz]','FontSize',16)
    ylabel('Amplitude [mm/N]','FontSize',16)
    set(gca,'YScale','log','FontSize',16)
    title([conditions(par.cond+1)],'FontWeight','normal','FontSize',16)
end