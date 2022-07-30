%% Initialization - Duffing-Van der Pol Symbolic
close all; clear all; clc;

%% Symbolic Fourier Expansion
syms m k c alpha beta f W phi
syms u(t) Uo Uc1 Us1 Uc3 Us3 Uc5 Us5

promptMessage = sprintf('How much harmonic fields do you want?');
button = questdlg(promptMessage, 'Question', '1', '2', '1');
switch button
    case '1'
        % Single Harmonic Field
        u(t) = Uo + Uc1 * cos(W*t) + Us1 * sin(W*t);
        du(t) = diff(u,t);
        ddu(t) = diff(du,t);
        r(t) = m*ddu(t) + c*du(t) + k*u(t) + alpha*u(t)^3 + beta*u(t)^2*du(t) - f*cos(W*t+phi);

        T = 2*pi/W;

        Ro(Uo,Uc1,Us1,W) = expand(1/T * int(r(t),0,T))
        Rc1(Uo,Uc1,Us1,W) = expand(2/T * int(r(t)*cos(W*t),0,T))
        Rs1(Uo,Uc1,Us1,W) = expand(2/T * int(r(t)*sin(W*t),0,T))

        % Derivation of Ro
        dRodUo = expand(diff(Ro,Uo))
        dRodUc1 = expand(diff(Ro,Uc1))
        dRodUs1 = expand(diff(Ro,Us1))
        dRodW = expand(diff(Ro,W))

        % Derivation of Rc
        dRcdUo = expand(diff(Rc1,Uo))
        dRcdUc1 = expand(diff(Rc1,Uc1))
        dRcdUs1 = expand(diff(Rc1,Us1))
        dRcdW = expand(diff(Rc1,W))

        % Derivation of Rs
        dRsdUo = expand(diff(Rs1,Uo))
        dRsdUc1 = expand(diff(Rs1,Uc1))
        dRsdUs1 = expand(diff(Rs1,Us1))
        dRsdW = expand(diff(Rs1,W))
        
    case '2'
        % Two Harmonic Field
        u(t) = Uo + Uc1 * cos(W*t) + Us1 * sin(W*t) + Uc3 * cos(3*W*t) + Us3 * sin(3*W*t);
        du(t) = diff(u,t);
        ddu(t) = diff(du,t);
        r(t) = m*ddu(t) + c*du(t) + k*u(t) + alpha*u(t)^3 + beta*u(t)^2*du(t) - f*cos(W*t+phi);

        T = 2*pi/W;

        Ro(Uo,Uc1,Us1,Uc3,Us3,W) = expand(1/T * int(r(t),0,T))
        Rc1(Uo,Uc1,Us1,Uc3,Us3,W) = expand(2/T * int(r(t)*cos(W*t),0,T))
        Rs1(Uo,Uc1,Us1,Uc3,Us3,W) = expand(2/T * int(r(t)*sin(W*t),0,T))
        Rc3(Uo,Uc1,Us1,Uc3,Us3,W) = expand(2/T * int(r(t)*cos(3*W*t),0,T))
        Rs3(Uo,Uc1,Us1,Uc3,Us3,W) = expand(2/T * int(r(t)*sin(3*W*t),0,T))

        % Derivation of Ro
        dRodUo = expand(diff(Ro,Uo))
        dRodUc1 = expand(diff(Ro,Uc1))
        dRodUs1 = expand(diff(Ro,Us1))
        dRodUc3 = expand(diff(Ro,Uc3))
        dRodUs3 = expand(diff(Ro,Us3))
        dRodW = expand(diff(Ro,W))

        % Derivation of Rc1
        dRc1dUo = expand(diff(Rc1,Uo))
        dRc1dUc1 = expand(diff(Rc1,Uc1))
        dRc1dUs1 = expand(diff(Rc1,Us1))
        dRc1dUc3 = expand(diff(Rc1,Uc3))
        dRc1dUs3 = expand(diff(Rc1,Us3))
        dRc1dW = expand(diff(Rc1,W))

        % Derivation of Rs1
        dRs1dUo = expand(diff(Rs1,Uo))
        dRs1dUc1 = expand(diff(Rs1,Uc1))
        dRs1dUs1 = expand(diff(Rs1,Us1))
        dRs1dUc3 = expand(diff(Rs1,Uc3))
        dRs1dUs3 = expand(diff(Rs1,Us3))
        dRs1dW = expand(diff(Rs1,W))

        % Derivation of Rc3
        dRc3dUo = expand(diff(Rc3,Uo))
        dRc3dUc1 = expand(diff(Rc3,Uc1))
        dRc3dUs1 = expand(diff(Rc3,Us1))
        dRc3dUc3 = expand(diff(Rc3,Uc3))
        dRc3dUs3 = expand(diff(Rc3,Us3))
        dRc3dW = expand(diff(Rc3,W))

        % Derivation of Rs3
        dRs3dUo = expand(diff(Rs3,Uo))
        dRs3dUc1 = expand(diff(Rs3,Uc1))
        dRs3dUs1 = expand(diff(Rs3,Us1))
        dRs3dUc3 = expand(diff(Rs3,Uc3))
        dRs3dUs3 = expand(diff(Rs3,Us3))
        dRs3dW = expand(diff(Rs3,W))

end