function [Xp,t] = Prediction(X0,Xold,TanRef,par,Comp,ds,dir)

Lo = length(X0);
J = HBM_Model(X0,par,Comp,'J');
C = ones(1,Lo);
%C = [0,0,0,1];
C = TanRef;
Z = zeros(Lo-1,1);

A = [ J ; C ];
B = [Z;1];
t = A\B;

t = t/norm(t);

if ((X0 - Xold)'*t*dir) < 0
    Xp = X0 - dir * ds * t;
    disp('chgt direction')
    pause
else
    Xp = X0 + dir * ds * t;
end

end
