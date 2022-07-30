function [Xc,iter,exitflag,fval] = Correction(Xp,X0,par,Comp,ds)

% arc-length
pc = @(X) ([X-X0].'*[X-X0] - ds^2);

% orthogonal
%pc = @(X) ([Xp-X0].'*[X-Xp] - ds^2);

R = @(X) HBM_Model(X,par,Comp,'R');
cor = @(X) [R(X); pc(X)];

options = optimset(optimset(@fsolve),'Display','off','Jacobian','off','MaxIter',5000,'TolFun', 1.0000e-30,'TolX', 1.0000e-30);
[Xc,fval,exitflag,output] = fsolve(cor,Xp,options);
iter = output.iterations;
end