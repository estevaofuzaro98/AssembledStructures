function [OUT] = HBM_Model(X,par,Comp,output)

H = Comp.H;

cond = par.cond;

m = par.m(cond);
c = par.c(cond);
k = par.k(cond);
alpha = par.alpha(cond);
beta = par.beta(cond);
f = par.f;
phi = par.phi;

Q = X(1:end-1,1);
W = X(end,1);

if H==1     % ONE HARMONIC FIELD
    Qo = Q(1,1);
    Qc1 = Q(2,1);
    Qs1 = Q(3,1);

    % Analytic Definition
    Ro = (3*alpha*Qc1^2*Qo)/2 + alpha*Qo^3 + (3*alpha*Qo*Qs1^2)/2 + k*Qo;
    Rc = (3*alpha*Qc1^3)/4 + (beta*Qc1^2*Qs1*W)/4 + 3*alpha*Qc1*Qo^2 + (3*alpha*Qc1*Qs1^2)/4 - m*Qc1*W^2 + k*Qc1 + beta*Qo^2*Qs1*W + (beta*Qs1^3*W)/4 + c*Qs1*W - f*cos(phi);
    Rs = - (beta*Qc1^3*W)/4 + (3*alpha*Qc1^2*Qs1)/4 - beta*Qc1*Qo^2*W - (beta*Qc1*Qs1^2*W)/4 - c*Qc1*W + 3*alpha*Qo^2*Qs1 + (3*alpha*Qs1^3)/4 - m*Qs1*W^2 + k*Qs1 + f*sin(phi);

    dRodQo = (3*alpha*Qc1^2)/2 + 3*alpha*Qo^2 + (3*alpha*Qs1^2)/2 + k;
    dRodQc = 3*Qc1*Qo*alpha;
    dRodQs = 3*Qo*Qs1*alpha;
    dRodW = 0;

    dRcdQo = 6*Qc1*Qo*alpha + 2*Qo*Qs1*W*beta;
    dRcdQc = (9*alpha*Qc1^2)/4 + (beta*Qc1*Qs1*W)/2 + 3*alpha*Qo^2 + (3*alpha*Qs1^2)/4 - m*W^2 + k;
    dRcdQs = (W*beta*Qc1^2)/4 + (3*alpha*Qc1*Qs1)/2 + W*beta*Qo^2 + (3*W*beta*Qs1^2)/4 + W*c;
    dRcdW = (beta*Qc1^2*Qs1)/4 - 2*W*m*Qc1 + beta*Qo^2*Qs1 + (beta*Qs1^3)/4 + c*Qs1;

    dRsdQo = 6*Qo*Qs1*alpha - 2*Qc1*Qo*W*beta;
    dRsdQc = - W*beta*Qo^2 - W*c - (3*Qc1^2*W*beta)/4 - (Qs1^2*W*beta)/4 + (3*Qc1*Qs1*alpha)/2;
    dRsdQs = (3*alpha*Qc1^2)/4 - (beta*Qc1*Qs1*W)/2 + 3*alpha*Qo^2 + (9*alpha*Qs1^2)/4 - m*W^2 + k;
    dRsdW = - Qc1*beta*Qo^2 - Qc1*c - (Qc1^3*beta)/4 - (Qc1*Qs1^2*beta)/4 - 2*Qs1*W*m;
    
    R = [Ro;Rc;Rs];

    dRdQ = [dRodQo,dRodQc,dRodQs;...
        dRcdQo,dRcdQc,dRcdQs;...
        dRsdQo,dRsdQc,dRsdQs];

    dRdW = [dRodW;dRcdW;dRsdW];

elseif H==2     % TWO HARMONIC FIELD
    Qo = Q(1,1);
    Qc1 = Q(2,1);
    Qs1 = Q(3,1);
    Qc3 = Q(4,1);
    Qs3 = Q(5,1);

    % Analytic Definition
    Ro = (3*alpha*Qc1^2*Qo)/2 + (3*alpha*Qc3^2*Qo)/2 + alpha*Qo^3 + (3*alpha*Qo*Qs1^2)/2 + (3*alpha*Qo*Qs3^2)/2 + k*Qo;
    Rc1 = (3*alpha*Qc1^3)/4 + (3*alpha*Qc1^2*Qc3)/4 + (beta*Qc1^2*Qs1*W)/4 + (beta*Qc1^2*Qs3*W)/4 + (3*alpha*Qc1*Qc3^2)/2 - (beta*Qc1*Qc3*Qs1*W)/2 + 3*alpha*Qc1*Qo^2 + (3*alpha*Qc1*Qs1^2)/4 + (3*alpha*Qc1*Qs1*Qs3)/2 + (3*alpha*Qc1*Qs3^2)/2 - m*Qc1*W^2 + k*Qc1 + (beta*Qc3^2*Qs1*W)/2 - (3*alpha*Qc3*Qs1^2)/4 + beta*Qo^2*Qs1*W + (beta*Qs1^3*W)/4 - (beta*Qs1^2*Qs3*W)/4 + (beta*Qs1*Qs3^2*W)/2 + c*Qs1*W - f*cos(phi);
    Rs1 = - (beta*Qc1^3*W)/4 - (beta*Qc1^2*Qc3*W)/4 + (3*alpha*Qc1^2*Qs1)/4 + (3*alpha*Qc1^2*Qs3)/4 - (beta*Qc1*Qc3^2*W)/2 - (3*alpha*Qc1*Qc3*Qs1)/2 - beta*Qc1*Qo^2*W - (beta*Qc1*Qs1^2*W)/4 - (beta*Qc1*Qs1*Qs3*W)/2 - (beta*Qc1*Qs3^2*W)/2 - c*Qc1*W + (3*alpha*Qc3^2*Qs1)/2 + (beta*Qc3*Qs1^2*W)/4 + 3*alpha*Qo^2*Qs1 + (3*alpha*Qs1^3)/4 - (3*alpha*Qs1^2*Qs3)/4 + (3*alpha*Qs1*Qs3^2)/2 - m*Qs1*W^2 + k*Qs1 + f*sin(phi);
    Rc3 = (alpha*Qc1^3)/4 + (3*alpha*Qc1^2*Qc3)/2 + (3*beta*Qc1^2*Qs1*W)/4 + (3*beta*Qc1^2*Qs3*W)/2 - (3*alpha*Qc1*Qs1^2)/4 + (3*alpha*Qc3^3)/4 + (3*beta*Qc3^2*Qs3*W)/4 + 3*alpha*Qc3*Qo^2 + (3*alpha*Qc3*Qs1^2)/2 + (3*alpha*Qc3*Qs3^2)/4 - 9*m*Qc3*W^2 + k*Qc3 + 3*beta*Qo^2*Qs3*W - (beta*Qs1^3*W)/4 + (3*beta*Qs1^2*Qs3*W)/2 + (3*beta*Qs3^3*W)/4 + 3*c*Qs3*W;
    Rs3 = - (beta*Qc1^3*W)/4 - (3*beta*Qc1^2*Qc3*W)/2 + (3*alpha*Qc1^2*Qs1)/4 + (3*alpha*Qc1^2*Qs3)/2 + (3*beta*Qc1*Qs1^2*W)/4 - (3*beta*Qc3^3*W)/4 + (3*alpha*Qc3^2*Qs3)/4 - 3*beta*Qc3*Qo^2*W - (3*beta*Qc3*Qs1^2*W)/2 - (3*beta*Qc3*Qs3^2*W)/4 - 3*c*Qc3*W + 3*alpha*Qo^2*Qs3 - (alpha*Qs1^3)/4 + (3*alpha*Qs1^2*Qs3)/2 + (3*alpha*Qs3^3)/4 - 9*m*Qs3*W^2 + k*Qs3;

    dRodQo = (3*alpha*Qc1^2)/2 + (3*alpha*Qc3^2)/2 + 3*alpha*Qo^2 + (3*alpha*Qs1^2)/2 + (3*alpha*Qs3^2)/2 + k;
    dRodQc1 = 3*Qc1*Qo*alpha;
    dRodQs1 = 3*Qo*Qs1*alpha;
    dRodQc3 = 3*Qc3*Qo*alpha;
    dRodQs3 = 3*Qo*Qs3*alpha;
    dRodW = 0;

    dRc1dQo = 6*Qc1*Qo*alpha + 2*Qo*Qs1*W*beta;
    dRc1dQc1 = (9*alpha*Qc1^2)/4 + (3*alpha*Qc1*Qc3)/2 + (beta*Qc1*Qs1*W)/2 + (beta*Qc1*Qs3*W)/2 + (3*alpha*Qc3^2)/2 - (beta*Qc3*Qs1*W)/2 + 3*alpha*Qo^2 + (3*alpha*Qs1^2)/4 + (3*alpha*Qs1*Qs3)/2 + (3*alpha*Qs3^2)/2 - m*W^2 + k;
    dRc1dQs1 = (W*beta*Qc1^2)/4 - (W*beta*Qc1*Qc3)/2 + (3*alpha*Qc1*Qs1)/2 + (3*alpha*Qc1*Qs3)/2 + (W*beta*Qc3^2)/2 - (3*alpha*Qc3*Qs1)/2 + W*beta*Qo^2 + (3*W*beta*Qs1^2)/4 - (W*beta*Qs1*Qs3)/2 + (W*beta*Qs3^2)/2 + W*c;
    dRc1dQc3 = (3*Qc1^2*alpha)/4 - (3*Qs1^2*alpha)/4 + 3*Qc1*Qc3*alpha - (Qc1*Qs1*W*beta)/2 + Qc3*Qs1*W*beta;
    dRc1dQs3 = (3*Qc1*Qs1*alpha)/2 + 3*Qc1*Qs3*alpha + (Qc1^2*W*beta)/4 - (Qs1^2*W*beta)/4 + Qs1*Qs3*W*beta;
    dRc1dW = (beta*Qc1^2*Qs1)/4 + (beta*Qc1^2*Qs3)/4 - (beta*Qc1*Qc3*Qs1)/2 - 2*W*m*Qc1 + (beta*Qc3^2*Qs1)/2 + beta*Qo^2*Qs1 + (beta*Qs1^3)/4 - (beta*Qs1^2*Qs3)/4 + (beta*Qs1*Qs3^2)/2 + c*Qs1;

    dRs1dQo = 6*Qo*Qs1*alpha - 2*Qc1*Qo*W*beta;
    dRs1dQc1 = - W*beta*Qo^2 - W*c - (3*Qc1^2*W*beta)/4 - (Qc3^2*W*beta)/2 - (Qs1^2*W*beta)/4 - (Qs3^2*W*beta)/2 + (3*Qc1*Qs1*alpha)/2 + (3*Qc1*Qs3*alpha)/2 - (3*Qc3*Qs1*alpha)/2 - (Qc1*Qc3*W*beta)/2 - (Qs1*Qs3*W*beta)/2;
    dRs1dQs1 = (3*alpha*Qc1^2)/4 - (3*alpha*Qc1*Qc3)/2 - (beta*Qc1*Qs1*W)/2 - (beta*Qc1*Qs3*W)/2 + (3*alpha*Qc3^2)/2 + (beta*Qc3*Qs1*W)/2 + 3*alpha*Qo^2 + (9*alpha*Qs1^2)/4 - (3*alpha*Qs1*Qs3)/2 + (3*alpha*Qs3^2)/2 - m*W^2 + k;
    dRs1dQc3 = 3*Qc3*Qs1*alpha - (3*Qc1*Qs1*alpha)/2 - (Qc1^2*W*beta)/4 + (Qs1^2*W*beta)/4 - Qc1*Qc3*W*beta;
    dRs1dQs3 = (3*Qc1^2*alpha)/4 - (3*Qs1^2*alpha)/4 + 3*Qs1*Qs3*alpha - (Qc1*Qs1*W*beta)/2 - Qc1*Qs3*W*beta;
    dRs1dW = (Qc3*Qs1^2*beta)/4 - (Qc1^3*beta)/4 - 2*Qs1*W*m - (Qc1*Qc3^2*beta)/2 - (Qc1^2*Qc3*beta)/4 - Qc1*Qo^2*beta - (Qc1*Qs1^2*beta)/4 - (Qc1*Qs3^2*beta)/2 - Qc1*c - (Qc1*Qs1*Qs3*beta)/2;

    dRc3dQo = 6*Qc3*Qo*alpha + 6*Qo*Qs3*W*beta;
    dRc3dQc1 = (3*Qc1^2*alpha)/4 - (3*Qs1^2*alpha)/4 + 3*Qc1*Qc3*alpha + (3*Qc1*Qs1*W*beta)/2 + 3*Qc1*Qs3*W*beta;
    dRc3dQs1 = (3*W*beta*Qc1^2)/4 - (3*Qs1*alpha*Qc1)/2 - (3*Qs1^2*W*beta)/4 + 3*Qc3*Qs1*alpha + 3*Qs1*Qs3*W*beta;
    dRc3dQc3 = (3*alpha*Qc1^2)/2 + (9*alpha*Qc3^2)/4 + (3*beta*Qc3*Qs3*W)/2 + 3*alpha*Qo^2 + (3*alpha*Qs1^2)/2 + (3*alpha*Qs3^2)/4 - 9*m*W^2 + k;
    dRc3dQs3 = (3*W*beta*Qc1^2)/2 + (3*W*beta*Qc3^2)/4 + (3*alpha*Qc3*Qs3)/2 + 3*W*beta*Qo^2 + (3*W*beta*Qs1^2)/2 + (9*W*beta*Qs3^2)/4 + 3*W*c;
    dRc3dW = (3*beta*Qc1^2*Qs1)/4 + (3*beta*Qc1^2*Qs3)/2 + (3*beta*Qc3^2*Qs3)/4 - 18*W*m*Qc3 + 3*beta*Qo^2*Qs3 - (beta*Qs1^3)/4 + (3*beta*Qs1^2*Qs3)/2 + (3*beta*Qs3^3)/4 + 3*c*Qs3;
    
    dRs3dQo = 6*Qo*Qs3*alpha - 6*Qc3*Qo*W*beta;
    dRs3dQc1 = (3*W*beta*Qs1^2)/4 + (3*Qc1*alpha*Qs1)/2 - (3*Qc1^2*W*beta)/4 + 3*Qc1*Qs3*alpha - 3*Qc1*Qc3*W*beta;
    dRs3dQs1 = (3*alpha*Qc1^2)/4 + (3*Qs1*W*beta*Qc1)/2 - (3*Qs1^2*alpha)/4 + 3*Qs1*Qs3*alpha - 3*Qc3*Qs1*W*beta;
    dRs3dQc3 = - 3*W*beta*Qo^2 - 3*W*c - (3*Qc1^2*W*beta)/2 - (9*Qc3^2*W*beta)/4 - (3*Qs1^2*W*beta)/2 - (3*Qs3^2*W*beta)/4 + (3*Qc3*Qs3*alpha)/2;
    dRs3dQs3 = (3*alpha*Qc1^2)/2 + (3*alpha*Qc3^2)/4 - (3*beta*Qc3*Qs3*W)/2 + 3*alpha*Qo^2 + (3*alpha*Qs1^2)/2 + (9*alpha*Qs3^2)/4 - 9*m*W^2 + k;
    dRs3dW = (3*Qc1*Qs1^2*beta)/4 - (Qc1^3*beta)/4 - (3*Qc3^3*beta)/4 - 18*Qs3*W*m - (3*Qc1^2*Qc3*beta)/2 - 3*Qc3*Qo^2*beta - 3*Qc3*c - (3*Qc3*Qs1^2*beta)/2 - (3*Qc3*Qs3^2*beta)/4;

    R = [Ro;Rc1;Rs1;Rc3;Rs3];

    dRdQ = [dRodQo,dRodQc1,dRodQs1,dRodQc3,dRodQs3;...
        dRc1dQo,dRc1dQc1,dRc1dQs1,dRc1dQc3,dRc1dQs3;...
        dRs1dQo,dRs1dQc1,dRs1dQs1,dRs1dQc3,dRs1dQs3;...
        dRc3dQo,dRc3dQc1,dRc3dQs1,dRc3dQc3,dRc3dQs3;...
        dRs3dQo,dRs3dQc1,dRs3dQs1,dRs3dQc3,dRs3dQs3];

    dRdW = [dRodW;dRc1dW;dRs1dW;dRc3dW;dRs3dW];

end

if output == 'R'
    OUT = R;
elseif output == 'J'
    OUT = [dRdQ,dRdW];
end

end