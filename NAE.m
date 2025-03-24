function [Dode] = NAE(I,D,eff,a,b,c,c_p,d)
%I- indepenedent variable
%D - dependent variable


% naming the ode values I want
N = D(1);
A = D(2);
E = D(3);

%set of odes
dNdt =a*eff-(c*A*N)/(N+1)-b*N*eff;
dAdt = (c_p*N*A)/(1 + N) - A;
dEdt = d*A - E;

% odes in vector form
Dode = [dNdt; dAdt; dEdt];
end