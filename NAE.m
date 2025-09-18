function dDdt = NAE(t, D, p0, p, phi, psi, nu, rho, xi, delta, eta, gamma)

N = D(1); A = D(2); E = D(3);

% dynamics
dNdt = phi * (p0/p) - (nu*N*A)/(N+gamma) - psi*N*p/p0;
dAdt = (xi*rho*N*A)/(N+gamma) - delta*A;
dEdt = rho*A - eta*E;

dDdt = [dNdt; dAdt; dEdt];
end
