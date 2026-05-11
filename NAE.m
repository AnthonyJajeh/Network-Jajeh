function dDdt = NAE(t, D, k_0,rho_EPS, phi, psi, nu, rho, xi, delta, eta, gamma)

N = D(1); A = D(2); E = D(3);

% dynamics
dNdt = phi * k_0*exp(-2*E/rho_EPS) - (nu*N*A)/(N+gamma) - psi*N*k_0*exp(-2*E/rho_EPS);
dAdt = (xi*rho*N*A)/(N+gamma) - delta*A;
dEdt = rho*A - eta*E;

dDdt = [dNdt; dAdt; dEdt];
end
