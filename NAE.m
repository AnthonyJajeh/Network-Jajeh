function Dode = NAE(t, D, perm_0, eff_perm, phi, psi, nu_1, nu_2, xi, delta, eta, gamma)
        N = D(1); 
        A = D(2); 
        E = D(3);
    
        dNdt = phi * (eff_perm/perm_0) - (nu_1 * N * A)/(N + gamma) - psi * N * (eff_perm/perm_0);
        dAdt = (xi * nu_1 * A * N)/(N + gamma) - delta * A;
        dEdt = nu_2 * A *(eff_perm/perm_0) - eta * E;
    
        Dode = [dNdt; dAdt; dEdt];
    end