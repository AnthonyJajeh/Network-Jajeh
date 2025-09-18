  function Dode = NAE_base(t, D, phi, psi, nu, rho, xi, delta, eta, gamma,vf0)
        N = D(1); 
        A = D(2); 
        E = D(3);
    
        dNdt = phi  - (nu * N * A)/(N + gamma) - psi * N ;
        dAdt = (xi * nu * A * N)/(N + gamma) - delta * A;
        dEdt = rho * A - eta * E;
    
        Dode = [dNdt; dAdt; dEdt];
    end 