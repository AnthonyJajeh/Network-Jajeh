    clear all; clc; close all;
    %**************************************************************************
    %
    %     Implementation of the network model to compute the effective properties 
    %     such as conductivity and permeability of a nonhomogeneous material.
    %
    %     The resulting linear system is solved using a multigrid algorithm
    %     (kikmul.m).
    %
    %     The random network coefficients are drawn from a log-normal
    %     distribution specified by parameters (rmu, rsig). To be more precise,
    %     we assume that cross sections are circular disks, and the areas have
    %     a log-normal distribution:
    %
    %             A = e^X, X \sim N(rmu, rsig^2)
    %
    %     We note that
    %             E[A] = e^(E[X]+0.5*Var(X))
    % 
    %     In this code, we only fix the average cross section area E[A] from the 
    %     target volumn fraction, and assume a particular rsig value (entered 
    %     as an input), which will impact the resulting average rmu.
    %
    %     The boundary conditions are set so that u=0 at the top
    %     (index i=1 in this code) and u=pdrop at the bottom
    %     (index i=nx+1).
    %
    %     Last modified: 3/23/2025
    %     This implimentation takes into consideration EPS coming from 
    %     an ODE model
    %*************************************************************************
    % Network model parameters
    nx = 256; %size of grid nxn power of 2
    ny = nx;
    pdrop = 1.d0; %pressure drop
    rsig = .1; %variance 
    trials = 10; %number  of trials
    rho_EPS = 1500; % EPS density [kg/m^3]
    max_iters = 15; 
    vf0 = .05;
    vf0_values = linspace(0.05, 0.25, max_iters);  % Sweep through various volume fractions
    
    %ODE model parameters 
    N_0 = .2;
    A_0 = .03;
    phi = .01;
    psi = .01;
    gamma = .01; 
    nu_1 = .2; 
    nu_2 = .75; 
    xi = .2;
    delta = .007; 
    eta = .01;
    
    % Allocate storage for all vf0 runs
    N_vec = zeros(1, max_iters+1);        % +1 since we fill EPS_conc_vec(i+1)
    A_vec = zeros(1, max_iters+1);
    EPS_conc_vec = zeros(1, max_iters+1); % [mg/L]
    eff_perm_vec = zeros(1, max_iters);   % Store permeability for each iteration
    A_sv_mat = zeros(nx+1, ny+1, max_iters);
    A_sh_mat = zeros(nx+1, ny+1, max_iters);

    N_vec(1)=N_0;
    A_vec(1)=A_0;
    EPS_conc_vec(1)=0;
    % Compute base am0 and h from initial vf0
am0 = pi*(7e-5 + 1.6e-4*vf0)^2;
h = sqrt((2 * am0) / vf0);
rmu = log(am0) - 0.5 * rsig^2;

% Generate pipe areas once
[A_sv_base, A_sh_base, ~, ~] = gen_pipes(rmu, rsig, nx, ny, h);
area_total = sum(A_sv_base(:)) + sum(A_sh_base(:));
d_v = A_sv_base / area_total;
d_h = A_sh_base / area_total;
A_sv = zeros(1,max_iters);
A_sh = zeros(1,max_iters);
% Initialize pipe states and EPS tracking
A_sv_current = A_sv_base;
A_sh_current = A_sh_base;
EPS_prev_kgm3 = 0;
EPS_comulative =0;
        for iter = 1:max_iters
            EPS_conc_mgL = EPS_conc_vec(iter);        % EPS in mg/L
            EPS_now_kgm3 = EPS_conc_mgL *1e-3;
            EPS_prev_kgm3 = EPS_now_kgm3;

            area_total = sum(A_sv_current(:))+sum(A_sh_current(:));
            volume_total = area_total*h;
            mass_EPS_total = volume_total*EPS_now_kgm3;
            volume_EPS_total = mass_EPS_total/rho_EPS;
       
            A_sv_current = max(0, A_sv_current - d_v*volume_EPS_total*h);
            A_sh_current = max(0, A_sh_current - d_h*volume_EPS_total*h);
             
            eff_perm = run_network_model(A_sv_current, A_sh_current, nx, ny, h, pdrop, trials);
            eff_perm_vec(iter)= eff_perm;
            A_sv_mat(:,:,iter) = A_sv_current;
            A_sh_mat(:,:,iter) = A_sh_current;
            mean_area_vec(iter)=mean([A_sv_current(:); A_sh_current(:)]);

            %Solve NAE model wtih current effective permeability
            tspan = [0 10];
            IC = [N_vec(iter); A_vec(iter); EPS_conc_vec(iter)];
            if EPS_conc_mgL ==0
                ode_fun = @(t,D) NAE_base(t, D, phi, psi, nu_1, nu_2, xi, delta, eta, gamma);
                perm_0 = eff_perm;
            else
                ode_fun = @(t,D) NAE(t, D,perm_0, eff_perm, phi, psi, nu_1, nu_2, xi, delta, eta, gamma);
            end
            [T, Y] = ode23s(ode_fun, tspan, IC);
    
            N_vec(iter+1)= Y(end, 1);
            A_vec(iter+1)= Y(end, 2);
            EPS_conc_vec(iter+1) = Y(end,3);
            fprintf("Iteration %d: EPS = %.4f mg/L, eff_perm = %.4e\n", ...
    iter, EPS_conc_mgL, eff_perm);


        end
    % ----------------------------------------------------
    % Plotting
    % ----------------------------------------------------
figure;
plot(EPS_conc_vec(2:end)*1e-3, eff_perm_vec, '-o', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('EPS Concentration (kg/m^3)');
ylabel('Effective Permeability (m^2)');
title('Effective Permeability vs EPS Concentration');
grid on;

figure;
plot(EPS_conc_vec(2:end)*1e-3, mean_area_vec, '-s', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('EPS Concentration (kg/m^3)');
ylabel('Mean Pipe Area (m^2)');
title('Mean Pipe Area vs EPS Concentration');
grid on;

    
    % -------------------------------------------------------------------------
    function eff_perm = run_network_model(A_sv, A_sh, nx, ny, h, pdrop, trials)
    eff_trial = zeros(1, trials);
    for t = 1:trials
        sv = (A_sv/pi).^2;
        sh = (A_sh/pi).^2;
        fin = zeros(nx + 1, ny + 1);
        [phi_grid, ~] = kikmul(fin, pdrop, sv, sh, nx, ny);
        effleft = sum((pdrop - phi_grid(nx,:)) .* sh(nx,:)) / h^2;
        effright = sum(phi_grid(2,:) .* sh(1,:)) / h^2;
        eff_trial(t) = 0.5 * (effleft + effright) * (pi / 8);
    end
    eff_perm = mean(eff_trial);
end
    
    % -------------------------------------------------------------------------
    function Dode = NAE(t, D, perm_0, eff_perm, phi, psi, nu_1, nu_2, xi, delta, eta, gamma)
        N = D(1); 
        A = D(2); 
        E = D(3);
    
        dNdt = phi * eff_perm/perm_0 - (nu_1 * N * A)/(N + gamma) - psi * N * eff_perm/perm_0;
        dAdt = (xi * nu_1 * A * N)/(N + gamma) - delta * A;
        dEdt = nu_2 * A - eta * E;
    
        Dode = [dNdt; dAdt; dEdt];
    end
    
    function Dode = NAE_base(t, D, phi, psi, nu_1, nu_2, xi, delta, eta, gamma)
        N = D(1); 
        A = D(2); 
        E = D(3);
    
        dNdt = phi  - (nu_1 * N * A)/(N + gamma) - psi * N ;
        dAdt = (xi * nu_1 * A * N)/(N + gamma) - delta * A;
        dEdt = nu_2 * A - eta * E;
    
        Dode = [dNdt; dAdt; dEdt];
    end 
    
    % -------------------------------------------------------------------------
    function [A_sv, A_sh, vfrac, amean] = gen_pipes(rmu, rsig, nx, ny, h)
        xmin = rmu - 4 * rsig;
        xmax = rmu + 4 * rsig;
    
        X_sv = max(min(rmu + rsig * randn(nx + 1, ny + 1), xmax), xmin); %Generate normal RV
        X_sh = max(min(rmu + rsig * randn(nx + 1, ny + 1), xmax), xmin); %Generate normal RV
    
        A_sv = exp(X_sv); %Generate lognormal RV for area
        A_sh = exp(X_sh); %Generate lognormal RV for area
        asum = sum(A_sv(:)) + sum(A_sh(:));
        amean = 0.5 * asum / ((nx + 1) * (ny + 1)); %mean area
        vfrac = asum / ((nx + 1) * (ny + 1)) / h^2; %volume fraction
    end