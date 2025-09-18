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
    trials = 5; %number  of trials
    rho_EPS = 1500; % EPS density [kg/m^3]
    max_iters = 15; 
    vf0_values = linspace(0.05, 0.25, trials);  % Sweep through various volume fractions
    
    %ODE model parameters 
    N_0 = .2* 1e-3;
    A_0 = .03* 1e-3;
    phi = .01*1e-3;
    psi = .01;
    gamma = .01*1e-3; 
    nu_1 = .2; 
    nu_2 = .75; 
    xi = .2;
    delta = .007; 
    eta = .03;
    
    % Allocate storage for all vf0 runs
    all_results = cell(length(vf0_values), 1);
    
parfor vf_index = 1:length(vf0_values)
    eps_mass_total_vec   = zeros(max_iters, 1);   % kg
    eps_volume_total_vec = zeros(max_iters, 1);   % m^3
    EPS_conc_vec = zeros(max_iters+1,1);         % mg/L, EPS concentration state from ODE
    eff_perm_vec = zeros(max_iters,1);           % m^2
    sv_mat  = zeros(nx+1, ny+1, max_iters);      % conductance r^4 (v)
    sh_mat  = zeros(nx+1, ny+1, max_iters);      % conductance r^4 (h)
    A_sv_mat = zeros(nx+1, ny+1, max_iters);     % areas (v)
    A_sh_mat = zeros(nx+1, ny+1, max_iters);     % areas (h)

    EPS_conc_vec(1) = 0;     % initial EPS conc (mg/L)
    N_vec = zeros(max_iters+1,1); A_vec = zeros(max_iters+1,1);
    N_vec(1) = N_0; A_vec(1) = A_0;

    vf0 = vf0_values(vf_index);

    % ----- Step 1a: initial network geometry from vf0 -------------------
    am0 = pi*(7e-5 + 1.6e-4*vf0)^2;
    rmu = log(am0) - 0.5 * rsig^2;
    h   = sqrt((2 * am0) / vf0);
    Lx = nx * h;
    Ly = ny * h;
    V_pore = Lx * Ly * h * vf0;   % total brine (pore) volume for this run

    [A_sv_init, A_sh_init, ~, ~] = gen_pipes(rmu, rsig, nx, ny, h);

    % set current geometry to initial, then build weights
    A_sv_current = A_sv_init;
    A_sh_current = A_sh_init;

    % ----- Step 1b: baseline permeability, EPS = 0 ----------------------
    dEPS_conc_kgm3 = 0; % incremental addition (kg/m^3)
    [eff_perm_0, sv, sh, A_sv_current, A_sh_current] = run_network_model( ...
    0, nx, ny, rho_EPS, trials, pdrop, ...
    A_sv_current, A_sh_current, [], [], h, V_pore, 0);

    % ----- Step 1c: ODE with no EPS to get next state -------------------
    IC0 = [N_0; A_0; 0];
    [~, Y0] = ode23s(@(t,D) NAE_base(t, D, phi, psi, nu_1, nu_2, xi, delta, eta, gamma, vf0), [0 250], IC0);
    N_vec(2)       = Y0(end,1);
    A_vec(2)       = Y0(end,2);
    EPS_conc_vec(2)= Y0(end,3);     % mg/L

    % store initial outputs
    sv_mat(:,:,1)   = sv;
    sh_mat(:,:,1)   = sh;
    A_sv_mat(:,:,1) = A_sv_current;
    A_sh_mat(:,:,1) = A_sh_current;
    perm_0 = eff_perm_0;
    eff_perm_current = perm_0;   % <-- use this for iter 2 ODE
    eff_perm_vec(1) = perm_0;    % <-- baseline point for plots (p/p0, drop%)
    alphaT = 0.7; 
    beta_mu = 0;
    Eref = 1e-6;   % T knobs (beta_mu=0 means no extra viscosity term)
    w_phi  = 0.7; 
    w_nu  = 0.6;             % how transport-limited phi, nu1 are
    theta0 = 0.25; 
    theta_slope = 0.6;      % deposition fraction base & extra when clogged
    mu_opt = .001*1e-3;                            % or set to your old 'mu' if you used exp(-E/mu)
         % ----- Iterations ----------------------------------------------------
    for iter = 2:max_iters
       tspan = [0 5];
       E_current = EPS_conc_vec(iter);
       IC=[N_vec(iter); A_vec(iter); E_current];

       
       ode_fun = @(t,D) NAE(t, D, perm_0, eff_perm_current, ...
            phi, psi, nu_1, nu_2, xi, delta, eta, gamma, ...
            alphaT, beta_mu, Eref, w_phi, w_nu, mu_opt,vf0);
       [T,Y] = ode23s(ode_fun,tspan,IC)

       N_next = Y(end,1);
       A_next = Y(end,2);
       E_next = Y(end,3);


       % ---- 2) Compute EPS PRODUCTION over this step (mg/L) ----
       % From dE/dt = rho*A - eta*E  =>  ∫rho A dt = (E_next - E_i) + ∫eta E dt
       prod_kgm3 = (E_next - E_current) + trapz(T, eta * Y(:,3));
       prod_kgm3 = max(prod_kgm3, 0);                  % never let production be negative

       % make deposition stickier as network clogs via T based on permeability
       clog_base = 1.0;
       T_now = min(1, max(0, (eff_perm_current / max(perm_0, eps))^alphaT ));
       theta_dep_eff = min(1, theta0 + theta_slope*(1 - T_now));
       fixed_clog = clog_base * (1 + 0.5*(1 - T_now));
       fixed_clog = min(max(fixed_clog, 0), 1e12);

     
         % Deposit a fraction of production
       dEPS_conc_kgm3 = theta_dep_eff* prod_kgm3        % kg/m^3


       [eff_perm_next, sv, sh, A_sv_new, A_sh_new] = run_network_model( ...
        dEPS_conc_kgm3, nx, ny, rho_EPS, trials, pdrop, ...
        A_sv_current, A_sh_current, [], [], h, V_pore, fixed_clog);

    eff_perm_vec(iter) = eff_perm_next;
    sv_mat(:,:,iter)   = sv;
    sh_mat(:,:,iter)   = sh;
    A_sv_mat(:,:,iter) = A_sv_new;
    A_sh_mat(:,:,iter) = A_sh_new;

    % Update geometry & ODE state for next step
    A_sv_current      = A_sv_new;
    A_sh_current      = A_sh_new;
    eff_perm_current  = eff_perm_next;           % <-- use updated p next step
    N_vec(iter+1)        = N_next;
    A_vec(iter+1)        = A_next;
    EPS_conc_vec(iter+1) = E_next;               % fluid EPS (can go up or down)
     % deposited concentration this step (kg/m^3) and deposited mass (kg)
E_dep_step_kg = dEPS_conc_kgm3 * V_pore;

fprintf(['vf0 %.3f | iter %2d | dEPS_dep = %.3e kg (%.3e kg/m^3) | ' ...
         'E_fluid = %.3f kg/m^3 | p = %.4e m^2 | p/p0 = %.3f\n'], ...
        vf0, iter, E_dep_step_kg, dEPS_conc_kgm3, E_next, eff_perm_next, eff_perm_next/perm_0);
    end

    all_results{vf_index} = struct( ...
        'vf0', vf0, ...
        'EPS_conc_vec', EPS_conc_vec, ...
        'eff_perm_vec', eff_perm_vec, ...
        'N_vec', N_vec, ...
        'A_vec', A_vec, ...
        'sv_mat', sv_mat, 'sh_mat', sh_mat, ...
        'A_sv_mat', A_sv_mat, 'A_sh_mat', A_sh_mat, ...
        'eps_mass_total_vec', eps_mass_total_vec, ...
        'eps_volume_total_vec', eps_volume_total_vec);
end

% Struct array for easier indexing
all_results_struct = [all_results{:}];

% ==== build EPS_matrix (kg/m^3) and perm_matrix (m^2) ====
EPS_matrix  = zeros(length(vf0_values), max_iters);
perm_matrix = zeros(length(vf0_values), max_iters);

for i = 1:length(vf0_values)
    EPS_matrix(i,:)  = all_results_struct(i).EPS_conc_vec(2:end);  % mg/L, length=max_iters
    perm_matrix(i,:) = all_results_struct(i).eff_perm_vec;         % m^2,   length=max_iters
end
EPS_matrix = EPS_matrix * 1e-3;   % mg/L -> kg/m^3

% Plot: perm vs EPS
figure; hold on; colors = lines(length(vf0_values));
for i = 1:length(vf0_values)
    scatter(EPS_matrix(i,:), perm_matrix(i,:), 80, 'filled', ...
        'MarkerEdgeColor','k', 'MarkerFaceColor', colors(i,:));
end
xlabel('EPS Concentration (kg/m^3)'); ylabel('Effective Permeability (m^2)');
title('Effective Permeability vs EPS Concentration (Colored by \phi_0)');
legend(arrayfun(@(v) sprintf('\\phi_0 = %.2f', v), vf0_values, 'UniformOutput', false), 'Location','bestoutside');
grid on;

% 2. Mean pipe area per iter
mean_area_matrix = zeros(length(vf0_values), max_iters);
for i = 1:length(vf0_values)
    result = all_results_struct(i);
    for j = 1:max_iters
        A_sv_j = result.A_sv_mat(:,:,j);
        A_sh_j = result.A_sh_mat(:,:,j);
        A_vals = [A_sv_j(:); A_sh_j(:)];
        mean_area_matrix(i,j) = mean(A_vals);
    end
end

% Plot: mean area vs EPS
figure; hold on;
for i = 1:length(vf0_values)
    scatter(EPS_matrix(i,:), mean_area_matrix(i,:), 80, 'filled', ...
        'MarkerEdgeColor','k', 'MarkerFaceColor', colors(i,:));
end
xlabel('EPS Concentration (kg/m^3)'); ylabel('Mean Pipe Area (m^2)');
title('Mean Pipe Area vs EPS Concentration (Colored by \phi_0)');
legend(arrayfun(@(v) sprintf('\\phi_0 = %.2f', v), vf0_values, 'UniformOutput', false), 'Location','bestoutside');
grid on;

% EPS vs k (lines)
m = length(vf0_values); colors = lines(m);
figure; hold on;
for j = 1:m
    plot(EPS_matrix(j,:), perm_matrix(j,:), '-o', 'LineWidth', 2, ...
        'Color', colors(j,:), 'DisplayName', sprintf('\\phi_0 = %.2f', vf0_values(j)));
end
xlabel('EPS Concentration (kg/m^3)'); ylabel('Effective Permeability (m^2)');
title('EPS Concentration vs Effective Permeability'); legend('show','Location','bestoutside'); grid on;

% Mean area vs k (lines)
figure; hold on;
for j = 1:m
    plot(mean_area_matrix(j,:), perm_matrix(j,:), '-s', 'LineWidth', 2, ...
        'Color', colors(j,:), 'DisplayName', sprintf('\\phi_0 = %.2f', vf0_values(j)));
end
xlabel('Mean Pipe Area (m^2)'); ylabel('Effective Permeability (m^2)');
title('Mean Pipe Area vs Effective Permeability'); legend('show','Location','bestoutside'); grid on;

% EPS vs mean area (lines)
figure; hold on;
for j = 1:m
    plot(EPS_matrix(j,:), mean_area_matrix(j,:), '-^', 'LineWidth', 2, ...
        'Color', colors(j,:), 'DisplayName', sprintf('\\phi_0 = %.2f', vf0_values(j)));
end
xlabel('EPS Concentration (kg/m^3)'); ylabel('Mean Pipe Area (m^2)');
title('EPS Concentration vs Mean Pipe Area'); legend('show','Location','bestoutside'); grid on;

% For each EPS level, plot vf0 vs k
[mR, nC] = size(perm_matrix); colors = lines(nC);
markerStyles = {'x','+','o','s','d','^','v','>','<','p','h','*','.','x','+'};
figure; hold on;
for i = 1:nC
    marker = markerStyles{mod(i-1, numel(markerStyles)) + 1};
    eff_col = perm_matrix(:, i);
    eps_label = EPS_matrix(1, i);
    plot(vf0_values, eff_col, 'LineStyle','none', 'Marker', marker, ...
        'MarkerSize', 8, 'LineWidth', 1.5, 'Color', colors(i,:), ...
        'DisplayName', sprintf('EPS = %.4f kg/m^3', eps_label));
end
set(gca,'YScale','log');
xlabel('Brine Volume Fraction, \phi_0'); ylabel('Effective Permeability, k (m^2)');
title('Effective Permeability vs Brine Volume Fraction'); legend('show','Location','bestoutside'); grid on;

% Percent drops
perm_drop_percent = zeros(length(vf0_values),1);
for i = 1:length(vf0_values)
    perm_start = all_results_struct(i).eff_perm_vec(1);
    perm_end   = all_results_struct(i).eff_perm_vec(end);
    perm_drop_percent(i) = 100 * (perm_start - perm_end) / max(perm_start, eps);
end
figure; plot(vf0_values, perm_drop_percent, '-o', 'LineWidth', 2);
xlabel('Brine Volume Fraction \phi_0'); ylabel('Percent Drop in Effective Permeability (%)');
title('EPS-Induced Drop in Permeability vs \phi_0'); grid on;

mean_area_drop_percent = zeros(length(vf0_values),1);
for i = 1:length(vf0_values)
    A_start = mean_area_matrix(i,1);
    A_end   = mean_area_matrix(i,end);
    mean_area_drop_percent(i) = 100 * (A_start - A_end) / max(A_start, eps);
end
figure; plot(vf0_values, mean_area_drop_percent, '-s', 'LineWidth', 2);
xlabel('Brine Volume Fraction \phi_0'); ylabel('Percent Drop in Mean Pipe Area (%)');
title('EPS-Induced Drop in Pipe Area vs \phi_0'); grid on;

figure; hold on; colors = lines(length(vf0_values));
for i = 1:length(vf0_values)
    plot(0:max_iters, all_results_struct(i).EPS_conc_vec, '-o', ...
        'LineWidth', 2, 'Color', colors(i,:), ...
        'DisplayName', sprintf('\\phi_0 = %.2f', vf0_values(i)));
end
xlabel('Iteration'); ylabel('EPS Concentration (mg/L)');
title('EPS Accumulation per Iteration for each \phi_0');
legend('show','Location','bestoutside'); grid on;

  
