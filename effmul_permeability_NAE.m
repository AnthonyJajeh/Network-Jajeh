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
    N_0 = .2;
    A_0 = .03;
    phi = .01;
    psi = .01;
    gamma = .01; 
    nu = .2; 
    rho = .75; 
    xi = .2;
    delta = .007; 
    eta = .03;
    
    % Allocate storage for all vf0 runs
    all_results = cell(length(vf0_values), 1);
    
% Allocate storage for all vf0 runs
all_results = cell(length(vf0_values), 1);

parfor vf_index = 1:length(vf0_values)
    eps_mass_total_vec   = zeros(max_iters, 1);   % kg
    eps_volume_total_vec = zeros(max_iters, 1);   % m^3
    EPS_conc_vec = zeros(max_iters+1,1);         % mg/L (fluid EPS)
    eff_perm_vec = zeros(max_iters,1);           % m^2
    sv_mat  = zeros(nx+1, ny+1, max_iters);      % conductance r^4 (v)
    sh_mat  = zeros(nx+1, ny+1, max_iters);      % conductance r^4 (h)
    A_sv_mat = zeros(nx+1, ny+1, max_iters);     % areas (v)
    A_sh_mat = zeros(nx+1, ny+1, max_iters);     % areas (h)

    EPS_conc_vec(1) = 0;     % initial fluid EPS (mg/L)
    N_vec = zeros(max_iters+1,1); 
    A_vec = zeros(max_iters+1,1);
    N_vec(1) = N_0; 
    A_vec(1) = A_0;

    vf0 = vf0_values(vf_index);

    % ----- Step 1a: initial network geometry from vf0 -------------------
    am0 = pi*(7e-5 + 1.6e-4*vf0)^2;
    rmu = log(am0) - 0.5 * rsig^2;
    h   = sqrt((2 * am0) / vf0);
    [A_sv, A_sh, ~, ~] = gen_pipes(rmu, rsig, nx, ny, h);

    % ----- Step 1b: baseline permeability, EPS = 0 ----------------------
    dEPS_conc_kgm3 = 0;  % (kg/m^3) deposit this call
    [eff_perm_0, sv, sh, A_sv, A_sh] = run_network_model( ...
        dEPS_conc_kgm3, nx, ny, rho_EPS, trials, pdrop, A_sv, A_sh, h);


      % ----- Step 1c: ODE with no EPS to get next state -------------------
    IC0 = [N_0; A_0; 0];
    [~, Y0] = ode23s(@(t,D) NAE_base(t, D, phi, psi, nu, rho, xi, delta, eta, gamma, vf0), [0 250], IC0);
    N_vec(2)        = Y0(end,1);
    A_vec(2)        = Y0(end,2);
    EPS_conc_vec(2) = Y0(end,3);   % mg/L

    % store initial outputs
    eff_perm_vec(1) = eff_perm_0;
    sv_mat(:,:,1)   = sv;
    sh_mat(:,:,1)   = sh;
    A_sv_mat(:,:,1) = A_sv;
    A_sh_mat(:,:,1) = A_sh;

         % ----- Iterations ----------------------------------------------------
   for iter = 2:max_iters
        % incremental EPS concentration (mg/L → kg/m^3)
        dEPS_conc_mgL  = max(EPS_conc_vec(iter) - EPS_conc_vec(iter-1), 0);
        dEPS_conc_kgm3 = dEPS_conc_mgL * 1e-3;

        % Run network model with incremental EPS
        [eff_perm, sv, sh, A_sv_new, A_sh_new] = run_network_model( ...
            dEPS_conc_kgm3, nx, ny, rho_EPS, trials, pdrop, A_sv, A_sh, h);

        % Diagnostics
        cumEPS_conc_kgm3 = EPS_conc_vec(iter)*1e-3;  % current fluid EPS (kg/m^3)
        fprintf(['vf0 %.3f | iter %2d | dEPS = %.3e kg/m^3 | ' ...
                 'cumEPS = %.3e kg/m^3 | p = %.4e m^2 | p/p0 = %.3f\n'], ...
                 vf0, iter, dEPS_conc_kgm3, cumEPS_conc_kgm3, eff_perm, eff_perm/eff_perm_0);

        % Save fields
        eff_perm_vec(iter) = eff_perm;
        sv_mat(:,:,iter)   = sv;
        sh_mat(:,:,iter)   = sh;
        A_sv_mat(:,:,iter) = A_sv_new;
        A_sh_mat(:,:,iter) = A_sh_new;

        % Update geometry for next iteration  **(bug fix here)**
        A_sv = A_sv_new;
        A_sh = A_sh_new;

        % Advance ODE with current geometry-induced permeability
        tspan = [0 5];
        IC = [N_vec(iter); A_vec(iter); EPS_conc_vec(iter)];
        % Make sure this matches your NAE signature
        % Option A (short signature):
        ode_fun = @(t,D) NAE(t, D, eff_perm_0, eff_perm, phi, psi, nu, rho, xi, delta, eta, gamma);
        % Option B (if you need vf0 as well):
        % ode_fun = @(t,D) NAE(t, D, eff_perm_0, eff_perm, phi, psi, nu, rho, xi, delta, eta, gamma, vf0);

        [~, Y] = ode23s(ode_fun, tspan, IC);
        N_vec(iter+1)        = Y(end, 1);
        A_vec(iter+1)        = Y(end, 2);
        EPS_conc_vec(iter+1) = Y(end, 3);
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

  
