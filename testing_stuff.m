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
    phi_perm = .01*1e-3;
    psi = .01;
    gamma = .01*1e-3; 
    nu = .2; 
    rho = .75; 
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

    [A_sv_init, A_sh_init, ~, ~] = gen_pipes(rmu, rsig, nx, ny, h);

    % set current geometry to initial, then build weights
    A_sv_current = A_sv_init;
    A_sh_current = A_sh_init;

    % ----- Step 1b: baseline permeability, EPS = 0 ----------------------
    dEPS_conc_kgm3 = 0; % incremental addition (kg/m^3)
     % Provisional conductances from current areas (before deposition)
            sv_temp = (A_sv_acc/pi).^2;
            sh_temp = (A_sh_acc/pi).^2;

            % % Flow-capacity weights: C ∝ sv_tmp / h  (constants cancel in normalization)
            W_v = A_sv_acc/sum(A_sv_acc(:));
            W_h = A_sh_acc/sum(A_sh_acc(:));
           

            %current brine volume [m^3]
            vol_v = sum(A_sv_acc(:))*h;
            vol_h = sum(A_sh_acc(:))*h;
            
            %Incremental EPS mass from a concentration step dEPS [kg/m^3]
            dmass_EPS_v = vol_v * dEPS;
            dmass_EPS_h = vol_h * dEPS;
              % Allocate EPS volume per edge: V = (mass / density) * weight
            Veps_v = (dmass_EPS_v / rho_EPS) * W_v;   % m^3 per edge
            Veps_h = (dmass_EPS_h / rho_EPS) * W_h;   % m^3 per edge

              % Convert deposited volume to area decrement: ΔA = V / h
            dA_v = Veps_v / h;   % m^2
            dA_h = Veps_h / h;   % m^2

            % Update areas (clip at 0)
            A_sv_acc = max(0, A_sv_acc - dA_v);
            A_sh_acc = max(0, A_sh_acc - dA_h);

            % Conductances after deposition (used in solver)
            sv = (A_sv_acc/pi).^2; 
            sh = (A_sh_acc/pi).^2;

        
            % Solve using multigrid (kikmul)
            fin = zeros(nx + 1, ny + 1);
            [phi, final_error] = kikmul(fin, pdrop, sv, sh, nx, ny);

            % Effective permeability computation
            effleft = sum((pdrop - phi(nx, :)) .* sh(nx, :)) / pdrop / h^2;
            effright = sum(phi(2, :) .* sh(1, :)) / pdrop / h^2;
            effcoe = (0.5 * (effleft + effright))*(pi/8);
            eff_perm_0 = effcoe;


    % ----- Step 1c: ODE with no EPS to get next state -------------------
    IC0 = [N_0; A_0; 0];
    [~, Y0] = ode23s(@(t,D) NAE(t, D,k_0,rho_EPS, phi_perm, psi, nu, rho, xi, delta, eta, gamma), [0 5], IC0);
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
    eff_perm_vec(1) = perm_0;    % <-- baseline point for plots (p/p0, drop%)                        % or set to your old 'mu' if you used exp(-E/mu)
         % ----- Iterations ----------------------------------------------------
    for iter = 2:max_iters
       tspan = [0 5];
       E_current = EPS_conc_vec(iter);
       IC=[N_vec(iter); A_vec(iter); E_current];

       
       ode_fun = @(t,D)  NAE(t, D,k_0,rho_EPS, phi_perm, psi, nu, rho, xi, delta, eta, gamma), [0 5], IC0);
       [T,Y] = ode23s(ode_fun,tspan,IC)

       N_next = Y(end,1);
       A_next = Y(end,2);
       E_next = Y(end,3);

            sv_temp = (A_sv_new/pi).^2;
            sh_temp = (A_sh_new/pi).^2;

            % % Flow-capacity weights: C ∝ sv_tmp / h  (constants cancel in normalization)
            W_v = A_sv_acc/sum(A_sv_new(:));
            W_h = A_sh_acc/sum(A_sh_new(:));
           

            %current brine volume [m^3]
            vol_v = sum(A_sv_new(:))*h;
            vol_h = sum(A_sh_new(:))*h;
            
            %Incremental EPS mass from a concentration step dEPS [kg/m^3]
            dmass_EPS_v = vol_v * dEPS;
            dmass_EPS_h = vol_h * dEPS;
              % Allocate EPS volume per edge: V = (mass / density) * weight
            Veps_v = (dmass_EPS_v / rho_EPS) * W_v;   % m^3 per edge
            Veps_h = (dmass_EPS_h / rho_EPS) * W_h;   % m^3 per edge

              % Convert deposited volume to area decrement: ΔA = V / h
            dA_v = Veps_v / h;   % m^2
            dA_h = Veps_h / h;   % m^2

            % Update areas (clip at 0)
            A_sv_acc = max(0, A_sv_acc - dA_v);
            A_sh_acc = max(0, A_sh_acc - dA_h);

            % Conductances after deposition (used in solver)
            sv = (A_sv_acc/pi).^2; 
            sh = (A_sh_acc/pi).^2;
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

  
