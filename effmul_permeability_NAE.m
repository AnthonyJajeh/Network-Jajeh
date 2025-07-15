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
    max_iters = 5; 
    vf0_values = linspace(0.05, 0.25, max_iters);  % Sweep through various volume fractions
    
    %ODE model parameters 
    N_0 = 5;
    A_0 = .03;
    phi = .1;
    psi = .05;
    gamma = .01; 
    nu_1 = .2; 
    nu_2 = .05; 
    xi = .2;
    delta = .007; 
    eta = .03;
    
    % Allocate storage for all vf0 runs
    all_results = cell(length(vf0_values), 1);
    for vf_index = 1:length(vf0_values)
        eps_mass_total_vec = zeros(max_iters, 1);   % EPS mass [kg]
        eps_volume_total_vec = zeros(max_iters, 1); % EPS volume [m^3]
        EPS_conc_vec = zeros(max_iters+1,1); %store space for EPS from ODE
        eff_perm_vec = zeros(max_iters,1); %store space for effective perm
        sv_mat = zeros(nx+1, ny+1, max_iters); %storage for R^4 v
        sh_mat = zeros(nx+1, ny+1, max_iters); %storage for R^4 h
        A_sv_mat = zeros(nx+1, ny+1, max_iters); %storage for areas v
        A_sh_mat = zeros(nx+1, ny+1, max_iters); %storage for areas h
        EPS_conc_vec(1) = 0;  % Initial EPS concentration
        N_vec = zeros(max_iters+1,1); %storage for initial condition of nutrients
        A_vec = zeros(max_iters+1,1); % storage for initial condition of algae
        N_vec(1) = N_0;
        A_vec(1) = A_0;
        vf0 = vf0_values(vf_index);  % Step 1: new vf0
        
        % Step 1a: Generate initial network geometry (depends on vf0)
        am0 = pi*(7e-5 + 1.6e-4*vf0)^2;
        rmu = log(am0) - 0.5 * rsig^2;
        h = sqrt((2 * am0) / vf0);
        [A_sv_init, A_sh_init, ~, ~] = gen_pipes(rmu, rsig, nx, ny, h);
        d_v =.5* A_sv_init / sum(A_sv_init(:));
        d_h = .5*A_sh_init / sum(A_sh_init(:));
%%%
% 
% vf0=0.05;
%      % Step 1a: Generate initial network geometry (depends on vf0)
%         am0 = pi*(7e-5 + 1.6e-4*vf0)^2
%         rmu = log(am0) - 0.5 * rsig^2
%         h = sqrt((2 * am0) / vf0)
%         [A_sv_init, A_sh_init, ~, ~] = gen_pipes(rmu, rsig, nx, ny, h);
%         d_v = A_sv_init / sum(A_sv_init(:));
%         d_h = A_sh_init / sum(A_sh_init(:));
% 
%         mean(d_h,'all')
% 
%         vf0=0.25;
%      % Step 1a: Generate initial network geometry (depends on vf0)
%         am0 = pi*(7e-5 + 1.6e-4*vf0)^2
%         rmu = log(am0) - 0.5 * rsig^2
%         h = sqrt((2 * am0) / vf0)
%         [A_sv_init, A_sh_init, ~, ~] = gen_pipes(rmu, rsig, nx, ny, h);
%         d_v = A_sv_init / sum(A_sv_init(:));
%         d_h = A_sh_init / sum(A_sh_init(:));
% 
%         mean(d_h,'all')
    
        % Step 1b: Compute perm_0 using network model with EPS = 0
        EPS_conc_kgm3 = 0;
        [eff_perm_0, sv, sh, A_sv_current, A_sh_current] = run_network_model(...
            EPS_conc_kgm3, nx, ny, rho_EPS, trials, pdrop, ...
            A_sv_init, A_sh_init, d_v, d_h, h);
        perm_0 = eff_perm_0;  % This is the baseline

        % Step 1c: Run ODE with no EPS using NAE_base
        IC = [N_0; A_0; 0];
        [~, Y] = ode45(@(t,D) NAE_base(t, D, phi, psi, nu_1, nu_2, ...
                        xi, delta, eta, gamma,vf0), [0 250], IC);
        N_vec(2) = Y(end,1);
        A_vec(2) = Y(end,2);
        EPS_conc_vec(2) = Y(end,3);
            
            % Store initial network outputs
            eff_perm_vec(1) = perm_0;
            sv_mat(:,:,1) = sv;
            sh_mat(:,:,1) = sh;
            A_sv_mat(:,:,1) = A_sv_current;
            A_sh_mat(:,:,1) = A_sh_current;
            
        for iter = 2:max_iters
           EPS_conc_mgL = EPS_conc_vec(iter); % This is the total EPS concentration so far
           EPS_conc_kgm3 = EPS_conc_mgL * 1e-3;

            % Run network model using cumulative EPS volume
            [eff_perm, sv, sh, A_sv_new, A_sh_new] = run_network_model(EPS_conc_kgm3, nx, ny, rho_EPS, trials,pdrop, A_sv_current, A_sh_current, d_v, d_h, h);
            fprintf('vf0 %g: iteration %d: cumulative EPS = %.4e kg/m^3, eff_perm = %.4e m^2, eff_perm/perm %g:\n', ...
        vf0, iter, EPS_conc_kgm3, eff_perm, eff_perm/perm_0);
            eff_perm_vec(iter) = eff_perm;
            sv_mat(:,:,iter) = sv;
            sh_mat(:,:,iter) = sh;
            A_sv_mat(:,:,iter) = A_sv_new;
            A_sh_mat(:,:,iter) = A_sh_new;

            % Update pipe areas for next iteration
            A_sv_current = A_sv_new;
            A_sh_current = A_sh_new;
             % After you get new A_sv_current and A_sh_current:
            d_v = A_sv_current / sum(A_sv_current(:));
            d_h = A_sh_current / sum(A_sh_current(:));
            display(perm_0)
            %Solve NAE model wtih current effective permeability
            tspan = [0 250];
            IC = [N_vec(iter); A_vec(iter); EPS_conc_mgL];
                ode_fun = @(t,D) NAE(t, D,perm_0, eff_perm, phi, psi, nu_1, nu_2, xi, delta, eta, gamma);

            [T, Y] = ode45(ode_fun, tspan, IC);

            N_vec(iter+1)= Y(end, 1);
            A_vec(iter+1)= Y(end, 2);
            EPS_conc_vec(iter+1) =  Y(end,3);  % If Y(end,3) is dE/dt

        end
           all_results{vf_index} = struct('vf0', vf0, 'EPS_conc_vec', EPS_conc_vec, 'eff_perm_vec', eff_perm_vec, 'N_vec', N_vec, 'A_vec', A_vec, 'sv_mat', sv_mat, 'sh_mat', sh_mat,'A_sv_mat', A_sv_mat, 'A_sh_mat', A_sh_mat,'eps_mass_total_vec', eps_mass_total_vec,'eps_volume_total_vec', eps_volume_total_vec);
    end
    % Convert cell array to struct array to fix dot-indexing issues
    all_results_struct = [all_results{:}];
    % ----------------------------------------------------
    % Plotting
    % ----------------------------------------------------

    figure;
hold on;
for i = 1:length(vf0_values)
    plot(1:max_iters, all_results_struct(i).eps_mass_total_vec, '-o', ...
        'LineWidth', 2, 'DisplayName', sprintf('\\phi_0 = %.2f', vf0_values(i)));
end
xlabel('Iteration');
ylabel('Cumulative EPS Mass (kg)');
title('Cumulative EPS Mass per Iteration');
legend('show', 'Location', 'bestoutside');
grid on;
    % 1. Extract EPS and permeability data (excluding the initial EPS=0 point)
EPS_matrix = zeros(length(vf0_values), max_iters);
perm_matrix = zeros(length(vf0_values), max_iters);

for i = 1:length(vf0_values)
    EPS_matrix(i,:) = all_results_struct(i).EPS_conc_vec(2:end);  % skip initial EPS=0
    perm_matrix(i,:) = all_results_struct(i).eff_perm_vec;       % same length
end

%Convert EPS to kg/m³
EPS_matrix = EPS_matrix * 1e-3;

% Plot 1: Effective Permeability vs EPS Concentration
figure;
hold on;
colors = lines(length(vf0_values));

for i = 1:length(vf0_values)
    scatter(EPS_matrix(i,:), perm_matrix(i,:), 80, ...
        'filled', 'MarkerEdgeColor', 'k', ...
        'MarkerFaceColor', colors(i,:));
end

xlabel('EPS Concentration (kg/m^3)');
ylabel('Effective Permeability (m^2)');
title('Effective Permeability vs EPS Concentration (Colored by \phi_0)');
legend(arrayfun(@(v) sprintf('\\phi_0 = %.2f', v), vf0_values, 'UniformOutput', false), ...
       'Location', 'bestoutside');
grid on;

% 2. Compute Mean Pipe Area (same dimensions as EPS and perm)
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

% Plot 2: Mean Pipe Area vs EPS Concentration
figure;
hold on;
for i = 1:length(vf0_values)
    scatter(EPS_matrix(i,:), mean_area_matrix(i,:), 80, ...
        'filled', 'MarkerEdgeColor', 'k', ...
        'MarkerFaceColor', colors(i,:));
end

xlabel('EPS Concentration (kg/m^3)');
ylabel('Mean Pipe Area (m^2)');
title('Mean Pipe Area vs EPS Concentration (Colored by \phi_0)');
legend(arrayfun(@(v) sprintf('\\phi_0 = %.2f', v), vf0_values, 'UniformOutput', false), ...
       'Location', 'bestoutside');
grid on;
 

    m =length(vf0_values);
    colors = lines(m); % Create a colormap for m lines
    figure;
    hold on;
    for j = 1:m
        % EPS concentrations and permeabilities for the j-th volume fraction
        eps_row = EPS_matrix(j, :);
        perm_row = perm_matrix(j, :);
    
        plot(eps_row, perm_row, '-o', ...
            'LineWidth', 2, ...
            'Color', colors(j,:), ...
            'DisplayName', sprintf('\\phi_0 = %.2f', vf0_values(j)));
    end
    
    xlabel('EPS Concentration (kg/m^3)');
    ylabel('Effective Permeability (m^2)');
    title('EPS Concentration vs Effective Permeability');
    legend('show', 'Location', 'bestoutside');
    grid on;
    
    figure;
    hold on;
    for j = 1:m
        mean_area_row = mean_area_matrix(j, :);
        perm_row = perm_matrix(j, :);
    
        plot(mean_area_row, perm_row, '-s', ...
            'LineWidth', 2, ...
            'Color', colors(j,:), ...
            'DisplayName', sprintf('\\phi_0 = %.2f', vf0_values(j)));
    end
    
    xlabel('Mean Pipe Area');
    ylabel('Effective Permeability (m^2)');
    title('Mean Pipe Area vs Effective Permeability');
    legend('show', 'Location', 'bestoutside');
    grid on;
    % 6. EPS Concentration vs Mean Area
    figure;
    hold on;
    
    for j = 1:m
        eps_row = EPS_matrix(j, :);
        mean_area_row = mean_area_matrix(j, :);
    
        plot(eps_row, mean_area_row, '-^', ...
            'LineWidth', 2, ...
            'Color', colors(j,:), ...
            'DisplayName', sprintf('\\phi_0 = %.2f', vf0_values(j)));
    end
    
    xlabel('EPS Concentration (kg/m^3)');
    ylabel('Mean Pipe Area');
    title('EPS Concentration vs Mean Pipe Area');
    legend('show', 'Location', 'bestoutside');
    grid on;
    
    % 7. For each EPS concentration, plot vf0 vs effective permeability
    % Prepare dimensions
    [m, n] = size(perm_matrix); % m = # vf0 values, n = # EPS levels (i.e., max_iters)
    colors = lines(n);          % n different EPS levels (columns)
    markerStyles = {'x', '+', 'o', 's', 'd', '^', 'v', '>', '<', 'p', 'h', '*', '.', 'x', '+'};
    
    figure;
    hold on;
    
    for i = 1:n
        marker = markerStyles{mod(i-1, length(markerStyles)) + 1};
        eff_col = perm_matrix(:, i);        % fixed EPS column
        eps_label = EPS_matrix(1, i);       % use EPS value from first row (they vary across rows)
    
        plot(vf0_values, eff_col, ...
            'LineStyle', 'none', ...
            'Marker', marker, ...
            'MarkerSize', 8, ...
            'LineWidth', 1.5, ...
            'Color', colors(i,:), ...
            'DisplayName', sprintf('EPS = %.4f kg/m^3', eps_label));
    end
    
    set(gca, 'YScale', 'log');
    xlabel('Brine Volume Fraction, \phi_0');
    ylabel('Effective Permeability, k (m^2)');
    title('Effective Permeability vs Brine Volume Fraction');
    legend('show', 'Location', 'bestoutside');
    grid on;
    hold off;

    % Compute percent drop in permeability for each vf0
perm_drop_percent = zeros(length(vf0_values),1);

for i = 1:length(vf0_values)
    perm_start = all_results_struct(i).eff_perm_vec(1);     % after first network call
    perm_end = all_results_struct(i).eff_perm_vec(end);     % after final iteration
    perm_drop_percent(i) = 100 * (perm_start - perm_end) / perm_start;
end

% Plot
figure;
plot(vf0_values, perm_drop_percent, '-o', 'LineWidth', 2);
xlabel('Brine Volume Fraction \phi_0');
ylabel('Percent Drop in Effective Permeability (%)');
title('EPS-Induced Drop in Permeability vs \phi_0');
grid on;

    % Compute percent drop in mean pipe area for each vf0
mean_area_drop_percent = zeros(length(vf0_values),1);

for i = 1:length(vf0_values)
    A_start = mean_area_matrix(i,1);   % first EPS
    A_end   = mean_area_matrix(i,end); % last EPS
    mean_area_drop_percent(i) = 100 * (A_start - A_end) / A_start;
end

% Plot
figure;
plot(vf0_values, mean_area_drop_percent, '-s', 'LineWidth', 2);
xlabel('Brine Volume Fraction \phi_0');
ylabel('Percent Drop in Mean Pipe Area (%)');
title('EPS-Induced Drop in Pipe Area vs \phi_0');
grid on;

figure;
hold on;
colors = lines(length(vf0_values));
for i = 1:length(vf0_values)
    plot(0:max_iters, all_results_struct(i).EPS_conc_vec, '-o', ...
        'LineWidth', 2, 'Color', colors(i,:), ...
        'DisplayName', sprintf('\\phi_0 = %.2f', vf0_values(i)));
end
xlabel('Iteration');
ylabel('EPS Concentration (mg/L)');
title('EPS Accumulation per Iteration for each \phi_0');
legend('show', 'Location', 'bestoutside');
grid on;
    
  
