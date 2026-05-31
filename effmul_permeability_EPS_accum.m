clear all; clc; close all;
t_all = tic;
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
%     This implimentation takes into consideration EPS as some base value
%     rather then coming from an ODE.
%*************************************************************************
baseSeed = 12345;

nx = 256; %size of grid nxn power of 2
ny = nx;
pdrop = 1.d0; %pressure drop
rsig = .25; %variance of lognormal distribution
trials = 5; %number  of trials
mu_fluid = 1.8e-3;   % viscosity [Pa*s]
n = 20;%number of deposit steps;
m = 5;

eff_mean_perm_plot = zeros(m, n); % mean effective permeability
eff_trial_all_perm = zeros(m, n, trials); % all permeability trials
eff_resistance_perm_plot = zeros(m, n);
eff_resistance_all_perm  = zeros(m, n, trials);
mean_area_plot = zeros(m, n); % mean pipe area
mean_area_all = zeros(m, n, trials); % all mean area trials
mass_EPS_total_v = zeros(m, n); % total mass of EPS
mass_EPS_total_h = zeros(m, n); % total mass of EPS
volume_total_v = zeros(m, n); % total fluid volume
volume_total_h = zeros(m, n); % total fluid volume
am0_EPS_mean = zeros(m, n);% Before loops

domain = [0 n];

% Parameter sweeps
vf0_values = linspace(0.05, .25, m); % sweep volume fraction from 0.05 to 0.2 [m^2]
EPS_con = linspace(0, 1500, n); % EPS concentration sweep [kg/m^3]
rho_EPS = 1500; % EPS density [kg/m^3]
for j = 1:m
    vf0_current = vf0_values(j);
    am0 = pi*(7.d-5 + 1.6d-4*vf0_current)^2; %[m^2]
    h = sqrt((2*am0)/vf0_current);

    rmu = log(am0)-.5*rsig^2;
    parfor t = 1:trials
 % Re-init per-trial accumulators
         rng(baseSeed + 1000*j + t, 'twister');
        eff_trial              = zeros(1, n);
        eff_resistance_trial = zeros(1, n);
        mean_area_single_trial = zeros(1, n);

        % --------- ONE base geometry per (vf0, trial) ----------
        [A_sv0, A_sh0, vfrac, amean] = gen_pipes(rmu, rsig, nx, ny, h);
        fprintf('Trial %d (base): vfrac=%.5f, amean=%.6e\n', t, vfrac, amean);

        % Base volumes (mass = conc * volume)
        A_sv_acc = A_sv0;
        A_sh_acc=A_sh0;
        EPS_prev=0;

        for i = 1:n
            EPS_target = EPS_con(i);

            % -----------------------------------------------------
            % LINEAR NON-ACCUMULATING EPS FILLING MODEL
            %
            % Every EPS value starts from the same base geometry.
            % EPS_target/rho_EPS is the fraction of the original
            % brine volume filled by EPS.
            % -----------------------------------------------------
        
            open_fraction = max(1 - EPS_target/rho_EPS, 0);
        
            A_sv_acc = A_sv0 .* open_fraction;
            A_sh_acc = A_sh0 .* open_fraction;
        
            % Conductances after EPS reduction
            sv = (A_sv_acc/pi).^2;
            sh = (A_sh_acc/pi).^2;
        
            % Solve using multigrid
            fin = zeros(nx + 1, ny + 1);
            [phi, final_error] = kikmul(fin, pdrop, sv, sh, nx, ny);
        
            % Effective permeability computation
            effleft = sum((pdrop - phi(nx, :)) .* sh(nx, :)) / pdrop / h^2;
            effright = sum(phi(2, :) .* sh(1, :)) / pdrop / h^2;
            effcoe = (0.5 * (effleft + effright))*(pi/8);
        
            eff_trial(i) = effcoe;
            % -----------------------------------------------------
            % Resistance-based permeability calculation
            %
            % Treat each vertical column of A_sv_acc as pipes in series.
            % Treat columns as parallel flow paths.
            % -----------------------------------------------------
            
            n_cols_res = nx;    % number of vertical columns
            n_pipes_res = ny-1;      % pipes per column
            
            D_total_res = n_pipes_res * h;        % total vertical height [m]
            A_sample_res = n_cols_res * h^2;      % bulk sample cross-sectional area [m^2]
            
            G_col = zeros(n_cols_res, 1);
            
            for c_res = 1:n_cols_res
            
                R_col = 0;
            
                for q_res = 1:n_pipes_res
            
                    A_pipe = A_sv_acc(q_res, c_res);
            
                    if A_pipe <= 0
                        R_pipe = Inf;
                    else
                        a_pipe = sqrt(A_pipe/pi);
                        R_pipe = 8*mu_fluid*h/(pi*a_pipe^4);
                    end
            
                    R_col = R_col + R_pipe;
            
                end
            
                if isinf(R_col) || R_col <= 0
                    G_col(c_res) = 0;
                else
                    G_col(c_res) = 1/R_col;
                end
            
            end
            
            G_total = sum(G_col);
            
            k_resistance = mu_fluid * D_total_res * G_total / A_sample_res;
            
            eff_resistance_trial(i) = k_resistance;
            mean_area_single_trial(i) = mean([A_sv_acc(:); A_sh_acc(:)]);
        
            % Store for later inspection/plots
            sv_all{j,i,t} = sv;
            sh_all{j,i,t} = sh;
        
            fprintf('vf0=%.3f | step %2d | EPS=%.3e kg/m^3 | k=%.6e m^2\n', ...
                vf0_current, i, EPS_target, effcoe);
        
        end

        % After finishing all EPS levels for this trial, aggregate:
       eff_trial_all_perm(j,:,t) = eff_trial;
       eff_resistance_all_perm(j,:,t) = eff_resistance_trial;
       mean_area_all(j,:,t) = mean_area_single_trial;
    end

    % Means across trials for plotting
    eff_mean_perm_plot(j,:) = mean(eff_trial_all_perm(j,:,:), 3);
    eff_resistance_perm_plot(j,:) = mean(eff_resistance_all_perm(j,:,:), 3);
    mean_area_plot(j,:) = mean(mean_area_all(j,:,:), 3);
end
% ----------------------------------------------------
% Plotting
% ----------------------------------------------------

% Example random sample sv, sh for histogram plots
random_j = randi(m);
random_i = randi(n);
random_t = randi(trials);
sv_final = sv_all{random_j, random_i, random_t};
sh_final = sh_all{random_j, random_i, random_t};

% =========================================================
% Normalized permeability decay for each brine volume fraction
% Compare network/multigrid, series-parallel resistance, and theory
% =========================================================

figure('Position',[100 100 1200 700]);

rows = ceil(length(vf0_values)/3);
cols = min(3, length(vf0_values));

for j = 1:length(vf0_values)

    subplot(rows, cols, j);
    hold on; grid on; box on

    % Raw permeability curves
    k_net = eff_mean_perm_plot(j,:);
    k_res = eff_resistance_perm_plot(j,:);

    % Normalize by the no-EPS permeability k(0)
    if k_net(1) > 0
        k_net_norm = k_net / k_net(1);
    else
        k_net_norm = nan(size(k_net));
    end

    if k_res(1) > 0
        k_res_norm = k_res / k_res(1);
    else
        k_res_norm = nan(size(k_res));
    end

    % Theoretical normalized decay
    K_theory = max(1 - EPS_con(:).'/rho_EPS, 0).^2;

    plot(EPS_con, k_net_norm, 'o-', ...
        'LineWidth', 1.5, ...
        'MarkerSize', 6, ...
        'DisplayName', 'network / multigrid');

    plot(EPS_con, k_res_norm, 's--', ...
        'LineWidth', 1.8, ...
        'MarkerSize', 5, ...
        'DisplayName', 'series-parallel resistance');

    plot(EPS_con, K_theory, 'k-', ...
        'LineWidth', 2, ...
        'DisplayName', 'theory: (1-E/\rho)^2');

    xlabel('EPS concentration, E [kg/m^3]');
    ylabel('Normalized permeability, k(E)/k(0)');
    title(sprintf('\\phi = %.2f', vf0_values(j)));

    ylim([0 1.05]);
    legend('Location','best');

    hold off

end

sgtitle('Normalized permeability decay for different brine volume fractions');
ratio = eff_resistance_perm_plot(j,:) ./ eff_mean_perm_plot(j,:);
figure;
plot(EPS_con, ratio, 'o-', 'LineWidth', 1.5)
grid on
xlabel('EPS concentration, E [kg/m^3]');
ylabel('k_{resistance} / k_{network}');
title(sprintf('Ratio for \\phi = %.2f', vf0_values(j)));

% =========================================================
% Compare multigrid permeability and resistance permeability
% =========================================================

figure('Position',[100 100 1200 700]);

rows = ceil(m/3);
cols = min(3,m);

for j = 1:m

    subplot(rows, cols, j);
    hold on; grid on; box on

    plot(EPS_con, eff_mean_perm_plot(j,:), 'o-', ...
        'LineWidth', 1.5, ...
        'MarkerSize', 6, ...
        'DisplayName', 'network / multigrid');

    plot(EPS_con, eff_resistance_perm_plot(j,:), 's--', ...
        'LineWidth', 1.8, ...
        'MarkerSize', 5, ...
        'DisplayName', 'series-parallel resistance');

    xlabel('EPS concentration, E [kg/m^3]');
    ylabel('Effective permeability, k [m^2]');
    title(sprintf('\\phi = %.2f', vf0_values(j)));

    set(gca,'YScale','linear');
    legend('Location','best');

    hold off

end

sgtitle('Effective permeability: network solver vs resistance calculation');

% ----------------------------------------------------
% Plotting
% =========================================================
% MODEL COMPARISON: EPS vs Effective Permeability
%
% Data:
%   k_data(E)
%
% Theoretical curve:
%   k_theory(E) = k(0) max(1 - E/rho_EPS, 0)^2
%
% Fitted curve:
%   k_fit(E) = K0_fit max(1 - E/rho_EPS, 0)^2
%
% The fitted curve uses the same theoretical shape, but allows
% the prefactor K0_fit to be determined by least squares.
% =========================================================

C = EPS_con(:);                 % EPS concentrations [n x 1]
K = eff_mean_perm_plot;         % mean permeability per vf0 row [m x n]

[mRows, ~] = size(K);

% Storage for goodness of fit
K0_theory = zeros(mRows,1);
K0_fit    = zeros(mRows,1);

RMSE_theory  = zeros(mRows,1);
NRMSE_theory = zeros(mRows,1);
R2_theory    = zeros(mRows,1);

RMSE_fit  = zeros(mRows,1);
NRMSE_fit = zeros(mRows,1);
R2_fit    = zeros(mRows,1);

figure;
rows = ceil(mRows/3);
cols = min(3,mRows);

for j = 1:mRows

    subplot(rows, cols, j);
    hold on; grid on; box on

    % -----------------------------
    % Observed data
    % -----------------------------
    k_obs = K(j,:).';

    valid = isfinite(C) & isfinite(k_obs) & (k_obs >= 0);

    Cj = C(valid);
    kj = k_obs(valid);

    % Basis function from theory
    f = max(1 - Cj/rho_EPS, 0).^2;

    % -----------------------------
    % Theoretical curve
    % Uses k(0) from the first data point
    % -----------------------------
    K0_theory(j) = k_obs(1);
    k_theory = K0_theory(j) * f;

    % -----------------------------
    % Fitted curve
    % Fit only the prefactor K0_fit in least-squares sense:
    % minimize ||K0*f - k_obs||^2
    % -----------------------------
    if sum(f.^2) > 0
        K0_fit(j) = sum(f .* kj) / sum(f.^2);
    else
        K0_fit(j) = NaN;
    end

    k_fit = K0_fit(j) * f;

    % -----------------------------
    % Goodness of fit: theory
    % -----------------------------
    err_theory = k_theory - kj;
    RMSE_theory(j) = sqrt(mean(err_theory.^2));

    data_range = max(kj) - min(kj);
    if data_range > 0
        NRMSE_theory(j) = RMSE_theory(j) / data_range;
    else
        NRMSE_theory(j) = NaN;
    end

    ss_res = sum((kj - k_theory).^2);
    ss_tot = sum((kj - mean(kj)).^2);

    if ss_tot > 0
        R2_theory(j) = 1 - ss_res/ss_tot;
    else
        R2_theory(j) = NaN;
    end

    % -----------------------------
    % Goodness of fit: fitted prefactor
    % -----------------------------
    err_fit = k_fit - kj;
    RMSE_fit(j) = sqrt(mean(err_fit.^2));

    if data_range > 0
        NRMSE_fit(j) = RMSE_fit(j) / data_range;
    else
        NRMSE_fit(j) = NaN;
    end

    ss_res = sum((kj - k_fit).^2);

    if ss_tot > 0
        R2_fit(j) = 1 - ss_res/ss_tot;
    else
        R2_fit(j) = NaN;
    end

    % -----------------------------
    % Smooth curves for plotting
    % -----------------------------
    Cfine = linspace(min(EPS_con), max(EPS_con), 300).';
    ffine = max(1 - Cfine/rho_EPS, 0).^2;

    K_theory_fine = K0_theory(j) * ffine;
    K_fit_fine    = K0_fit(j)    * ffine;

    % -----------------------------
    % Plot
    % -----------------------------
    plot(EPS_con, k_obs, 'o-', ...
        'LineWidth', 1.5, ...
        'MarkerSize', 6, ...
        'DisplayName', 'data');

    plot(Cfine, K_theory_fine, 'k--', ...
        'LineWidth', 2, ...
        'DisplayName', 'theory: k(0)(1-E/\rho)^2');

    plot(Cfine, K_fit_fine, 'r-', ...
        'LineWidth', 2, ...
        'DisplayName', 'fit: K_0(1-E/\rho)^2');

    xlabel('EPS concentration, E [kg/m^3]');
    ylabel('Effective permeability, k [m^2]');
    title(sprintf('\\phi = %.2f', vf0_values(j)));

    legend('Location','best');

    hold off

end

sgtitle('EPS vs Effective Permeability: Data, Theory, and Fit');

% =========================================================
% Goodness-of-fit table
% =========================================================

GOF_Table = table( ...
    vf0_values(:), ...
    K0_theory(:), ...
    K0_fit(:), ...
    RMSE_theory(:), ...
    NRMSE_theory(:), ...
    R2_theory(:), ...
    RMSE_fit(:), ...
    NRMSE_fit(:), ...
    R2_fit(:), ...
    'VariableNames', { ...
        'vf0', ...
        'K0_theory', ...
        'K0_fit', ...
        'RMSE_theory', ...
        'NRMSE_theory', ...
        'R2_theory', ...
        'RMSE_fit', ...
        'NRMSE_fit', ...
        'R2_fit'} ...
);

disp(GOF_Table);
% --- Subplots: EPS concentration vs Expected Area for each vf0 ---
figure;
rows = ceil(m/3);  % 2 rows if m=5
cols = 3;          % 3 columns

for j = 1:m
    subplot(rows, cols, j);
    plot(EPS_con, mean_area_plot(j,:), '-^', 'LineWidth', 2);
    xlabel('EPS Concentration (kg/m^3)');
    ylabel('E[A] (m^2)');
    title(sprintf('vf0 = %.2f', vf0_values(j)));
    grid on;

end

sgtitle('EPS Concentration vs Expected Area for Different vf0');
% 3. Histograms: log areas after EPS accumulation
figure;

% Remove zero or negative values before taking the log
sv_positive = sv_final(sv_final > 0);
sh_positive = sh_final(sh_final > 0);

subplot(1,2,1);
histogram(log(sv_positive(:)), 50, 'FaceColor', 'b');
xlabel('ln(A) of Vertical Pipes');
ylabel('Frequency');
title(sprintf('Histogram of ln(sv)\n(vf0=%.2f, EPS=%.2f)', vf0_values(random_j), EPS_con(random_i)));

subplot(1,2,2);
histogram(log(sh_positive(:)), 50, 'FaceColor', 'r');
xlabel('ln(A) of Horizontal Pipes');
ylabel('Frequency');
title(sprintf('Histogram of ln(sh)\n(vf0=%.2f, EPS=%.2f)', vf0_values(random_j), EPS_con(random_i)));

sgtitle('Distribution of ln(A) for Pipe Cross-Sections After EPS Accumulation');
% 7. For each EPS concentration, plot vf0 vs effective permeability
% =========================================================
% Effective permeability vs brine volume fraction
% Show only 4 EPS concentrations
% =========================================================

figure;
hold on; grid on; box on

% Pick 4 representative EPS values:
% first, one-third, two-thirds, last
eps_indices = [1, round(n/4), round(n/2), n];

markerStyles = {'x', '+', 'o', 's', 'd', '^', 'v'};


for rr = 1:length(eps_indices)

    i = eps_indices(rr);

    plot(vf0_values, eff_mean_perm_plot(:,i), ...
        'LineStyle', 'none', ...
        'Marker', markerStyles{rr}, ...
        'MarkerSize', 9, ...
        'LineWidth', 1.8, ...
        'DisplayName', sprintf('EPS = %.1f kg/m^3', EPS_con(i)));

end

set(gca, 'YScale', 'log');
xlabel('Brine volume fraction, \phi');
ylabel('Effective permeability, k [m^2]');
title('Effective permeability vs brine volume fraction (network)');

legend('Location','northwest');
ylim([1e-13 1e-8]);   % adjust if needed
xlim([min(vf0_values)*0.8, max(vf0_values)*1.1]);
hold off

figure;
hold on; grid on; box on

% Pick 4 representative EPS values:
% first, one-third, two-thirds, last
eps_indices = [1, round(n/4), round(n/2), n];

markerStyles = {'x', '+', 'o', 's', 'd', '^', 'v'};


for rr = 1:length(eps_indices)

    i = eps_indices(rr);

    plot(vf0_values, eff_resistance_perm_plot(:,i), ...
        'LineStyle', 'none', ...
        'Marker', markerStyles{rr}, ...
        'MarkerSize', 9, ...
        'LineWidth', 1.8, ...
        'DisplayName', sprintf('EPS = %.1f kg/m^3', EPS_con(i)));

end

set(gca, 'YScale', 'log');


xlabel('Brine volume fraction, \phi');
ylabel('Effective permeability, k [m^2]');
title('Effective permeability vs brine volume fraction (column)');

legend('Location','northwest');
ylim([1e-13 1e-8]);   % adjust if needed
xlim([min(vf0_values)*0.8, max(vf0_values)*1.1]);

hold off
% Subplots: each EPS gets its own panel of k vs vf0
n_eps = numel(EPS_con);

% Near-square layout
rows = ceil(sqrt(n_eps));
cols = ceil(n_eps/rows);

% Consistent y-limits across all panels (log scale)
kmin = max( min(eff_mean_perm_plot(:)), 1e-16 );  % avoid nonpositive
kmax = max(eff_mean_perm_plot(:));
% If you want fixed paper-like limits, uncomment:
% kmin = 1e-13; kmax = 1e-8;

figure;
tiledlayout(rows, cols, "Padding","compact","TileSpacing","compact");

for i = 1:n_eps
    nexttile;
    plot(vf0_values, eff_mean_perm_plot(:, i), '-o', 'LineWidth', 1.6, 'MarkerSize', 5);
    set(gca, 'YScale', 'linear');
    ylim([kmin, kmax]);
    if numel(vf0_values) == 1
        xlim(vf0_values(1) + [-0.01, 0.01]);
    else
        xlim([min(vf0_values), max(vf0_values)]);
    end
    grid on;

    % Labels (sparse to reduce clutter)
    if i > (rows-1)*cols
        xlabel('Brine volume fraction, \phi');
    end
    if mod(i-1, cols) == 0
        ylabel('k (m^2)');
    end

    title(sprintf('EPS = %.3g kg m^{-3}', EPS_con(i)));
end

% --- Subplots: Expected Area vs Effective Permeability (one panel per vf0) ---
figure;
rows = ceil(m/3);
cols = 3;

for j = 1:m
    subplot(rows, cols, j);
    hold on; grid on;
    
    % plot mean area vs permeability for this φ
    plot(mean_area_plot(j,:), eff_mean_perm_plot(j,:), 'o-', 'LineWidth', 1.5);
    
    set(gca, 'YScale','linear'); % permeability usually log scale
    xlabel('E[A] (m^2)');
    ylabel('k (m^2)');
    title(sprintf('\\phi = %.2f', vf0_values(j)));
end

sgtitle('Expected Area vs Effective Permeability for Each Volume Fraction');
function [A_sv, A_sh, vfrac, amean] = gen_pipes(rmu, rsig, nx, ny, h)
    xmin = rmu - 4 * rsig;
    xmax = rmu + 4 * rsig;

    % Sample log-normal variables with truncation
    X_sv = max(min(rmu + rsig * randn(nx + 1, ny + 1), xmax), xmin);
    X_sh = max(min(rmu + rsig * randn(nx + 1, ny + 1), xmax), xmin);

    A_sv = exp(X_sv);
    A_sh = exp(X_sh);

    % % Compute conductance coefficients
    % sv = (A_sv / pi).^2;
    % sh = (A_sh / pi).^2;

    % % Apply random pipe failures (broken ducts)
    % broken_h = rand(nx + 1, ny + 1) <= ph;  % horizontal breakage
    % broken_v = rand(nx + 1, ny + 1) <= pv;  % vertical breakage

    % sh(broken_h) = 0;
    % sv(broken_v) = 0;

    % Compute volume fraction and average area (like Fortran netstr.f)
    asum = sum(A_sv(:)) + sum(A_sh(:));
    amean = 0.5 * asum / ((nx + 1) * (ny + 1));
    vfrac = asum / ((nx + 1) * (ny + 1)) / h^2;
end
fprintf('Total time elapsed: %.2f s\n', toc(t_all));