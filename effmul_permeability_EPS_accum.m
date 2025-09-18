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

nx = 1024; %size of grid nxn power of 2
ny = nx;
pdrop = 1.d0; %pressure drop
vf0 = .05; %targeted volume fraction
rsig = .5; %variance of lognormal distribution
trials = 5; %number  of trials
ph = 0; %probability of horizontal pipe breaking
pv = 0; %probability of vertical pipe breaking
n = 15;%number of deposit steps;
m = 10;

eff_mean_perm_plot = zeros(m, n); % mean effective permeability
eff_trial_all_perm = zeros(m, n, trials); % all permeability trials
mean_area_plot = zeros(m, n); % mean pipe area
mean_area_all = zeros(m, n, trials); % all mean area trials
mass_EPS_total_v = zeros(m, n); % total mass of EPS
mass_EPS_total_h = zeros(m, n); % total mass of EPS
volume_total_v = zeros(m, n); % total fluid volume
volume_total_h = zeros(m, n); % total fluid volume
am0_EPS_mean = zeros(m, n);% Before loops

domain = [0 n];

% Parameter sweeps
vf0_values = linspace(0.05, 0.25, m); % sweep volume fraction from 0.05 to 0.2 [m^2]
EPS_con = linspace(0, 3000, n); % EPS concentration sweep [kg/m^3]
rho_EPS = 1500; % EPS density [kg/m^3]
for j = 1:m
    vf0_current = vf0_values(j);
    am0 = pi*(7.d-5 + 1.6d-4*vf0_current)^2; %[m^2]
    h = sqrt((2*am0)/vf0_current);

    rmu = log(am0)-.5*rsig^2;
    parfor t = 1:trials
 % Re-init per-trial accumulators
        eff_trial              = zeros(1, n);
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
            dEPS = max(0,EPS_target-EPS_prev);
            EPS_prev = EPS_target;
               % ---- Start from the SAME base geometry each EPS level ----
          
     % ---------------- Conductance-weighted EPS deposition ----------------
            % Provisional conductances from current areas (before deposition)
            sv_temp = (A_sv_acc/pi).^2;
            sh_temp = (A_sh_acc/pi).^2;

            % % Flow-capacity weights: C ∝ sv_tmp / h  (constants cancel in normalization)
            % C_v = sv_temp/h;
            % C_h = sh_temp/h;
            % W_v = C_v / (sum(C_v(:)) + eps);
            % W_h = C_h / (sum(C_h(:)) + eps);
            %W_v = A_sv_acc/sum(A_sv_acc(:));
            %W_h = A_sh_acc/sum(A_sh_acc(:));
            W_v = (pi/(8)*eta)

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
            eff_trial(i) = effcoe;
            mean_area_single_trial(i) = mean([A_sv_acc(:); A_sh_acc(:)]);
                % Store for later inspection/plots
            sv_all{j,i,t} = sv;
            sh_all{j,i,t} = sh;

            % Debug output
           fprintf('vf0=%.3f | step %2d | dEPS=%.3e kg/m^3 | cumEPS=%.3e | k=%.6e m^2\n', ...
            vf0_current, i, dEPS, EPS_target, effcoe);
        end

        % After finishing all EPS levels for this trial, aggregate:
        eff_trial_all_perm(j,:,t) = eff_trial;                    % 1×n
        mean_area_all(j,:,t) = mean_area_single_trial;       % 1×n
    end

    % Means across trials for plotting
    eff_mean_perm_plot(j,:) = mean(eff_trial_all_perm(j,:,:), 3);
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

% ----------------------------------------------------
% Plotting
% ----------------------------------------------------
% --- Separate plots: EPS concentration vs effective permeability (one per vf0) ---

% --- Subplots: EPS vs Effective Permeability for each vf0 ---
figure;
rows = ceil(m/3);  % 2 rows if m=5
cols = 3;          % 3 columns

for j = 1:m
    subplot(rows, cols, j);
    plot(EPS_con, eff_mean_perm_plot(j,:), '-o', 'LineWidth', 2);
    xlabel('EPS Concentration (kg/m^3)');
    ylabel('k (m^2)');
    title(sprintf('vf0 = %.2f', vf0_values(j)));
    grid on;
end

sgtitle('EPS Concentration vs Effective Permeability for Different vf0');

% --- Subplots: E[A] vs Effective Permeability for each vf0 ---
figure;
rows = ceil(m/3);  % make it 2 rows if m=5
cols = 3;          % 3 columns
for j = 1:m
    subplot(rows, cols, j);
    plot(mean_area_plot(j,:), eff_mean_perm_plot(j,:), '-o', 'LineWidth', 2);
    xlabel('E[A] (m^2)');
    ylabel('k (m^2)');
    title(sprintf('vf0 = %.2f', vf0_values(j)));
    grid on;
end

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
colors = lines(n); % Different color for each EPS concentration
figure;
hold on;

% Define marker styles (to mimic paper’s different symbols)
markerStyles = {'x', '+', 'o', 's', 'd', '^', 'v', '>', '<', 'p', 'h', '*', '.', 'x', '+'};

colors = lines(n); % Different color for each EPS concentration


for i = 1:n
    % Select marker style (cycle if you have more EPS points than styles)
    marker = markerStyles{mod(i-1, length(markerStyles)) + 1};
    
    % Plot with distinct marker and color
    plot(vf0_values, eff_mean_perm_plot(:,i), ...
        'LineStyle', 'none', ...
        'Marker', marker, ...
        'MarkerSize', 8, ...
        'LineWidth', 1.5, ...
        'Color', colors(i,:), ...
        'DisplayName', sprintf('EPS=%.4f kg/m^3', EPS_con(i)));
end

set(gca, 'YScale', 'log'); % Set y-axis to log scale
xlabel('Brine Volume Fraction, \phi');
ylabel('Effective Permeability, k (m^2)');
title('Effective Permeability vs Brine Volume Fraction');
legend('show', 'Location', 'bestoutside');
grid on;
hold off;

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
    set(gca, 'YScale', 'log');
    ylim([kmin, kmax]);
    xlim([min(vf0_values), max(vf0_values)]);
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

sgtitle('Effective permeability k vs brine volume fraction \phi (one subplot per EPS)');
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