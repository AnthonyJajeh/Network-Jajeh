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
%     This implimentation takes into consideration EPS as some base value
%     rather then coming from an ODE.
%*************************************************************************

nx = 256; %size of grid nxn power of 2
ny = nx;
pdrop = 1.d0; %pressure drop
vf0 = .05; %targeted volume fraction
rsig = .1; %variance of lognormal distribution
trials = 10; %number  of trials
ph = 0; %probability of horizontal pipe breaking
pv = 0; %probability of vertical pipe breaking
n = 5;
m = 5;

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

sv_all = cell(n, trials);
sh_all = cell(n, trials);

% Parameter sweeps
vf0_values = linspace(0.05, 0.25, m); % sweep volume fraction from 0.05 to 0.2 [m^2]
EPS_con = linspace(0, 25, n); % EPS concentration sweep [kg/m^3]
rho_EPS = 1500; % EPS density [kg/m^3]
for j = 1:m
    vf0_current = vf0_values(j);
    for i = 1:n
        EPS_concentration = EPS_con(i); %[kg/m^3]
        eff_trial = zeros(1, trials);
        mean_area_single_trial = zeros(1, trials);
        am0_EPS_sum = 0;
        for t = 1:trials
            % Compute modified am0_EPS based on current vf0 and EPS
            am0 = pi*(7.d-5 + 1.6d-4*vf0_current)^2; %[m^2]
            am0_EPS = am0 * (1 - ((EPS_concentration)/rho_EPS));
                 fprintf("vf0 = %.3f, EPS = %.4f kg/m^3, am0 = %.2e, am0_EPS = %.2e\n", ...
        vf0, EPS_concentration, am0, am0_EPS);

            am0_EPS_sum = am0_EPS_sum + am0_EPS; %[m^2]
            
            rmu = log(am0_EPS) - 0.5*rsig^2;
            h = sqrt((2*am0_EPS)/vf0_current);

           [A_sv, A_sh, vfrac, amean] = gen_pipes(rmu, rsig, nx, ny, h);
           fprintf('Trial %d: vfrac=%.5f, amean=%.6e\n', t, vfrac, amean);

            % Volume and mass for EPS accumulation
            volume_total_v(j,i) = (sum(A_sv(:))) * h;
            volume_total_h(j,i)= (sum(A_sh(:))) * h;

            mass_EPS_total_v(j,i) = volume_total_v(j,i) * EPS_concentration;
            mass_EPS_total_h(j,i) = volume_total_h(j,i) * EPS_concentration;

            % Distribute EPS mass proportionally
            d_v = A_sv / sum(A_sv(:));
            d_h = A_sh / sum(A_sh(:));
            volume_EPS_v = (d_v .* mass_EPS_total_v(j,i)) / rho_EPS;
            volume_EPS_h = (d_h .* mass_EPS_total_h(j,i)) / rho_EPS;

            % Update areas after EPS accumulation
            A_sv = max(0, A_sv - volume_EPS_v / h);
            A_sh = max(0, A_sh - volume_EPS_h / h);
            % Compute average EPS-adjusted pipe area
            sv = (A_sv/pi).^2;
            sh = (A_sh/pi).^2;

            sv_all{j,i,t} = sv;
            sh_all{j,i,t} = sh;

            % Solve using multigrid (kikmul)
            fin = zeros(nx + 1, ny + 1);
            [phi, final_error] = kikmul(fin, pdrop, sv, sh, nx, ny);

            % Effective permeability computation
            effleft = sum((pdrop - phi(nx, :)) .* sh(nx, :)) / pdrop / h^2;
            effright = sum(phi(2, :) .* sh(1, :)) / pdrop / h^2;
            effcoe = (0.5 * (effleft + effright))*(pi/8);
            eff_trial(t) = effcoe;
            mean_area_single_trial(t) = mean([A_sv(:); A_sh(:)]);
            % Debug output
            fprintf('vf0=%.3e EPS=%.2e Trial=%d effcoe=%.6e\n', vf0_current, EPS_concentration, t, effcoe);
        end

        % Store results for this (vf0, EPS) pair
        eff_trial_all_perm(j,i,:) = eff_trial;
        eff_mean_perm_plot(j,i) = mean(eff_trial);
        mean_area_all(j,i,:) = mean_area_single_trial;
        mean_area_plot(j,i) = mean(mean_area_single_trial);
        am0_EPS_mean(j,i) = am0_EPS_sum / trials;

    end
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
% 1. Heatmap: Effective permeability
figure;
imagesc(EPS_con, vf0_values, eff_mean_perm_plot);
colorbar;
xlabel('EPS Concentration');
ylabel('Volume Fraction (vf0)');
title('Effective Permeability vs EPS Concentration and Volume Fraction');
set(gca, 'YDir', 'normal');
grid on;

% 2. Heatmap: Mean Pipe Area
figure;
imagesc(EPS_con, vf0_values, mean_area_plot);
colorbar;
xlabel('EPS Concentration');
ylabel('Volume Fraction (vf0)');
title('Mean Pipe Area vs EPS Concentration and Volume Fraction');
set(gca, 'YDir', 'normal');
grid on;

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

% 4. EPS Concentration vs Effective Permeability
colors = lines(m); % Create a colormap for m lines
figure;
hold on;
for j = 1:m
    plot(EPS_con, eff_mean_perm_plot(j,:), '-o','LineWidth',2, 'Color', colors(j,:), 'DisplayName', sprintf('vf0=%.2f', vf0_values(j)));
end
xlabel('EPS Concentration');
ylabel('Effective Permeability');
title('EPS Concentration vs Effective Permeability');
legend('show', 'Location', 'bestoutside');
grid on;

% 5. Mean Area vs Effective Permeability
figure;
hold on;
for j = 1:m
    plot(mean_area_plot(j,:), eff_mean_perm_plot(j,:), '-s','LineWidth',2, 'Color', colors(j,:), 'DisplayName', sprintf('vf0=%.2f', vf0_values(j)));
end
xlabel('Mean Pipe Area');
ylabel('Effective Permeability');
title('Mean Pipe Area vs Effective Permeability');
legend('show', 'Location', 'bestoutside');
grid on;

% 6. EPS Concentration vs Mean Area
figure;
hold on;
for j = 1:m
    plot(EPS_con, mean_area_plot(j,:), '-^','LineWidth',2, 'Color', colors(j,:), 'DisplayName', sprintf('vf0=%.2f', vf0_values(j)));
end
xlabel('EPS Concentration');
ylabel('Mean Pipe Area');
title('EPS Concentration vs Mean Pipe Area');
legend('show', 'Location', 'bestoutside');
grid on;

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
    % 
    % sh(broken_h) = 0;
    % sv(broken_v) = 0;

    % Compute volume fraction and average area (like Fortran netstr.f)
    asum = sum(A_sv(:)) + sum(A_sh(:));
    amean = 0.5 * asum / ((nx + 1) * (ny + 1));
    vfrac = asum / ((nx + 1) * (ny + 1)) / h^2;
end