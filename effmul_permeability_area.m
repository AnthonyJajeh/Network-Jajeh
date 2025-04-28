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
%
%*************************************************************************

nx = 256; %size of grid nxn power of 2
ny = nx;
pdrop = 1.d0; %pressure drop
vf0 = .1; %targeted volume fraction
rsig = .5; %mean of lognormal distribution
%rmu = 1;
%h=.9;
trials = 3; %number  of trials
n = 5;
m = 5;

eff_mean_perm_plot = zeros(m, n); % mean effective permeability
eff_trial_all_perm = zeros(m, n, trials); % all permeability trials
mean_area_plot = zeros(m, n); % mean pipe area
mean_area_all = zeros(m, n, trials); % all mean area trials
mass_EPS_total = zeros(m, n); % total mass of EPS
volume_total = zeros(m, n); % total fluid volume

domain = [0 n];

sv_all = cell(n, trials);
sh_all = cell(n, trials);

% Parameter sweeps
vf0_values = linspace(0.05, 0.2, m); % sweep volume fraction from 0.05 to 0.2
EPS_con = linspace(0, 100, n); % EPS concentration sweep
rho_EPS = 1500; % EPS density
parfor j = 1:m
    vf0_current = vf0_values(j);
    for i = 1:n
        EPS_concentration = EPS_con(i);

        eff_trial = zeros(1, trials);
        mean_area_single_trial = zeros(1, trials);

        for t = 1:trials
            % Compute modified am0_EPS based on current vf0 and EPS
            am0 = pi*(7.d-2 + 1.6d-1*vf0_current)^2;
            am0_EPS = am0 * (1 - (EPS_concentration/rho_EPS));
            rmu = log(am0_EPS) - 0.5*rsig^2;
            h = sqrt((2*am0_EPS)/vf0_current);

            % Initialize sv, sh
            sv = exp(rmu + rsig * randn(nx + 1, ny + 1));
            sh = exp(rmu + rsig * randn(nx + 1, ny + 1));

            % Volume and mass for EPS accumulation
            volume_total(j,i) = (sum(sv(:)) + sum(sh(:))) * h;
            mass_EPS_total(j,i) = volume_total(j,i) * EPS_concentration;

            % Distribute EPS mass proportionally
            d_v = sv / sum(sv(:));
            d_h = sh / sum(sh(:));
            volume_EPS_v = (d_v .* mass_EPS_total(j,i)) / rho_EPS;
            volume_EPS_h = (d_h .* mass_EPS_total(j,i)) / rho_EPS;

            % Update areas after EPS accumulation
            sv = max(0, sv - volume_EPS_v / h);
            sh = max(0, sh - volume_EPS_h / h);

            sv_all{j,i,t} = sv;
            sh_all{j,i,t} = sh;

            % Solve using multigrid (kikmul)
            fin = zeros(nx + 1, ny + 1);
            [phi, final_error] = kikmul(fin, pdrop, sv, sh, nx, ny);

            % Effective permeability computation
            effleft = sum((pdrop - phi(nx, :)) .* sh(nx, :)) / pdrop / h^2;
            effright = sum(phi(2, :) .* sh(1, :)) / pdrop / h^2;
            effcoe = 0.5 * (effleft + effright);
            eff_trial(t) = effcoe;

            % Mean area for this trial
            mean_area_single_trial(t) = mean([sv(:); sh(:)]);

            % Debug output
            fprintf('vf0=%.3f EPS=%.2f Trial=%d effcoe=%.6f\n', vf0_current, EPS_concentration, t, effcoe);
        end

        % Store results for this (vf0, EPS) pair
        eff_trial_all_perm(j,i,:) = eff_trial;
        eff_mean_perm_plot(j,i) = mean(eff_trial);
        mean_area_all(j,i,:) = mean_area_single_trial;
        mean_area_plot(j,i) = mean(mean_area_single_trial);

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
subplot(1,2,1);
histogram(log(sv_final(:)), 50, 'FaceColor', 'b');
xlabel('ln(A) of Vertical Pipes');
ylabel('Frequency');
title('Histogram of ln(sv)');

subplot(1,2,2);
histogram(log(sh_final(:)), 50, 'FaceColor', 'r');
xlabel('ln(A) of Horizontal Pipes');
ylabel('Frequency');
title('Histogram of ln(sh)');

sgtitle('Distribution of ln(A) for Pipe Cross-Sections After EPS Accumulation');

% 1. EPS Concentration vs Effective Permeability
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

% 2. Mean Area vs Effective Permeability
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

% 3. EPS Concentration vs Mean Area
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
