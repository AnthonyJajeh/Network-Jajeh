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
vf0 = .25; %targeted volume fraction
rsig = .5; %mean of lognormal distribution
trials = 5; %number  of trials
n = 10;

eff_mean_plot = zeros(1, n); %allocated spot for effective permeability
eff_trial_all=zeros(n,trials);
EPS_con = linspace(0, 50, n); % EPS concentration values
eff_mean = zeros(1, n); %allocated spot for effective permaebility
domain = [0 n];

     %The following is to specify the distribution parameters for the
%     connecting tube bond coefficients, and generate samples according 
%     to the target distribution.
%     
%     vf0 is the targeted volume fraction, and am0 is the average of the
%     cross-sectional areas, which is assumed to be a function of the target 
%     volumn fraction:
%      am0 = 3.d-2+0.13d0*(vf0-2.d-2)
      am0 = pi*(7.d-2 + 1.6d-1*vf0)^2;
%     so the average cross section area can lead to the average of cross
%     sectional radius:
      rmu = log(am0) - 0.5*rsig^2;
%     which is the mean of log area and
      amu = exp(2.0*rmu)*(exp(2.0*rsig^2)-exp(rsig^2));
%     which is the variance of the cross section areas, and cell size
      h = sqrt(2*am0/vf0);
      fprintf('std of cross section areas is %8.5f\n',sqrt(amu));

mass_EPS_total = zeros(1, n); %total mass of eps allocation
volume_total = zeros(1, n); %total volume allocation
rho_EPS = 1500;

parfor i = 1:n
    eff_trial = zeros(1,trials);
    EPS_concentration = EPS_con(i)
    for t = 1:trials
        % Initialize cross-sectional areas for the first iteration
        sv = exp(rmu + rsig * randn(nx + 1, ny + 1));
        sh = exp(rmu + rsig * randn(nx + 1, ny + 1));
        % Update volume and mass for each EPS concentration
        volume_total(i) = (sum(sv(:)) + sum(sh(:))) * h;  % Total volume of pipes
        mass_EPS_total(i) = volume_total(i) * EPS_concentration;  % Total mass of EPS
        
        % Distribution factor for each pipe (proportional)
        d_v = sv / sum(sv(:));  % Vertical distribution
        d_h = sh / sum(sh(:));  % Horizontal distribution
        
        % Compute EPS volume based on mass and current density
        volume_EPS_v = (d_v .* mass_EPS_total(i)) / rho_EPS;
        volume_EPS_h = (d_h .* mass_EPS_total(i)) / rho_EPS;
        
        % Calculate total EPS volume and update rho_EPS_eff dynamically

        
        % Update cross-sectional areas after EPS accumulation
        sv = max(0, sv - volume_EPS_v / h);  % Update vertical pipes (sv)
        sh = max(0, sh - volume_EPS_h / h);  % Update horizontal pipes (sh)

    % Solve the permeability problem using the multigrid method (kikmul)
    fin = zeros(nx + 1, ny + 1);
    [phi, final_error] = kikmul(fin, pdrop, sv, sh, nx, ny);
    
    % Compute the effective permeability
    effleft = sum((pdrop - phi(nx, :)) .* sh(nx, :)) / pdrop / h^2;
    effright = sum(phi(2, :) .* sh(1, :)) / pdrop / h^2;
    effcoe = 0.5 * (effleft + effright);
    eff_trial(t) = effcoe;
    % Debugging output
    fprintf('Trial %d, EPS concentration = %8.2f, effcoe = %8.6f\n', t, EPS_concentration, effcoe);
    end
       if i == n-5
            sv_final = sv;
            sh_final = sh;
       end
    %store all trial results
    eff_trial_all(i,:)=eff_trial;
    % compute thee maen permaebility across all trials for specific EPS
    % concentriation
    eff_mean_plot(i) = mean(eff_trial);
    % Debugging output
      fprintf('effcoe = %8.6f\n',effcoe);
      fprintf('rmu and rsig: %12.6f %12.6f\n',rmu,rsig);
end

% Plot permeability vs EPS concentration
% Plot EPS concentration vs Effective Permeability for all trials
figure;
hold on;

% Plot the individual trial data points
for i = 1:n
    plot(repmat(EPS_con(i), 1, trials), eff_trial_all(i, :), 'bo'); % Trial data points
end
% Calculate the standard deviation for each EPS concentration
eff_std = std(eff_trial_all, 0, 2); % Standard deviation along the second dimension (across trials)

% Plot the mean effective permeability with error bars
errorbar(EPS_con, eff_mean_plot, eff_std, 'r-', 'LineWidth', 2, 'MarkerSize', 8); % Mean with error bars

xlabel('EPS concentration');
ylabel('Effective Permeability');
title('Effect of EPS on Permeability for Multiple Trials');
grid on;
hold off;
 
figure;
subplot(1,2,1);
histogram(log(sv_final(:)));
xlabel('ln(A) of Vertical Pipes');
ylabel('Frequency');
title('Histogram of ln(sv)');

subplot(1,2,2);
histogram(log(sh_final(:)));
xlabel('ln(A) of Horizontal Pipes');
ylabel('Frequency');
title('Histogram of ln(sh)');

sgtitle('Distribution of ln(A) for Pipe Cross-Sections After EPS Accumulation');

