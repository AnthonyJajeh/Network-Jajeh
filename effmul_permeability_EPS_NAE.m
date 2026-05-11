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

nx = 256; %size of grid nxn power of 2
ny = nx;
pdrop = 1.d0; %pressure drop
vf0 = .05; %targeted volume fraction
rsig = .5; %variance of lognormal distribution
trials = 5; %number  of trials
ph = 0; %probability of horizontal pipe breaking
pv = 0; %probability of vertical pipe breaking
n = 100;%number of deposit steps;
m = 10;

T = 30;      % ODE horizon (tune)
t_grid = linspace(0, T, n+1);

%ode def 
phi_par = .01;
psi = .01;
mu = .001;
gamma = .01; 
nu = .2; 
rho = .75; 
xi = .2;
delta = .007; 
eta = .03;

N0 = 5.0;
A0 = 0.05;
E0 = 0.0;

eff_mean_perm_plot = zeros(m, n); % mean effective permeability
eff_trial_all_perm = zeros(m, n, trials); % all permeability trials
mean_area_plot = zeros(m, n); % mean pipe area
mean_area_all = zeros(m, n, trials); % all mean area trials
mass_EPS_total_v = zeros(m, n); % total mass of EPS
mass_EPS_total_h = zeros(m, n); % total mass of EPS
volume_total_v = zeros(m, n); % total fluid volume
volume_total_h = zeros(m, n); % total fluid volume
am0_EPS_mean = zeros(m, n);% Before loops
k_base = zeros(m, trials);   % baseline k per (vf0, trial)

% A lightweight, parfor-safe dump of one snapshot + time series per trial
trial_dump(m, trials) = struct('sv', [], 'sh', [], 'E', [], 'k', [], 'Amean', []);
domain = [0 n];

% Parameter sweeps
vf0_values = linspace(0.05, 0.25, m); % sweep volume fraction from 0.05 to 0.2 [m^2]
rho_EPS = 1500; % EPS density [kg/m^3]
for j = 1:m
    vf0_current = vf0_values(j);
    am0 = pi*(7.d-5 + 1.6d-4*vf0_current)^2; %[m^2]
    h = sqrt((2*am0)/vf0_current);
    rmu = log(am0)-.5*rsig^2;
    
    parfor t = 1:trials
 % Re-init per-trial accumulators
   
        % --------- ONE base geometry per (vf0, trial) ----------
        [A_sv0, A_sh0, vfrac, amean] = gen_pipes(rmu, rsig, nx, ny, h);
        fprintf('Trial %d (base): vfrac=%.5f, amean=%.6e\n', t, vfrac, amean);

        % Base volumes (mass = conc * volume)
        A_sv_acc = A_sv0;
        A_sh_acc=A_sh0;
        sv = (A_sv_acc/pi).^2;
    sh = (A_sh_acc/pi).^2;

    fin = zeros(nx + 1, ny + 1);
    [phi_field, final_error] = kikmul(fin, pdrop, sv, sh, nx, ny);

    effleft  = sum((pdrop - phi_field(nx, :)) .* sh(nx, :)) / pdrop / h^2;
    effright = sum(phi_field(2, :) .* sh(1, :)) / pdrop / h^2;
    k_base_trial = (0.5 * (effleft + effright))*(pi/8);   % baseline k for THIS geometry

    % keep the baseline (we'll average across trials later)
    k_base(j, t) = k_base_trial;

    fprintf('vf0=%.3f | baseline (no EPS) | k=%.6e m^2\n', vf0_current, k_base_trial);

    % Optionally: use baseline k as reference scale in ODE
    k_0 = k_base_trial;   % gives sharper coupling than a fixed k0

    % =================== (B) ODE-driven EPS accumulation ===================
    % Reset accumulators to the SAME base geometry to start EPS deposition
    A_sv_acc = A_sv0;
    A_sh_acc = A_sh0;

    k_0 = k_base_trial;
    
    % ODE state over bins
    IC = [N0; A0; E0];      % mg/L
    E_prev_kgm3 = 0;
    
    E_plot = zeros(1,n);
    k_plot = zeros(1,n);
    A_mean_plot = zeros(1,n);
    sv_snap = []; 
    sh_snap = [];
        for i = 1:n
             % --- ODE evolve over [t_{i-1}, t_i] ---
        
        tspan_bin = [t_grid(i), t_grid(i+1)];
        
        opts = odeset('RelTol',1e-6,'AbsTol',1e-9);

        % unit-safe RHS (see function at end)
        [~, Y_sol] = ode15s(@(t,y) NAE(t, y, k_0, rho_EPS, ...
                                   phi_par, psi, nu, rho, xi, ...
                                   delta, eta, gamma), ...
                          tspan_bin, IC, opts);
        IC = Y_sol(end,:).';   % mg/L

      % EPS from ODE (mg/L) -> kg/m^3
        E_mgL  = max(IC(3), 0);
        E_kgm3 = 1e-3 * E_mgL;            % 1 mg/L = 1e-3 kg/m^3

        % incremental deposition (only if EPS increased)
        dEPS = max(0, E_kgm3 - E_prev_kgm3);
        E_prev_kgm3 = E_kgm3;
               % ---- Start from the SAME base geometry each EPS level ----
          
     % ---------------- Conductance-weighted EPS deposition ----------------
        % ---- Deposition (conductance/area-weighted) ----
            sum_v = sum(A_sv_acc(:)); sum_h = sum(A_sh_acc(:));
            if sum_v == 0 || sum_h == 0
                % fully clogged: stop and set remaining ks to ~0
                eff_val = 0;
                eff_trial_local = zeros(1, n);
                eff_trial_local(1:i-1) = k_plot(1:i-1); % keep the ones we already have
                eff_trial_local(i:end) = eff_val;
                % write to outputs after loop (below)
                % fill series
                k_plot(i:end) = eff_val;
                A_mean_plot(i:end) = mean([A_sv_acc(:); A_sh_acc(:)]);
                E_plot(i:end) = E_kgm3;
                sv_snap = sv; sh_snap = sh;
                break;
            end

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
            
                % Store for later inspection/plots
            E_plot(i) = E_kgm3;
            k_plot(i)=effcoe;
            A_mean_plot(i) = mean([A_sv_acc(:); A_sh_acc(:)]);
            if i == n
                sv_snap = sv;  
                sh_snap = sh;  
            end

            % Debug output
           fprintf('vf0=%.3f | Trial %d | step %3d/%3d | dEPS=%.3e kg/m^3 | cumEPS=%.3e | k=%.6e m^2\n', ...
                    vf0_current, t, i, n, dEPS, E_kgm3, effcoe);
        end

       
        % write per-trial series to global arrays (parfor-safe)
        eff_trial_all_perm(j,:,t) = k_plot;
        mean_area_all(j,:,t)      = A_mean_plot;

        trial_dump(j,t).sv    = sv_snap;
        trial_dump(j,t).sh    = sh_snap;
        trial_dump(j,t).E     = E_plot;
        trial_dump(j,t).k     = k_plot;
        trial_dump(j,t).Amean = A_mean_plot;
    end % parfor trials

    % Averages across trials (for this vf0)
    eff_mean_perm_plot(j,:) = mean(eff_trial_all_perm(j,:,:), 3);
    mean_area_plot(j,:)     = mean(mean_area_all(j,:,:), 3);
end % for j


fprintf('Total time elapsed: %.2f s\n', toc(t_all));


