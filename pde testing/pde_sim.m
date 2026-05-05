clear; clc; close all

% Parameters

D_B = 1e-9;
D_E = 1e-10;
lambda = 1;
mu_fluid = 1.8e-3;   % viscosity [Pa*s]

N_B = 10;     % number of brine radial nodes
N_E = 10;     % number of EPS radial nodes
N_x = 10;     % number of axial nodes

vf0_current = .25;
sigma_A = .1;
am0 = pi*(7.d-5 + 1.6e-4*vf0_current)^2; %[m^2]
mu_A = log(am0)-.5*sigma_A^2;

n_cols = 6; %number of column
n_pipes_per_col=10; %number of pipes per column
n_pipes = n_cols*n_pipes_per_col; %total number of pipes
A = (exp(mu_A+sigma_A*randn(1,n_pipes)))'; %logN cross sectional area
a = sqrt(A/pi);  %lognormal distribution of radius of pipes
R = a + 2e-5;
L   = sqrt((2 * am0) / vf0_current);

% =========================================================
% BUILD RADIAL GRIDS AND PRESSURE-DRIVEN VELOCITY PROFILES
% FOR STACKED VERTICAL COLUMNS
% =========================================================

r_B = zeros(N_B, n_pipes);
r_E = zeros(N_E, n_pipes);
u   = zeros(N_B, n_pipes);

Delta_r_B = zeros(1,n_pipes);
Delta_r_E = zeros(1,n_pipes);

Delta_x = L/(N_x-1);
x_local = linspace(0,L,N_x);

% Total pressure drop from bottom to top of each vertical column
p_drop = 1;  % [Pa]

u_max_pipe = zeros(n_pipes,1);
Q_col = zeros(n_cols,1);

% ---------------------------------------------------------
% Each column is a stack of pipes in series.
% Same Q through all pipes in one column.
% Same total pressure drop p_drop across each column.
% ---------------------------------------------------------
for c = 1:n_cols

    R_col = 0;

    % Total hydraulic resistance of this column
    for q = 1:n_pipes_per_col

        p = (c-1)*n_pipes_per_col + q;

        R_pipe = 8*mu_fluid*L/(pi*a(p)^4);

        R_col = R_col + R_pipe;

    end
    R_col_all(c) = R_col;
    % Flow rate through this column
    Q_col(c) = p_drop / R_col;

    % Convert column flow rate into pipe-specific u_max
    for q = 1:n_pipes_per_col

        p = (c-1)*n_pipes_per_col + q;

        A_pipe = pi*a(p)^2;

        U_avg_pipe = Q_col(c) / A_pipe;

        % Poiseuille profile has U_avg = u_max/2
        u_max_pipe(p) = 2*U_avg_pipe;

    end

end

% ---------------------------------------------------------
% Build radial grids and velocity profile for each pipe
% ---------------------------------------------------------
for p = 1:n_pipes

    Delta_r_B(p) = a(p)/(N_B-1);
    Delta_r_E(p) = (R(p)-a(p))/(N_E-1);

    r_B(:,p) = linspace(0,a(p),N_B)';
    r_E(:,p) = linspace(a(p),R(p),N_E)';

    % Pressure-driven Poiseuille velocity profile
    u(:,p) = u_max_pipe(p)*(1 - (r_B(:,p)/a(p)).^2);

end

%Initial condition concentrations 
C_in  = .01;

params.a= a;
params.R = R;
params.D_B = D_B;
params.D_E = D_E;
params.lambda = lambda;
params.L = L;

params.C_in = C_in;
params.JE_out =  1e-11;
params.JB_out =  1e-10;
params.N_B = N_B;
params.N_E = N_E;
params.N_x = N_x;
params.n_cols = n_cols;
params.n_pipes_per_col = n_pipes_per_col;
params.n_pipes = n_pipes;

params.Delta_r_B = Delta_r_B;
params.Delta_r_E = Delta_r_E;
params.Delta_x   = Delta_x;

params.r_B = r_B;
params.r_E = r_E;
params.u   = u;

% Initial conditions
n_B = (N_B-1)*(N_x-2);
n_E = (N_E-1)*(N_x-2);


% Stack into one long vector
y0 = zeros(n_pipes*(n_B+n_E),1);

tspan = linspace(0,5000);
steady_tol = 1e-6;
opts = odeset('RelTol',1e-4,'AbsTol',1e-6,'Events', @(t,y) steady_state_event(t,y,params,steady_tol));
tic
[t,y,te,ye,ie] = ode15s(@(t,y) pde_rhs(t,y,params), tspan, y0, opts);
toc
% Final time solution
% Final time solution

if ~isempty(te)
    fprintf('Steady state reached at t = %.6g\n', te(end));
else
    fprintf('Steady state not reached by final time t = %.6g\n', t(end));
end

dydt_final = pde_rhs(t(end), y(end,:)', params);
fprintf('Final ||dy/dt||_inf = %.6e\n', norm(dydt_final, inf));
y_m = y(end,:)';

% =========================================================
% PERMEABILITY CHANGE USING EPS-LAYER CONCENTRATION PROXY
% =========================================================
n_trials = 5;

rho_EPS  = 1500;     % EPS density [kg/m^3]

D_total  = n_pipes_per_col * L;  % total vertical height of one column [m]
A_sample = n_cols * L^2;         % sample cross-sectional area [m^2]

n_EPS = 200;
EPS_con = linspace(0,rho_EPS,n_EPS)';
vf0_values = [0.05 0.10 0.15 0.20 0.25];
m = length(vf0_values);
% Store mean curves for each vf0
eff_mean_perm_plot = zeros(m,n_EPS);
eff_std_perm_plot  = zeros(m,n_EPS);

mean_area_plot = zeros(m,n_EPS);
std_area_plot  = zeros(m,n_EPS);

alpha_hat_all = zeros(m,1);
k0_hat_all = zeros(m,1);
R2_log_all = zeros(m,1);
for j = 1:m
     vf0_current = vf0_values(j);

    % Recompute geometry statistics for this vf0
    am0 = pi*(7e-5 + 1.6e-4*vf0_current)^2;  % mean area [m^2]
    mu_A = log(am0) - 0.5*sigma_A^2;

    L = sqrt((2*am0)/vf0_current);

    D_total  = n_pipes_per_col * L;
    A_sample = n_cols * L^2;

    k_eff_all = zeros(n_EPS,n_trials);
    A_mean_all = zeros(n_EPS,n_trials);
    for trial = 1:n_trials
         % -----------------------------------------------------
        % Generate one base geometry for this trial
        % -----------------------------------------------------
        A_trial = (exp(mu_A + sigma_A*randn(1,n_pipes)))';
    
        % Start from the same base geometry for this EPS sweep
        A_acc = A_trial;
    
        EPS_prev = 0;
    
        for i = 1:n_EPS
            EPS_target = EPS_con(i);
        
            dEPS = max(0,EPS_target-EPS_prev);
            EPS_prev = EPS_target;
        
            % -----------------------------------------------------
            % 2. Deposit EPS volume into pipe cross-sectional areas
            % -----------------------------------------------------
            if dEPS > 0 && sum(A_acc(:)) > 0
        
                % Larger open areas receive more EPS
                W_v = A_acc / sum(A_acc(:));
        
                % Current total brine volume across all pipes
                vol_v = sum(A_acc(:)) * L;
        
                % Incremental EPS mass
                dmass_EPS_v = vol_v * dEPS;
        
                % Convert EPS mass to EPS volume and distribute
                Veps_v = (dmass_EPS_v / rho_EPS) * W_v;
        
                % Convert deposited volume to area loss
                dA_v = Veps_v / L;
        
                % Update open brine area
                A_acc = max(0, A_acc - dA_v);
        
            end
        
            A_mean_all(i,trial) = mean(A_acc);
        
            % -----------------------------------------------------
            % 3. Compute effective permeability
            % -----------------------------------------------------
            G_col = zeros(n_cols,1);
        
            for c = 1:n_cols
        
                R_col = 0;
        
                for q = 1:n_pipes_per_col
        
                    p = (c-1)*n_pipes_per_col + q;
        
                    a_eff = sqrt(A_acc(p)/pi);
                    a_eff = max(a_eff, 1e-12);
        
                    R_pipe = 8*mu_fluid*L/(pi*a_eff^4);
        
                    R_col = R_col + R_pipe;
        
                end
        
                G_col(c) = 1/R_col;
        
            end
        
            G_total = sum(G_col); %columns in parallel
        
            k_eff_all(i,trial) = mu_fluid * D_total * G_total / A_sample; %darcy equiv
        
        end
              fprintf('vf0 = %.2f | Finished permeability trial %d of %d\n', ...
            vf0_current, trial, n_trials);
    end
    
   % Average over trials for this vf0
    k_eff_mean = mean(k_eff_all, 2);
    k_eff_std  = std(k_eff_all, 0, 2);

    A_mean = mean(A_mean_all, 2);
    A_std  = std(A_mean_all, 0, 2);

    eff_mean_perm_plot(j,:) = k_eff_mean.';
    eff_std_perm_plot(j,:)  = k_eff_std.';

    mean_area_plot(j,:) = A_mean.';
    std_area_plot(j,:)  = A_std.';

    % -----------------------------------------------------
    % Exponential fit for this vf0
    % k(E) = k0 exp(-alpha E)
    % -----------------------------------------------------
    valid = isfinite(k_eff_mean) & k_eff_mean > 0;

    pfit = polyfit(EPS_con(valid), log(k_eff_mean(valid)), 1);

    alpha_hat = -pfit(1);
    k0_hat = exp(pfit(2));

    logK = log(k_eff_mean(valid));
    logK_hat = polyval(pfit, EPS_con(valid));

    SS_res = sum((logK - logK_hat).^2);
    SS_tot = sum((logK - mean(logK)).^2);
    R2_log = 1 - SS_res/SS_tot;

    alpha_hat_all(j) = alpha_hat;
    k0_hat_all(j) = k0_hat;
    R2_log_all(j) = R2_log;

end
% =========================================================
% EFFECTIVE PERMEABILITY VS BRINE VOLUME FRACTION
% Similar style to Zhu/Golden plot
% =========================================================

figure('Position',[100 100 850 650]);
hold on; grid on; box on

% Choose a few EPS concentrations to show
eps_indices = [1, round(n_EPS/4), round(n_EPS/2), n_EPS];

markerStyles = {'x', '+', 'o', 's'};
lineStyles = {'none', 'none', 'none', 'none'};

for r = 1:length(eps_indices)

    i = eps_indices(r);

    plot(vf0_values, eff_mean_perm_plot(:,i), ...
        'LineStyle', lineStyles{r}, ...
        'Marker', markerStyles{r}, ...
        'MarkerSize', 9, ...
        'LineWidth', 1.8, ...
        'DisplayName', sprintf('EPS = %.1f kg/m^3', EPS_con(i)));

end

set(gca, 'YScale', 'log');

xlabel('Brine volume fraction, \phi');
ylabel('Effective permeability, k_{eff} [m^2]');
title('Effective permeability vs brine volume fraction');

legend('Location','northwest');
ylim([1e-13 1e-8]);   % adjust if needed
xlim([min(vf0_values)*0.8, max(vf0_values)*1.1]);

hold off
% =========================================================
% EPS VS EFFECTIVE PERMEABILITY WITH FITTED CURVE
% ONE PANEL PER vf0
% =========================================================

figure('Position',[100 100 1200 700]);
rows = ceil(m/3);
cols = 3;
for j = 1:m
    subplot(rows, cols, j);
    hold on; grid on;
    
    % Actual mean permeability data
    plot(EPS_con, eff_mean_perm_plot(j,:), 'o-', ...
        'LineWidth', 1.5, ...
        'DisplayName','data');
    
    % Fitted exponential: k(E) = k0 exp(-alpha E)
    Cfine = linspace(min(EPS_con), max(EPS_con), 200);
    Kfit  = k0_hat_all(j) * exp(-alpha_hat_all(j) * Cfine);

    plot(Cfine, Kfit, 'r--', ...
        'LineWidth', 2, ...
        'DisplayName', sprintf('fit: k(E)=%.2e e^{-%.3gE}', ...
        k0_hat_all(j), alpha_hat_all(j)));

    % Theoretical exponential: k(E) = k(0) exp(-2E/rho_EPS)
    Ktheory = eff_mean_perm_plot(j,1) * exp(-2*Cfine/rho_EPS);

    plot(Cfine, Ktheory, 'k:', ...
        'LineWidth', 2.5, ...
        'DisplayName', 'theory: k(0)e^{-2E/\rho}');
    
    xlabel('EPS Concentration, E (kg/m^3)');
    ylabel('k(E) (m^2)');
    title(sprintf('\\phi = %.2f', vf0_values(j)));

    legend('Location','best');
    hold off;
end

sgtitle('EPS vs Effective Permeability with Exponential Fit');
% =========================================================
% EXPECTED AREA VS EFFECTIVE PERMEABILITY
% ONE PANEL PER vf0
% =========================================================

figure('Position',[100 100 1200 700]);

rows = ceil(m/3);
cols = 3;
for j = 1:m

    subplot(rows, cols, j);
    hold on; grid on;

    plot(mean_area_plot(j,:), eff_mean_perm_plot(j,:), ...
        'o-', ...
        'LineWidth', 1.5, ...
        'MarkerSize', 5);

    set(gca, 'YScale','log');

    xlabel('E[A] [m^2]');
    ylabel('k_{eff} [m^2]');
    title(sprintf('\\phi = %.2f', vf0_values(j)));

end

sgtitle('Expected area vs effective permeability for each volume fraction');
height_stretch=1;
% =========================================================
% EPS AVERAGE PLOT BY COLUMN
% =========================================================
figure('Position',[100 100 900 600]);
hold on

for c = 1:n_cols

    x_all = [];
    E_all = [];

    for q = 1:n_pipes_per_col

        % Global pipe index
        p = (c-1)*n_pipes_per_col + q;

        offset = (p-1)*(n_B+n_E);

        % Extract stored solution for this pipe
        B_stored = reshape(y_m(offset + (1:n_B)), N_B-1, N_x-2);
        E_stored = reshape(y_m(offset + n_B + (1:n_E)), N_E-1, N_x-2);

        % Vertical coordinate inside this column
        x_plot = x_local(2:N_x-1) + (q-1)*(L - 2*Delta_x);

        % Average nutrient in EPS across EPS thickness
        E_avg_x = mean(E_stored, 1);

        x_all = [x_all, x_plot];
        E_all = [E_all, E_avg_x];

    end

    plot(x_all, E_all, 'LineWidth', 2, ...
        'DisplayName', sprintf('Column %d', c));

end

xlabel('Vertical distance, x');
ylabel('Average nutrient in EPS');
title('Average nutrient concentration in EPS by column');
legend('Location','best');
grid on
hold off

% =========================================================
% EPS AVERAGE: EACH COLUMN + AVERAGE OVER ALL COLUMNS
% =========================================================

figure('Position',[100 100 900 600]);
hold on

E_by_col_full = zeros(n_cols, n_pipes_per_col*(N_x-2));
x_full = [];

% Build common vertical coordinate
for q = 1:n_pipes_per_col
    x_plot = x_local(2:N_x-1) + (q-1)*(L - 2*Delta_x);
    x_full = [x_full, x_plot];
end

for c = 1:n_cols

    E_col_full = [];

    for q = 1:n_pipes_per_col

        p = (c-1)*n_pipes_per_col + q;

        offset = (p-1)*(n_B+n_E);

        E_stored = reshape(y_m(offset + n_B + (1:n_E)), N_E-1, N_x-2);

        E_avg_x = mean(E_stored, 1);

        E_col_full = [E_col_full, E_avg_x];

    end

    E_by_col_full(c,:) = E_col_full;

    plot(x_full, E_col_full, '--', 'LineWidth', 1.2, ...
        'DisplayName', sprintf('Column %d', c));

end

E_mean_all_cols = mean(E_by_col_full, 1);

plot(x_full, E_mean_all_cols, 'k-', 'LineWidth', 3, ...
    'DisplayName', 'Average over all columns');

xlabel('Vertical distance, x');
ylabel('Average nutrient in EPS');
title('Average EPS nutrient concentration over all columns');
legend('Location','best');
grid on
hold off
% =========================================================
% MULTIPLE COLUMN FINAL PLOT
% =========================================================
figure('Position',[100 100 650 1000]);
hold on

cmin = inf;
cmax = -inf;

Rmax = max(R);
Nr_plot = 60;
r_plot = linspace(-Rmax, Rmax, Nr_plot);

% Horizontal spacing between columns
col_spacing = 3.0*Rmax;

% ---------------------------------------------------------
% Get common color limits over all pipes
% ---------------------------------------------------------
for p = 1:n_pipes

    offset = (p-1)*(n_B+n_E);

    B_stored = reshape(y_m(offset + (1:n_B)), N_B-1, N_x-2);
    E_stored = reshape(y_m(offset + n_B + (1:n_E)), N_E-1, N_x-2);

    drB = Delta_r_B(p);
    drE = Delta_r_E(p);

    B_interface = (D_B*drE*B_stored(end,:) + D_E*drB*E_stored(1,:)) ...
                / (D_B*drE + D_E*drB);

    F = build_display_field(B_stored,E_stored,B_interface, ...
        r_B(:,p),r_E(:,p),a(p),R(p),r_plot);

    cmin = min(cmin, min(F(:), [], 'omitnan'));
    cmax = max(cmax, max(F(:), [], 'omitnan'));

end

% ---------------------------------------------------------
% Plot each column side-by-side
% ---------------------------------------------------------
for c = 1:n_cols

    x_shift = (c-1)*col_spacing;

    for q = 1:n_pipes_per_col

        % Global pipe index
        p = (c-1)*n_pipes_per_col + q;

        offset = (p-1)*(n_B+n_E);

        B_stored = reshape(y_m(offset + (1:n_B)), N_B-1, N_x-2);
        E_stored = reshape(y_m(offset + n_B + (1:n_E)), N_E-1, N_x-2);

        drB = Delta_r_B(p);
        drE = Delta_r_E(p);

        B_interface = (D_B*drE*B_stored(end,:) + D_E*drB*E_stored(1,:)) ...
                    / (D_B*drE + D_E*drB);

        F = build_display_field(B_stored,E_stored,B_interface, ...
            r_B(:,p),r_E(:,p),a(p),R(p),r_plot);

        % Vertical position inside the column
        y_plot = height_stretch*(x_local(2:N_x-1) + (q-1)*(L - 2*Delta_x));

        % Horizontal shift for this column
        r_plot_shifted = r_plot + x_shift;

        imagesc(r_plot_shifted, y_plot, F);

        % Interface lines
        plot(x_shift + [ a(p)  a(p)], [y_plot(1) y_plot(end)], ...
            'w--', 'LineWidth', 1.5);

        plot(x_shift + [-a(p) -a(p)], [y_plot(1) y_plot(end)], ...
            'w--', 'LineWidth', 1.5);

        % Outer EPS wall lines
        plot(x_shift + [ R(p)  R(p)], [y_plot(1) y_plot(end)], ...
            'k-', 'LineWidth', 1.5);

        plot(x_shift + [-R(p) -R(p)], [y_plot(1) y_plot(end)], ...
            'k-', 'LineWidth', 1.5);

    end

    % Label each column
    text(x_shift, -0.08*n_pipes_per_col*L, sprintf('Column %d', c), ...
        'HorizontalAlignment','center', ...
        'FontWeight','bold');

end

set(gca,'YDir','normal');
axis tight

xlim([-Rmax, (n_cols-1)*col_spacing + Rmax]);
ylim(height_stretch*[x_local(2), x_local(end-1) + (n_pipes_per_col-1)*(L - 2*Delta_x)]);

if isfinite(cmin) && isfinite(cmax)
    if cmax > cmin
        caxis([cmin cmax]);
    else
        caxis([cmin cmin + 1e-12]);
    end
end

xlabel('column position / radius');
ylabel('vertical distance');
title('Nutrient concentration across multiple vertical columns', 'FontSize', 12);
colorbar
hold off

% =========================================================
% STACKED TIME ANIMATION ONLY
% =========================================================
figure('Position',[100 100 650 1000]);
v = VideoWriter('nutrient_stacked_pipes.mp4', 'MPEG-4');
v.FrameRate = 10;
open(v);
% =========================================================
% MULTIPLE COLUMN TIME ANIMATION
% =========================================================
cmin_t = inf;
cmax_t = -inf;

Rmax = max(R);
Nr_plot = 60;
r_plot = linspace(-Rmax, Rmax, Nr_plot);

col_spacing = 3.0*Rmax;

% ---------------------------------------------------------
% Get common color limits over all times and pipes
% ---------------------------------------------------------
for n = 1:length(t)

    y_n = y(n,:)';

    for p = 1:n_pipes

        offset = (p-1)*(n_B+n_E);

        B_stored = reshape(y_n(offset + (1:n_B)), N_B-1, N_x-2);
        E_stored = reshape(y_n(offset + n_B + (1:n_E)), N_E-1, N_x-2);

        drB = Delta_r_B(p);
        drE = Delta_r_E(p);

        B_interface = (D_B*drE*B_stored(end,:) + D_E*drB*E_stored(1,:)) ...
                    / (D_B*drE + D_E*drB);

        F = build_display_field(B_stored, E_stored, B_interface, ...
            r_B(:,p), r_E(:,p), a(p), R(p), r_plot);

        cmin_t = min(cmin_t, min(F(:), [], 'omitnan'));
        cmax_t = max(cmax_t, max(F(:), [], 'omitnan'));

    end
end

% ---------------------------------------------------------
% Make animation
% ---------------------------------------------------------
for n = 1:length(t)

    clf
    hold on

    y_n = y(n,:)';

    for c = 1:n_cols

        x_shift = (c-1)*col_spacing;

        for q = 1:n_pipes_per_col

            % Global pipe index
            p = (c-1)*n_pipes_per_col + q;

            offset = (p-1)*(n_B+n_E);

            B_stored = reshape(y_n(offset + (1:n_B)), N_B-1, N_x-2);
            E_stored = reshape(y_n(offset + n_B + (1:n_E)), N_E-1, N_x-2);

            drB = Delta_r_B(p);
            drE = Delta_r_E(p);

            B_interface = (D_B*drE*B_stored(end,:) + D_E*drB*E_stored(1,:)) ...
                        / (D_B*drE + D_E*drB);

            F = build_display_field(B_stored, E_stored, B_interface, ...
                r_B(:,p), r_E(:,p), a(p), R(p), r_plot);

            % Vertical position inside each column
            y_plot = height_stretch*(x_local(2:N_x-1) + (q-1)*(L - 2*Delta_x));

            % Horizontal shift for each column
            r_plot_shifted = r_plot + x_shift;

            imagesc(r_plot_shifted, y_plot, F);

            % Interface lines
            plot(x_shift + [ a(p)  a(p)], [y_plot(1) y_plot(end)], ...
                'w--', 'LineWidth', 1.5);

            plot(x_shift + [-a(p) -a(p)], [y_plot(1) y_plot(end)], ...
                'w--', 'LineWidth', 1.5);

            % Outer EPS wall lines
            plot(x_shift + [ R(p)  R(p)], [y_plot(1) y_plot(end)], ...
                'k-', 'LineWidth', 1.5);

            plot(x_shift + [-R(p) -R(p)], [y_plot(1) y_plot(end)], ...
                'k-', 'LineWidth', 1.5);

        end

        text(x_shift, -0.08*n_pipes_per_col*L, sprintf('Column %d', c), ...
            'HorizontalAlignment','center', ...
            'FontWeight','bold');

    end

    set(gca,'YDir','normal');
    axis tight

    xlim([-Rmax, (n_cols-1)*col_spacing + Rmax]);
    ylim(height_stretch*[x_local(2), x_local(end-1) + (n_pipes_per_col-1)*(L - 2*Delta_x)]);

    if isfinite(cmin_t) && isfinite(cmax_t)
        if cmax_t > cmin_t
            caxis([cmin_t cmax_t]);
        else
            caxis([cmin_t cmin_t + 1e-12]);
        end
    end

    xlabel('column position / radius');
    ylabel('vertical distance');

    h = title(sprintf('Nutrient concentration across multiple columns, t = %.2f', t(n)), ...
        'FontSize', 12);

    h.Units = 'normalized';
    h.Position(2) = 1.05;

    colorbar
    hold off
    drawnow

end

% =========================================================
% Local function
% =========================================================
function F = build_display_field(B_stored,E_stored,B_interface,rB,rE,a_val,R_val,r_plot)

    Nx_int = size(B_stored,2);
    Nr_plot = numel(r_plot);

    B_prof = [B_stored; B_interface];
    E_prof = [B_interface; E_stored];

    F = nan(Nx_int, Nr_plot);

    rr = abs(r_plot);
    brine_mask = rr <= a_val;
    eps_mask   = (rr > a_val) & (rr <= R_val);

    rr_brine = rr(brine_mask);
    rr_eps   = rr(eps_mask);

    for k = 1:Nx_int
        if ~isempty(rr_brine)
            F(k,brine_mask) = interp1(rB, B_prof(:,k), rr_brine, 'linear', 'extrap');
        end
        if ~isempty(rr_eps)
            % for debugging: only use EPS data itself
            F(k,eps_mask) = interp1(rE(2:end), E_stored(:,k), rr_eps, 'nearest', 'extrap');
        end
    end
end
function [value,isterminal,direction] = steady_state_event(t,y,params,steady_tol)

    dydt = pde_rhs(t,y,params);
    value = norm(dydt, inf) - steady_tol;

    isterminal = 1;
    direction = -1;
end