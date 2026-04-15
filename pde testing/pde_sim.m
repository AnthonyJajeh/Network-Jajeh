clear; clc; close all

% Parameters
D_B = 1;
D_E = .02;
lambda = 0.1;
L=.5;

N_B = 10;     % number of brine radial nodes
N_E = 10;     % number of EPS radial nodes
N_x = 50;     % number of axial nodes


% Pipe radii: one entry per pipe
a = [.05, .075, .1, .125];     % brine/interface radii
R = [.25, .25, .25, .25];     % outer EPS radii
n_pipes = numel(a);

r_B = zeros(N_B, n_pipes);
r_E = zeros(N_E, n_pipes);
u   = zeros(N_B, n_pipes);

Delta_r_B = zeros(1,n_pipes);
Delta_r_E = zeros(1,n_pipes);
Delta_x =L/(N_x-1);
x_local   = linspace(0,L,N_x);
for p = 1:n_pipes
    Delta_r_B(p) = a(p)/(N_B-1);
    Delta_r_E(p) = (R(p)-a(p))/(N_E-1);

    r_B(:,p) = linspace(0,a(p),N_B)';
    r_E(:,p) = linspace(a(p),R(p),N_E)';

    % Example velocity profile in brine region
    u(:,p) = 1 - (r_B(:,p)/a(p)).^2;
end

%Initial condition concentrations 
C_in  = 1.0;

params.a= a;
params.R = R;
params.D_B = D_B;
params.D_E = D_E;
params.lambda = lambda;
params.L = L;

params.C_in = C_in;
params.C_out_E = 0;
params.C_out_B = 0;
params.N_B = N_B;
params.N_E = N_E;
params.N_x = N_x;
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

tspan = [0 100];
opts = odeset('RelTol',1e-5,'AbsTol',1e-7);
tic
[t,y] = ode15s(@(t,y) pde_rhs(t,y,params), tspan, y0, opts);
toc
% Final time solution
% Final time solution
y_m = y(end,:)';

% Common display grid in radius
Rmax = max(R);
Nr_plot = 250;
r_plot = linspace(-Rmax, Rmax, Nr_plot);

% ========= final-time single-pipe plots =========
for p = 1:n_pipes

    offset = (p-1)*(n_B+n_E);

    B_stored = reshape(y_m(offset + (1:n_B)), N_B-1, N_x-2);
    E_stored = reshape(y_m(offset + n_B + (1:n_E)), N_E-1, N_x-2);

    drB = Delta_r_B(p);
    drE = Delta_r_E(p);

    B_interface = (D_B*drE*B_stored(end,:) + D_E*drB*E_stored(1,:)) ...
                / (D_B*drE + D_E*drB);

    F = build_display_field(B_stored,E_stored,B_interface,r_B(:,p),r_E(:,p),a(p),R(p),r_plot);

    x_plot = x_local(2:N_x-1) + (p-1)*L;

    figure;
    imagesc(r_plot, x_plot, F);
    set(gca,'YDir','normal');
    xlim([-R(p) R(p)]);
    ylim([x_plot(1) x_plot(end)]);
    xlabel('r');
    ylabel('x');
    title(sprintf('Pipe %d nutrient everywhere', p));
    colorbar
    hold on
    plot([ a(p)  a(p)], [x_plot(1) x_plot(end)], 'w--', 'LineWidth', 1.5);
    plot([-a(p) -a(p)], [x_plot(1) x_plot(end)], 'w--', 'LineWidth', 1.5);
    plot([ R(p)  R(p)], [x_plot(1) x_plot(end)], 'k-', 'LineWidth', 1.5);
    plot([-R(p) -R(p)], [x_plot(1) x_plot(end)], 'k-', 'LineWidth', 1.5);
    hold off
end

% ========= stacked plot =========
figure;
hold on

cmin = inf;
cmax = -inf;

% get common color limits
for p = 1:n_pipes
    offset = (p-1)*(n_B+n_E);

    B_stored = reshape(y_m(offset + (1:n_B)), N_B-1, N_x-2);
    E_stored = reshape(y_m(offset + n_B + (1:n_E)), N_E-1, N_x-2);

    drB = Delta_r_B(p);
    drE = Delta_r_E(p);

    B_interface = (D_B*drE*B_stored(end,:) + D_E*drB*E_stored(1,:)) ...
                / (D_B*drE + D_E*drB);

    F = build_display_field(B_stored,E_stored,B_interface,r_B(:,p),r_E(:,p),a(p),R(p),r_plot);

    cmin = min(cmin, min(F(:), [], 'omitnan'));
    cmax = max(cmax, max(F(:), [], 'omitnan'));
end

for p = 1:n_pipes
    offset = (p-1)*(n_B+n_E);

    B_stored = reshape(y_m(offset + (1:n_B)), N_B-1, N_x-2);
    E_stored = reshape(y_m(offset + n_B + (1:n_E)), N_E-1, N_x-2);

    drB = Delta_r_B(p);
    drE = Delta_r_E(p);

    B_interface = (D_B*drE*B_stored(end,:) + D_E*drB*E_stored(1,:)) ...
                / (D_B*drE + D_E*drB);

    F = build_display_field(B_stored,E_stored,B_interface,r_B(:,p),r_E(:,p),a(p),R(p),r_plot);

    x_plot = x_local(2:N_x-1) + (p-1)*L;

    imagesc(r_plot, x_plot, F);

    plot([ a(p)  a(p)], [x_plot(1) x_plot(end)], 'w--', 'LineWidth', 1.5);
    plot([-a(p) -a(p)], [x_plot(1) x_plot(end)], 'w--', 'LineWidth', 1.5);
    plot([ R(p)  R(p)], [x_plot(1) x_plot(end)], 'k-', 'LineWidth', 1.5);
    plot([-R(p) -R(p)], [x_plot(1) x_plot(end)], 'k-', 'LineWidth', 1.5);
end

set(gca,'YDir','normal');
xlim([-Rmax Rmax]);
ylim([x_local(2), x_local(end-1) + (n_pipes-1)*L]);

if isfinite(cmin) && isfinite(cmax)
    if cmax > cmin
        caxis([cmin cmax]);
    else
        caxis([cmin cmin + 1e-12]);
    end
end

xlabel('r');
ylabel('global x');
title('Nutrient everywhere across stacked pipes');
colorbar
hold off

% =========================================================
% STACKED TIME ANIMATION ONLY
% =========================================================
figure;

cmin_t = inf;
cmax_t = -inf;

% common color limits over all times and all pipes
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

for n = 1:length(t)
    clf
    hold on

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

        x_plot = x_local(2:N_x-1) + (p-1)*L;

        imagesc(r_plot, x_plot, F);

        plot([ a(p)  a(p)], [x_plot(1) x_plot(end)], 'w--', 'LineWidth', 1.5);
        plot([-a(p) -a(p)], [x_plot(1) x_plot(end)], 'w--', 'LineWidth', 1.5);
        plot([ R(p)  R(p)], [x_plot(1) x_plot(end)], 'k-', 'LineWidth', 1.5);
        plot([-R(p) -R(p)], [x_plot(1) x_plot(end)], 'k-', 'LineWidth', 1.5);
    end

    set(gca,'YDir','normal');
    xlim([-Rmax Rmax]);
    ylim([x_local(2), x_local(end-1) + (n_pipes-1)*L]);

    if isfinite(cmin_t) && isfinite(cmax_t)
        if cmax_t > cmin_t
            caxis([cmin_t cmax_t]);
        else
            caxis([cmin_t cmin_t + 1e-12]);
        end
    end

    xlabel('r');
    ylabel('global x');
    title(sprintf('Nutrient everywhere across stacked pipes, t = %.2f', t(n)));
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
            F(k,eps_mask) = interp1(rE, E_prof(:,k), rr_eps, 'linear', 'extrap');
        end
    end
end