clear; clc; close all

% Parameters
D_B = 1;
D_E = .02;
lambda = 0.1;
L=5;

N_B = 50;     % number of brine radial nodes
N_E = 50;     % number of EPS radial nodes
N_x = 50;     % number of axial nodes


% Pipe radii: one entry per pipe
a = [1.0, 1.3];     % brine/interface radii
R = [2.0, 1.7];     % outer EPS radii
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

tspan = [0 10];
opts = odeset('RelTol',1e-5,'AbsTol',1e-7);
tic
[t,y] = ode15s(@(t,y) pde_rhs(t,y,params), tspan, y0, opts);
toc
% Final time solution
y_m = y(end,:)';

for p = 1:n_pipes

    offset = (p-1)*(n_B+n_E);

    B_stored = reshape(y_m(offset + (1:n_B)), N_B-1, N_x-2);
    E_stored = reshape(y_m(offset + n_B + (1:n_E)), N_E-1, N_x-2);

    drB = Delta_r_B(p);
    drE = Delta_r_E(p);

    % Interface nutrient value at r = a(p)
    B_interface = (D_B*drE*B_stored(end,:) + D_E*drB*E_stored(1,:)) ...
                / (D_B*drE + D_E*drB);

    % Right half of the pipe: brine + interface + EPS
    N_right = [B_stored; B_interface; E_stored(2:end,:)];
    r_right = [r_B(1:N_B-1,p); a(p); r_E(2:end,p)];

    % Mirror to full pipe
    N_full = [flipud(N_right(2:end,:)); N_right];
    r_full = [-flipud(r_right(2:end)); r_right];

    % Axial location of this pipe in the global domain
    x_plot = x_local(2:N_x-1) + (p-1)*L;

    figure;
   imagesc(r_full, x_plot, N_full');
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

% -------- Combined stacked plot for all pipes --------
figure;
hold on

% global radial extent across all pipes
Rmax = max(R);

% optional: common color scale from final-time data
cmin = inf;
cmax = -inf;

% first pass: get common color limits
for p = 1:n_pipes
    offset = (p-1)*(n_B+n_E);

    B_stored = reshape(y_m(offset + (1:n_B)), N_B-1, N_x-2);
    E_stored = reshape(y_m(offset + n_B + (1:n_E)), N_E-1, N_x-2);

    drB = Delta_r_B(p);
    drE = Delta_r_E(p);

    B_interface = (D_B*drE*B_stored(end,:) + D_E*drB*E_stored(1,:)) ...
                / (D_B*drE + D_E*drB);

    N_right = [B_stored; B_interface; E_stored(2:end,:)];
    N_full  = [flipud(N_right(2:end,:)); N_right];

    cmin = min(cmin, min(N_full(:)));
    cmax = max(cmax, max(N_full(:)));
end

% second pass: draw each pipe on its own x-interval
for p = 1:n_pipes
    offset = (p-1)*(n_B+n_E);

    B_stored = reshape(y_m(offset + (1:n_B)), N_B-1, N_x-2);
    E_stored = reshape(y_m(offset + n_B + (1:n_E)), N_E-1, N_x-2);

    drB = Delta_r_B(p);
    drE = Delta_r_E(p);

    B_interface = (D_B*drE*B_stored(end,:) + D_E*drB*E_stored(1,:)) ...
                / (D_B*drE + D_E*drB);

    % right half: brine core + interface + EPS layer
    N_right = [B_stored; B_interface; E_stored(2:end,:)];
    r_right = [r_B(1:N_B-1,p); a(p); r_E(2:end,p)];

    % mirror to full pipe
    N_full = [flipud(N_right(2:end,:)); N_right];
    r_full = [-flipud(r_right(2:end)); r_right];

    % global x-location for this pipe
    x_plot = x_local(2:N_x-1) + (p-1)*L;

    % draw the data
    imagesc(r_full, x_plot, N_full');
    
    % draw interface and wall lines
    plot([ a(p)  a(p)], [x_plot(1) x_plot(end)], 'w--', 'LineWidth', 1.5);
    plot([-a(p) -a(p)], [x_plot(1) x_plot(end)], 'w--', 'LineWidth', 1.5);
    plot([ R(p)  R(p)], [x_plot(1) x_plot(end)], 'k-', 'LineWidth', 1.5);
    plot([-R(p) -R(p)], [x_plot(1) x_plot(end)], 'k-', 'LineWidth', 1.5);
end

set(gca,'YDir','normal');
axis tight
xlim([-Rmax Rmax]);
ylim([x_local(2) x_local(end-1) + (n_pipes-1)*L]);
caxis([cmin cmax]);

xlabel('r');
ylabel('global x');
title('Nutrient everywhere across stacked pipes');
colorbar
hold off

% -------- Time animation for nutrient everywhere --------
for p = 1:n_pipes   % choose which pipe to animate
figure;

% optional fixed color scale
cmin = inf;
cmax = -inf;
for n = 1:length(t)
    y_n = y(n,:)';

    offset = (p-1)*(n_B+n_E);
    B_stored = reshape(y_n(offset + (1:n_B)), N_B-1, N_x-2);
    E_stored = reshape(y_n(offset + n_B + (1:n_E)), N_E-1, N_x-2);

    drB = Delta_r_B(p);
    drE = Delta_r_E(p);

    B_interface = (D_B*drE*B_stored(end,:) + D_E*drB*E_stored(1,:)) ...
                / (D_B*drE + D_E*drB);

    N_right = [B_stored; B_interface; E_stored(2:end,:)];
    N_full = [flipud(N_right(2:end,:)); N_right];

    cmin = min(cmin, min(N_full(:)));
    cmax = max(cmax, max(N_full(:)));
end

for n = 1:length(t)
    y_n = y(n,:)';

    offset = (p-1)*(n_B+n_E);

    B_stored = reshape(y_n(offset + (1:n_B)), N_B-1, N_x-2);
    E_stored = reshape(y_n(offset + n_B + (1:n_E)), N_E-1, N_x-2);

    drB = Delta_r_B(p);
    drE = Delta_r_E(p);

    B_interface = (D_B*drE*B_stored(end,:) + D_E*drB*E_stored(1,:)) ...
                / (D_B*drE + D_E*drB);

    N_right = [B_stored; B_interface; E_stored(2:end,:)];
    r_right = [r_B(1:N_B-1,p); a(p); r_E(2:end,p)];

    N_full = [flipud(N_right(2:end,:)); N_right];
    r_full = [-flipud(r_right(2:end)); r_right];

    x_plot = x_local(2:N_x-1) + (p-1)*L;

    clf
    imagesc(r_full, x_plot, N_full');
    set(gca,'YDir','normal');
    caxis([cmin cmax]);
    xlabel('r');
    ylabel('global x');
    title(sprintf('Pipe %d nutrient everywhere, t = %.2f', p, t(n)));
    colorbar
    hold on
    plot([ a(p)  a(p)], [x_plot(1) x_plot(end)], 'w--', 'LineWidth', 1.5);
    plot([-a(p) -a(p)], [x_plot(1) x_plot(end)], 'w--', 'LineWidth', 1.5);
    plot([ R(p)  R(p)], [x_plot(1) x_plot(end)], 'k-', 'LineWidth', 1.5);
    plot([-R(p) -R(p)], [x_plot(1) x_plot(end)], 'k-', 'LineWidth', 1.5);
    hold off
    drawnow
end
end