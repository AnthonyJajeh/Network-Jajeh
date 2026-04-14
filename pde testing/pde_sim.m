clear; clc; close all

% Parameters
D_B = 1;
D_E = .02;
lambda = 0.1;

% Pipe 1 radii
a_1 = 1.0;
R_1 = 2.0;

% Pipe 2 radii
a_2 = 0.7;
R_2 = 1.6;
L=5;

C_in  = 1.0;
C_con = 0.5;
C_out = 0.0;

N_B = 10;     % number of brine radial nodes
N_E = 10;     % number of EPS radial nodes
N_x = 20;     % number of axial nodes

Delta_r_B1 = a_1/(N_B-1);
Delta_r_E1 = (R_1-a_1)/(N_E-1);
Delta_r_B2 = a_2/(N_B-1);
Delta_r_E2 = (R_2-a_2)/(N_E-1);
Delta_x   = L/(N_x-1);

r_B1 = linspace(0,a_1,N_B)';     % r_{B,j}
r_E1 = linspace(a_1,R_1,N_E)';     % r_{E,j}
r_B2 = linspace(0,a_2,N_B)';     % r_{B,j}
r_E1 = linspace(a_2,R_2,N_B)';     % r_{B,j}

x   = linspace(0,L,N_x);      % x_k

% Example velocity u(r) in the brine region
% Velocity profiles
u1 = 1 - (r_B1/a1).^2;
u2 = 1 - (r_B2/a2).^2;

params.D_B = D_B;
params.D_E = D_E;
params.lambda = lambda;
params.a_1 = a_1;
params.a_2= a_2;
params.R_1=R_1;
params.R_2 = R_2;
params.L = L;

params.C_in = C_in;
params.C_con = C_con;
params.C_out = C_out;

params.N_B = N_B;
params.N_E = N_E;
params.N_x = N_x;

params.Delta_r_B1 = Delta_r_B1;
params.Delta_r_E1 = Delta_r_E1;

params.Delta_r_B2 = Delta_r_B2;
params.Delta_r_E2 = Delta_r_E2;
params.Delta_x = Delta_x;

params.r_B1 = r_B1;
params.r_E1 = r_E1;

params.r_B2 = r_B2;
params.r_E2 = r_E2;
params.u = u;

% Initial conditions for stored unknowns
% B^1_{j,k}, j=1,...,N_B-1, k=2,...,N_x-1
B1_initial = zeros(N_B-1, N_x-2);

% E^1_{j,k}, j=2,...,N_E, k=2,...,N_x-1
E1_initial = zeros(N_E-1, N_x-2);

% B^2_{j,k}
B2_initial = zeros(N_B-1, N_x-2);

% E^2_{j,k}
E2_initial = zeros(N_E-1, N_x-2);

% Stack into one long vector
y0 = [B1_initial(:); E1_initial(:); B2_initial(:); E2_initial(:)];

tspan = [0 10];
opts = odeset('RelTol',1e-6,'AbsTol',1e-8);

[t,y] = ode23s(@(t,y) pde_rhs(t,y,params), tspan, y0, opts);

m = length(t);
y_m = y(m,:)';

n_B = (N_B-1)*(N_x-2);
n_E = (N_E-1)*(N_x-2);

B1_stored = reshape(y_m(1:n_B), N_B-1, N_x-2);
E1_stored = reshape(y_m(n_B+(1:n_E)), N_E-1, N_x-2);
B2_stored = reshape(y_m(n_B+n_E+(1:n_B)), N_B-1, N_x-2);
E2_stored = reshape(y_m(n_B+n_E+n_B+(1:n_E)), N_E-1, N_x-2);

B1_interface = (D_B*Delta_r_E*B1_stored(end,:) + D_E*Delta_r_B*E1_stored(1,:)) ...
               /(D_B*Delta_r_B + D_E*Delta_r_E);

E1_interface = B1_interface;

B2_interface = (D_B*Delta_r_E*B2_stored(end,:) + D_E*Delta_r_B*E2_stored(1,:)) ...
               /(D_B*Delta_r_B + D_E*Delta_r_E);

E2_interface = B2_interface;

x_int = x(2:N_x-1);
B1_plot = [B1_stored; B1_interface];
B2_plot = [B2_stored; B2_interface];

rB_plot = [r_B(1:N_B-1); a];

E1_plot = [E1_interface; E1_stored];
E2_plot = [E2_interface; E2_stored];

rE_plot = [a; r_E(2:N_E)];

figure;
imagesc(rB_plot,x_int, B1_plot);
set(gca,'YDir','normal');
xlabel('r_B');
ylabel('x');
title('Pipe 1 brine including interface');
colorbar;

figure;
imagesc(rB_plot,x_int, B2_plot);
set(gca,'YDir','normal');
xlabel('r_B');
ylabel('x');
title('Pipe 2 brine including interface');
colorbar;

figure;
imagesc(rE_plot,x_int, E1_plot);
set(gca,'YDir','normal');
xlabel('r_E');
ylabel('x');
title('Pipe 1 EPS including interface');
colorbar;

figure;
imagesc(rE_plot,x_int, E2_plot);
set(gca,'YDir','normal');
xlabel('r_E');
ylabel('x');
title('Pipe 2 EPS including interface');
colorbar;