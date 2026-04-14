function dydt = pde_rhs(~,y,param)
N_B= param.N_B;
N_E = param.N_E;
N_x = param.N_x;
Delta_r_B = param.Delta_r_B;
Delta_r_E = param.Delta_r_E;
Delta_x   = param.Delta_x;

r_B = param.r_B;     % r_{B,j}
r_E = param.r_E;     % r_{E,j}
u   = param.u;       % u(r_{B,j})

D_B = param.D_B;
D_E = param.D_E;
lambda = param.lambda;

%Number of stored unknowns per pipe
n_B = (N_B-1)*(N_x-2);
n_E = (N_E-1)*(N_x-2);

% Block 1: B^1_{j,k}, j=1,...,N_B-1, k=2,...,N_x-1
% Block 2: E^1_{j,k}, j=2,...,N_E,   k=2,...,N_x-1
% Block 3: B^2_{j,k}, j=1,...,N_B-1, k=2,...,N_x-1
% Block 4: E^2_{j,k}, j=2,...,N_E,   k=2,...,N_x-1

I_B1 = 1:n_B;
I_E1 = n_B + (1:n_E);
I_B2 = n_B + n_E + (1:n_B);
I_E2 = n_B + n_E + n_B + (1:n_E);

% Unpack y into matrices of stored unknowns

B1_stored = reshape(y(I_B1), N_B-1, N_x-2);
E1_stored = reshape(y(I_E1), N_E-1, N_x-2);
B2_stored = reshape(y(I_B2), N_B-1, N_x-2);
E2_stored = reshape(y(I_E2), N_E-1, N_x-2);

B1 = zeros(N_B,N_x);
B2 = zeros(N_B,N_x);
E1 = zeros(N_E,N_x);
E2 = zeros(N_E,N_x);

B1(1:N_B-1,2:N_x-1) = B1_stored;
B2(1:N_B-1,2:N_x-1) = B2_stored;
E1(2:N_E,2:N_x-1) = E1_stored;
E2(2:N_E,2:N_x-1) = E2_stored;

%BC for brine
B1(:,1) = param.C_in;
B1(:,N_x)=param.C_con;

B2(:,1) = param.C_con;
B2(:,N_x) = param.C_out;

%BC for EPS
E1(:,1)= E1(:,2);
E1(:,N_x) = E1(:,N_x-1);

E2(:,1)= E2(:,2);
E2(:,N_x) = E2(:,N_x-1);

%at r=a matching 
% Continuity, E_{1,k} = B_{N_B,k}
% Derivative matching  = (D_BDelta_r_E B^i_{N_B-1,k} + D_EDelta_r_B E^i_{2,k})
%   /(D_BDelta_r_B + D_EDelta_r_E)

for k = 2:N_x-1
    B1(N_B,k) = (D_B*Delta_r_E*B1(N_B-1,k)+D_E*Delta_r_B*E1(2,k))/(D_B*Delta_r_B-D_E*Delta_r_E);
    E1(1,k)=B1(N_B,k);
    B2(N_B,k) = (D_B*Delta_r_E*B2(N_B-1,k)+D_E*Delta_r_B*E2(2,k))/(D_B*Delta_r_B-D_E*Delta_r_E);
    E2(1,k) = B2(N_B,k);
end
%Interface values at axial boundaries
E1(1,1)=B1(N_B,1);
E1(1,N_x)=B1(N_B,N_x);
E2(1,1)=B2(N_B,1);
E2(1,N_x)=B2(N_B,N_x);

%time derivative allocation
dB1_dt = zeros(N_B-1, N_x-2);
dE1_dt = zeros(N_E-1, N_x-2);
dB2_dt = zeros(N_B-1, N_x-2);
dE2_dt = zeros(N_E-1, N_x-2);

%pipe 1 brine 
for k = 2:N_x-1
    k_interior = k - 1;

    % ----------------------------
    % j = 1 : centerline r = 0
    % ----------------------------
    j = 1;

    radial_term = 4*(B1(2,k) - B1(1,k))/Delta_r_B^2;

    axial_diffusion = (B1(j,k+1) - 2*B1(j,k) + B1(j,k-1))/Delta_x^2;

    axial_advection = u(j)*(B1(j,k) - B1(j,k-1))/Delta_x;

    dB1_dt(j,k_interior) = D_B*(radial_term + axial_diffusion) - axial_advection;

    % ----------------------------
    % j = 2,...,N_B-1 : interior brine nodes
    % ----------------------------
    for j = 2:N_B-1
        radial_term = (B1(j+1,k) - 2*B1(j,k) + B1(j-1,k))/Delta_r_B^2 ...
                    + (1/r_B(j))*(B1(j+1,k) - B1(j-1,k))/(2*Delta_r_B);

        axial_diffusion = (B1(j,k+1) - 2*B1(j,k) + B1(j,k-1))/Delta_x^2;

        axial_advection = u(j)*(B1(j,k) - B1(j,k-1))/Delta_x;

        dB1_dt(j,k_interior) = D_B*(radial_term + axial_diffusion)- axial_advection;
                               
    end
end

%pipe 1 EPS
for k = 2:N_x-1
    k_interior = k - 1;

    % ----------------------------
    % j = 2,...,N_E-1 : interior EPS nodes
    % ----------------------------
    for j = 2:N_E-1
        radial_term = (E1(j+1,k) - 2*E1(j,k) + E1(j-1,k))/Delta_r_E^2 ...
                    + (1/r_E(j))*(E1(j+1,k) - E1(j-1,k))/(2*Delta_r_E);

        axial_diffusion = (E1(j,k+1) - 2*E1(j,k) + E1(j,k-1))/Delta_x^2;

        dE1_dt(j-1,k_interior) = D_E*(radial_term + axial_diffusion) ...
                                 - lambda*E1(j,k);
    end

    % ----------------------------
    % j = N_E : outer wall r = R
    % Use no-flux condition dE/dr = 0
    % ----------------------------
    j = N_E;

    radial_term = 2*(E1(N_E-1,k) - E1(N_E,k))/Delta_r_E^2;

    axial_diffusion = (E1(j,k+1) - 2*E1(j,k) + E1(j,k-1))/Delta_x^2;

    dE1_dt(N_E-1,k_interior) = D_E*(radial_term + axial_diffusion) ...
                               - lambda*E1(j,k);
end

%pipe 2 brine
for k = 2:N_x-1
    k_interior = k - 1;

    % j = 1
    j = 1;

    radial_term = 4*(B2(2,k) - B2(1,k))/Delta_r_B^2;

    axial_diffusion = (B2(j,k+1) - 2*B2(j,k) + B2(j,k-1))/Delta_x^2;

    axial_advection = u(j)*(B2(j,k) - B2(j,k-1))/Delta_x;

    dB2_dt(j,k_interior) = D_B*(radial_term + axial_diffusion) ...
                           - axial_advection;

    % j = 2,...,N_B-1
    for j = 2:N_B-1
        radial_term = (B2(j+1,k) - 2*B2(j,k) + B2(j-1,k))/Delta_r_B^2 ...
                    + (1/r_B(j))*(B2(j+1,k) - B2(j-1,k))/(2*Delta_r_B);

        axial_diffusion = (B2(j,k+1) - 2*B2(j,k) + B2(j,k-1))/Delta_x^2;

        axial_advection = u(j)*(B2(j,k) - B2(j,k-1))/Delta_x;

        dB2_dt(j,k_interior) = D_B*(radial_term + axial_diffusion) ...
                               - axial_advection;
    end
end

%pipe 2 EPS
for k = 2:N_x-1
    k_interior = k - 1;

    for j = 2:N_E-1
        radial_term = (E2(j+1,k) - 2*E2(j,k) + E2(j-1,k))/Delta_r_E^2 ...
                    + (1/r_E(j))*(E2(j+1,k) - E2(j-1,k))/(2*Delta_r_E);

        axial_diffusion = (E2(j,k+1) - 2*E2(j,k) + E2(j,k-1))/Delta_x^2;

        dE2_dt(j-1,k_interior) = D_E*(radial_term + axial_diffusion) ...
                                 - lambda*E2(j,k);
    end

    j = N_E;

    radial_term = 2*(E2(N_E-1,k) - E2(N_E,k))/Delta_r_E^2;

    axial_diffusion = (E2(j,k+1) - 2*E2(j,k) + E2(j,k-1))/Delta_x^2;

    dE2_dt(N_E-1,k_interior) = D_E*(radial_term + axial_diffusion) ...
                               - lambda*E2(j,k);
end
dydt = [dB1_dt(:); dE1_dt(:); dB2_dt(:); dE2_dt(:)];
end