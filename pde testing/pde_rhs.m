function dydt = pde_rhs(~,y,param)
N_B = param.N_B;
N_E = param.N_E;
N_x = param.N_x;
n_pipes = param.n_pipes;

D_B = param.D_B;
D_E = param.D_E;
lambda = param.lambda;

dx = param.Delta_x;

dr_B = param.Delta_r_B;
dr_E = param.Delta_r_E;

r_B_all = param.r_B;
r_E_all = param.r_E;
u_all   = param.u;

n_B = (N_B-1)*(N_x-2);
n_E = (N_E-1)*(N_x-2);

dydt_blocks = cell(2*n_pipes,1);

offset=0;

for p =1:n_pipes
    drB = dr_B(p);
    drE = dr_E(p);

    rB = r_B_all(:,p);
    rE = r_E_all(:,p);
    u  = u_all(:,p);
      % Unpack this pipe
    B_stored = reshape(y(offset + (1:n_B)), N_B-1, N_x-2);
    offset = offset + n_B;

    E_stored = reshape(y(offset + (1:n_E)), N_E-1, N_x-2);
    offset = offset + n_E;

    % Reconstruct full arrays
    B = zeros(N_B,N_x);
    E = zeros(N_E,N_x);

    B(1:N_B-1,2:N_x-1) = B_stored;
    E(2:N_E,2:N_x-1)   = E_stored;
  % Axial BCs for brine
    if p == 1
        B(:,1) = param.C_in;
    else
        B(:,1) = prev_pipe_outlet;
    end

    % Top outflow for brine interior rows only
    B(1:N_B-1,N_x) = B(1:N_B-1,N_x-1);

    % Axial BCs for EPS: zero x-derivative
    E(:,1)   = E(:,2);
    E(:,N_x) = E(:,N_x-1);

    % Interface reconstruction at r = a(p)
    for k = 2:N_x-1
        B(N_B,k) = (D_B*drE*B(N_B-1,k) + D_E*drB*E(2,k)) ...
                 / (D_B*drE + D_E*drB);
        E(1,k) = B(N_B,k);
    end

 % Boundary interface values
    B(N_B,1)   = B(N_B,2);
    B(N_B,N_x) = B(N_B,N_x-1);

    E(1,1)   = B(N_B,1);
    E(1,N_x) = B(N_B,N_x);

    % Save outlet of this pipe for next pipe
    prev_pipe_outlet = B(:,N_x);

    % Time derivatives
    dB_dt = zeros(N_B-1, N_x-2);
    dE_dt = zeros(N_E-1, N_x-2);

 %Brine equation
  for k = 2:N_x-1
        kk = k - 1;

        % centerline j = 1
        radial_term = 4*(B(2,k) - B(1,k))/drB^2;
        axial_diff  = (B(1,k+1) - 2*B(1,k) + B(1,k-1))/dx^2;
        axial_adv   = u(1)*(B(1,k) - B(1,k-1))/dx;

        dB_dt(1,kk) = D_B*(radial_term + axial_diff) - axial_adv;


        % interior j = 2,...,N_B-1
        for j = 2:N_B-1
            radial_term = (B(j+1,k) - 2*B(j,k) + B(j-1,k))/drB^2 ...
                        + (1/rB(j))*(B(j+1,k) - B(j-1,k))/(2*drB);

            axial_diff = (B(j,k+1) - 2*B(j,k) + B(j,k-1))/dx^2;
            axial_adv  = u(j)*(B(j,k) - B(j,k-1))/dx;

            dB_dt(j,kk) = D_B*(radial_term + axial_diff) - axial_adv;
        end
  end
 % -------------------------
    % EPS equation
    % -------------------------
    for k = 2:N_x-1
        kk = k - 1;

        % interior j = 2,...,N_E-1
        for j = 2:N_E-1
            radial_term = (E(j+1,k) - 2*E(j,k) + E(j-1,k))/drE^2 ...
                        + (1/rE(j))*(E(j+1,k) - E(j-1,k))/(2*drE);

            axial_diff = (E(j,k+1) - 2*E(j,k) + E(j,k-1))/dx^2;

            dE_dt(j-1,kk) = D_E*(radial_term + axial_diff) - lambda*E(j,k);
        end

        % outer wall j = N_E, no-flux
        j = N_E;
        radial_term = 2*(E(N_E-1,k) - E(N_E,k))/drE^2;
        axial_diff  = (E(j,k+1) - 2*E(j,k) + E(j,k-1))/dx^2;

        dE_dt(N_E-1,kk) = D_E*(radial_term + axial_diff) - lambda*E(j,k);
    end
    dydt_blocks{2*p-1} = dB_dt(:);
    dydt_blocks{2*p}   = dE_dt(:);
end

dydt = vertcat(dydt_blocks{:});
end