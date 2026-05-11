%Boundary layer 1D PDE 
%Brine: 
% x=0: Dirichlet C_B(1)=C_in
% x=L: Neumenn dC_B/dx=0 -> ghost node  C_B_{N+1}= C_B_{N-1}
% EPS:
% x=0 and x=L: Neumenn, dC_E/dx=0 
% Using finite difference 

function dydt = pipe_PDE(t,y,params)

Nx = params.Nx;
x = params.x;
dx= params.dx;
a = params.a;
R= params.R;
D_B = params.D_B;
D_E = params.D_E;
U = params.U;
lambda = params.lambda;
k = params.k;
C_in_node = params.C_in_node;
C_node = params.C_node;

C_B = y(1:Nx);
C_E = y(Nx+1:2*Nx);

C_B = max(C_B,0);
C_E = max(C_E,0);
dC_Bdt = zeros(Nx,1);
dC_Edt = zeros(Nx,1);
beta_eps = 2*a/(R^2-a^2);
q = interface_flux(t,C_B,C_E,params);
    % ghost points for brine Robin BCs
    CB_left_ghost  = C_B(2) - (2*dx/D_B)*k*(C_in_node - C_B(1));
    CB_right_ghost = C_B(Nx-1) + (2*dx/D_B)*k*(C_node - C_B(Nx));

    % ghost points for EPS Neumann BCs
    CE_left_ghost  = C_E(2);
    CE_right_ghost = C_E(Nx-1);
% Interior PDE %
for i = 1:Nx
    %brine derivatives 
    if i== 1
        dC_Bdx = (C_B(2)-CB_left_ghost)/(2*dx);
        d2C_Bdx2 = (C_B(2)-2*C_B(1)+CB_left_ghost)/(dx^2);
    elseif i==Nx 
        dC_Bdx = (CB_right_ghost-C_B(Nx-1))/(2*dx);
        d2C_Bdx2 = (CB_right_ghost-2*C_B(Nx)+C_B(Nx-1))/(dx^2);
    else
        dC_Bdx = (C_B(i+1)-C_B(i-1))/(2*dx);
        d2C_Bdx2 = (C_B(i+1)-2*C_B(i)+C_B(i-1))/(dx^2);
    end
    %EPS derivatives
    if i==1 
        d2C_Edx2 = (C_E(2)-2*C_E(1)+CE_left_ghost)/(dx^2);
    elseif i==Nx
        d2C_Edx2= (CE_right_ghost - 2*C_E(Nx)+C_E(Nx-1))/(dx^2);
    else 
        d2C_Edx2 = (C_E(i+1)-2*C_E(i)+C_E(i-1))/(dx^2);
    end
dC_Bdt(i) = -U*dC_Bdx +D_B*d2C_Bdx2-(2/a)*q(i);
dC_Edt(i) = D_E*d2C_Edx2 -lambda*C_E(i)+beta_eps*q(i);

end
dydt=  [dC_Bdt;dC_Edt]