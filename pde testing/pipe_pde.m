function [c,f,s] = pipe_pde(x,t,u,dudx,params)

C_B = u(1);
C_E = u(2);

C_Bx = dudx(1);
C_Ex = dudx(2);

D_B = params.D_B;
D_E = params.D_E;
U = params.U;
a = params.a;
R=  params.R;
lambda = params.lambda;
gamma = params.gamma;

q = gamma*(C_B-C_E);

beta_EPS = 2*a/(R^2-a^2);

%pde of form : c*u_t = d/dx(f)+s

c = [1;1];

f=[D_B*C_Bx;
    D_E*C_Ex];

s = [-U*C_Bx-(2/a)*q;
    -lambda*C_E+beta_EPS*q];
end
