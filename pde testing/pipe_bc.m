function [pl,ql,pr,qr] = pipe_bc(xl,ul, xr,ur,t,params)
D_B = params.D_B;
k = params.k;
C_in_node = params.C_in_node;
C_node = params.C_node;

% pdepe boundary form:
% p + q*f = 0
%
% Since f1 = DB*CBx, to impose:
% DB*CBx(0,t) = k*(C_in_node - CB(0,t))
% we write:
% k*(C_in_node - CB_left) - f_left = 0
%
% Likewise at x=L:
% DB*CBx(L,t) = k*(C_node - CB(L,t))
% so
% k*(C_node - CB_right) - f_right = 0
 
%left BC
pl = [k*(C_in_node-ul(1));
    0];
ql = [1;
    1];
%Right BC

pr = [k*(C_node-ur(1));
    0];
qr = [-1;
    1];

end
