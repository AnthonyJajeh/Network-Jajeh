% Solving the coarest (3X3) grid exactly by inverting the matrix.

function u = dirsol(f,sv,sh)

u = zeros(3,3);

r1 = sv(2,2) + sh(1,3) + sh(2,3);
r2 = sv(2,1) + sh(1,1) + sh(2,1);
u(2,2) = (f(2,2) + sv(2,2)*f(2,3)/r1 + sv(2,1)*f(2,1)/r2)/...
         (sv(2,2)^2/r1 + sv(2,1)^2/r2 - sv(2,2) - sv(2,1) - sh(1,2) - sh(2,2));
u(2,3) = (sv(2,2)*u(2,2) - f(2,3))/r1;
u(2,1) = (sv(2,1)*u(2,2) - f(2,1))/r2;