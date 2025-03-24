%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the residual of the equations
% 
% Last modified: 5/15/2020

function [res, el2] = residu(nx,ny,u,f,sv,sh)

res(2:nx,1) = f(2:nx,1) + (sh(1:nx-1,1)+sh(2:nx,1)+sv(2:nx,1)).*u(2:nx,1) -...
                           sh(2:nx,1).*u(3:nx+1,1) - sh(1:nx-1,1).*u(1:nx-1,1) - ...
                           sv(2:nx,1).*u(2:nx,2);

res(2:nx,2:ny) = f(2:nx,2:ny) + (sh(1:nx-1,2:ny)+sh(2:nx,2:ny)+sv(2:nx,1:ny-1)+sv(2:nx,2:ny)).*u(2:nx,2:ny) - ...
                                 sh(2:nx,2:ny).*u(3:nx+1,2:ny) - sh(1:nx-1,2:ny).*u(1:nx-1,2:ny) - ...
                                 sv(2:nx,2:ny).*u(2:nx,3:ny+1) - sv(2:nx,1:ny-1).*u(2:nx,1:ny-1);

res(2:nx,ny+1) = f(2:nx,ny+1) + (sh(1:nx-1,ny+1)+sh(2:nx,ny+1)+sv(2:nx,ny)).*u(2:nx,ny+1) -...
                                 sh(2:nx,ny+1).*u(3:nx+1,ny+1) - sh(1:nx-1,ny+1).*u(1:nx-1,ny+1) - ...
                                 sv(2:nx,ny).*u(2:nx,ny);
 
el2 = sum(sum(res(2:nx,1:ny+1).^2,1));
return
