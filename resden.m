%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Averaging the edge conductivity dvf and dhf in 4 neighboring
% cells in the fine grid to obtain the conductivity dvc and dhc
% in the coarse grid.
%
% lx = nx/2, ly = ny/2.
%
% Last modified: 5/15/2020

function [dvc,dhc] = resden(dvf,dhf)

[m,n] = size(dvf);
nx = m-1;
ny = n-1;
lx = nx/2;
ly = ny/2;

dhc = zeros(lx+1,ly+1);
dvc = zeros(lx+1,ly+1);

dhc(1:lx,1) = (dhf(1:2:nx-1,1) + dhf(2:2:nx,1) + dvf(2:2:nx,1))/3;
dhc(1:lx,2:ly) = .25*(dhf(1:2:nx-1,3:2:ny-1) + dhf(2:2:nx,3:2:ny-1) + ...
                      dvf(2:2:nx,3:2:ny-1) + dvf(2:2:nx,2:2:ny-2));
dhc(1:lx,ly+1) = (dhf(1:2:nx-1,ny+1) + dhf(2:2:nx,ny+1) + dvf(2:2:nx,ny))/3;


dvc(2:lx,1:ly) = .25*(dvf(3:2:nx-1,1:2:ny-1) + dvf(3:2:nx-1,2:2:ny) + ...
                      dhf(2:2:nx-2,2:2:ny) + dhf(3:2:nx-1,2:2:ny));
dvc(lx+1,1:ly) = (dvf(nx+1,1:2:ny-1) + dvf(nx+1,2:2:ny) + dhf(nx,2:2:ny))/3;

return