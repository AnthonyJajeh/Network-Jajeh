%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Restrict the residual from a fine grid to the next coarse grid
%  lx = nx/2, ly = ny/2.
%
%  Last modified: 5/15/2020

function rc = restrt(nx,ny,rf)

lx = nx/2;
ly = ny/2;

rc = zeros(lx+1,ly+1);

rc(2:lx,1) = rf(3:2:nx-1,1) + .5*(rf(4:2:nx,1) + rf(2:2:nx-2,1) + rf(3:2:nx-1,2)) +...
                              .25*(rf(2:2:nx-2,2) + rf(4:2:nx,2));
                          
rc(2:lx,2:ly) = rf(3:2:nx-1,3:2:ny-1) + .5*(rf(4:2:nx,3:2:ny-1) + rf(2:2:nx-2,3:2:ny-1) + ...
                                            rf(3:2:nx-1,4:2:ny) + rf(3:2:nx-1,2:2:ny-2)) + ...
                                       .25*(rf(4:2:nx,4:2:ny) + rf(2:2:nx-2,4:2:ny) + ...
                                            rf(4:2:nx,2:2:ny-2) + rf(2:2:nx-2,2:2:ny-2));

rc(2:lx,ly+1) = rf(3:2:nx-1,ny+1) + .5*(rf(4:2:nx,ny+1) + rf(2:2:nx-2,ny+1) + rf(3:2:nx-1,ny)) +...
                                   .25*(rf(2:2:nx-2,ny) + rf(4:2:nx,ny));
return
