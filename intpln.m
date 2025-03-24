%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Interpolate u from the coarse grid to the fine grid
%  mx = 2*nx, my = 2*ny
%  
%  Last modified: 5/15/2020

function uf = intpln(nx,ny,uc,uf)

mx = 2*nx;
my = 2*ny;

%Bottom boundary points:
uadd = uc(2:nx,1);
uf(3:2:mx-1,1) = uf(3:2:mx-1,1) + uadd;
uf(4:2:mx,1) = uf(4:2:mx,1) + 0.5*uadd;
uf(2:2:mx-2,1) = uf(2:2:mx-2,1) + 0.5*uadd;
uf(3:2:mx-1,2) = uf(3:2:mx-1,2) + 0.5*uadd;
uf(2:2:mx-2,2) = uf(2:2:mx-2,2) + 0.25*uadd;
uf(4:2:mx,2) = uf(4:2:mx,2) + 0.25*uadd;

%Interior points:

uadd = uc(2:nx,2:ny);

uf(3:2:mx-1,3:2:my-1) = uf(3:2:mx-1,3:2:my-1) + uadd;
uf(4:2:mx,3:2:my-1) = uf(4:2:mx,3:2:my-1) + .5*uadd;
uf(2:2:mx-2,3:2:my-1) = uf(2:2:mx-2,3:2:my-1) + .5*uadd;

uf(3:2:mx-1,4:2:my) = uf(3:2:mx-1,4:2:my) + .5*uadd;
uf(4:2:mx,4:2:my) = uf(4:2:mx,4:2:my) + .25*uadd;
uf(2:2:mx-2,4:2:my) = uf(2:2:mx-2,4:2:my) + .25*uadd;


uf(3:2:mx-1,2:2:my-2) = uf(3:2:mx-1,2:2:my-2) + .5*uadd;
uf(4:2:mx,2:2:my-2) = uf(4:2:mx,2:2:my-2)+.25*uadd;
uf(2:2:mx-2,2:2:my-2) = uf(2:2:mx-2,2:2:my-2) + .25*uadd;



%Top boundary points:

uadd = uc(2:nx,ny+1);

uf(3:2:mx-1,my+1) = uf(3:2:mx-1,my+1) + uadd;
uf(4:2:mx,my+1) = uf(4:2:mx,my+1) + .5*uadd;
uf(2:2:mx-2,my+1) = uf(2:2:mx-2,my+1) + .5*uadd;

uf(3:2:mx-1,my) = uf(3:2:mx-1,my) + .5*uadd;
uf(4:2:mx,my) = uf(4:2:mx,my) + .25*uadd;
uf(2:2:mx-2,my) = uf(2:2:mx-2,my) + .25*uadd;

return
