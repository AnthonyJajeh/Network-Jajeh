%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  V-cycle multigrid algorithm to solve the linear system from the network model
%
%  Inputs:
%   fin     External sources into the nodes for the finest grid,
%           zero for this problem.  For the coarse grids, it is
%           used for the residual carried from fine to the current
%           grid.
%   pdrop   Potential drop = u_{right end} - u_{left end}
%   sv      Random conductivity of the vertical edges
%   sh      Random conductivity of the horizontal edges
%  (nx,ny)  Size of the finest grid
%
%  Outputs:
%   uf             solution at the finest level
%   final_error    final residual
%
%  Crucial parameters fixed in this code:
%   eps     Tolerance for L-2 residual
%   nv      Max number of V-cycles
%   kit     Number of iterations for each level

%  Variable used within this code:
%  (nx,ny)  current grid size
%  (mx,my)  next fine grid size
%  (lx,ly)  next coarse grid size
%
%  Grid structure:
%
%               x--------x--------x
%               |                 |
%               |                 |
%               |                 |
%               x <--sv(i,j)      x
%               |                 |
%               |                 |
%               |                 |
%    u(i,j) --> x--------x--------x
%                        ^
%                        |
%                      sh(i,j)
%
%   Last modified: 5/15/2020 by Jingyi Zhu

function [uf, final_error] = kikmul(fin,pdrop,sv,sh,nx,ny)

kit = 8;           % number of iterations for each GS step
nv = 10;           % number of V-cycles
eps = 1e-8;        % tolerance for convergence
nlevel = log2(ny); % number of levels in the multigrid algorithm

%Storage: all contained in the structure array um

um{1}.u = zeros(nx+1,ny+1);
um{1}.u(nx+1,1:ny+1) = pdrop;    % Set BC for u
um{1}.f = fin;
um{1}.sv = sv;
um{1}.sh = sh;
um{1}.nx = nx;
um{1}.ny = ny;

for j=2:nlevel
   um{j}.nx = um{j-1}.nx/2;
   um{j}.ny = um{j-1}.ny/2;
   um{j}.u = zeros(um{j}.nx+1,um{j}.ny+1);
   um{j}.f = zeros(um{j}.nx+1,um{j}.ny+1);    
   [um{j}.sv, um{j}.sh] = resden(um{j-1}.sv, um{j-1}.sh);
end

%Main V-Cycle:

kv = 1;
el2 = 1;

while (kv <= nv) && (el2 >= eps)

%   Descending steps:
    for j = 1:nlevel-1
        nx = um{j}.nx;
        ny = um{j}.ny;
        um{j}.u = gsnetw(nx,ny,um{j}.u,um{j}.f,um{j}.sv,um{j}.sh,kit);
        [res,el2] = residu(nx,ny,um{j}.u,um{j}.f,um{j}.sv,um{j}.sh);
        um{j+1}.f = restrt(nx,ny,res);
    end
    
%       Solve the coarse grid equations exactly by interting the matrix:
        um{nlevel}.u = dirsol(um{nlevel}.f,um{nlevel}.sv,um{nlevel}.sh);

%   Ascending steps:
    for j = nlevel:-1:2
        nx = um{j}.nx;
        ny = um{j}.ny;
        um{j-1}.u = intpln(nx,ny,um{j}.u,um{j-1}.u);
        um{j}.u = zeros(nx+1,ny+1);
        nnx = 2*nx;
        nny = 2*ny;
        um{j-1}.u = gsnetw(nnx,nny,um{j-1}.u,um{j-1}.f,um{j-1}.sv,um{j-1}.sh,kit);
    end


%Compute the final residue:
    [res,el2] = residu(um{1}.nx,um{1}.ny,um{1}.u,um{1}.f,um{1}.sv,um{1}.sh);
    el2 = sqrt(el2)*ny;
    %fprintf('residual at V-cycle number %3.0f is %e\n',kv,el2);
    kv = kv+1;
end

uf = um{1}.u;
surf(uf);

final_error = el2;