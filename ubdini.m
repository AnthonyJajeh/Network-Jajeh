%Initialize the field and prescribe boundary conditions for the refined
%grid.

function u = ubdini(nx,ny,nxy1,pdrop)

u = zeros(nxy1,1);

for j = 0:ny
    u((nx+1)*(j+1)) = pdrop;
end