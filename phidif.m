% compute the effective properties by pdrop/sum (currents) 

function [eff1 eff2 effcoe] = phidif(u,sh,pdrop,n,h)

cleft = 0;
cright = 0;
for j=0:n-1
    cleft = cleft + (pdrop-u(end-j))*sh;
    cright = cright + u(j+1)*sh;
end

eff1 = cleft/pdrop*pi/h/h/8;
eff2 = cright/pdrop*pi/h/h/8;
effcoe = 0.5*(eff1+eff2);