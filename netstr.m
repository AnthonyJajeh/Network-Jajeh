%obtain a sample of bonds with conductances from a log-normal distribution, namely:
%   ln(A) has distribution N(rmu,rsig*rsig)
%   where A is the area of the brine cross section and rmu and rsig
%   are prescribed.

%Once the sample is obtained, the conductances of the bonds
%and the volume fraction of the brines are readily computed.

function [sv, sh, vfrac] = netstr(nx,ny,rmu,rsig,h)

%These are cutoffs to avoid extremely large or small samples.

%xmax = rmu + 4.0*rsig;
%xmin = rmu - 4.0*rsig;

% sv = exp(rmu+rsig*randn(nx+1,ny+1));
% sh = exp(rmu+rsig*randn(nx+1,ny+1));
sv = ones(nx+1,ny+1);
sh = ones(nx+1,ny+1);
asum = 0;


