clear all;
%**************************************************************************
%
%     Implementation of the network model to compute the effective properties 
%     such as conductivity and permeability of a nonhomogeneous material.
%
%     The resulting linear system is solved using a multigrid algorithm
%     (kikmul.m).
%
%     The random network coefficients are drawn from a log-normal
%     distribution specified by parameters (rmu, rsig). To be more precise,
%     we assume that cross sections are circular disks, and the areas have
%     a log-normal distribution:
%
%             A = e^X, X \sim N(rmu, rsig^2)
%
%     We note that
%             E[A] = e^(E[X]+0.5*Var(X))
% 
%     In this code, we only fix the average cross section area E[A] from the 
%     target volumn fraction, and assume a particular rsig value (entered 
%     as an input), which will impact the resulting average rmu.
%
%     The boundary conditions are set so that u=0 at the top
%     (index i=1 in this code) and u=pdrop at the bottom
%     (index i=nx+1).
%
%     Last modified: 3/23/2025
%
%*************************************************************************

      nx = input('enter nx (must be a power of 2):');
      pdrop = input('enter pressure drop:');
      vf0 = input('enter the desired volumn fraction:');
      rsig = input('enter the std of the log-normal dist:');
      ny = nx;

%     *******************************
%     The following is a reference sample for all parameters:      
%      n=256
%      pdrop = 1.d0
%      rmu = 0
%      rsig = 0.5
%      nv = 10    % number of V-cycles in the multigrid package.
%     *******************************

%     *******************************
%     The following is to specify the distribution parameters for the
%     connecting tube bond coefficients, and generate samples according 
%     to the target distribution.
%     
%     vf0 is the targeted volume fraction, and am0 is the average of the
%     cross-sectional areas, which is assumed to be a function of the target 
%     volumn fraction:
%      am0 = 3.d-2+0.13d0*(vf0-2.d-2)
      am0 = pi*(7.d-2 + 1.6d-1*vf0)^2;
%     so the average cross section area can lead to the average of cross
%     sectional radius:
      rmu = log(am0) - 0.5*rsig^2;
%     which is the mean of log area and
      amu = exp(2.0*rmu)*(exp(2.0*rsig^2)-exp(rsig^2));
%     which is the variance of the cross section areas, and cell size
      h = sqrt(2*am0/vf0);
      fprintf('std of cross section areas is %8.5f\n',sqrt(amu));

%     draw a sample of bond coefficients:
      sv = exp(rmu+rsig*randn(nx+1,ny+1));
      sh = exp(rmu+rsig*randn(nx+1,ny+1));
%    *************************************

%     multigrid algorithm to solve the elliptic system:

     [phi, final_error] = kikmul(pdrop,sv,sh,nx,ny);

%     compute the effective coefficients based on the solution:

      effleft = sum((pdrop-phi(nx,:)).*sh(nx,:))/pdrop/h^2;
      effright = sum(phi(2,:).*sh(1,:))/pdrop/h^2;
      effcoe = 0.5*(effleft+effright);


%     Output:

      fprintf('effcoe = %8.6f\n',effcoe);
      fprintf('rmu and rsig: %12.6f %12.6f\n',rmu,rsig);



