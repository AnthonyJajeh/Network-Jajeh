clear all;clc;close all;
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

      %nx = input('enter nx (must be a power of 2):');
      %pdrop = input('enter pressure drop:');
      %vf0 = input('enter the desired volumn fraction:');
      %rsig = input('enter the std of the log-normal dist:');
      nx = 256;
      ny = nx;
      pdrop = 1.d0;
      vf0 = .05;
      rsig = .5;
      trials = 2;
      rho_EPS = 100;

n=25;
EPS_con = linspace(1,50,n);
eff_plot = zeros(1,n);
domain = [0 n];




% a =11; %infow of nutrients
% b = .1; %outflow of nutrients
% c = .8; %Nutrient uptake by algae 
% c_p = 1.3; %algal growth rate
% d = .5; %EPS growth rate due to algae 

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
%    *************************************

mass_EPS_total = zeros(1,n);
volume_total = zeros(1,n);
rho_EPS_eff=zeros(1,n);
rho_EPS_eff(1) = 1000;
    for i=1:n
        %First cross sectional area logN
        if i ==1
        sv = exp(rmu+rsig*randn(nx+1,ny+1));
        sh = exp(rmu+rsig*randn(nx+1,ny+1));
        
%     multigrid algorithm to solve the elliptic system:
       %Updating cross seectional area with EPS addition
        else
            % if i==2 
            %     IC_N = 25;
            %     IC_A = 25;
            %     IC_E = 0;
            %     IC = [IC_N IC_A IC_E];
            % else
            %     IC_N = N_sol(end) ;
            %     IC_A = A_sol(end) ;
            %     IC_E = E_sol(end);
            %     E_plot(i) = E_sol(end);
            %     IC = [IC_N IC_A IC_E];
            % end
            % [IV_sol,DVsol] = ode23(@(t, y) NAE(t, y, effcoe, a,b,c,c_p,d), domain, IC);
            %     N_sol = DVsol(:, 1);
            %     A_sol = DVsol(:, 2);
            %     E_sol = DVsol(:, 3);

                volume_total(i) = (sum(sv(:)) + sum(sh(:)))*h; %calc total volume
                mass_EPS_total(i) = volume_total(i) * (EPS_con(i)); %calc total mass of eps
                d_v = sv / sum(sv(:));  % Distribution factor for each pipe (proportional) vertical
                d_h = sh / sum(sh(:));  % Distribution factor for each pipe (proportional) horizontal
             
        % Compute EPS volume based on mass and an initial assumption for rho_EPS
             
             volume_EPS_v = (d_v .* mass_EPS_total(i)) / rho_EPS_eff(i);  
             volume_EPS_h = (d_h .* mass_EPS_total(i)) / rho_EPS_eff(i);  
        
        % Compute total EPS volume and update rho_EPS dynamically
                total_EPS_volume = sum(volume_EPS_v(:)) + sum(volume_EPS_h(:)); % m^3
                if total_EPS_volume>0
                     rho_EPS_eff(i) = mass_EPS_total(i) / total_EPS_volume; % kg/m^3
                else 
                    rho_EPS_eff(i) = rho_EPS_eff(i-1);
                end
                sv = max(0, sv - volume_EPS_v/h);  % Update vertical pipes (sv)
                sh = max(0, sh - volume_EPS_h/h);  % Update horizontal pipes (sh)

%                 figure;
%                 subplot(1,2,1);
%                 histogram(log(sv(:)));
%                 xlabel('ln(A) of Vertical Pipes');
%                 ylabel('Frequency');
%                 title('Histogram of ln(sv)');
% 
%                 subplot(1,2,2);
%                 histogram(log(sh(:)));
%                 xlabel('ln(A) of Horizontal Pipes');
%                 ylabel('Frequency');
%                 title('Histogram of ln(sh)');
% 
% sgtitle('Distribution of ln(A) for Pipe Cross-Sections After EPS Accumulation');
                
        end
     fin = zeros(nx+1, ny+1);
     [phi, final_error] = kikmul(fin,pdrop,sv,sh,nx,ny);

%     compute the effective coefficients based on the solution:

      effleft = sum((pdrop-phi(nx,:)).*sh(nx,:))/pdrop/h^2;
      effright = sum(phi(2,:).*sh(1,:))/pdrop/h^2;
      effcoe = 0.5*(effleft+effright);
      eff_plot(i) = effcoe;

%     Output:

      fprintf('effcoe = %8.6f\n',effcoe);
      fprintf('rmu and rsig: %12.6f %12.6f\n',rmu,rsig);

    end
    figure;
plot(EPS_con,eff_plot, 'LineWidth',2)
xlabel('EPS concentration')
ylabel('Permeability')
