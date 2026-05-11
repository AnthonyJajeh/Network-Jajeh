clear all; clc; close all;
params.L = 1; %pipe length
params.a = .1; %brine core radius
params.R = .2; %outer pipe radius
params.D_B = 1; %brine diffusion coefficient
params.D_E = .002; %EPS diffusion coefficient
params.U = .2;
params.lambda = .2; %EPS loss rate
params.k = .5;

params.C_in_node = 1;
params.C_node = .2;

params.gamma=.15;

params.Nx = 120;
params.x  = linspace(0,params.L,params.Nx)';
params.dx = params.x(2)-params.x(1);

% Initial conditions
C_B_0 = zeros(params.Nx,1);
C_E_0 = zeros(params.Nx,1);
y_0  = [C_B_0; C_E_0];

t_span = [0 100]; %time interval

opts = odeset('RelTol',1e-6,'AbsTol',1e-8);
[t,y] = ode15s(@(t,y) pipe_PDE(t,y,params),t_span,y_0,opts);

C_B = y(:,1:params.Nx);
C_E = y(:,params.Nx+1:end);

% Plot
figure;
plot(params.x,C_B(end,:),'LineWidth',2); hold on;
plot(params.x,C_E(end,:),'LineWidth',2);
xlabel('x');
ylabel('Concentration');
legend('\bar C_B','\bar C_E','Location','best');
title('Final concentration profiles');
grid on;

figure;
imagesc(params.x,t,C_B);
set(gca,'YDir','normal');
colorbar;
xlabel('x'); ylabel('t');
title('Brine concentration');

figure;
imagesc(params.x,t,C_E);
set(gca,'YDir','normal');
colorbar;
xlabel('x'); ylabel('t');
title('EPS concentration');

% 3. Animation
figure;
for n = 1:length(t)
    plot(params.x, C_B(n,:), 'LineWidth', 2); hold on;
    plot(params.x, C_E(n,:), 'LineWidth', 2);
    hold off;
    xlabel('x');
    ylabel('Concentration');
    title(sprintf('Concentration in pipe at t = %.2f', t(n)));
    legend('Brine','EPS','Location','best');
    grid on;
    ylim([0, max([C_B(:); C_E(:); 1])]);
    drawnow;
end