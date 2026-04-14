clear; clc; close all;

% Parameters
params.L = 1;
params.a = 0.1;
params.R = 0.2;
params.D_B = 1;
params.D_E = 0.002;
params.U = 0.2;
params.lambda = 0.2;
params.k = 0.5;
params.gamma = 0.15;

params.C_in_node = 1;
params.C_node = 0.2;

% Grids
x = linspace(0, params.L, 120);
t = linspace(0, 10, 10);

m = 0;   % slab symmetry for 1D pipe axis

% Solve with pdepe
sol = pdepe(m, ...
    @(x,t,u,dudx) pipe_pde(x,t,u,dudx,params), ...
    @(x) pipe_ic(x,params), ...
    @(xl,ul,xr,ur,t) pipe_bc(xl,ul,xr,ur,t,params), ...
    x, t);

C_B = sol(:,:,1);
C_E = sol(:,:,2);

% Diagnostics
disp(['Minimum brine concentration: ', num2str(min(C_B(:)))]);
disp(['Minimum EPS concentration:   ', num2str(min(C_E(:)))]);

% Final profiles
figure;
plot(x, C_B(end,:), 'LineWidth', 2); hold on;
plot(x, C_E(end,:), 'LineWidth', 2);
xlabel('x');
ylabel('Concentration');
legend('C_B','C_E','Location','best');
title(['Final profiles at t = ', num2str(t(end))]);
grid on;

% Heatmaps
figure;
imagesc(x,t,C_B);
set(gca,'YDir','normal');
xlabel('x');
ylabel('t');
title('Brine concentration C_B(x,t)');
colorbar;

figure;
imagesc(x,t,C_E);
set(gca,'YDir','normal');
xlabel('x');
ylabel('t');
title('EPS concentration C_E(x,t)');
colorbar;

% Surface plots
figure;
surf(x,t,C_B,'EdgeColor','none');
xlabel('x'); ylabel('t'); zlabel('C_B');
title('Brine concentration surface');
view(45,30); colorbar;

figure;
surf(x,t,C_E,'EdgeColor','none');
xlabel('x'); ylabel('t'); zlabel('C_E');
title('EPS concentration surface');
view(45,30); colorbar;

% Selected time slices
figure;
idx = round(linspace(1,length(t),6));
for j = 1:length(idx)
    plot(x, C_B(idx(j),:), 'LineWidth', 1.5); hold on;
end
xlabel('x');
ylabel('C_B');
title('Brine concentration at selected times');
legend(compose('t = %.1f', t(idx)), 'Location', 'best');
grid on;

figure;
for j = 1:length(idx)
    plot(x, C_E(idx(j),:), 'LineWidth', 1.5); hold on;
end
xlabel('x');
ylabel('C_E');
title('EPS concentration at selected times');
legend(compose('t = %.1f', t(idx)), 'Location', 'best');
grid on;

% Animation
figure;
for n = 1:length(t)
    plot(x, C_B(n,:), 'LineWidth', 2); hold on;
    plot(x, C_E(n,:), 'LineWidth', 2);
    hold off;
    xlabel('x');
    ylabel('Concentration');
    title(sprintf('Pipe concentration at t = %.2f', t(n)));
    legend('Brine','EPS','Location','best');
    grid on;
    ylim([0, max([C_B(:); C_E(:); 1])]);
    drawnow;
end