   
    % -------------------------------------------------------------------------
    function [A_sv, A_sh, vfrac, amean] = gen_pipes(rmu, rsig, nx, ny, h)
        xmin = rmu - 4 * rsig;
        xmax = rmu + 4 * rsig;
    
        X_sv = max(min(rmu + rsig * randn(nx + 1, ny + 1), xmax), xmin); %Generate normal RV
        X_sh = max(min(rmu + rsig * randn(nx + 1, ny + 1), xmax), xmin); %Generate normal RV
    
        A_sv = exp(X_sv); %Generate lognormal RV for area
        A_sh = exp(X_sh); %Generate lognormal RV for area
        asum = sum(A_sv(:)) + sum(A_sh(:));
        amean = 0.5 * asum / ((nx + 1) * (ny + 1)); %mean area
        vfrac = asum / ((nx + 1) * (ny + 1)) / h^2; %volume fraction
    end