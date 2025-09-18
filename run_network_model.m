function [eff_perm, sv_avg, sh_avg, A_sv_avg, A_sh_avg] = run_network_model( ...
        dEPS_conc_kgm3, nx, ny, rho_EPS, trials, pdrop, A_sv, A_sh, h)

    % dEPS_conc_kgm3 : EPS concentration to deposit this call [kg/m^3]
    % nx,ny          : grid size
    % rho_EPS        : EPS density [kg/m^3]
    % trials         : # pressure-solve trials to average (geometry fixed)
    % pdrop          : pressure drop
    % A_sv, A_sh     : vertical/horizontal edge AREAS [m^2], same shapes kikmul expects
    % h              : edge length [m]

    % ----------------- knobs (for visibility & calibration) -----------------
    alpha_weight  = 6.0;   % deposit weight exponent (>=2 concentrates on high-G edges)
    S_deposit     = 1.0;   % global scale on deposition; try 10–1e4 to see sensitivity
    A_floor       = 0.0;   % hard floor on area (m^2); set >0 to avoid zeroing
    % -----------------------------------------------------------------------

    % ---- 1) Estimate pore volume from current geometry ----
    % Each edge volume ~ A_edge * h. Average v/h to avoid double counting.
    V_v = sum(A_sv(:)) * h;
    V_h = sum(A_sh(:)) * h;
    V_pore_est = 0.5 * (V_v + V_h);

    % ---- 2) Total EPS mass/volume to deposit this call ----
    mass_EPS_total = max(0, dEPS_conc_kgm3) * V_pore_est * S_deposit;  % [kg]
    vol_EPS_total  = mass_EPS_total / rho_EPS;                          % [m^3]

    if vol_EPS_total <= 0 || ~isfinite(vol_EPS_total)
        % nothing to deposit; just compute permeability with current areas
        sv = (A_sv / pi).^2;
        sh = (A_sh / pi).^2;
        [eff_perm, sv_avg, sh_avg] = local_effective_perm(nx, ny, pdrop, sv, sh, h, trials);
        A_sv_avg = A_sv;
        A_sh_avg = A_sh;
        return
    end

    % ---- 3) Build deposition weights (conductance^alpha) ------------------
    % Conductance proxy: r^4 ~ (A/pi)^2. Use power alpha_weight to concentrate.
    Gv = (A_sv / pi).^2;
    Gh = (A_sh / pi).^2;

    % small epsilon so zero-area edges don't kill normalization
    epsW = 1e-30;
    Wv = max(Gv, epsW) .^ alpha_weight;
    Wh = max(Gh, epsW) .^ alpha_weight;

    Wsum = sum(Wv(:)) + sum(Wh(:));
    if Wsum <= 0 || ~isfinite(Wsum)
        % fallback to area weights
        Wv = max(A_sv, epsW);
        Wh = max(A_sh, epsW);
        Wsum = sum(Wv(:)) + sum(Wh(:));
    end
    Wv = Wv / Wsum;
    Wh = Wh / Wsum;

    % ---- 4) Assign EPS volume to edges and update areas -------------------
    vol_v = vol_EPS_total * Wv;     % [m^3] per vertical edge
    vol_h = vol_EPS_total * Wh;     % [m^3] per horizontal edge

    % ΔA = volume / h
    A_sv_upd = max(A_floor, A_sv - vol_v / h);
    A_sh_upd = max(A_floor, A_sh - vol_h / h);

    % ---- 5) Conductances (Hagen–Poiseuille): r^4 ~ (A/pi)^2 --------------
    sv = (A_sv_upd / pi).^2;
    sh = (A_sh_upd / pi).^2;

    % ---- 6) Multigrid pressure solves & k estimation ----------------------
    [eff_perm, sv_avg, sh_avg] = local_effective_perm(nx, ny, pdrop, sv, sh, h, trials);

    % Return averaged updated geometry (geometry doesn’t vary across trials here)
    A_sv_avg = A_sv_upd;
    A_sh_avg = A_sh_upd;
end

% ========================= helper: permeability ============================
function [eff_perm, sv_avg, sh_avg] = local_effective_perm(nx, ny, pdrop, sv, sh, h, trials)
    eff_trial = zeros(1, trials);
    sv_sum = zeros(size(sv));
    sh_sum = zeros(size(sh));
    for t = 1:trials
        fin = zeros(nx + 1, ny + 1);
        [phi_grid, ~] = kikmul(fin, pdrop, sv, sh, nx, ny);

        % Effective permeability from left/right flux (symmetric average)
        effleft  = sum((pdrop - phi_grid(nx,:)) .* sh(nx,:)) / h^2;
        effright = sum(phi_grid(2,:) .* sh(1,:)) / h^2;
        eff_trial(t) = 0.5 * (effleft + effright) * (pi / 8);

        sv_sum = sv_sum + sv;
        sh_sum = sh_sum + sh;
    end
    eff_perm = mean(eff_trial);
    sv_avg   = sv_sum / trials;
    sh_avg   = sh_sum / trials;
end
