%% run_lightcone_mass_shell_epsilon.m
% Quantum SOC complex mass–shell + light–cone families in w-space,
% and 2D deviation from the weak (classical) mass shell.
%
% Figure 1:
%   - Real classical light cone and mass shell w^2 = c^2  (reference).
%   - Complex mass–shell + cone families (weak closure, σ>0).
%   - SOC electron/positron drifts w(t) overlaid.
%
% Figure 2:
%   - Deviation (w^T eta w - c^2)/c^2 vs time, for e^- and e^+,
%     averaged over trajectories (weak mass-shell residual).

clear;
rng(7);   % reproducible bundle

%% Simulation parameters
n      = 0;
s_spin = -1;
B      = 0.1;                 % Tesla (physical field for the simulation)
ky     = 5e11;                % transverse momentum (guiding center)
pz     = 0;                   % choose pz=0 -> wz=0 Landau slice
Ntraj  = 50;                  % number of trajectories in bundle
Nt     = 10000;               % time steps
save_figs = false;            % set true to export PDFs
lambda_gauge = 'weak_mass_sell';

% ---- constants & derived scalars ----
C   = phys_constants();
wc  = abs(C.e_q)*B/C.m;
Tc  = 2*pi/wc;
T   = 5*Tc;                   % simulate several cyclotron periods
dt  = T / Nt;

fprintf('Running %d trajectories: B=%.2g T, Nt=%d, T=%.3e s (Tc=%.3e s)\n', ...
        Ntraj, B, Nt, T, Tc);

% Electron branch: ε=+1
sim_elec = simulate_bundle('electron', +1, n, s_spin,  ky,  pz, B, Ntraj, Nt, dt, []);
% Positron branch: ε=-1, with spin and ky flipped (charge conjugation)
sim_pos  = simulate_bundle('positron', -1, n, -s_spin, -ky, pz, B, Ntraj, Nt, dt, []);

% Use drifts (w_store) from the simulator; real parts for plotting scale
w_e = real(sim_elec.w_store);
w_p = real(sim_pos.w_store);
w_e0_all = reshape(w_e(1,:,:), Nt*Ntraj, 1);
w_ex_all = reshape(w_e(2,:,:), Nt*Ntraj, 1);
w_ey_all = reshape(w_e(3,:,:), Nt*Ntraj, 1);
w_p0_all = reshape(w_p(1,:,:), Nt*Ntraj, 1);
w_px_all = reshape(w_p(2,:,:), Nt*Ntraj, 1);
w_py_all = reshape(w_p(3,:,:), Nt*Ntraj, 1);

all_sp_rad = [vecnorm([w_ex_all w_ey_all], 2, 2); vecnorm([w_px_all w_py_all], 2, 2)];
w_max = max(all_sp_rad) / C.c;
wr_limit = max([1.2*w_max, 1e-2]);  % dimensionless span in w/c

%% ============================================================
%  FIGURE 1: Classical shell + complex families + SOC drifts
% =============================================================

% --- XY-only spatial drift (consistent with plotted wx, wy)
wsp_e_xy = sim_elec.w_store(2:3,:,:);
wsp_p_xy = sim_pos.w_store(2:3,:,:);

u_e = real(wsp_e_xy);  v_e = imag(wsp_e_xy);
u_p = real(wsp_p_xy);  v_p = imag(wsp_p_xy);

r_e = squeeze(vecnorm(u_e,2,1)) / C.c;          % Nt x Ntraj, r=|Re w_xy|/c
r_p = squeeze(vecnorm(u_p,2,1)) / C.c;

sigma_e = squeeze(vecnorm(v_e,2,1)) / C.c;      % Nt x Ntraj, σ=|Im w_xy|/c
sigma_p = squeeze(vecnorm(v_p,2,1)) / C.c;

cosang_e = squeeze( sum(u_e.*v_e,1) ./ (vecnorm(u_e,2,1).*vecnorm(v_e,2,1) + 1e-30) );
cosang_p = squeeze( sum(u_p.*v_p,1) ./ (vecnorm(u_p,2,1).*vecnorm(v_p,2,1) + 1e-30) );

sigma_all  = [sigma_e(:); sigma_p(:)];
cosang_all = [cosang_e(:); cosang_p(:)];

sigma_all  = sigma_all(isfinite(sigma_all));
cosang_all = cosang_all(isfinite(cosang_all));

% --- Choose σ-levels from data; always include 0 for the real (classical) shell
sig_levels = unique([0, quantile(sigma_all, [0.1 0.25 0.5 0.75 0.9])]);
sig_levels = sig_levels(sig_levels < 0.95);

% --- For each σ-level, estimate typical Re/Im angle from data
alpha_levels = zeros(size(sig_levels));
delta_sig = 0.02;   % σ-bin half-width

for k = 1:numel(sig_levels)
    s0 = sig_levels(k);
    mask = abs(sigma_all - s0) < delta_sig;

    if any(mask)
        c0 = median(cosang_all(mask), 'omitnan');
        c0 = min(max(c0,-1),1);
        alpha_levels(k) = acos(c0);
    else
        alpha_levels(k) = pi/2;  % fallback if no data
    end
end

% --- Build mesh in PHYSICAL units then normalize by c in plotting
rho_max = 1.2 * wr_limit * C.c;
rho = linspace(0, rho_max, 140);
th  = linspace(0, 2*pi, 200);
[RHO, TH] = meshgrid(rho, th);

ux = RHO .* cos(TH);
uy = RHO .* sin(TH);

fC = figure; hold on; box on; grid on;
title('Classical and complex mass–shell geometry with SOC drifts');
xlabel('Re w_x / c'); ylabel('Re w_y / c'); zlabel('Re w_0 / c');

% --- Classical real light cone and real mass shell (reference)
wr = linspace(0, wr_limit, 150);
th2 = linspace(0, 2*pi, 220);
[WRc, THc] = meshgrid(wr, th2);

WXc = WRc .* cos(THc);
WYc = WRc .* sin(THc);

cone_real  = WRc;                 % (w0/c) = |w_sp|/c
shell_real = sqrt(1 + WRc.^2);    % (w0/c) = sqrt(1 + |w_sp|^2/c^2)

% Real light cone (upper + lower)
surf(WXc, WYc,  cone_real, ...
    'FaceAlpha',0.25, ...
    'EdgeAlpha',0.18, ...
    'EdgeColor',[0.2 0.2 0.2], ...
    'FaceColor',[0.3 0.3 0.3], ...
    'DisplayName','light cone');
surf(WXc, WYc, -cone_real, ...
    'FaceAlpha',0.25, ...
    'EdgeAlpha',0.18, ...
    'EdgeColor',[0.2 0.2 0.2], ...
    'FaceColor',[0.3 0.3 0.3], ...
    'HandleVisibility','off');

% Real classical mass shell (upper + lower)
h_class_shell = surf(WXc, WYc,  shell_real, ...
    'FaceAlpha',0.20, 'EdgeColor','none', ...
    'FaceColor',[0 0.55 0.95], ...
    'DisplayName','w^2 - c^2 = 0 (+)');
surf(WXc, WYc, -shell_real, ...
    'FaceAlpha',0.20, 'EdgeColor','none', ...
    'FaceColor',[0.95 0.45 0.2], ...
    'DisplayName','w^2 - c^2 = 0 (-)');

% --- Draw complex mass–shell & cone families for σ > 0 (weak closure)
for k = 1:numel(sig_levels)
    sig   = sig_levels(k);
    alpha = alpha_levels(k);

    if sig <= 0
        % σ=0 is exactly the real shell already drawn above
        continue;
    end

    % Imaginary spatial part with empirical angle alpha
    vx = sig*C.c .* cos(TH + alpha);
    vy = sig*C.c .* sin(TH + alpha);

    wx = ux + 1i*vx;
    wy = uy + 1i*vy;
    wz = 0;

    % ---- Complex MASS SHELL via weak closure: w^2 = c^2 (complex)
    w0s = sqrt(C.c^2 + wx.^2 + wy.^2 + wz.^2);

    surf(real(wx)/C.c, real(wy)/C.c,  real(w0s)/C.c, ...
        'FaceAlpha',0.06, 'EdgeColor','none', ...
        'FaceColor',[0 0.55 0.95], 'HandleVisibility','off');
    surf(real(wx)/C.c, real(wy)/C.c, -real(w0s)/C.c, ...
        'FaceAlpha',0.06, 'EdgeColor','none', ...
        'FaceColor',[0.95 0.45 0.2], 'HandleVisibility','off');

    % ---- Complex LIGHT CONE via weak closure: w^2 = 0 (complex)
    w0c = sqrt(wx.^2 + wy.^2 + wz.^2);

    % Mask tiny Re(w0) so cone families are visually recognizable
    w0c_re = real(w0c);
    w0c_re(abs(w0c_re) < 1e-6) = NaN;

    surf(real(wx)/C.c, real(wy)/C.c,  w0c_re/C.c, ...
        'FaceAlpha',0.04, 'EdgeAlpha',0.15, ...
        'EdgeColor',[0.2 0.2 0.2], 'FaceColor',[0.8 0.8 0.8], ...
        'HandleVisibility','off');
    surf(real(wx)/C.c, real(wy)/C.c, -w0c_re/C.c, ...
        'FaceAlpha',0.04, 'EdgeAlpha',0.15, ...
        'EdgeColor',[0.2 0.2 0.2], 'FaceColor',[0.8 0.8 0.8], ...
        'HandleVisibility','off');
end

% --- Overlay a few simulated SOC drifts (electron blue, positron red)
n_show = 5;
for j = 1:n_show
    we = sim_elec.w_store(:,:,j);
    plot3(real(we(2,:))/C.c, real(we(3,:))/C.c, real(we(1,:))/C.c, ...
        'LineWidth',1.3, 'Color',[0 0.45 0.74]);

    wp = sim_pos.w_store(:,:,j);
    plot3(real(wp(2,:))/C.c, real(wp(3,:))/C.c, real(wp(1,:))/C.c, ...
        'LineWidth',1.3, 'Color',[0.85 0.33 0.1]);
end

view(28,24); axis vis3d;
xlim([-1 1]*wr_limit);
ylim([-1 1]*wr_limit);
zlim([-1.2*max(shell_real(:)), 1.2*max(shell_real(:))]);

legend('Location','northwest');

if save_figs
    exportgraphics(fC, fullfile('figures','fig_complex_shell_and_drifts.pdf'), ...
        'ContentType','vector', 'BackgroundColor','none', 'Resolution',600);
end

%% ============================================================
%  FIGURE 2: Deviation from weak mass shell vs time
% ============================================================

t_vec = (0:Nt-1).' * dt;

[mean_rel_e, max_abs_e] = mass_shell_error(sim_elec.w_store, C);
[mean_rel_p, max_abs_p] = mass_shell_error(sim_pos.w_store,  C);
max_abs_rel = max(max_abs_e, max_abs_p);

f_res = figure; hold on; box on; grid on;
plot(t_vec, mean_rel_e, 'LineWidth', 1.2, 'Color',[0 0.45 0.74], ...
    'DisplayName','electron (mean over trajectories)');
plot(t_vec, mean_rel_p, 'LineWidth', 1.2, 'Color',[0.85 0.33 0.1], ...
    'DisplayName','positron (mean over trajectories)');
yline(0, 'k--', 'LineWidth', 1.0, 'DisplayName','on-shell');

xlabel('t [s]');
ylabel('(w^{T}\eta w - c^{2})/c^{2}');
title('Deviation from weak mass shell w^{T}\eta w = c^{2}');
ylim(1.2 * [-1 1] * max(1e-18, max_abs_rel));
legend('Location','best');

if save_figs
    exportgraphics(f_res, fullfile('figures','fig_mass_shell_residual.pdf'), ...
        'ContentType','vector', 'BackgroundColor','none', 'Resolution',600);
end

%% Helper: bundle simulator
function out = simulate_bundle(mode_label, epsilon, n, s_spin, ky, pz, B, Ntraj, Nt, dt, pdev)
    % epsilon: SOC charge sign (+1 electron, -1 positron)

    z_store = complex(zeros(4, Nt+1, Ntraj));
    w_store = complex(zeros(4, Nt,   Ntraj));
    z_init = complex(zeros(4,1));

    for j = 1:Ntraj
        rng(137 + j*1000 + double(mode_label(1)), 'twister'); % independent noise per trajectory
        [~,~,~,~,~, w_j, z_j] = dirac_landau_action(n, s_spin, B, dt, Nt, z_init, ky, pz, pdev, epsilon, lambda_gauge);
        w_store(:,:,j) = w_j;
        z_store(:,:,j) = z_j;
    end

    out.z_store = z_store;
    out.w_store = w_store;
end

%% Helper: weak mass-shell residual
function [mean_rel_err, max_abs_rel] = mass_shell_error(w_store, C)
    % w_store: 4 x Nt x Ntraj (complex)
    % SOC weak shell: w^T eta w = c^2
    eta      = C.eta;
    Nt       = size(w_store, 2);
    Ntraj    = size(w_store, 3);
    rel_err  = zeros(Nt, Ntraj);

    for j = 1:Ntraj
        w   = w_store(:,:,j);                 % 4 x Nt
        s_v = sum(w .* (eta*w), 1);          % 1 x Nt, bilinear (no conj)
        w2_err = s_v - C.c^2;
        rel_err(:,j) = w2_err.'/C.c^2;       % keep complex if you like
    end

    % For plotting: real part of the mean (imag part is usually tiny too)
    mean_rel_err = mean(real(rel_err), 2);   % Nt x 1
    max_abs_rel  = max(abs(mean_rel_err));
end