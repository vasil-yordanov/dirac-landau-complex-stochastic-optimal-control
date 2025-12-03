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
n      = 1;
s_spin = -1;
B      = 0.1;                 % Tesla (physical field for the simulation)
ky     = 1e12;                % transverse momentum (guiding center)
pz     = 0;                   % choose pz=0 -> wz=0 Landau slice
Ntraj  = 1;                  % number of trajectories in bundle
Nt     = 1000;               % time steps

% ---- constants & derived scalars ----
C   = phys_constants();
wc  = abs(C.e_q)*B/C.m;
Tc  = 2*pi/wc;
T   = 1*Tc;                   % simulate several cyclotron periods
dt  = T / Nt;
epsilon = 1; % 1- electron ; 2 -positron

fprintf('Running %d trajectories: B=%.2g T, Nt=%d, T=%.3e s (Tc=%.3e s)\n', ...
        Ntraj, B, Nt, T, Tc);


z_init = complex(zeros(4,1));
rng(137 + 1000, 'twister');
[~,~,~,~,~, w_e, z_z] = dirac_landau_action(n, s_spin, B, dt, Nt, z_init, ky, pz, [], epsilon);


% Use drifts (w_store) from the simulator; real parts for plotting scale
w_e0_all = reshape(real(w_e(1,:)), Nt*Ntraj, 1);
w_ex_all = reshape(real(w_e(2,:)), Nt*Ntraj, 1);
w_ey_all = reshape(real(w_e(3,:)), Nt*Ntraj, 1);


all_sp_rad = [vecnorm([w_ex_all w_ey_all], 2, 2)];
w_max = max(all_sp_rad) / C.c;
wr_limit = max([1.2*w_max, 1e-2]);  % dimensionless span in w/c

%% ----------------- Styling knobs -----------------
FS_AX_3D  = 12;   % tick labels
FS_LAB_3D = 16;   % axis labels
FS_LEG_3D = 12;
LW        = 1.6;

%% ============================================================
%  FIGURE 1: Classical shell + complex families + SOC drifts
% =============================================================

% --- XY-only spatial drift (consistent with plotted wx, wy)
wsp_e_xy = w_e(2:3);

u_e = real(wsp_e_xy);  v_e = imag(wsp_e_xy);

r_e = squeeze(vecnorm(u_e,2,1)) / C.c;          % Nt x Ntraj, r=|Re w_xy|/c

sigma_e = squeeze(vecnorm(v_e,2,1)) / C.c;      % Nt x Ntraj, σ=|Im w_xy|/c

cosang_e = squeeze( sum(u_e.*v_e,1) ./ (vecnorm(u_e,2,1).*vecnorm(v_e,2,1) + 1e-30) );

sigma_all  = [sigma_e(:)];
cosang_all = [cosang_e(:)];

sigma_all  = sigma_all(isfinite(sigma_all));
cosang_all = cosang_all(isfinite(cosang_all));

% Imaginary spatial drift (xy-plane) along the path
wsp_e_xy = w_e(2:3,:);      % 2xNt: wx, wy
v_e      = imag(wsp_e_xy);  % 2xNt

% --- Choose σ-levels from data;
%sig_levels = unique([quantile(sigma_all, [0.1 0.25 0.5 0.75 0.9])]);
%sig_levels = sig_levels(sig_levels < 0.95 & sig_levels >0);

rho_tau  = vecnorm(v_e, 2, 1) / C.c;   % 1xNt, this is ρ(τ)
sig_levels  = sqrt(mean(rho_tau.^2));   % effective ρ

% --- For each  σ-level, estimate typical Re/Im angle from data
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
rho_max = 4 * wr_limit * C.c;
rho = linspace(0, rho_max, 140);
th  = linspace(0, 2*pi, 200);
[RHO, TH] = meshgrid(rho, th);

ux = RHO .* cos(TH);
uy = RHO .* sin(TH);

fC = figure('Color','w', 'Units','inches', 'Position',[1 1 7.0 6.0]);
axC = axes('Parent',fC); hold(axC,'on'); box(axC,'on'); grid(axC,'on');
xlabel(axC,'$\Re \pi_x / mc$', 'Interpreter', 'latex','FontSize',FS_LAB_3D); 
ylabel(axC,'$\Re \pi_y / mc$', 'Interpreter', 'latex','FontSize',FS_LAB_3D);
zlabel(axC,'$\Re \pi_0 / mc$', 'Interpreter', 'latex','FontSize',FS_LAB_3D);

% --- Classical real light cone and real mass shell (reference)
wr = linspace(0, 3.7*wr_limit, 150);
th2 = linspace(0, 2*pi, 220);
[WRc, THc] = meshgrid(wr, th2);

WXc = WRc .* cos(THc);
WYc = WRc .* sin(THc);

cone_real  = WRc;                 % (w0/c) = |w_sp|/c
shell_real = sqrt(1 + WRc.^2);    % (w0/c) = sqrt(1 + |w_sp|^2/c^2)

% Real light cone (upper + lower)
surf(WXc, WYc,  cone_real, ...
    'FaceAlpha',0.2, ...
    'EdgeAlpha',0.15, ...
    'EdgeColor',[0.8 0.8 0.8], ...
    'FaceColor',[0.8 0.8 0.8], ...
    'DisplayName','light cone');
surf(WXc, WYc, -cone_real, ...
    'FaceAlpha',0.2, ...
    'EdgeAlpha',0.15, ...
    'EdgeColor',[0.8 0.8 0.8], ...
    'FaceColor',[0.8 0.8 0.8], ...
    'HandleVisibility','off');

% Real classical mass shell (upper + lower)
h_class_shell = surf(WXc, WYc,  shell_real, ...
    'FaceAlpha',0.20, 'EdgeColor','none', ...
    'FaceColor',[0 0.55 0.95], ...
    'DisplayName','$S_0$ - classical mass shell reference');
surf(WXc, WYc, -shell_real, ...
    'FaceAlpha',0.20, 'EdgeColor','none', ...
    'FaceColor',[0 0.55 0.95], ... 
    'HandleVisibility','off'); 
  

% --- Draw complex mass–shell for varrho > 0 (weak closure)
for k = 1:numel(sig_levels)
    sig   = sig_levels(k);
    alpha = alpha_levels(k);

    if sig <= 0
        % varrho=0 is exactly the real shell already drawn above
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
        'FaceAlpha',0.1, ...
        'EdgeColor','none', ...
        'FaceColor',[0 0.55 0.95], ...
        'DisplayName','$\widetilde S_\varrho$ - representative weak mass shell');
    surf(real(wx)/C.c, real(wy)/C.c, -real(w0s)/C.c, ...
        'FaceAlpha',0.1, ...
        'EdgeColor','none', ...
        'FaceColor',[0 0.55 0.95], ... %'FaceColor',[0.95 0.45 0.2], ...
        'HandleVisibility','off');
end

% --- Overlay a few simulated SOC drifts (electron blue, positron red)
plot3(real(w_e(2,:))/C.c, real(w_e(3,:))/C.c, real(w_e(1,:))/C.c, ...
    'LineWidth',LW, ...
    'Color',[0 0.45 0.74], ...
     'DisplayName','$\Re \pi(\tau) / mc$ - projected momentum trajectory');

set(axC,'FontSize',FS_AX_3D);

view(axC,60, 15);

% Camera position:
campos(axC,[18, -32, 10]);

% Camera view angle:
camva(axC,10);

%view(28,24); 
axis(axC,'vis3d');
xlim(axC,[-5 5]*wr_limit);
ylim(axC,[-5 5]*wr_limit);
zlim(axC,[-1.5*max(shell_real(:)), 1.5*max(shell_real(:))]);

legC = legend(axC,'Location','northeast', 'Interpreter','latex');
legC.FontSize = FS_LEG_3D;

exportgraphics(fC, fullfile('figures','figure_complex_shell_and_drifts.pdf'), ...
    'ContentType','image', 'BackgroundColor','white', 'Resolution',600);
