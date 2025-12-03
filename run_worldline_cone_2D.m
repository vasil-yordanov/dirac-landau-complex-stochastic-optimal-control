% run_worldline_cone_2D.m
% Plot the two 2D light-cone diagnostics figures.

clear; close all; clc;

n       = 1;
s_spin  = -1;
B       = 0.1;
ky      = 1e12;
pz      = 0;
epsilon = 1;
Nt      = 1200;

%% ----------------- Timing -----------------
C  = phys_constants();
wc = abs(C.e_q)*B/C.m;
Tc = 2*pi/wc;
T  = 1*Tc;
dt = T/Nt;

t   = (0:Nt)*dt;
tau_norm = t / Tc;

%% ----------------- Simulate -----------------
z_init = complex(zeros(4,1));
rng(137+1000,'twister');

[~,~,~,~,~, w_store, z_store] = dirac_landau_action( ...
    n, s_spin, B, dt, Nt, z_init, ky, pz, [], epsilon);

zR = real(z_store);    % 4 x (Nt+1)
wR = real(w_store);    % 4 x Nt

% Shift to origin
z0 = zR(1,:) - zR(1,1);
x  = zR(2,:) - zR(2,1);
y  = zR(3,:) - zR(3,1);

% Drift-only reconstruction
z_drift = zeros(4, Nt+1);
for k = 1:Nt
    z_drift(:,k+1) = z_drift(:,k) + wR(:,k)*dt;
end
z0_d = z_drift(1,:);
x_d  = z_drift(2,:);
y_d  = z_drift(3,:);

% Causal diagnostics (xy projection)
s2_perp      = (z0.^2)   - (x.^2   + y.^2);
s2_perp_drift= (z0_d.^2) - (x_d.^2 + y_d.^2);
ds2_perp     = s2_perp - s2_perp_drift;

% Normalization
L0 = C.c * T;
x_hat  = x  / L0;  y_hat  = y  / L0;  z0_hat = z0 / L0;
x_d_hat= x_d/ L0;  y_d_hat= y_d/ L0;  z0_d_hat = z0_d / L0;

ds2_hat = ds2_perp / (L0^2);

%% Styling
FS_AX   = 16;   % tick labels for 2D figures
FS_LAB  = 20;   % axis labels
FS_LEG  = 15;
FS_TIT  = 20;
LW      = 1.6;

%% 2D projection r_perp vs z0
rh_all = sqrt(x_hat.^2 + y_hat.^2);
mask = z0_hat >= 0;
z0h = z0_hat(mask);
rh  = rh_all(mask);

zmax = max([max(z0h), max(rh), 1e-9]);

f_r = figure('Color','w', 'Units','inches', 'Position',[1 1 6.2 5.2]);
ax_r = axes('Parent',f_r); hold(ax_r,'on'); grid(ax_r,'on'); box(ax_r,'on');

% Cone boundary (not in legend)
plot(ax_r, [0 zmax], [0 zmax], 'k:', 'LineWidth',1.3, 'DisplayName','light-cone boundary');

% Worldline only
plot(ax_r, z0h, rh, 'LineWidth',LW, 'DisplayName','worldline');

xlabel(ax_r,'$c t/(c T_C)$','Interpreter','latex','FontSize',FS_LAB);
ylabel(ax_r,'$r_\perp/(c T_C)$','Interpreter','latex','FontSize',FS_LAB);
title(ax_r,'Worldline projection vs. light-cone bound', ...
    'Interpreter','latex','FontSize',FS_TIT);

set(ax_r,'FontSize',FS_AX);
axis(ax_r,'equal');
xlim(ax_r,[0 zmax]); ylim(ax_r,[0 zmax]);

leg_r = legend(ax_r,'Location','northwest','Interpreter','latex');
leg_r.FontSize = FS_LEG;

tightAxesLocal(ax_r,0);

exportgraphics(f_r, fullfile('figures','figure_cone_2D_rperp_vs_z0.pdf'), ...
    'ContentType','image', 'BackgroundColor','none', 'Resolution',600);

%% Causal residual about drift
f_res = figure('Color','w', 'Units','inches', 'Position',[1 1 6.2 5.2]);
ax_res = axes('Parent',f_res); hold(ax_res,'on'); grid(ax_res,'on'); box(ax_res,'on');

plot(ax_res, tau_norm, ds2_hat, 'LineWidth',LW);

yline(ax_res, 0, '-', 'Color',[0.2 0.2 0.2], 'LineWidth',1.0, 'HandleVisibility','off');

xlabel(ax_res,'$\tau/T_C$','Interpreter','latex','FontSize',FS_LAB);
ylabel(ax_res,'$\delta s_\perp^2/(c T_C)^2$','Interpreter','latex','FontSize',FS_LAB);
title(ax_res,'Causal residual about drift', ...
    'Interpreter','latex','FontSize',FS_TIT);

set(ax_res,'FontSize',FS_AX);

tightAxesLocal(ax_res,0);

exportgraphics(f_res, fullfile('figures','figure_causal_residual_about_drift.pdf'), ...
    'ContentType','image', 'BackgroundColor','none', 'Resolution',600);

disp('Generated figure_cone_2D_rperp_vs_z0.pdf and figure_causal_residual_about_drift.pdf');

function tightAxesLocal(ax, pad)
    % Reduce white space around axes by padding a tight layout (normalized units).
    ax.Units = 'normalized';
    outerpos = ax.OuterPosition;
    ti = ax.TightInset;
    ax.Position = [outerpos(1)+ti(1)+pad, outerpos(2)+ti(2)+pad, ...
                   outerpos(3)-ti(1)-ti(3)-2*pad, outerpos(4)-ti(2)-ti(4)-2*pad];
end
