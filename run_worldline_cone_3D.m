% run_worldline_cone_3D.m
% Plot only the 3D light-cone/worldline figure.

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
tau = t / T;

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
FS_AX   = 15;   % tick labels for 3D figures
FS_LAB  = 19;   % axis labels
FS_LEG  = 15;
LW      = 1.6;

%% Cone surface mesh
z_top = max(z0_hat); 
z_top = max(z_top, 1e-6);
theta = linspace(0, 2*pi, 180);
zg    = linspace(0, z_top, 140);
[TH, ZG] = meshgrid(theta, zg);
R  = ZG;
Xc = R .* cos(TH);
Yc = R .* sin(TH);
Zc = ZG;

%% Figure
f = figure('Color','w', 'Units','inches', 'Position',[1 1 6.8 5.2]);
ax = axes('Parent',f); hold(ax,'on'); grid(ax,'on'); box(ax,'on');

xlabel(ax,'$x/(cT_c)$','Interpreter','latex','FontSize',FS_LAB);
ylabel(ax,'$y/(cT_c)$','Interpreter','latex','FontSize',FS_LAB);
zlabel(ax,'$c t/(cT_c)$','Interpreter','latex','FontSize',FS_LAB);

surf(ax, Xc, Yc, Zc, ...
    'FaceAlpha',0.12,'EdgeAlpha',0.06, ...
    'EdgeColor',[0.4 0.4 0.4],'FaceColor',[0.4 0.4 0.4], ...
    'DisplayName','light cone');

plot3(ax, x_hat, y_hat, z0_hat, ...
    'LineWidth',LW, 'DisplayName','worldline');

set(ax,'FontSize',FS_AX);

daspect(ax,[1 1 1]);
axis(ax,'vis3d');
view(ax,69,20);

leg = legend(ax,'Location','east','Interpreter','latex');
leg.FontSize = FS_LEG;

pbaspect(ax,[1 1 1]); 
zlim(ax,[0 1]);
xlim(ax,[-1 1]);
ylim(ax,[-1 1]);

exportgraphics(f, fullfile('figures','figure_worldline_cone_3D.pdf'), ...
    'ContentType','image', 'BackgroundColor','none', 'Resolution',600);

disp('Generated figure_worldline_cone_3D.pdf');
