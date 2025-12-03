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
ky     = 5e7;                % transverse momentum (guiding center)
pz     = 0;                   % choose pz=0 -> wz=0 Landau slice
Ntraj  = 1;                  % number of trajectories in bundle
Nt     = 100000;               % time steps

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

w_e0_all = reshape(real(w_e(1,:)), Nt*Ntraj, 1);
w_ex_all = reshape(real(w_e(2,:)), Nt*Ntraj, 1);
w_ey_all = reshape(real(w_e(3,:)), Nt*Ntraj, 1);


all_sp_rad = [vecnorm([w_ex_all w_ey_all], 2, 2)];
w_max = max(all_sp_rad) / C.c;
wr_limit = max([1.2*w_max, 1e-2]);  % dimensionless span in w/c

%% ----------------- Styling knobs -----------------
FS_AX  = 15;   % tick labels
FS_LAB = 19;   % axis labels
FS_LEG = 14;
FS_TIT = 19;
LW     = 1.6;

%% ============================================================
%  FIGURE: Deviation from weak mass shell vs time
% ============================================================

 t_vec = (0:Nt-1).' * dt/Tc;
 
 [rel_e] = mass_shell_deviation(w_e, C);
 
 f_res = figure;
 ax_res = axes('Parent',f_res); hold(ax_res,'on'); box(ax_res,'on'); grid(ax_res,'on');
 plot(ax_res, t_vec, real(rel_e), 'LineWidth', 1.4, 'Color',[0.85 0.33 0.1], ...
     'DisplayName','$\Re\,\delta(\tau)$ (Real part)');
 plot(ax_res, t_vec, imag(rel_e), 'LineWidth', 1.4, 'Color',[0 0.45 0.74], ...
     'DisplayName','$\Im\,\delta(\tau)$ (Imaginary part)');
 yline(ax_res, 0, 'k--', 'LineWidth', 1.0, 'DisplayName','$\delta=0$ (weak closure)');
 
 set(ax_res, 'FontSize', FS_AX);

 %title(ax_res,'Deviation from approximate quantum mass shell', ...
 %    'Interpreter','latex','FontSize',FS_TIT);
 xlabel(ax_res,'$\tau/Tc$', 'Interpreter', 'latex','FontSize',FS_LAB);
 ylabel(ax_res,'$(\pi_\mu \pi^\mu - m^{2} c^2)/(m^2 c^{2})$', 'Interpreter', 'latex','FontSize',FS_LAB);
 lg = legend(ax_res,'Location','northeast','Interpreter', 'latex');
 lg.FontSize = FS_LEG;

 ylim(ax_res,[-3.2e-12 3.2e-12])
 exportgraphics(f_res, fullfile('figures','figure_mass_shell_deviation.pdf'), ...
     'ContentType','image', 'BackgroundColor','none', 'Resolution',600);
 


%%
function [rel_err] = mass_shell_deviation(w, C)
    wnorm   = w / C.c;
    rel_err = sum(wnorm .* (C.eta*wnorm), 1) - 1;
end
