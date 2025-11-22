% run_demo_local_optimality.m
% Scatter of ΔRe(S) vs R with per-step tangent spacelike deviations.
% Optimal policy (no deviation) is shown at (0,0) with a star.

clear;
rng(7);  % fix seed so the baseline is reproducible

%% Simulation params 
n         = 0;            % Landau level index
s_spin    = -1;           % spin ±1
B         = 0.1;          % magnetic field [Tesla]
z_init    = complex(zeros(4,1));  % starting 4-position/phase-space state
ky        = 0;
pz        = 0;
N_tc      = 800;          % number of time steps
N_pol     = 200;          % number of scattered policies
N_opt_avg = 200;          % number of optimal poliies to average

lambda_gauge = 'weak_mass_sell';
% deterministic-mode control: deterministic directions, among several fixed comtrol directions
mode      = 'deterministic';
A_max     = 2e-5;         % max deviation amplitude
%N_avg     = 1;          % number of policies to average
N_avg     = 100;

% random-mode control: fully random direction & amplitude
% in this mode M != 1 makes no sense since we wll not averaging
% under the same determinicstic controls according the SOC theory
% mode = 'random';
% A_max = 2e-3;           % max deviation amplitude
% N_avg    = 1;           % number of policies to average

C   = phys_constants();
wc  = abs(C.e_q)*B/C.m;
Tc  = 2*pi/wc;        % Cyclotron period
T   = 5*Tc;
dt  = T / N_tc;
Nt  = N_tc;
epsilon = 1; % 1- electron ; 2 -positron

%% Optimal control, no deviation
    
S0_accum = 0;
for m = 1:N_opt_avg
    rng(7 + 1000*m, 'twister');   % independent noise per replicate
    [Sm, ~, ~, ~, ~] = dirac_landau_action(n, s_spin, B, dt, Nt, z_init, ky, pz, [], epsilon, lambda_gauge);
    S0_accum = S0_accum + Sm;
end
S0 = S0_accum / N_opt_avg;

R = zeros(1,N_pol);
S = complex(zeros(1,N_pol));

p_dev = struct('a_max', NaN, 'mode', mode);

t_T = tic;
for i = 1:N_pol
    p_dev.a_max = A_max * (i-1) / max(1, N_pol-1);

    % Average S and R over M independent Brownian paths for this *same* control
    S_accum = 0;
    R_accum = 0;
    for m = 1:N_avg
        rng(137 + 1000*i + m, 'twister');   % independent noise per replicate
        [S_im, ~, ~, ~, R_im] = dirac_landau_action(n, s_spin, B, dt, Nt, z_init, ky, pz, p_dev, epsilon, lambda_gauge);
        S_accum = S_accum + S_im;
        R_accum = R_accum + R_im;
    end

    S(i) = S_accum / N_avg;
    R(i) = R_accum / N_avg;

    if mod(i, max(1,round(N_pol/10)))==0 || i == N_pol
            fprintf('   %3.0f%% (%d/%d) | elapsed %.1fs\r', 100*i/N_pol, i, N_pol, toc(t_T));
    end
end

%% --- Theoretical curvature: ΔS ≈ κ_th * R^2 ---
dS   = real(S) - real(S0);
R2   = R.^2;

xfit = linspace(0, max(R)*1.05, 200);

kappa_fit = (R2 * dS.') / (R2 * R2.');    % scalar slope
kappa_th = kappa_theory(Nt, dt, mode);

%5 -------- Plot --------------
f = figure; hold on; box on; grid on;

scatter(R, dS, 10, 'filled', 'MarkerFaceAlpha', 0.75, 'DisplayName','policies');
plot(0, 0, 'k*', 'MarkerSize', 7, 'DisplayName','optimal policy');

if N_avg >1
    plot(xfit, kappa_th * (xfit.^2), 'k--', 'LineWidth', 1.4, ...
         'DisplayName', 'theory: \Re \langle\Delta S\rangle = 1/2 m c^2 T R^2');
    plot(xfit, kappa_fit * (xfit.^2), 'r-', 'LineWidth', 1.2, ...
         'DisplayName', sprintf('fit: \\Re \\langle\\Delta S\\rangle = k R^2'));
else
    plot(xfit, kappa_th * (xfit.^2), 'k--', 'LineWidth', 1.4, ...
         'DisplayName', 'theory: \Re \Delta S = 1/2 m c^2 T R^2');
    plot(xfit, kappa_fit * (xfit.^2), 'r-', 'LineWidth', 1.2, ...
         'DisplayName', sprintf('fit: \\Re \\Delta S = k R^2'));

end

set(gca, 'FontSize', 14);

if N_avg > 1
    t = title(sprintf('Averaged over %d samples per control', N_avg), 'Interpreter','latex');
else
    t = title('Single sample per control', 'Interpreter','latex');
end
t.FontSize = 19;


xt = xlabel('$R_{\delta}/c$', 'Interpreter', 'latex');
xt.FontSize = 19;
if N_avg >1
    xt = ylabel('$\Re \, \langle S\rangle - \Re \, \langle S^{*}\rangle  \, [\mathrm{J \cdot s}]$', 'Interpreter','latex');
else
    xt = ylabel('$\Re \, S - \Re \, S^{*} \, [\mathrm{J \cdot s}]$', 'Interpreter','latex');
end
xt.FontSize = 19;
lt = legend('Location','northwest');
lt.FontSize = 14;

ylim([-1.5e-32 3.5e-32]);

filename = sprintf('local_optimality_M%d.pdf', N_avg);

exportgraphics(f, filename, ...
    'ContentType', 'vector', ...
    'BackgroundColor', 'none', ...
    'Resolution', 600);

fprintf('Exported: %s\n', filename);

%% Functions

function kappa_th = kappa_theory(Nt, dt, mode)
    %KAPPA_THEORY  Theoretical curvature coefficient for ΔS ≈ kappa_th * R^2.
    %
    % Assumptions:
    %   • Small Deviations δ are tangent & spacelike (SR mass-shell geometry).
    %   • If mode == 'random' 
    %       - Per-step amplitude a = sqrt(|δ^T η δ|)/c is drawn Uniform(0, a_max)  → χ = E[a^2]/E[a]^2 = 4/3.
    %       - Complex spacelike sampler uses pure rejection in the left half-plane,
    %         with phase θ ~ Uniform([π/2, 3π/2]) → γ = E[-cos θ] = 2/π.
    %
    % Inputs:
    %   Nt   : number of time steps
    %   dt   : time step size
    %   mode : 'random' or 'deterministic'
    %
    % Output:
    %   kappa_th : theoretical curvature coefficient

    C = phys_constants();
    T = Nt * dt;

    if strcmpi(mode, 'random')
        chi   = 4/3;      % Uniform(0, a_max)
        gamma = 2/pi;     % rejection sampler, uniform phase in left half-plane
        kappa_th = 0.5 * C.m * C.c^2 * T * chi * gamma;
        return
    end

    kappa_th = 0.5 * C.m * C.c^2 * T;
end