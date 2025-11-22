%% run_optimal_action_report.m
% Report action components for all (n,s) pairs

clear; 

B         = 0.1;         % Tesla
Ntraj     = 10000;
Nt        = 10000;       % steps per T
ky  = 0;
pz  = 0;
epsilon = 1; % 1- electron ; 2 -positron
lambda_gauge = 'weak_mass_sell';

% ---- constants & scalars ----
C   = phys_constants();
wc  = abs(C.e_q)*B/C.m;
Tc  = 2*pi/wc;

T   = 5*Tc;
dt  = T / Nt;

% ---- initial x ----
lB     = sqrt(C.hbar/(abs(C.e_q)*B));
X0     = - C.hbar * ky / (C.e_q * B);
x_init = X0 + lB;

s_norm = (C.hbar*wc)*T;
ito_covar_corr = 0.5;

% ---- (n,s) pairs to run ----
pairs = [0 -1; 0 +1; 1 -1; 1 +1];
np    = size(pairs,1);

fprintf('run_optimal_action_report: B=%.3g T, Ntraj=%d, Nt=%d, T=%.3e s (Tc=%.3e s)\n', ...
        B, Ntraj, Nt, T, Tc);

fprintf('All numbers below are normalized by (ħ ω_c T).\n\n');

z_init = complex([0; x_init; 0; 0]);

for i = 1:np
    n      = pairs(i,1);
    s_spin = pairs(i,2);

    S_tot_vec = zeros(1, Ntraj);
    S_kin_vec = zeros(1, Ntraj);
    S_EM_vec  = zeros(1, Ntraj);
    S_SP_vec  = zeros(1, Ntraj);

    parfor r = 1:Ntraj
        rng(137 + r*1000 + i*5000);  % reproducible
        [St, Sk, Se, Ss] = dirac_landau_action(n, s_spin, B, dt, Nt, z_init, ky, pz, [], epsilon, lambda_gauge);
        S_tot_vec(r) = real(St);
        S_kin_vec(r) = real(Sk);
        S_EM_vec(r)  = real(Se);
        S_SP_vec(r)  = real(Ss);
    end

    S_tot = mean(S_tot_vec);
    S_kin = mean(S_kin_vec);
    S_EM  = mean(S_EM_vec);
    S_SP  = mean(S_SP_vec);

    S_th_tot = dirac_landau_action_theory(n, s_spin, B, T);
    S_th_EM  = -(n + 0.5) * (C.hbar*wc) * T;              % EM piece
    S_th_SP  = -(C.g_factor/4) * (C.hbar*wc) * s_spin * T;  % spin piece
    S_th_kin = S_th_tot - S_th_EM - S_th_SP;              % remainder

    vEM     = real(S_EM) / s_norm + ito_covar_corr;
    vEM_th  = real(S_th_EM) / s_norm;
    vSP     = real(S_SP) / s_norm;
    vSP_th  = real(S_th_SP) / s_norm;
    vKIN    = real(S_kin);
    vKIN_th = real(S_th_kin) / s_norm;
    vTOT    = real(S_tot) / s_norm + ito_covar_corr;
    vTOT_th = real(S_th_tot) / s_norm;

    % ---- 95% percentile bootstrap CIs (normalized) ----
    nBoot = 5000;
    rng(42);
    boot_mean = @(x) mean(x / s_norm);

    bmTOT = zeros(nBoot,1);
    bmKIN = zeros(nBoot,1);
    bmEM = zeros(nBoot,1);
    bmSP = zeros(nBoot,1);
    N = numel(S_tot_vec);

    for bb = 1:nBoot
        idx = randi(N, [N,1]);         % resample with replacement
        bmTOT(bb) = boot_mean(S_tot_vec(idx));
        bmKIN(bb) = boot_mean(S_kin_vec(idx));
        bmEM(bb)  = boot_mean(S_EM_vec(idx));
        bmSP(bb)  = boot_mean(S_SP_vec(idx));
    end
    vTOT_ci = quantile(bmTOT, [0.025 0.975]).' + ito_covar_corr;
    vKIN_ci = quantile(bmKIN, [0.025 0.975]).';
    vEM_ci  = quantile(bmEM,  [0.025 0.975]).' + ito_covar_corr;
    vSP_ci  = quantile(bmSP,  [0.025 0.975]).';

    fprintf('n=%d, s=%+d\n', n, s_spin);
    fprintf('  S_EM/(ħω_c T)   : %.4f  [%.4f, %.4f]   (theory %.4f)   Δ=%.4f\n', ...
            vEM,  vEM_ci(1), vEM_ci(2),  vEM_th,  vEM -  vEM_th);
    fprintf('  S_SP/(ħω_c T)   : %.4f  [%.4f, %.4f]   (theory %.4f)   Δ=%.4f\n', ...
            vSP,  vSP_ci(1), vSP_ci(2),  vSP_th,  vSP -  vSP_th);
    fprintf('  S_kin/(ħω_c T)  : %.4f  [%.4f, %.4f]   (theory %.4f)   Δ=%.4f\n', ...
            vKIN, vKIN_ci(1), vKIN_ci(2), vKIN_th, vKIN - vKIN_th);
    fprintf('  S_tot/(ħω_c T)  : %.4f  [%.4f, %.4f]   (theory %.4f)   Δ=%.4f\n\n', ...
            vTOT, vTOT_ci(1), vTOT_ci(2), vTOT_th, vTOT - vTOT_th);
end