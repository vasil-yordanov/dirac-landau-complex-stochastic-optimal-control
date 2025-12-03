% run_sweep_T_report_S_statistics.m

clear; clc;

% ---- start/ensure parallel pool (uses default local profile) ----
if isempty(gcp('nocreate'))
    parpool('local');   % or parpool to let MATLAB choose the default
end

%% ----- fixed simulation params -----
n      = 0;
s_spin = 1;
B      = 0.1;                % Tesla
ky     = 0;
pz     = 0;
epsilon = 1; % 1- electron ; 2 -positron

C   = phys_constants();
wc  = abs(C.e_q)*B/C.m;
Tc  = 2*pi/wc;

% initial condition (x two magnetic length off the guiding center)
lB     = sqrt(C.hbar/(abs(C.e_q)*B));
X0     = - C.hbar * ky / (C.e_q * B);
z_init = complex(zeros(4,1));
z_init(2) = X0 + 2*lB;

% No control deviations
p_dev = [];

%% ----- sweep settings -----
Ntc_per_Tc = 1000;   % steps per cyclotron period
m_list     = [0.001 0.005 0.007 0.01 0.02 0.1 0.2 0.3 0.4 0.5 1 1.5 2 5 10 20 50 70 100 150 200 250 300, 400, 500, 1000];

if Ntc_per_Tc == 100
   m_list = m_list(4:15);
elseif Ntc_per_Tc == 1000
   m_list = m_list(1:20);
end

T_list  = m_list * Tc;
Nt_list = m_list * Ntc_per_Tc;

Nrun      = 1000;   % seeds per T
base_seed = 4242;

%% ----- storage (normalized only) -----
mean_S_norm = nan(size(T_list));
std_S_norm  = nan(size(T_list));
S_norm_by_T = cell(size(T_list));   % each cell: 1xNrun S_norm samples

%% ----- sweep -----
t_all = tic;
for it = 1:numel(T_list)
    m  = m_list(it);
    T  = T_list(it);
    Nt = Nt_list(it);
    dt = T / Nt;

    fprintf('==> T = %g * Tc | Nt = %d | %d seeds...\n', m, Nt, Nrun);
    t_T = tic;

    % per-T buffer (sliced variable for parfor)
    S_vec = nan(1, Nrun);

    ito_covar_corr = 0.5;

    % ---- parallel seed loop ----
    parfor ir = 1:Nrun
        % deterministic per-seed RNG (order-independent)
        rng(base_seed + ir, 'twister');

        [S, S_Lk, S_LEM, S_Lsp] = dirac_landau_action(n, s_spin, B, dt, Nt, z_init, ky, pz, [], epsilon);
        S_vec(ir) = real(S_LEM + S_Lsp) / (C.hbar * wc * T)  + ito_covar_corr;   % See Appendix A.5
    end

    % aggregate
    mean_S_norm(it) = mean(S_vec);
    std_S_norm(it)  = std(S_vec);
    S_norm_by_T{it} = S_vec;

    fprintf('   done: mean(S_norm)=%.6g, std=%.3g  | T=%.3g s (m=%g) | elapsed %.1fs\n', ...
            mean_S_norm(it), std_S_norm(it), T, m, toc(t_T));
end
fprintf('All T done in %.1fs\n', toc(t_all));

%% ----- save for post-processing -----
res = struct();

% minimal metadata
res.meta = struct('n', n, 's_spin', s_spin, 'B', B, 'ky', ky, 'pz', pz, ...
                  'Ntc_per_Tc', Ntc_per_Tc, 'Nrun', Nrun);

% sweep axes
res.m_list  = m_list;
res.T_list  = T_list;
res.Nt_list = Nt_list;

% raw normalized totals (one 1xNrun vector per T)
res.S_norm_by_T = S_norm_by_T;

% descriptive filename: include Nrun, Ntc_per_Tc, and m range
fname = sprintf('S_sweep_T_raw_Nrun%d_NtcPerTc%d.mat', Nrun, Ntc_per_Tc);

save(fname, 'res', '-v7.3');
fprintf('Saved raw results to %s\n', fname);