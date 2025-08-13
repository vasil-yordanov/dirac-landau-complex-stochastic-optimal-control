clear;

files = {
   'data/S_sweep_T_raw_Nrun1000_NtcPerTc100.mat'
   'data/S_sweep_T_raw_Nrun1000_NtcPerTc1000.mat',
   'data/S_sweep_T_raw_Nrun1000_NtcPerTc10000.mat'
};

C = phys_constants();   % must define: C.hbar, C.e_q, C.m

%% --- plot ---
f = figure; hold on; grid on; box on;
set(gca,'XScale','log');   % always log-x

leg = cell(1, numel(files));

for k = 1:numel(files)
    S = load(files{k});
    R = S.res;                         
    B = R.meta.B;

    wc  = abs(C.e_q)*B/C.m;
    Tc  = 2*pi/wc;

    T  = R.T_list(:).'/Tc;             % 1 x NT
    Sn = R.S_norm_by_T;                % 1 x NT cells

    mean_Sn  = cellfun(@mean, Sn);     % 1 x NT

    plot(T, mean_Sn, '-o', 'LineWidth', 1.6, 'MarkerSize', 5);

    leg{k} = sprintf('\\Delta t = T_c/%d', R.meta.Ntc_per_Tc);
end

ylim([-5 0]);

set(gca, 'FontSize', 14);

t.FontSize = 19;

xt = xlabel('$T/T_c$', 'Interpreter', 'latex');
xt.FontSize = 19;
xt = ylabel('$\Re\, \langle S \rangle /(\hbar \omega_c T)$', 'Interpreter','latex');
xt.FontSize = 19;
lt = legend(leg, 'Location','northwest');
lt.FontSize = 14;

filename = sprintf('S_sweep_T.pdf');

exportgraphics(f, filename, ...
    'ContentType', 'vector', ...
    'BackgroundColor', 'none', ...
    'Resolution', 600);

fprintf('Exported: %s\n', filename);
