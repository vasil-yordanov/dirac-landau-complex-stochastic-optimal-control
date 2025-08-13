function [S_th_tot, S_th_EM, S_th_SP, S_th_kin] = dirac_landau_action_theory(n, s_spin, B, T)

    C   = phys_constants();
    wc  = abs(C.e_q)*B/C.m;                   % cyclotron frequency
    alpha = (2*n + 1 + s_spin) * (C.hbar*wc) / (C.m*C.c^2);
    Erel  = C.m*C.c^2 * sqrt(1 + alpha);      % relativistic Landau energy (pz=0)

    % Total theoretical action (Dirac):
    S_th_tot = -(Erel - C.m*C.c^2) * T;

    % EM and spin pieces (closed forms):
    S_th_EM  = - (n + 0.5)   * (C.hbar*wc) * T;
    S_th_SP  = - (C.g_factor/4)* (C.hbar*wc) * s_spin * T;

    % Kinetic by consistency with the exact total:
    S_th_kin = -S_th_tot + S_th_EM + S_th_SP;
end