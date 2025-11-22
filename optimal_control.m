function [w] = optimal_control(n, s_spin, z, ky, pz, B, epsilon, lambda_gauge)
    C = phys_constants();

    A  = complex([0; 0; B*z(2); 0]);    
    g = grad_ln_phi(n, z(2), ky, pz, B, lambda_gauge);

    w = epsilon*(1i*C.hbar/C.m) * g - (epsilon * C.e_q/C.m) * A;

    if strcmpi(lambda_gauge, 'weak_mass_sell')
        sp2 = w(2)^2 + w(3)^2 + w(4)^2;
        w0_shell = sqrt(C.c^2 + sp2);
        if real(w0_shell) < 0, w0_shell = -w0_shell; end
        w(1) = w0_shell;
    elseif strcmpi(lambda_gauge, 'dirac')
        E = dirac_landau_energy(n, s_spin, B, pz, C);
        w(1) = epsilon * (E / (C.m * C.c));
    else
        error('optimal_control: only lambda_gauge=dirac or lambda_gauge=weak_mass_sell are implemented.');
    end
end


function E = dirac_landau_energy(n, s_spin, B, pz, C)
%DIRAC_LANDAU_ENERGY  Relativistic Landau spectrum E_{n,s}, Eq. (A.3).
%
%   E^2 = m^2 c^4 + p_z^2 c^2 + (2n + 1 + s) ħ ω_c m c^2,
%   with ω_c = |e| B / m.  (See Appendix A.1 of the paper.)
%

    wc = abs(C.e_q) * B / C.m;  % cyclotron frequency ω_c

    term_landau = (2*n + 1 + s_spin) * C.hbar * wc * C.m * C.c^2;

    E2 = (C.m^2 * C.c^4) + (pz^2) * (C.c^2) + term_landau;

    E = sqrt(E2);
end