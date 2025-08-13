% grad_ln_phi: contravariant 4-gradient of ln φ_n for Landau gauge A=(0,0,Bx,0)
% g = [ NaN;  -∂x lnφ;  -∂y lnφ;  -∂z lnφ ]
% Supports n=0,1; x may be complex.
function g = grad_ln_phi(n, x, ky, pz, B)
    C  = phys_constants();
    lB = sqrt(C.hbar/(abs(C.e_q)*B));

    % Guiding-center shift (Landau gauge): X0 = -ħ ky / (e B)
    X0 = - C.hbar * ky / (C.e_q * B);
    xi = x - X0;

    % ∂x ln φ_n(ξ)
    switch n
        case 0
            % φ0 ∝ exp(-ξ^2/(2 lB^2)) => d/dx lnφ0 = -ξ/lB^2
            dlogphi_dx = -xi / (lB^2);

        case 1
            % φ1 ∝ ξ exp(-ξ^2/(2 lB^2)) => d/dx lnφ1 = 1/ξ - ξ/lB^2
            dlogphi_dx = (1/xi) - xi/(lB^2);
        otherwise
            error('grad_ln_phi: only n=0 or n=1 are implemented.');
    end

    g_x =  dlogphi_dx;
    g_y =  1i * ky;
    g_z =  1i * pz / C.hbar;

    g = complex([NaN; -g_x; -g_y; -g_z]);
end