% grad_ln_phi: contravariant 4-gradient of ln φ_n for Landau gauge A=(0,0,Bx,0)
% g = [ NaN;  -∂x lnφ;  -∂y lnφ;  -∂z lnφ ]
% Supports n=0,1
% x may be complex
% Supports lambda_gauge=dirac or lambda_gauge=weak_mass_sell

function g = grad_ln_phi(n, x, ky, pz, B, lambda_gauge)
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
            if strcmpi(lambda_gauge, 'weak_mass_sell')
                dlogphi_dx = (1./xi) - xi/(lB^2);
            elseif strcmpi(lambda_gauge, 'dirac')
                xi_cut = 1e-3 * lB;  % small fraction of magnetic length
                xi_reg = xi;
                if abs(xi_reg) < xi_cut
                    % keep the sign of Re(xi) so we don't get stuck at exactly 0
                    if real(xi_reg) >= 0
                        xi_reg = xi_cut + 1i*imag(xi_reg);
                    else
                        xi_reg = -xi_cut + 1i*imag(xi_reg);
                    end
                end
                dlogphi_dx = (1./xi_reg) - xi/(lB^2);
            else
                error('grad_ln_phi: only lambda_gauge=dirac or lambda_gauge=weak_mass_sell are implemented.');
            end
        otherwise
            error('grad_ln_phi: only n=0 or n=1 are implemented.');
    end

    g_x =  dlogphi_dx;
    g_y =  1i * ky;
    g_z =  1i * pz / C.hbar;

    g = complex([NaN; -g_x; -g_y; -g_z]);
end