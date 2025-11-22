function [S, S_Lk, S_LEM, S_Lsp, R, w_store, z_store] = dirac_landau_action(
        n, s_spin, B, dt, Nt, z_init, ky, pz, p_dev, epsilon, lambda_gauge
    )

    C = phys_constants();
    wc = abs(C.e_q)*B/C.m;

    sigma =  sqrt(C.hbar/C.m);

    S = 0+0i; S_Lk = 0+0i; S_LEM = 0+0i; S_Lsp = 0+0i;
    R_accum = 0;
    z = z_init;

    % Optional trajectory recording
    if nargout >= 6
        w_store = complex(zeros(4, Nt));
        z_store = complex(zeros(4, Nt+1));
        z_store(:,1) = z;
    end

    for k = 1:Nt
        A  = complex([0; 0; B*z(2); 0]);    

        [w] = optimal_control(n, s_spin, z, ky, pz, B, epsilon, lambda_gauge);

        if ~isempty(p_dev)               
            delta = tangent_spacelike_delta(w, p_dev, k);  
            R_accum = R_accum + sqrt(abs(delta.'*(C.eta*delta)))/C.c;

            w = w + delta;
        end

        s_val = w.'*(C.eta*w);

        Lk  = -C.m*C.c*sqrt(s_val) + C.m*C.c^2;
        LEM =  - epsilon * C.e_q * (A.'*(C.eta*w));
        Lsp = -0.25 * C.g_factor * C.hbar * wc * s_spin;

        S_Lk  = S_Lk  + Lk*dt;
        S_LEM = S_LEM + LEM*dt;
        S_Lsp = S_Lsp + Lsp*dt;

        S     = S + (Lk + LEM + Lsp)*dt;

        z = z + w*dt + sigma*((eye(4) + 1i*epsilon*C.eta)*randn(4,1))*sqrt(dt);


        if nargout >= 6
            w_store(:,k) = w;
        end

        if nargout >= 6
            z_store(:,k+1) = z;
        end
    end

    R = R_accum / Nt;
end
