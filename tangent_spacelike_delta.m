function delta = tangent_spacelike_delta(wopt, p_dev, k)
% tangent_spacelike_delta: Generates a tangent spacelike perturbation δ for wopt.
%
% Parameters:
%   p_dev.a_max  - maximum amplitude scaling (dimensionless, R/c)
%   p_dev.mode   - 'random' -> fully random directions (not applicable for SOC minimality test)
%   p_dev.mode   - 'deterministic'  -> deterministic direction, among M fixed tangent directions

    C = phys_constants();
    den = wopt.'*(C.eta*wopt);  % bilinear (no conjugation)

    if abs(den) < 1e-18
        delta = complex(zeros(4,1));
        return;
    end

    % Default values
    if ~isfield(p_dev, 'max_trials'), p_dev.max_trials = 100; end
    if ~isfield(p_dev, 'a_max'), p_dev.a_max = 0.02; end
    if ~isfield(p_dev, 'mode'), p_dev.mode = 'random'; end

    tol = 1e-22;  % avoid nearly-null q
    max_trials = 100; % max number of rejection iterations

    % =====================================================
    % MODE 1: Fully random directions
    % =====================================================
    if strcmpi(p_dev.mode, 'random')
        for t = 1:max_trials
            v = (randn(4,1) + 1i*randn(4,1)) / sqrt(2);
            v = v - wopt * ((wopt.'*(C.eta*v)) / den);  % project to tangent

            q = v.'*(C.eta*v);
            if abs(q) < tol, continue; end
            if real(q) >= 0, continue; end

            a = p_dev.a_max * rand();
            delta = (a*C.c)/sqrt(abs(q)) * v;
            return
        end
        error('Failed to generate δ after %d trials in RANDOM mode', p_dev.max_trials);
    end

    % =====================================================
    % MODE 2: Deterministic fixed directions
    % =====================================================
    
    if strcmpi(p_dev.mode, 'deterministic')    
        Uraw = [ [ones(1,6), zeros(1,6)] ;  
            [eye(3), 1i*eye(3), eye(3), 1i*eye(3)] ];

        % Require at least some spatial component in the seed
        assert( norm(Uraw(2:4)) > 1e-12, ...
            '[ASSERT FAIL] Pure-time seed detected: u_raw = [%g %g %g %g]', ...
            Uraw(1), Uraw(2), Uraw(3), Uraw(4) );

        M = size(Uraw, 2);
    
        % Pick direction deterministically from timestep k
        j = 1 + mod(k - 1, M);
        
        u_raw = Uraw(:, j);
    
        % Project to tangent: wopt.'*(eta*u) = 0
        u = u_raw - wopt * ((wopt.'*(C.eta*u_raw)) / den);
    
        % Fallback if degenerate after projection
        if norm(u) < 1e-18
            u_raw = [0; 1; 0; 1i];
            u = u_raw - wopt * ((wopt.'*(C.eta*u_raw)) / den);
        end
    
        % Enforce spacelike branch: rotate phase so q is negative real
        q = u.'*(C.eta*u);
        if ~(real(q) < 0 && imag(q) < 0)
            % Add a tiny negative phase bias to push q strictly into lower half-plane
            %phi_star = -3*pi/4;
            phi_star = pi;
            theta = 0.5*(phi_star - angle(q));  % tiny bias ensures imag(q)<0
            u = exp(1i*theta) * u;
            q = u.'*(C.eta*u);
        end
    
        assert( real(q) < 0, '[ASSERT FAIL] real(q)<0' );

        % Normalise Minkowski magnitude so u.'*eta*u ≈ -1
        u = u / sqrt(abs(q));
    
        % Deterministic amplitude: always use a_max
        a = p_dev.a_max;
    
        % Scale to target Minkowski amplitude: sqrt(|u.'*eta*u|)/c = a  (here |…|=1)
        delta = (a*C.c) * u;
        return
    end

    error('Unknown p_dev.mode = %s', p_dev.mode);
end