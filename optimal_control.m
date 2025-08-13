function [w] = optimal_control(n, z, ky, pz, B)
    C = phys_constants();

    A  = complex([0; 0; B*z(2); 0]);    
    g = grad_ln_phi(n, z(2), ky, pz, B);

    w = (1i*C.hbar/C.m) * g - (C.e_q/C.m) * A;

    sp2 = w(2)^2 + w(3)^2 + w(4)^2;
    w0_shell = sqrt(C.c^2 + sp2);
    if real(w0_shell) < 0, w0_shell = -w0_shell; end

    w(1) = w0_shell;
    
end

