
function C = phys_constants()
    C.hbar = 1.0545718e-34;
    C.m    = 9.10938356e-31;
    C.e_q  = -1.60217662e-19;      % electron (q = -e)
    C.c    = 2.99792458e8;
    C.eta  = diag([1 -1 -1 -1]);   % Minkowski (+,-,-,-)
    C.g_factor   = 2.0;
end