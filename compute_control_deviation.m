function R = compute_control_deviation(dw)
    C  = phys_constants();
    Nt = size(dw,2);
    R_delta_sum = 0;
    for k = 1:Nt
        delta = dw(:,k);
        R_delta_sum = R_delta_sum + sqrt(abs(delta.'*(C.eta*delta)))/C.c;
    end
    R = R_delta_sum / Nt;
end