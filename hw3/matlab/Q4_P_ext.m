function P_e = Q4_P_ext(tau, xi, eta, F_ref, edge)
    N_ip = N(xi,eta)
    if edge == 1
        N_t = N_ip(:,1:2);
    elseif edge == 2
        N_t = N_ip(:,2:3);
    elseif edge == 3
        N_t = N_ip(:,3:4);
    else
        N_t = [N_ip(:,4) N_ip(:,1)];
    end
    
    P_e = N_t' * tau * det(F_ref);
endfunction
