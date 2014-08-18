function F_ref = Q4_F_ref(xi, eta, X_coords)
    F_ref = Q4_dN_dxi(xi, eta) * X_coords;
endfunction
