function dN_dref = Q4_dN_dref(F_ref, dN_dxi)
    dN_dref = inv(F_ref)' * dN_dxi;
endfunction
