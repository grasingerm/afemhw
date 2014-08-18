function dN_dxi = Q4_dN_dxi(xi, eta)

    dN_dxi(1,1) = Q4_dN1_dxi(eta);
    dN_dxi(1,2) = Q4_dN2_dxi(eta);
    dN_dxi(1,3) = Q4_dN3_dxi(eta);
    dN_dxi(1,4) = Q4_dN4_dxi(eta);
    
    dN_dxi(2,1) = Q4_dN1_deta(xi);
    dN_dxi(2,2) = Q4_dN2_deta(xi);
    dN_dxi(2,3) = Q4_dN3_deta(xi);
    dN_dxi(2,4) = Q4_dN4_deta(xi);
    
endfunction
