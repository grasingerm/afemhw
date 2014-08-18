function N = Q4_N(xi, eta)
    N = zeros(2,8);
    
    N(1,1) = Q4_N1(xi, eta);
    N(1,3) = Q4_N2(xi, eta);
    N(1,5) = Q4_N3(xi, eta);
    N(1,7) = Q4_N4(xi, eta);
    
    for i = 1:4
        N(2,2*i) = N(1,2*i-1);
    end
    
endfunction
