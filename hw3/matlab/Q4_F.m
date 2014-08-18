function F = Q4_F(dN_dref, u)
    F = eye(2);
    
    for i = 1:2
        for j = 1:2
            for k = 1:4
                F(i,j) += dN_dref(i,k) * u(2*k+j-2)
            end
        end
    end
endfunction
