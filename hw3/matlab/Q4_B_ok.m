function B_o = Q4_B_ok(F, dN_dref, node_number)
    
    for i = 1:2
        for j = 1:2
            B(i,j) = F(i,j)*dN_dref(i,node_number);
        end
    end
    
    B(3,1) = F(1,1)*dN_dref(2,node_number) + F(2,1)*dN_dref(1,node_number);
    B(3,2) = F(2,2)*dN_dref(1,node_number) + F(1,2)*dN_dref(2,node_number);
endfunction
