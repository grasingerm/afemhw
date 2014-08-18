function sb = S_bar(S_vgt)
    sb = zeros(4,4);
    sb(1,1) = S_vgt(1);
    sb(1,2) = S_vgt(3);
    sb(2,1) = S_vgt(3);
    sb(2,2) = S_vgt(2);
    
    for i = 1:2
        for j = 1:2
            sb(i+2,j+2) = sb(i,j);
        end
    end
endfunction
