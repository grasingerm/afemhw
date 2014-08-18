function s = S_vgt(C_SE, F)
    em = E(F);
    evgt = [em(1,1); em(2,2); em(1,2)];
    s = C_SE * evgt;
endfunction
