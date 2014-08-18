function green_strain = E(F)
    green_strain = F'F - eye(size(F))
endfunction
