require("mat.jl")

Ey = 100;
nu = 0.45;

Ey_star = Ey / (1-nu^2);
nu_star = nu / (1-nu);

C_SE = tan_stiff_moduli_plane_strain(Ey_star, nu_star);

coords = [[0 0; 5 2; 5 5; 0 8]]; # global coords
# no connectivity needed as there is only one element ...
tau = [1000; 500]; # traction

