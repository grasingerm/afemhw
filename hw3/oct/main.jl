require("debug.jl");
require("fem.jl");
require("integration.jl");
require("mat.jl");
require("nonlinear.jl");

const Ey = 100;
const nu = 0.45;

const Ey_star = Ey / (1-nu^2);
const nu_star = nu / (1-nu);

const C_SE = tan_stiff_moduli_plane_strain(Ey_star, nu_star);

const X = [[0 0; 5 2; 5 5; 0 8]]; # global coords
# no connectivity needed as there is only one element ...
const tau = [1000; 500]; # traction

const P = ext_force_vec();

function main(load_steps::Integer, max_iters::Intger, eps::FloatingPoint)
  # Maps with respect to referential coords only need to be computed once
  # F_xi_o_t =
  # F_xi_o =
  # J_xi_o =

  const dlambda = 1/load_steps;

  for i = 1:load_steps
    lambda = i*dlambda;
    P_ext = P*lambda;
  end
end
