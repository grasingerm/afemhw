require("debug.jl");
require("fem.jl");
require("integration.jl");
require("mat.jl");
require("nonlinear.jl");

# remind user that debugging may be disabled
if !@ndebug
  write(STDERR, "Warning! Debugging is on. To improve performance turn off debugging\n\n");
end

const Ey = 100;
const nu = 0.45;

const Ey_star = Ey / (1-nu^2);
const nu_star = nu / (1-nu);

const C_SE = tan_stiff_moduli_plane_strain(Ey_star, nu_star);

const X = Array{Float64,2}[Float64[0 0; 5 2; 5 5; 0 8]]; # global coords
# no connectivity needed as there is only one element ...
const tau = Float64[1000; 500]; # traction

# Calculate external force vector
xi = 1.0; #!< surface that traction is applied to
eta = gp1[1]; #!< only need one integration point
unit_normal = Float64[1; 0];
const P = gw1[1] * traction_to_extforcevec(xi, eta, X[1], tau, unit_normal)
println("P_ext = \n", P);

# Check that only one integration point is needed
if !@ndebug
  P2 = zeros(8,1);
  xi = 1.0;
  for i in 1:2
    eta = gp2[i];
    P2 += gw2[i] * traction_to_extforcevec(xi, eta, X[1], tau, unit_normal);
  end
  @dassert(norm(P - P2) < 1e-6, 
    "only one integration pt should be needed for exact int!");

  # what if element was same geom as parent element?
  # x = [-1.0 -1.0; 1.0 -1.0; 1.0 1.0; -1.0 1.0];
  # println("P_par = ", gw1[1] * traction_to_extforcevec(1.0, 0.0, x, tau, unit_normal));
end

function main(load_steps::Int, max_iters::Int, eps::FloatingPoint)

  # Maps with respect to referential coords only need to be computed once
  # F_xi_o_t =
  # F_xi_o =
  # J_xi_o =

  #=
  const dlambda = 1/load_steps;

  for i = 1:load_steps
    lambda = i*dlambda;
    P_ext = P*lambda;
  end
  =#
end
