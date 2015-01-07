require("debug.jl");

macro N1_Q4(xi, eta) 
  return :((1-($xi))*(1-($eta))/4.);
end
macro N2_Q4(xi, eta)
  return :((1+($xi))*(1-($eta))/4.);
end
macro N3_Q4(xi, eta)
  return :((1+($xi))*(1+($eta))/4.);
end
macro N4_Q4(xi, eta)
  return :((1-($xi))*(1+($eta))/4.);
end

macro dN1dxi_Q4(eta) 
  return :(-(1-($eta))/4.);
end
macro dN2dxi_Q4(eta)
  return :((1-($eta))/4.);
end
macro dN3dxi_Q4(eta)
  return :((1+($eta))/4.);
end
macro dN4dxi_Q4(eta)
  return :(-(1+($eta))/4.);
end

macro dN1deta_Q4(xi)
  return :(-(1-($xi))/4.);
end
macro dN2deta_Q4(xi)
  return :((1-($xi))/4.);
end
macro dN3deta_Q4(xi)
  return :((1+($xi))/4.);
end
macro dN4deta_Q4(xi)
  return :(-(1+($xi))/4.);
end

#! \macro Check parent coordinates to make sure they fall within parent element
#!
#! \param xi Parent coordinate in x-direction
#! \param eta Parent coordinate in y-direction
#! \return dassertions to check parent coordinates for correctness
macro check_par_coords_2D(xi, eta)
  return quote
    @dassert(xi >= -1.0 && xi <= 1.0, "-1.0 <= xi <= 1.0 must be true");
    @dassert(eta >= -1.0 && eta <= 1.0, "-1.0 <= eta <= 1.0 must be true");
  end
end

#! N matrix for a Q4 element
#!
#! \param xi Parent coordinate in x-direction
#! \param eta Parent coordinate in y-direction
#! \return N matrix
function N_Q4(xi::Real, eta::Real)
  @check_par_coords_2D(xi, eta);

  return (Float64[
    @N1_Q4(xi, eta) 0 @N2_Q4(xi, eta) 0 @N3_Q4(xi, eta) 0 @N4_Q4(xi, eta) 0;
    0 @N1_Q4(xi, eta) 0 @N2_Q4(xi, eta) 0 @N3_Q4(xi, eta) 0 @N4_Q4(xi, eta);
  ]);
end

#! Tranpose of N matrix for a Q4 element
#!
#! \param xi Parent coordinate in x-direction
#! \param eta Parent coordinate in y-direction
#! \return N matrix
function N_t_Q4(xi::Real, eta::Real)
  @check_par_coords_2D(xi, eta);

  N_t = (Float64[
    @N1_Q4(xi, eta) 0;
    0 @N1_Q4(xi, eta);
    @N2_Q4(xi, eta) 0;
    0 @N2_Q4(xi, eta);
    @N3_Q4(xi, eta) 0;
    0 @N3_Q4(xi, eta);
    @N4_Q4(xi, eta) 0;
    0 @N4_Q4(xi, eta);
  ]);

  @dassert(N_Q4(xi, eta)' == N_t, "(N')' should always equal N");

  return N_t;
end

#! Map global referential coordinates to parent referential coordinates
#!
#! \param N N matrix
#! \param X Nodal global referential coordinates
#! \return Referential global coords mapped to referential parent coords
function ref_glob_to_par(N::Array{Float64,2}, X::Array{Float64,2})
  return N*X;
end

#! Derivative of shape functions with respect to parent coords, Q4
#!
#! \param xi Parent coordinate in x-direction
#! \param eta Parent coordinate in y-direction
#! \return D X N matrix, dN_j / dx_i
macro dN_dxi_Q4(xi, eta)
  return :(Float64[
    @dN1dxi_Q4($eta) @dN2dxi_Q4($eta) @dN3dxi_Q4($eta) @dN4dxi_Q4($eta);
    @dN1deta_Q4($xi) @dN2deta_Q4($xi) @dN3deta_Q4($xi) @dN4deta_Q4($xi);
    ]);
end

#! Transpose derivative of shape functions with respect to parent coords, Q4
#!
#! \param xi Parent coordinate in x-direction
#! \param eta Parent coordinate in y-direction
#! \return N X D matrix, dN_i / dx_j
macro dN_dxi_t_Q4(xi, eta)
  return :(Float64[
    @dN1dxi_Q4($eta) @dN1deta_Q4($xi);
    @dN2dxi_Q4($eta) @dN2deta_Q4($xi);
    @dN3dxi_Q4($eta) @dN3deta_Q4($xi);
    @dN4dxi_Q4($eta) @dN4deta_Q4($xi);
    ]);
end

#! Transpose of deformation gradient to map parent coords to deformed coords
#!
#! \param xi Parent coordinate in x-direction
#! \param eta Parent coordinate in y-direction
#! \param x Matrix of deformed, spatial, nodal coordinates in table notation
#! \return Deformation gradient
function defgrad_par_to_deformed_t(xi::Real, eta::Real, x::Array{Float64,2})
  @check_par_coords_2D(xi, eta);

  return @dN_dxi_Q4(xi, eta) * x;
end

#! \macro Jacobian determinant that relates diff areas from parrent to current
#!
#! \param F_xi_t Transpose deformation gradient mapping parent to spatial
#! \param N Unit vector normal to referential area
#! \return @expr Jacobian determinant for mapping area
macro jac_area(F_xi_t, N)
  return :(det($F_xi_t) * norm(($F_xi_t) \ ($N))); #!< A \ B is a polyalgorithm solver
end

#! Jacobian determinant that relates differential areas from parent to current
#!
#! \param xi Parent coordinate in x-direction
#! \param eta Parent coordinate in y-direction
#! \param X Matrix of referential global coordinates in table notation
#! \param N Unit vector normal to the surface of the parent element
#! \return Jacobian determinant for mapping area
function jac_area_par_to_deformed(xi::Real, eta::Real, X::Array{Float64,2},
  N::Array{Float64,1})

  @check_par_coords_2D(xi, eta);

  F_xi = defgrad_par_to_deformed(xi, eta, X);
  return @jac_area(F_xi, N);
end

#!!! there are two ways to calculate deformation gradient ...
# F = F_xi * F_o_xi^-1 !< probably more computationally expensive
# grad(u) + eye(2), or outer_product(partial{N}{X}, u) + eye(2)

#! Gradient of shape functions with respect to [spatial or referential] coords
#!
#! \param F_xi Deformation gradient mapping parent to spatial
#! \param xi Parent coordinate in x-direction
#! \param eta Parent coordinate in y-direction
#! \param x Matrix of global coordinates in table notation
#! \return Gradient of shape functions, N X D
function grad_Q4(F_xi::Array{Float64,2}, xi::Real, eta::Real, 
  x::Array{Float64,2})

  @check_par_coords_2D(xi, eta);
  dN_dxi_t = @dN_dxi_t_Q4(xi, eta);
  @dassert(@dN_dxi_Q4(xi, eta) == dN_dxi_t', 
    "Tranpose of a tranpose should equal itself?");

  @dassert(size(dN_dx_t) == (4,2), "Size should be N X D, or 4 X 2");
  @dassert(size(F_xi) == (2,2), "Def grad should be N X N");
  return dN_dxi_t * inv(F_xi);
end

#! Take gradient of deformation, velocity, etc.
#!
#! \param u Variable of interest in table notation, D X N
#! \param grad Gradient of shape functions, N X D (call "grad_" function)
macro grad_u(u, grad)
  return :(($u) * ($grad));
end

#! Traction's contribution to external force vector
function traction_to_extforcevec(xi::Real, eta::Real, x::Array{Float64,2},
  tau::Array{Float64,1}, unit_normal::Array{Float64,1})
  F_xi_o_t = defgrad_par_to_deformed_t(xi, eta, x); #!< elem 1
  J_xi_o = det(F_xi_o_t);
  J_xi_o_A = J_xi_o * norm(inv(F_xi_o_t) * unit_normal);
  return zeros(8,1) + N_t_Q4(xi, eta) * tau * J_xi_o_A;
end