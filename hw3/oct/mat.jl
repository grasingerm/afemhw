require("debug.jl");

#! Tanget stiffness moduli for plan strain condition using Kirchoff material
#! constitutive model
#!
#! \param Ey Material elastic modulus
#! \param nu Poisson's ratio
#! \return Tangent stiffness moduli matrix
function tan_stiff_moduli_plane_strain(Ey::Float64, nu::Float64)
  return Ey / (1-nu^2) * [1 nu 0; nu 1 0; 0 0 (1-nu)/2.0];
end
