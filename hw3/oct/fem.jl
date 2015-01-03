@N1_Q4(xi, eta) ((1-(xi))*(1-(eta))/4.)
@N2_Q4(xi, eta) ((1+(xi))*(1-(eta))/4.)
@N3_Q4(xi, eta) ((1+(xi))*(1+(eta))/4.)
@N4_Q4(xi, eta) ((1-(xi))*(1+(eta))/4.)

@dN1dxi_Q4(eta) (-(eta)/4.)
@dN2dxi_Q4(eta) (-(eta)/4.)
@dN3dxi_Q4(eta) ((eta)/4.)
@dN4dxi_Q4(eta) ((eta)/4.)

@dN1deta_Q4(xi) (-(xi)/4.)
@dN2deta_Q4(xi) ((xi)/4.)
@dN3deta_Q4(xi) ((xi)/4.)
@dN4deta_Q4(xi) (-(xi)/4.)

#! N matrix for a Q4 element
#!
#! \param xi Parent coordinate in x-direction
#! \param eta Parent coordinate in y-direction
#! \return N matrix
function N_Q4(xi::Real, eta::Real)
  Float64[
    N1_Q4(xi, eta) 0 N2_Q4(xi, eta) 0 N3_Q4(xi, eta) 0 N4_Q4(xi, eta) 0;
    0 N1_Q4(xi, eta) 0 N2_Q4(xi, eta) 0 N3_Q4(xi, eta) 0 N4_Q4(xi, eta);
  ];
end

#! Map global referential coordinates to parent referential coordinates
#!
#! \param N N matrix
#! \param X Nodal global referential coordinates
#! \return Referential global coords mapped to referential parent coords
function ref_glob_to_par(N::Array{Float64,2}, X::Array{Float64,2})
  return N*X;
end

#! Deformation gradient to map parent spatial coords to deformed spatial coords
#!
#! \param xi Parent coordinate in x-direction
#! \param eta Parent coordinate in y-direction
#! \param x Two dimensional matrix of deformed, spatial, nodal coordinates
#! \return Deformation gradient
function defgrad_par_to_deformed(xi::Real, eta::Real, x::Array{Float64,2})
  dN_dxi = Float64[
    dN1_dxi(eta) dN2_dxi(eta) dN3_dxi(eta) dN4_dxi(eta);
    dN1_deta(xi) dN2_deta(xi) dN3_deta(xi) dN4_deta(xi);
    ];

  F_t = dN_dxi * x;
  return F_t';
end

