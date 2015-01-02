# TODO: make sure res and jac return same type
# TODO: make init_guess type more general

#! Newton-Raphson solver
#!
#! \param res Function to compute residual
#! \param jac Function to compute jacobian
#! \param init_guess Initial guess
#! \param max_iters
function newton_raphson(res::Function, jac::Function, init_guess::Array,
  max_iters::Integer, eps::FloatingPoint, alpha::FloatingPoint=1.0)

  u = init_guess;
  k = 1;
  eps = 1e9;

  while k < max_iters && err > eps
    r = res();
    j = jac();
    delta_u = -r / j;
    u_star = u + delta_u;
    u = alpha*u_star + (1-alpha)*u;
    eps = norm(delta_u)/norm(u);
    k+=1;
  end

  return (u, eps, r, k)
end