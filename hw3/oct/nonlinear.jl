# TODO: make sure res and jac return same type
# TODO: make init_guess type more general

#! Newton-Raphson solver
#!
#! \param res Function to compute residual
#! \param jac Function to compute jacobian
#! \param fin "Finally" function to clean up before next iteration
#! \param init_guess Initial guess
#! \param max_iters Maximum number of iterations to perform
#! \param eps Epsilon, error tolerance convergence criteria
#! \param alpha Relaxation factor
#! \return tuple(answer, converged?, error, residual, number of iterations)
function newton_raphson(res::Function, jac::Function, fin::Function,
  init_guess::Array, max_iters::Integer, eps::FloatingPoint, 
  alpha::FloatingPoint=1.0)

  @assert(max_iters > 0, "Maximum number of iterations should be positive.");
  @assert(eps >= 0, "Error tolerance must be greater than or equal to zero.");

  u = init_guess;
  k = 1;
  err = 1e9;

  while k < max_iters && err > eps
    r = res();
    j = jac();
    fin();
    delta_u = -r / j;
    u_star = u + delta_u;
    u = alpha*u_star + (1-alpha)*u;
    err = norm(delta_u)/norm(u);
    k+=1;
  end

  converged = err <= eps;

  return (u, converged, err, r, k);
end
