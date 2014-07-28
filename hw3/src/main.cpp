#include <iostream>
#include <armadillo>
#include "element.hpp"
#include "point.hpp"

using namespace arma;
using namespace std;

const static double E = 100.;
const static double v = 0.45;

mat material_tangent_stiffness(const double, const double);

int main(int argc, *char argv[])
{
    
    
    return 0;
}

/**
 * Calculate the material tangent moduli matrix
 *
 * @param youngs_mod Young's modulus
 * @param poissons_ratio Poisson's ratio
 * @return Matrix
 */
mat material_tangent_moduli(const double youngs_mod, 
    const double poissons_ratio)
{
    mat C_SE;
    
    double E_star = youngs_mod/(1-poissons_ratio*poissons_ratio);
    double v_star = poissons_ratio/(1-poissons_ratio);
    
    C_SE << 1.      << v_star   << 0                << endr
         << v_star  << 1.       << 0                << endr
         << 0       << 0        << (1-v_star)/2.    << endr;
    C_SE *= E_star/(1-v_star*v_star);
    
    return C_SE;
}
