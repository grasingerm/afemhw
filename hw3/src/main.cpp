#include <iostream>
#include <armadillo>
#include <list>
#include "element.hpp"
#include "point.hpp"

using namespace arma;
using namespace std;
using namespace pt;

/*
 * global program specific data
 */
const static double E = 100.;
const static double v = 0.45;

/* TODO: encapsulate this in an object */
mat material_tangent_stiffness(const double, const double);

int main(int argc, *char argv[])
{
    mat C_SE = material_tangent_moduli(E, v);
    
    /* surface traction */
    vec::fixed<2> tau;
    tau << 1000. << 500.;
    
    /* mesh data */
    unsigned int num_nodes;
    double x_coords { {0., 5., 5., 0.} };
    double y_coords { {0., 2., 5., 8.} };
    
    list<Node2D> node_list;
    for (unsigned int i = 0; i < num_nodes; i++)
        node_list.emplace_back(Point2D(x_coords[i], y_coords[i]), 2*i, 2*i+1);
    
    cout << "Nodes: " << endl;
    for (auto& node : node_list) cout << node.pt << endl;
    
    list<Q4> element_list;
    element_list.emplace_back(Q4(
    
    return 0;
}

/**
 * Calculate the material tangent moduli matrix
 *
 * @param youngs_mod Young's modulus
 * @param poissons_ratio Poisson's ratio
 * @return Matrix
 */
mat material_tangent_moduli
    (const double youngs_mod, const double poissons_ratio)
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
