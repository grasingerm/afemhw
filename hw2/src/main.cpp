#include <iostream>
#include <armadillo>
#include <array>
#include "element.h"

using namespace std;

/**
 * Essential boundary condition
 */
class EBC
{
public:
    EBC() {}
    EBC(const unsigned int g_dof, const double val) : g_dof(g_dof), val(val) {}
    unsigned int g_dof;
    double val;
};

int main(int argc, char* argv[])
{
    const array<double, 2> x_domain { {0.0, 1.0} };
    array<EBC, 2> u_ebcs { {EBC(), EBC()} };
}

