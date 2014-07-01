#include <armadillo>
#include "integration.h"

GaussQuad1DLinear::integrate(Element1D &elem, double (*func)(const Element1D&,
    const double))
{
    arma::colvec u(elem.nodes.size());
    return weight * func(elem, x);
}

