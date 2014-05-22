//#define ARMA_DONT_USE_WRAPPER // need because the wrapper caused linking errors
#include <armadillo>
#include "integration.h"

GaussQuad1D::integrate(Element1D &elem, double (*func)(const Element1D &elem))
{
    arma::colvec u(elem.nodes.size());
    for (auto node : elem.nodes)
}

