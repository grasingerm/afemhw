#ifndef _INTEGRATION_H
#define _INTEGRATION_H 1

//#define ARMA_DONT_USE_WRAPPER // need because the wrapper caused linking errors
#include <armadillo>

/**
 * Handler for numerical integration
 */
class Integrator
{
public:
    virtual ~Integrator() {}
    virtual arma::colvec integrate(Element1D&, double (*)(const Element1D&)) = 0;
};

/**
 * Handler for integration by Gauss Quadrature
 */
class GaussQuad1D : Integrator
{
private:
    const int weight = 2;
    const int int_pt = 0;
public:
    GaussQuad1D() {};
    arma::colvec integrate(Element1D&, double (*)(const Element1D&));
};

#endif /* _INTEGRATION_H */

