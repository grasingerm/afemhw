#ifndef _ELEMENT_H
#define _ELEMENT_H 1

//#define ARMA_DONT_USE_WRAPPER // need because the wrapper caused linking errors
#include <armadillo>
#include <array>

/* linear shape functions */
inline double linear_N1(const double xi) { return (1 - xi)/2.0; }
inline double linear_N2(const double xi) { return (1 + xi)/2.0; }

/* derivatives of linear shape functions */
const double linear_dN1dxi = -1.0/2.0;
const double linear_dN2dxi = 1.0/2.0;

/* Node classes */
class Node1D
{
public:
    double x_coord;
    unsigned int g_dof;
    Node1D(const double, const unsigned int, double (*)(const double),
        const double);
    Node1D() {}
    ~Node1D();
    double (*N)(const double);
    double dNdxi;
};

/* Element classes */
/* 1 dimensional element */
class Element1D
{
public:
    virtual ~Element1D() {}
    virtual arma::rowvec N(const double) = 0;
    virtual arma::rowvec B() = 0; // TODO: generalize this interface better
    virtual double J() = 0;
};

/* 2 node bar element */
class Bar2 : Element1D
{
public:
    Bar2(const double, const double, const unsigned int, const unsigned int);
    ~Bar2();
    std::array<Node1D, 2> nodes { {Node1D(), Node1D()} };
    arma::rowvec N(const double);
    arma::rowvec B();
    inline double J()
    {
        return nodes[0].dNdxi * nodes[0].x_coord +
            nodes[1].dNdxi * nodes[1].x_coord;
    }
};

#endif /* _ELEMENT_H */

