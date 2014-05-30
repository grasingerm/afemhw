//#define ARMA_DONT_USE_WRAPPER // need because the wrapper caused linking errors
#include <armadillo>
#include "element.h"

//#define NDEBUG 1
#include <cassert>

/**
 * Constructor for 1 dimensional node
 *
 * @param x_coord X coordinate of the node
 * @param g_dof Global degree of freedom of the nodal displacement
 * @param N Address of shape function
 * @param dNdxi Derivative of shape function with respect to xi
 * @return Node1D
 */
Node1D::Node1D(const double x_coord, const unsigned int g_dof,
    double (*N)(const double), const double dNdxi) : x_coord(x_coord),
    g_dof(g_dof), N(N), dNdxi(dNdxi) {}

/**
 * Destructor for 1 dimensional node
 */
Node1D::~Node1D() {}

/**
 * Constructor for 2 node bar element
 *
 * @param x1 X coordinate of first node
 * @param x2 X coordinate of second node
 * @param g_dof1 Global degree of freedom of nodal displacement 1
 * @param g_dof2 Global degree of freedom of nodal displacement 2
 * @return Bar2
 */
Bar2::Bar2(const double x1, const double x2, const unsigned int g_dof1,
    const unsigned int g_dof2)
{
    nodes[0].x_coord = x1;
    nodes[0].g_dof = g_dof1;
    nodes[0].N = &linear_N1;
    nodes[0].dNdxi = linear_dN1dxi;

    nodes[1].x_coord = x2;
    nodes[1].g_dof = g_dof2;
    nodes[1].N = &linear_N2;
    nodes[1].dNdxi = linear_dN2dxi;
}

/**
 * Destructor for 2 node bar
 */
Bar2::~Bar2() {}

/**
 * Matrix of shape functions
 *
 * @param xi Parent element xi coordinate
 * @return N matrix
 */
arma::rowvec Bar2::N(const double xi)
{
    arma::rowvec::fixed<2> N;
    N << nodes[0].N(xi) << nodes[1].N(xi);
    return N;
}

/**
 * Matrix of shape function derivatives
 *
 * @param xi Parent element xi coordinate
 * @return B matrix
 */
arma::rowvec Bar2::B()
{
    arma::rowvec::fixed<2> B;
    B << nodes[0].dNdxi << nodes[1].dNdxi;
    return B/J();
}

