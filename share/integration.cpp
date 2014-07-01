#include "integration.hpp"
#include <numeric>

using namespace std;

/**
 * Gaussian quadrature in 1D
 *
 * @param f Lambda function to integrate
 * @param num_pts Number of points in one direction
 * @return Function integrated
 */
template <class T> gquad::GaussQuad1D
    ::integrate(function<T(const double)> f, const unsigned int num_pts)
{
    array *int_pts = gquad::int_pts[num_pts];
    array *weights = gquad::weights[num_pts];
    
    T sum = *weights[0] * f(*int_pts[0]);
    for (unsigned int i = 1; i < num_pts; i++)
        sum += *weights[i] * f(*int_pts[i]);
    
    return sum;
}

/**
 * Gaussian quadrature in 2D
 *
 * @param f Lambda function to integrate
 * @param num_pts Number of points in one direction
 * @return Function integrated
 */
template <class T> gquad::GaussQuad2D
    ::integrate(function<T(const double, const double)> f, 
        const unsigned int num_pts)
{
    array *int_pts = gquad::int_pts[num_pts];
    array *weights = gquad::weights[num_pts];
    
    T sum = *weights[0] * *weights[0] * f(*int_pts[0], *int_pts[0]);
    for (unsigned int i = 1; i < num_pts; i++)
        for (unsigned int j = 0; j < num_pts; j++)
            sum += *weights[i] * *weights[j] * f(*int_pts[i], *int_pts[j]);
    
    return sum;
}

/**
 * Gaussian quadrature in 3D
 *
 * @param f Lambda function to integrate
 * @param num_pts Number of points in one direction
 * @return Function integrated
 */
template <class T> gquad::GaussQuad3D
    ::integrate(function<T(const double, const double, const double)> f, 
        const unsigned int num_pts)
{
    array *int_pts = gquad::int_pts[num_pts];
    array *weights = gquad::weights[num_pts];
    
    T sum = *weights[0] * *weights[0] * *weights[0] * 
        f(*int_pts[0], *int_pts[0], *int_pts[0]);
    for (unsigned int i = 1; i < num_pts; i++)
        for (unsigned int j = 0; j < num_pts; j++)
            for (unsigned int k = 0; k < num_pts; k++)
                sum += *weights[i] * *weights[j] * *weights[k] *
                    f(*int_pts[i], *int_pts[j], *int_pts[k]);
    
    return sum;
}
