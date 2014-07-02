#ifndef __INTEGRATION_HPP__
#define __INTEGRATION_HPP__

#include <vector>
#include <functional>

using namespace std;

//TODO: are these interfaces necessary?
/* interfaces for 1D, 2D, and 3D integration */
template <class T> class Integrator1D
{
public:
    virtual T integrate(function<T(const double)>, const unsigned int) = 0;
};

template <class T> class Integrator2D
{
public:
    virtual T integrate(function<T(const double, const double)>, 
        const unsigned int) = 0;
};

template <class T> class Integrator3D
{
public:
    virtual T integrate(function<T(const double, const double, const double)>, 
        const unsigned int) = 0;
};

/* gauss quadrature */
namespace gquad
{
    /* set up integration points and weights */
    static double int_pts_1[] = { 0.0 };
    static double weights_1[] = { 2.0 };
    static double int_pts_2[] = 
        { -0.5773502691896257, 0.5773502691896257 };
    static double weights_2[] = { 1.0, 1.0 };
    static double int_pts_3[] =
        { 0.0, -0.7745966692414834, 0.7745966692414834 };
    static double weights_3[] = 
        { 0.8888888888888888, 0.5555555555555556, 0.5555555555555556 };
    static double int_pts_4[] =
        { -0.3399810435848563, 0.3399810435848563, -0.8611363115940526, 
            0.8611363115940526 };
    static double weights_4[] =
        { 0.6521451548625461, 0.6521451548625461, 0.3478548451374538,
            0.3478548451374538 };
    static double int_pts_5[] =
        { 0.0, -0.5384693101056831, 0.5384693101056831, -0.9061798459386640,
            0.9061798459386640 };
    static double weights_5[] =
        { 0.5688888888888889, 0.4786286704993665, 0.4786286704993665,
            0.2369268850561891, 0.2369268850561891 };
    static double int_pts_6[] = { 0.6612093864662645, -0.6612093864662645,
        -0.2386191860831969, 0.2386191860831969, -0.9324695142031521,
        0.9324695142031521 };
    static double weights_6[] = { 0.3607615730481386, 0.3607615730481386,
        0.4679139345726910, 0.4679139345726910, 0.1713244923791704,
        0.1713244923791704 };
            
    static double *int_pts[] = { int_pts_1, int_pts_2, int_pts_3,
        int_pts_4, int_pts_5, int_pts_6 };
    static double *weights[] = { weights_1, weights_2, weights_3,
        weights_4, weights_5, weights_6 };

    /* write integration templates */
    template <class T> class GaussQuad1D : public Integrator1D<T>
    {
    public:
        GaussQuad1D() {};
        T integrate(function<T(const double)>, const unsigned int);
    };

    template <class T> class GaussQuad2D : public Integrator2D<T>
    {
    public:
        GaussQuad2D() {};
        T integrate(function<T(const double, const double)>, 
            const unsigned int);
    };

    template <class T> class GaussQuad3D : public Integrator3D<T>
    {
    public:
        GaussQuad3D() {};
        T integrate(function<T(const double, const double, const double)>, 
            const unsigned int);
    };
    
    
    /**
     * Gaussian quadrature in 1D
     *
     * @param f Lambda function to integrate
     * @param num_pts Number of points in one direction
     * @return Function integrated
     */
    template <class T> T gquad::GaussQuad1D<T>
        ::integrate(function<T(const double)> f, const unsigned int num_pts)
    {
        double *int_pts = gquad::int_pts[num_pts-1];
        double *weights = gquad::weights[num_pts-1];
        
        T sum = *weights * f(*int_pts);
        for (unsigned int i = 1; i < num_pts; i++)
            sum += *(weights+i) * f(*(int_pts+i));
        
        return sum;
    }

    /**
     * Gaussian quadrature in 2D
     *
     * @param f Lambda function to integrate
     * @param num_pts Number of points in one direction
     * @return Function integrated
     */
    template <class T> T gquad::GaussQuad2D<T>
        ::integrate(function<T(const double, const double)> f, 
            const unsigned int num_pts)
    {
        double *int_pts = gquad::int_pts[num_pts-1];
        double *weights = gquad::weights[num_pts-1];
        
        // TODO: is there a cleaner way to do this?
        // TODO: will this translate well to finite element implementations?
        // init by computing first row
        T sum = *weights * *weights * f(*int_pts, *int_pts);
        for (unsigned int j = 1; j < num_pts; j++)
            sum += *(weights) * *(weights+j) * f(*(int_pts), *(int_pts+j));
        
        // finish summation
        for (unsigned int i = 1; i < num_pts; i++)
            for (unsigned int j = 0; j < num_pts; j++)
                sum += *(weights+i) * *(weights+j) * 
                    f(*(int_pts+i), *(int_pts+j));
        
        return sum;
    }

    /**
     * Gaussian quadrature in 3D
     *
     * @param f Lambda function to integrate
     * @param num_pts Number of points in one direction
     * @return Function integrated
     */
    template <class T> T gquad::GaussQuad3D<T>
        ::integrate(function<T(const double, const double, const double)> f, 
            const unsigned int num_pts)
    {
        double *int_pts = gquad::int_pts[num_pts-1];
        double *weights = gquad::weights[num_pts-1];
        
        // init by computing first row
        // TODO: can this be cleaner? more efficient?
        T sum = *weights * *weights * *weights * 
            f(*int_pts, *int_pts, *int_pts);
        for (unsigned int k = 1; k < num_pts; k++)
            sum += *(weights) * *(weights) * *(weights+k) *
                        f(*(int_pts), *(int_pts), *(int_pts+k));
        for (unsigned int j = 1; j < num_pts; j++)
            for (unsigned int k = 0; k < num_pts; k++)
                sum += *(weights) * *(weights+j) * *(weights+k) *
                        f(*(int_pts), *(int_pts+j), *(int_pts+k));
        
        // finish summation
        for (unsigned int i = 1; i < num_pts; i++)
            for (unsigned int j = 0; j < num_pts; j++)
                for (unsigned int k = 0; k < num_pts; k++)
                    sum += *(weights+i) * *(weights+j) * *(weights+k) *
                        f(*(int_pts+i), *(int_pts+j), *(int_pts+k));
        
        return sum;
    }

}

#endif /* __INTEGRATION_HPP__ */
