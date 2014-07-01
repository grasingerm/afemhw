#ifndef __INTEGRATION_HPP__
#define __INTEGRATION_HPP__

#include <vector>
#include <functional>

using namespace std;

template <class T> class Integrator
{
public:
    virtual T integrate(function<T(const double)>, const double) = 0;
};

namespace gquad
{
    static const array<double,1> int_pts_1 = { 0.0 };
    static const array<double,1> weights_1 = { 2.0 };
    static const array<double,2> int_pts_2 = 
        { -0.5773502691896257, 0.5773502691896257 };
    static const array<double,2> weights_2 = { 1.0, 1.0 };
    static const array<double,3> int_pts_3 =
        { 0.0, -0.7745966692414834, 0.7745966692414834 };
    static const array<double,3> weights_3 = 
        { 0.8888888888888888, 0.5555555555555556, 0.5555555555555556 };
    static const array<double,4> int_pts_4 =
        { -0.3399810435848563, 0.3399810435848563, -0.8611363115940526, 
            0.8611363115940526 };
    static const array<double,4> weights_4 =
        { 0.6521451548625461, 0.6521451548625461, 0.3478548451374538,
            0.3478548451374538 };
    static const array<double,5> int_pts_5 =
        { 0.0, -0.5384693101056831, 0.5384693101056831, -0.9061798459386640,
            0.9061798459386640 };
    static const array<double,5> weights_5 =
        { 0.5688888888888889, 0.4786286704993665, 0.4786286704993665,
            0.2369268850561891, 0.2369268850561891 };
            
    //TODO: can't have a pointer to an array? use a bunch of vectors?
    static const vector<array*> int_pts = { &int_pts_1, &int_pts_2, &int_pts_3,
        &int_pts_4, &int_pts_5 };
    static const vector<array*> weights = { &weights_1, &weights_2, &weights_3,
        &weights_4, &weights_5 };

    template <class T> class GaussQuad1D : Integrator
    {
    public:
        GaussQuad1D() {};
        T integrate(function<T(const double)>, const unsigned int);
    };

    template <class T> class GaussQuad2D : Integrator
    {
    public:
        GaussQuad2D() {};
        T integrate(function<T(const double, const double)>, 
            const unsigned int);
    };

    template <class T> class GaussQuad3D : Integrator
    {
    public:
        GaussQuad3D() {};
        T integrate(function<T(const double, const double, const double)>, 
            const unsigned int);
    };

}

#endif /* __INTEGRATION_HPP__ */
