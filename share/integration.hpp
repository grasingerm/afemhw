#ifndef __INTEGRATION_HPP__
#define __INTEGRATION_HPP__

#include <vector>

using namespace std;

template <class T> class Integrator
{
public:
    virtual Integrator() = 0;
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
        { 0.0, 0.5555555555555556, 0.5555555555555556 };
    
    template <class T> class GaussQuad : Integrator
    {
    public:
        virtual GaussQuad() = 0;
        virtual T integrate(function<T(const double)>, const double) = 0;
    protected:
        array *int_pts;
        array *weights;
    }

    /* 1D integration */
    template <class T> class GaussQuad1D : Integrator
    {
    public:
        virtual GaussQuad1D() = 0;
        T integrate(function<T(const double)>, const double);
    };

    class GaussQuad1D1Pt : GaussQuad1D
    {
    public:
        GaussQuad1D1Pt() : int_pts(&int_pts_1), 
            weights(&weights_1) {};
    };

    class GaussQuad1D2Pts : GaussQuad1D
    {
    public:
        GaussQuad1D2Pts() : int_pts(&int_pts_2), 
            weights(&weights_2) {};
    };

    class GaussQuad1D3Pts : GaussQuad1D
    {
    public:
        GaussQuad1D3Pts() : int_pts(&int_pts_3), 
            weights(&weights_3) {};
    };
    /* end of 1D integration */

    /* 2D integration */
    template <class T> class GaussQuad2D : Integrator
    {
    public:
        virtual GaussQuad2D() = 0;
        T integrate(function<T(const double)>, const double);
    };

    class GaussQuad2D1Pt : GaussQuad2D
    {
    public:
        GaussQuad2D1Pt() : int_pts(&int_pts_1), 
            weights(&weights_1) {};
    };

    class GaussQuad2D2Pts : GaussQuad2D
    {
    public:
        GaussQuad2D2Pts() : int_pts(&int_pts_2), 
            weights(&weights_2) {};
    };

    class GaussQuad2D3Pts : GaussQuad2D
    {
    public:
        GaussQuad2D3Pts() : int_pts(&int_pts_3), 
            weights(&weights_3) {};
    };
    /* end of 2D integration */
    
    /* 3D integration */
    template <class T> class GaussQuad3D : Integrator
    {
    public:
        virtual GaussQuad3D() = 0;
        T integrate(function<T(const double)>, const double);
    };

    class GaussQuad3D1Pt : GaussQuad3D
    {
    public:
        GaussQuad3D1Pt() : int_pts(&int_pts_1), 
            weights(&weights_1) {};
    };

    class GaussQuad3D2Pts : GaussQuad3D
    {
    public:
        GaussQuad3D2Pts() : int_pts(&int_pts_2), 
            weights(&weights_2) {};
    };

    class GaussQuad3D3Pts : GaussQuad3D
    {
    public:
        GaussQuad3D3Pts() : int_pts(&int_pts_3), 
            weights(&weights_3) {};
    };
    /* end of 3D integration */

}

#endif /* __INTEGRATION_HPP__ */
