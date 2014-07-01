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

}

#endif /* __INTEGRATION_HPP__ */
