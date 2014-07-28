#ifndef __ELEMENT_HPP__
#define __ELEMENT_HPP__

#include <array>
#include <tuple>
#include <memory>
#include <armadillo>
#include "point.hpp"

/* shape function definitions */
namespace shape_fncts
{
    inline double linear_N1(const double x) { return (1-x)/2.; }
    inline double linear_N2(const double x) { return (1+x)/2.; }
    
    const double linear_dN1 = -0.5;
    const double linear_dN2 = 0.5;
    
    inline double linear_dN1N1(const double x) { return linear_dN1*linear_N1(x); }
    inline double linear_dN2N1(const double x) { return linear_dN2*linear_N1(x); }
    inline double linear_dN1N2(const double x) { return linear_dN1*linear_N2(x); }
    inline double linear_dN2N2(const double x) { return linear_dN2*linear_N2(x); }

    class ShapeFunction2D
    {
    public:
        ShapeFunction2D(double(*const)(const double), 
            double(*const)(const double));
        inline double operator()(const double xi, const double eta) const
            { return ns[0](xi) * ns[1](eta); }
    private:
        std::array<double(*)(const double),2> ns;
    };

}

class Q4
{
public:
    Q4(pt::Point2D&, pt::Point2D&, pt::Point2D&, pt::Point2D&);

    arma::mat N(const double, const double);
    arma::mat B(const double, const double);
    arma::mat J(const double, const double);
private:
    static const std::array<shape_fncts::ShapeFunction2D,4> shape_fncts;
    static const std::array<double(*)(const double),4> dxi_shape_fncts;
    static const std::array<double(*)(const double),4> deta_shape_fncts;
    static const arma::Mat<int>::fixed<3,4> e;
    std::array<std::shared_ptr<pt::Point2D>,4> nodes;
    std::array<double,8> g_dof;
};

#endif /* __ELEMENT_HPP__ */
