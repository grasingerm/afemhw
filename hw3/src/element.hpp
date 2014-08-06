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

/* TODO: try to minimize/simply interface as much as possible */
class Q4
{
public:
    Q4(Node2D&, Node2D&, Node2D&, Node2D&);

    arma::mat::fixed<2,8> N(const double, const double);
    arma::mat::fixed<2,2> F_o_xi(const double, const double);
    arma::mat::fixed<2,4> dN(const double, const double);
    arma::mat::fixed<2,2> F(const double, const double, const arma::vec&);
    
    /* TODO: probably a more efficient way to write this */
    inline arma::mat::fixed<2,4> dN_dX(const double xi, const double eta)
        { return F_o_xi(xi, eta).i().t() * dN(xi, eta); } /* F^-t * dN */
    std::array<arma::mat,4> B_o
        (const double, const double, const arma::vec&);
    
        
private:
    static const std::array<Edges,4> edges;
    std::array<std::shared_ptr<pt::Point2D>,4> nodes;
};

#endif /* __ELEMENT_HPP__ */
