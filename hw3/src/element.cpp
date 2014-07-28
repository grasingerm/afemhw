#include <armadillo>
#include <array>
#include "element.hpp"

using namespace shape_fncts;

ShapeFunction2D::ShapeFunction2D(double (* const nxi)(const double),
    double (* const neta)(const double))
{
    ns[0] = nxi;
    ns[1] = neta;
}

const std::array<ShapeFunction2D,4> Q4::shape_fncts =
{
    ShapeFunction2D(&linear_N2, &linear_N2),
    ShapeFunction2D(&linear_N1, &linear_N2),
    ShapeFunction2D(&linear_N1, &linear_N1),
    ShapeFunction2D(&linear_N2, &linear_N1)
};

const std::array<double(*)(const double),4> Q4::dxi_shape_fncts =
{
    &linear_dN2N2,
    &linear_dN2N1,
    &linear_dN1N1,
    &linear_dN1N2
};

const std::array<double(*)(const double),4> Q4::deta_shape_fncts =
{
    &linear_dN2N2,
    &linear_dN1N2,
    &linear_dN1N1,
    &linear_dN2N1
};

static const int e_aux_data[12] =
{
    1, 0, 0,
    0, 0, 1,
    0, 0, 1,
    0, 1, 0
};

const arma::Mat<int>::fixed<3,4> Q4::e(e_aux_data);

Q4::Q4(pt::Point2D& p, pt::Point2D& q, pt::Point2D& r, pt::Point2D& s)
{
    nodes[0] = std::shared_ptr<pt::Point2D>(&p);
    nodes[1] = std::shared_ptr<pt::Point2D>(&q);
    nodes[2] = std::shared_ptr<pt::Point2D>(&r);
    nodes[3] = std::shared_ptr<pt::Point2D>(&s);
}

arma::mat Q4::N(const double xi, const double eta)
{
    arma::mat N(2,8);
    for (unsigned int i = 0; i < 4; i++)
    {
        N(0,2*i) = shape_fncts[i](xi,eta);
        N(0,2*i+1) = 0;
        N(1,2*i) = 0;
        N(1,2*i+1) = shape_fncts[i](xi,eta); 
    }
    
    return N;
}

arma::mat Q4::J(const double xi, const double eta)
{
    arma::mat J(2,2);
    J.zeros();
    
    for (unsigned int j = 0; j < 4; j++)
        for (unsigned int k = 0; k < 2; k++)
            J(0,k) += dxi_shape_fncts[j](eta) * (*nodes[j])[k];
            
    for (unsigned int j = 0; j < 4; j++)
        for (unsigned int k = 0; k < 2; k++)
            J(1,k) += deta_shape_fncts[j](xi) * (*nodes[j])[k];
    
    return J;
}
