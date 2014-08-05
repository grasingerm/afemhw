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

arma::mat::fixed<2,8> Q4::N(const double xi, const double eta)
{
    arma::mat::fixed<2,8> N;
    for (unsigned int i = 0; i < 4; i++)
    {
        N(0,2*i) = shape_fncts[i](xi,eta);
        N(0,2*i+1) = 0;
        N(1,2*i) = 0;
        N(1,2*i+1) = shape_fncts[i](xi,eta); 
    }
    
    return N;
}

/**
 * Map undeformed coordinates to parent coordinates
 */
arma::mat::fixed<2,2> Q4::F_o_xi(const double xi, const double eta)
{
    arma::mat::fixed<2,2> J;
    J.zeros();
    
    for (unsigned int j = 0; j < 4; j++)
        for (unsigned int k = 0; k < 2; k++)
            J(0,k) += dxi_shape_fncts[j](eta) * (*nodes[j])[k];
            
    for (unsigned int j = 0; j < 4; j++)
        for (unsigned int k = 0; k < 2; k++)
            J(1,k) += deta_shape_fncts[j](xi) * (*nodes[j])[k];
    
    return J;
}

/**
 * Map deformed coordinates to undeformed coordinates
 */
arma::mat::fixed<2,2> Q4::F
    (const double xi, const double eta, const arma::vec& u)
{
    arma::mat::fixed<2,2> F;
    F.zeros();
    
    for (unsigned int j = 0; j < 4; j++)
        for (unsigned int k = 0; k < 2; k++)
            F(0,k) += dxi_shape_fncts[j](eta) * u(2*j+k);
            
    for (unsigned int j = 0; j < 4; j++)
        for (unsigned int k = 0; k < 2; k++)
            F(1,k) += deta_shape_fncts[j](xi) * u(2*j+k);
    
    /* F = grad(u) + I */
    for (unsigned int i = 0; i < 2; i++) F(i,i) += 1;
    
    return F;
}

/**
 * dN matrix
 */
arma::mat::fixed<2,4> Q4::dN(const double xi, const double eta)
{
    arma::mat::fixed<2,4> dN;
    
    for (unsigned int i = 0; i < 4; i++)
        dN << dxi_shape_fncts[i](eta);
    dN << arma::endr;
    
    for (unsigned int i = 0; i < 4; i++)
        dN << deta_shape_fncts[i](xi);
    dN << arma::endr;
    
    return dN;
}

/**
 * B_o
 */
std::array<arma::mat,4> B_o
    (const double xi, const double eta, const arma::vec& u)
{
    /* set up data */
    std::array<arma::mat,4> B_o;
    arma::mat B_k;
    
    arma::mat::fixed<2,2> s_F = F(xi, eta, u);
    arma::mat::fixed<2,4> s_dN_dX = dN_dX(xi, eta);
    
    for (unsigned int i = 0; i < 4; i++)
    {
        B_k = arma::mat(3,2);
        for (unsigned int j = 0; j < 2; j++)
            for (unsigned int k = 0; k < 2; k++)
                B_k(j,k) = s_F(k,j) * s_dN_dX(j,i);
        for (unsigned int j = 0; j < 2; j++)
            B_k(2,j) = s_F(j,0)*s_dN_dX(1,i) + s_F(j,1)*s_dN_dX(0,i);
            
        B_o[i] = B_k;
    }
    
    return B_o;
}
