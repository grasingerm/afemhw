#include <armadillo>
#include <array>
#include <memory>
#include <tuple>
#include "element.hpp"

/*
 * TODO: write inputs and outputs for each "piece"/"unit"/function. 
            highlight redundancies and similarities to minimize calculating 
            the same value twice
 */

using namespace shape_fncts;

ShapeFunction2D::ShapeFunction2D(double (* const nxi)(const double),
    double (* const neta)(const double))
{
    ns[0] = nxi;
    ns[1] = neta;
}

const std::array<ShapeFunction2D,4> Q4::shape_fncts =
{
    {
        ShapeFunction2D(&linear_N2, &linear_N2),
        ShapeFunction2D(&linear_N1, &linear_N2),
        ShapeFunction2D(&linear_N1, &linear_N1),
        ShapeFunction2D(&linear_N2, &linear_N1)
    }
};

const std::array<double(*)(const double),4> Q4::dxi_shape_fncts =
{
    {
        &linear_dN2N2,
        &linear_dN2N1,
        &linear_dN1N1,
        &linear_dN1N2
    }
};

const std::array<double(*)(const double),4> Q4::deta_shape_fncts =
{
    {
        &linear_dN2N2,
        &linear_dN1N2,
        &linear_dN1N1,
        &linear_dN2N1
    }
};

/*
static const int e_aux_data[12] =
{
    1, 0, 0,
    0, 0, 1,
    0, 0, 1,
    0, 1, 0
};

const arma::Mat<int>::fixed<3,4> Q4::e(e_aux_data);
*/

Q4::Q4(pt::Node2D& p, pt::Node2D& q, pt::Node2D& r, pt::Node2D& s)
{
    nodes[0] = std::shared_ptr<pt::Node2D>(&p);
    nodes[1] = std::shared_ptr<pt::Node2D>(&q);
    nodes[2] = std::shared_ptr<pt::Node2D>(&r);
    nodes[3] = std::shared_ptr<pt::Node2D>(&s);
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
            J(0,k) += dxi_shape_fncts[j](eta) * (nodes[j]->pt)[k];
            
    for (unsigned int j = 0; j < 4; j++)
        for (unsigned int k = 0; k < 2; k++)
            J(1,k) += deta_shape_fncts[j](xi) * (nodes[j]->pt)[k];
    
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
    {
        dN(0,i) = dxi_shape_fncts[i](eta);
        dN(1,i) = deta_shape_fncts[i](xi);
    }
    
    return dN;
}

/**
 * TODO: consider rewrite to take s_F and s_dN_dX as input
 * B_o
 */
std::array<arma::mat,4> Q4::B_o
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

/*
 * TODO: consider taking dN_dX as input
 * B_geom
 *
 * if move semantics aren't used here this will be a resource nightmare
 */
std::array<arma::mat,4> Q4::B_geom (const double xi, const double eta)
{
    /* set up data */
    std::array<arma::mat,4> B_geom;
    arma::mat B_k;
    
    arma::mat::fixed<2,4> s_dN_dX = dN_dX(xi, eta);
    
    for (unsigned int j = 0; j < s_dN_dX.n_cols; j++)
    {
        B_k = arma::mat(4,2);
        B_k.zeros();
        
        for (unsigned int i = 0; i < s_dN_dX.n_rows; i++)
            for (unsigned int k = 0; k < s_dN_dX.n_rows; k++)
                B_k(i+2*k,0) = s_dN_dX(i,j);
                
        B_geom[j] = B_k;
    }
    
    return B_geom;
}

/**
 * Gets nodes from edge number
 */
std::array<std::shared_ptr<pt::Node2D>,2>
    Q4::nodes_from_edge(const unsigned int e)
{
    std::array<std::shared_ptr<pt::Node2D>,2> e_nodes;
    
    e_nodes[0] = nodes[e];
    
    if (e < 3)
        e_nodes[1] = nodes[e+1];
    else
        e_nodes[1] = nodes[0];
        
    return e_nodes;
}

/**
 * Given an edge nodes return unit normal to edge
 */
inline arma::vec::fixed<2> unit_normal_from_edge_nodes
    (const std::array<std::shared_ptr<pt::Node2D>,2> e_nodes)
{
    arma::vec::fixed<2> u = pt::vec_from_pts(e_nodes[0]->pt, e_nodes[1]->pt);
    return pt::unit_normal(u);
}

/* TODO: anything that isn't a mapping and is a mechanical calculation should 
 *          be contained within a 'physics'/'model'/'mechanics'/etc class */
 // THIS DOESN'T BELONG HERE. MOVE IT
/**
 * Given an edge number and traction vector return nodes and P
 */
std::tuple<
    arma::vec::fixed<4>,
    std::array<std::shared_ptr<pt::Node2D>,2>
    > 
    Q4::P_ext_from_traction
    (const arma::vec& t, const unsigned int e, const arma::mat& s_F_o_xi,
    const arma::mat& s_Nt)
{
    std::array<std::shared_ptr<pt::Node2D>,2> e_nodes = nodes_from_edge(e);
    arma::vec::fixed<2> N_un = ::unit_normal_from_edge_nodes(e_nodes);
    double J_oA_xi;
    
    J_oA_xi = arma::det(s_F_o_xi) * pt::norm_2D(s_F_o_xi.i().t() * N_un);
    arma::vec::fixed<4> integrand = s_Nt * t * J_oA_xi;
    
    return std::tuple<
        arma::vec::fixed<4>,
        std::array<std::shared_ptr<pt::Node2D>,2>
        >
        (integrand,e_nodes);
}
