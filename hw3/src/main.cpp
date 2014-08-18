#include <cassert>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <fstream>
#include <armadillo>
#include <list>
#include <vector>
#include "element.hpp"
#include "point.hpp"
#include "integration.hpp"

using namespace arma;
using namespace std;
using namespace pt;

/* TODO: implement something like this........
 * class FEMSolver
   {
   public:
        FEMSolver(list<Node> node_list, list<Elem> elem_list,
            list<ConstitutiveRelationship> cr_list, list<ExtLoad> ext_load_list)
            : blah blah blah {}
        vec solve(vec initial_guess, const double tolerance, const blah blah)
   };
 */

/*
 * global program specific data
 */
const static double E = 100.;
const static double v = 0.45;

/* TODO: encapsulate this in an object */
mat material_tangent_moduli(const double, const double);
vec solve_nonlinear_static_by_TL
    (const list<Q4> elem_list, vec u_guess, const double delta_load_factor, 
    const double tolerance, const unsigned int max_iter);
double norm_L2(const vec&);    
    

int main(int argc, char* argv[])
{
    ofstream ofile;
    ofile.open("hw3.txt");

    mat C_SE = material_tangent_moduli(E, v);
    
    /* mesh data */
    unsigned int num_nodes = 4;
    double x_coords[] = {0., 5., 5., 0.};
    double y_coords[] = {0., 2., 5., 8.};
    
    const double tol = 1e-4;
    const unsigned int max_iter = 1e5;
    const double delta_load_factor = 0.05;
    
    vector<Node2D> node_list;
    for (unsigned int i = 0; i < num_nodes; i++)
        node_list.emplace_back(Point2D(x_coords[i], y_coords[i]), 2*i, 2*i+1);
    
    ofile << "Nodes: " << endl;
    for (auto& node : node_list) ofile << node.pt << endl;
   
    list<Q4> element_list;
    element_list.emplace_back(Q4(node_list[0], node_list[1], node_list[2],
        node_list[3]));
        
    /* surface traction */
    vec::fixed<2> tau;
    tau << 1000. << 500.;

    /* code to solve static problem */
    const static array<double,2> int_pts 
        { { -0.5773502691896257, 0.5773502691896257 } };
    const static array<double,2> weights { { 1., 1. } };
    
    /* TODO:pass by reference, or stay static? */
    // auto fn_ptr = &compute_residual_and_jacobian;
    
    vec u(8), delta_u(8);
    u.zeros();

    unsigned int ndofs = u.n_rows;
    double lf = 0.;
    
    vec P_o(ndofs);
    P_o.zeros();
    
    vec I_f(ndofs);
    I_f.zeros();
    
    mat A(ndofs, ndofs);
    A.zeros();
    
    double xi, eta;
    
    array<mat,4> B_o;
    array<mat,4> B_geom;
    
    mat::fixed<2,2> s_F;
    mat::fixed<2,2> s_F_o_xi;
    mat::fixed<2,8> s_N;
    
    tuple<vec::fixed<2>, array<shared_ptr<Node2D>,2>> tresult;
    mat::fixed<4,2> s_Nt;
    
    ofile << "calculate P_o ... ";
    /* calculate P_o */
    P_o.zeros();
    for (auto& elem : element_list)
    {
        double weight = 2;
        double xi = 0.;
        double eta = 1.;
        
        /* calculate mappings */
        s_F_o_xi = elem.F_o_xi(xi, eta);
        s_N = elem.N(xi, eta);
        
        /* TODO: fix heuristic, slice by "edge number" */
        s_Nt = s_N.cols(0,3).t();
        
        /* assemble P_ext */
        vec::fixed<4> p_e;
        array<shared_ptr<Node2D>,2> nodes;
        tie(p_e, nodes) = elem.P_ext_from_traction
                    (tau, 1, s_F_o_xi, s_Nt);
        
        ofile << "p_e = " << p_e << endl;
        ofile << "xi = " << xi << endl;
        ofile << "eta = " << eta << endl;
        
        for (unsigned int i = 0; i < nodes.size(); i++)
        {
            auto& node_ptr = nodes[i];
            for (unsigned int j = 0; j < 2; j++)
                P_o(node_ptr->gdofs[j]) += 
                    weight*p_e(2*i+j);
        }
    }
    
    ofile << " ... done" << endl;
    ofile << "P_o = " << endl << P_o << endl;
    
    vec I_k(8);
    I_k.zeros();
    
    ofile << "starting iterative solution ..." << endl;
    unsigned int total_iterations = 0;
    while ( (lf += delta_load_factor) <= 1. )
    {
        ofile << "load factor: " << lf << endl;
        unsigned int iter = 0;
        double err = 1e9;
        
        vec R_o, P_ext = P_o*lf;
        while (++iter < max_iter && err > tol)
        {
            total_iterations++;
            ofile << "iter: " << iter << " total_iter: " << total_iterations << endl;
        
            /* integrate over elements */
            /* compute residual, and jacobian */
            for (auto& elem : element_list)
            {
                for (unsigned int i = 0; i < 2; i++)
                {
                    for (unsigned int j = 0; j < 2; j++)
                    {
                        xi = int_pts[i];
                        eta = int_pts[j];
                    
                        /* calculate mappings */
                        s_F = elem.F(xi, eta, u);
                        s_F_o_xi = elem.F_o_xi(xi, eta);
                        double s_J_o_xi = det(s_F_o_xi);
                        s_N = elem.N(xi, eta);
                        
                        B_o = elem.B_o(xi, eta, u);
                        B_geom = elem.B_geom(xi, eta);
                        
                        /* TODO: start to think about how this can be
                         *          better encapsulated/generalized */
                        // TODO: clean this up
                        /* compute second PK stress tensor */
                        mat I(2,2);
                        I.eye();
                        mat E_green = s_F * trans(s_F) - I;
                        vec Ev_green(3);
                        Ev_green << E_green(0,0) << E_green(1,1) << E_green(1,0);
                        vec S = C_SE * Ev_green;
                        mat S_bar(4,4);
                        S_bar.zeros();
                        
                        S_bar(0,0) = S(0);
                        S_bar(0,1) = S(2);
                        S_bar(1,0) = S(2);
                        S_bar(1,1) = S(1);
                        S_bar(2,2) = S(0);
                        S_bar(2,3) = S(2);
                        S_bar(3,2) = S(2);
                        S_bar(3,3) = S(1);
                        
                        ofile << "F " << endl << s_F << endl;
                        ofile << "S " << endl << S << endl;
                        ofile << "S_bar " << endl << S_bar << endl;
                        ofile << "E " << endl << E_green << endl;
                        ofile << "Evgt " << endl << Ev_green << endl;
                        
                        /* material and geometric tangent stiffness */
                        for (unsigned int k = 0; k < elem.n_nodes(); k++)
                        {
                            mat A_e = B_o[k].t() * C_SE * B_o[k];
                            A_e += B_geom[k].t() * S_bar * B_geom[k];
                            
                            /* assemble */
                            for (unsigned int i = 0; i < 2; i++)
                                for (unsigned int j = 0; j < 2; j++)
                                    A(elem.gdof(k,i), elem.gdof(k,j)) 
                                        += weights[i]*weights[j] * A_e(i,j);
                                        
                            vec I_e = B_o[k].t() * S * s_J_o_xi;
                            ofile << "I_e " << endl << I_e << endl;
                            
                            /* assemble */
                            for (unsigned int i = 0; i < 2; i++)
                                I_k(elem.gdof(k,i)) += 
                                    weights[i]*weights[j]*I_e(i);
                        }
                    }
                }
            } /* end loop over elements/assembly */
            
            /* compute residual */
            vec R = I_k - P_ext;
            ofile << "I = " << endl << I_k << endl;
            ofile << "R = " << endl << R << endl;
            
            /* enforce homogeneous essential boundary conditions */
            //TODO: fix heuristic/generalize
            u(1) = 0;
            u(7) = 0;

            /* compute correction */
            assert(solve(delta_u, A, -R));
            
            /* update solution */
            u += delta_u;
            
            ofile << "A = " << endl << A << endl
                 << "delta_u = " << endl << delta_u << endl;
            
            /* compute error */
            if (iter == 1)
                R_o = R;
                
            err = norm_L2(R) / norm_L2(R_o);
            
            ofile << "err = " << err << endl;
        }
        ofile << "===========================================" << endl << endl;
        
        /* make updates */
        I_k = P_ext;
        
    }
    
    ofile << "solution:" << endl
         << "u = " << endl << u << endl
         << "in " << total_iterations << " iterations." << endl;
    
    ofile.close();
    return 0;
}

/**
 * Calculate the material tangent moduli matrix
 *
 * @param youngs_mod Young's modulus
 * @param poissons_ratio Poisson's ratio
 * @return Matrix
 */
mat material_tangent_moduli
    (const double youngs_mod, const double poissons_ratio)
{
    mat C_SE;
    
    double E_star = youngs_mod/(1-poissons_ratio*poissons_ratio);
    double v_star = poissons_ratio/(1-poissons_ratio);
    
    C_SE << 1.      << v_star   << 0                << endr
         << v_star  << 1.       << 0                << endr
         << 0       << 0        << (1-v_star)/2.    << endr;
    C_SE *= E_star/(1-v_star*v_star);
    
    return C_SE;
}

/*
 *
 */
double norm_L2(const vec& u)
{
    double sum = 0;
    for (unsigned int i = 0; i < u.n_rows; i++) sum += u(i)*u(i);
    return sqrt(sum);
}
