#include <iostream>
//#define ARMA_DONT_USE_WRAPPER // need because the wrapper caused linking errors
#include <armadillo>
#include <array>
#include <forward_list>
#include <cmath>
#include "element.h"

#include <cassert>

using namespace std;
using namespace arma;

int main(int argc, char* argv[])
{
    /* setting up model data structures */
    const array<double, 2> x_domain { {0.0, 1.0} };
    const double u__0 = 0;
    const double u__1 = 0.5;
    forward_list<Bar2> elem_list;

    cout << "Problem definition..." << endl;
    cout << "PDE: u * d^2u / dx^2 + x^2 = 0" << endl;
    cout << "Domain: (" << x_domain[0] << ", " << x_domain[1] << ")" << endl;
    cout << "EBCs: u(" << x_domain[0] << ")=" << u__0 << ", u(" << x_domain[1]
         << ")=" << u__1 << endl;
    cout << endl;

    const int max_iterations = 100;
    cout << "Maximum iterations: " << max_iterations << endl;

    /* Newton-Raphson */
    /* TODO: encapsulate NR logic in an object */

    /* mesh data */
    int g_dof, num_elem, total_dof, dof_per_elem;
    double elem_size, x1, x2, bar_length = x_domain[1] - x_domain[0];

    /* solution data */
    mat dR, dR_e(2, 2);     // TODO: fix these heuristics, define dof_per_elem
    vec u_g, delta_u_g, u_e(2), x_e(2), R, R_e(2);
    vec x(1), u(1), du_dx(1);
    double J;
    rowvec N(2), B(2);

    /* iteration data */
    const double eps = 1.e-5;
    double rel_error;
    int iteration;
    double prev_norm = eps;
    double R_norm;

    /* integration data */
    const array<double, 1> weights { {2.0} };
    const array<double, 1> int_pts { {0.0} };

    for (elem_size = 1.0; elem_size >= 0.125; elem_size /= 2)
    {
        num_elem = ceil(bar_length / elem_size);
        cout << "Element size: " << elem_size << " | Number of elements: "
             << num_elem << endl;

        /* clear element list */
        elem_list.clear();

        /* mesh the domain */
        for (x1 = x_domain[0], g_dof = 0; x1 < x_domain[1]; x1 += elem_size,
            g_dof++)
        {
            x2 = (x1 + elem_size > x_domain[1]) ? x_domain[1] : x1 + elem_size;
            elem_list.emplace_front(x1, x2, g_dof, g_dof+1);
        }
        dof_per_elem = elem_list.front().nodes.size();
        
        // subtract shared nodes then add an end node
        total_dof = num_elem * (dof_per_elem - 1) + 1;

        cout << "Total dof: " << total_dof << endl;

        /* resize global solution vectors to total degrees of freedom */
        u_g.set_size(total_dof);
        delta_u_g.set_size(total_dof);
        R.set_size(total_dof);
        dR.set_size(total_dof, total_dof);

        /* enforce EBCs, once before iteration */
        u_g(0) = 0;
        u_g(total_dof-1) = 0.5;

        /* reset iteration data */
        rel_error = 1e9;
        iteration = 0;
        
        /* Integrate, assmeble, NR */ /* TODO: primary function of NR object */
        while ( rel_error > eps && (iteration++) < max_iterations)
        {
            /* init iteration */
            iteration++;
            
            R.zeros();
            dR.zeros();

            cout << "... calculating B matrix ...";
            /* **Because elements are 1D-linear, B-matrix will not change.
             * Calculate once. */
            B = elem_list.front().B();
            cout << " ...calculated." << endl;

            for (auto &curr_elem : elem_list)
            {
                cout << "... calculating Jacobian ...";
                /* **Only need to calc J once per element** */
                J = curr_elem.J();
                assert(J > 0);              // sanity check
                cout << " ...calculated." << endl;

                u_e(0) = u_g(curr_elem.nodes[0].g_dof);
                u_e(1) = u_g(curr_elem.nodes[1].g_dof);

                /* col vector of element x coordinates */
                x_e(0) = curr_elem.nodes[0].x_coord;
                x_e(1) = curr_elem.nodes[1].x_coord;
                
                /* zero element solution vars */
                R_e.zeros();
                dR_e.zeros();
                
                /* TODO: encapsulate integration in an object */
                for (auto weight : weights)
                {
                    for (auto int_pt : int_pts)
                    {
                        /* cache repeated linear algebra */
                        N = curr_elem.N(int_pt);
                        x = N * x_e;
                        u = N * u_e;
                        du_dx = B * u_e; // B has been previously calculated
                        
                        /* debugging information
                        cout << "N" << endl << N << endl;
                        cout << "B" << endl << B << endl;
                        cout << "x" << endl << x_e << endl;
                        cout << "u" << endl << u_e << endl;
                        cout << "Nu" << endl << u << endl;
                        cout << "Bu" << endl << x << endl;
                        cout << "J" << endl << J << endl;
                        */
                        
                        cout << "... calculating R_e ...";
                        /* calculate R */
                        /* TODO: encapsulate "equation" in an object/lambda */
                        R_e += weight * 
                            (   trans(N) * (x*x)                   // n x 1
                              - trans(N) * du_dx * du_dx           // n x 1
                              - trans(B) * u * du_dx     ) * J;    // n x 1
                                
                        cout << " ...calculated." << endl;

                        cout << "... calculating dR ...";
                        /* calculate dR */ /* fix matrix mult, dR cannot be 0 */
                        dR_e += weight *
                                ( - 2 * trans(N) * B * du_dx(0)        // n x n
                                  - trans(B) * N * du_dx(0)            // n x n
                                  - trans(B) * B * u(0)         ) * J; // n x n
                        cout << " ...calculated." << endl;
                    }
                } /* end of integration */

                /* assemble global R */
                cout << "... assemble global R and dR ...";
                /* debugging information
                cout << "dR" << endl << dR << endl;
                cout << "R" << endl << R << endl;
                cout << "dR_e" << endl << dR_e << endl;
                cout << "R_e" << endl << R_e << endl;
                */
                for (int i = 0; i < dof_per_elem; i++)
                {
                    R(curr_elem.nodes[i].g_dof) += R_e(i);
                    for (int j = 0; j < dof_per_elem; j++)
                        dR(curr_elem.nodes[i].g_dof, curr_elem.nodes[j].g_dof)
                            += dR_e(i, j);
                }
                /* debugging information
                cout << "dR" << endl << dR << endl;
                cout << "R" << endl << R << endl;
                cout << " ...assembled." << endl;
                */

            } /* looped over all elements */

            /* update displacements */
            cout << "... updating displacements ...";
            arma::solve(delta_u_g, dR, -R);     // solve is faster for square
            //delta_u_g = inv(dR) * -R;
            u_g += delta_u_g;
            cout << " ...updated." << endl;
            cout << "R = " << endl << R << endl;
            cout << "dR = " << endl << dR << endl;
            cout << "delta_u_g = " << endl << delta_u_g << endl;

            /* check error */
            cout << "... calculate error ...";
            R_norm = norm(R, 2);
            rel_error = R_norm/prev_norm;
            prev_norm = R_norm;
            cout << " ...calculated." << endl;
            cout << "R_norm = " << R_norm <<endl;

            /* update user with iteration details */
            cout << "iteration " << iteration << ", rel_error " << rel_error
                 << endl << "u_g = " << endl << u_g << endl;

        } /* NR solution loop */
    }

    /* clean up and exit */
    return 0;
}

