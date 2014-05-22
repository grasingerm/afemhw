#include <iostream>
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

    const int max_iterations = 1000;
    cout << "Maximum iterations: " << max_iterations << endl;

    /* Newton-Raphson */

    /* mesh data */
    int g_dof, num_elem, total_dof;
    double elem_size, x1, x2, bar_length = x_domain[1] - x_domain[0];

    /* solution data */
    vec u_g, delta_u_g, u_e(2), x(1);
    double J;
    rowvec N(2), x_e(2), B(2), R, R_e(2), dR(1);

    /* iteration data */
    const double eps = 1.e-5;
    double rel_error = 1.e9;
    int iteration = 0;
    double prev_norm = eps;
    double R_norm;

    /* integration data */
    const array<double, 1> weights { {2.0} };
    const array<double, 1> int_pts { {0.0} };

    for (elem_size = 0.5; elem_size >= 0.5; elem_size /= 2)
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
        total_dof = num_elem * elem_list.front().nodes.size();

        cout << "Total dof: " << total_dof << endl;

        /* resize global solution vectors to total degrees of freedom */
        u_g.set_size(total_dof);
        delta_u_g.set_size(total_dof);
        R.set_size(total_dof);

        /* enforce EBCs, once before iteration */
        u_g(0) = 0;
        u_g(total_dof-1) = 0.5;

        /* Integrate, assmeble, NR */
        while ( rel_error > eps && (iteration++) < max_iterations)
        {
            R.zeros();
            dR.zeros();

            /* **Because elements are 1D-linear, B-matrix will not change.
             * Calculate once. */
            B = elem_list.front().B();

            for (auto &curr_elem : elem_list)
            {
                /* **Only need to calc J once per element** */
                J = curr_elem.J();
                assert(J > 0);              // sanity check

                u_e(0) = u_g(curr_elem.nodes[0].g_dof);
                u_e(1) = u_g(curr_elem.nodes[1].g_dof);

                /* col vector of element x coordinates */
                x_e(0) = curr_elem.nodes[0].x_coord;
                x_e(1) = curr_elem.nodes[1].x_coord;

                R_e.zeros();
                /* TODO: encapsulate integration in an object */
                for (auto weight : weights)
                {
                    for (auto int_pt : int_pts)
                    {
                        N = curr_elem.N(int_pt);
                        x = N * x_e;

                        /* TODO: simplify, a lot of div and mult of J */
                        /* calculate R */
                        R_e += weight *
                                (   trans(N)/J * (x*x)
                                  - trans(N)/J * (B/J*u_e)
                                  - N/J*u_e * trans(B)/J * B/J*u_e  ) * J;

                        /* calculate dR */
                        dR += weight *
                                ( - trans(N)/J * B/J
                                  - trans(N)/J * B/J * B/J*u_e
                                  - N/J*u_e * trans(B)/J * B/J      ) * J;
                    }
                } /* end of integration */

                /* assemble global R */
                R(0) += R_e(curr_elem.nodes[0].g_dof);
                R(1) += R_e(curr_elem.nodes[1].g_dof);

            } /* looped over all elements */

            /* update displacements */
            delta_u_g = -R/dR;
            u_g += delta_u_g;

            /* check error */
            /* TODO: use armadillo "norm" function, solve linking errors? */
            R_norm = norm(R, 2);
            rel_error = R_norm/prev_norm;
            prev_norm = R_norm;

            /* update user with iteration details */
            cout << "iteration " << iteration << ", rel_error " << rel_error
                 << endl << "u_g = " << u_g << endl;

        } /* NR solution loop */
    }

    /* clean up and exit */
    return 0;
}

