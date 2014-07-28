#include <cassert>
#include <iostream>
#include <cstdlib>
#include <armadillo>
#include "element.hpp"

using namespace shape_fncts;

double norm1(arma::mat& a)
{
    double sum = 0;
    for (unsigned int i = 0; i < a.n_rows; i++)
        for (unsigned int j = 0; j < a.n_cols; j++)
            sum += fabs(a(i,j));
            
    return sum;
}

int main()
{
    int r, s;
    double sum;
    arma::mat N;
    pt::Point2D o(0,0);
    Q4 elem(o, o, o, o);
    
    std::cout << elem.N(-1,1) << std::endl;
    
    for (unsigned int iter = 0; iter < 1000; iter++)
    {
        r = rand() % 15 - 30;
        s = rand() % 15 - 30;
        sum = 0;
        
        N = elem.N(r,s);
        
        for (unsigned int i = 0; i < 8; i++)
            sum += N(0,i);
        assert(fabs(sum - 1) < 1e-3);
        
        sum = 0;
        for (unsigned int i = 0; i < 8; i++)
            sum += N(1,i);
        assert(fabs(sum - 1) < 1e-3);
        
        std::cout << iter << std::endl;
    }
    
    assert(fabs(linear_N1(1)) < 1e-5);
    assert(fabs(linear_N1(-1)-1) < 1e-5);
    assert(fabs(linear_N2(1)-1) < 1e-5);
    assert(fabs(linear_N2(-1)) < 1e-5);

    return 0;
}
