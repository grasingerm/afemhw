#include "../integration.hpp"
#include "assertion_helper.h"
#include <iostream>
#include <cassert>
#include <cmath>
#include <functional>

using namespace std;
using namespace gquad;

int main()
{
    auto f1 = [](const double x) -> double { return x*x; };
    auto f2 = [](const double x) -> double { return 3*x+5; };
    auto f3 = [](const double x) -> double { return x; };
    auto f4 = [](const double x, const double y) -> double { return x*x + y; };
    auto f5 = [](const double x, const double y, const double z) -> double
        {
            return pow(x,5) - pow(y,3) + x*y*z - 4;
        };
    
    GaussQuad1D<double> Int1D;
    GaussQuad2D<double> Int2D;
    GaussQuad3D<double> Int3D;
    
    cout << "integration of x^2: " << Int1D.integrate(f1, 2) << " = 2/3?" << 
            endl;
    cout << "integration of 3x+5: " << Int1D.integrate(f2, 1) << " = 10?" << 
            endl;
    cout << "integration of x: " << Int1D.integrate(f3, 1) << " = 0?" << endl;
    cout << "integration of x^2+y: " << Int2D.integrate(f4, 2) << " = 4/3?" << 
            endl;
    cout << "integration of x^5 - y^3 + xyz -4: " << Int3D.integrate(f5, 3) << 
            " = -32?" << endl;
    
    ASSERT_NEAR(2.0/3.0, Int1D.integrate(f1, 2), 1e-5);
    ASSERT_NEAR(10.0, Int1D.integrate(f2, 1), 1e-5);
    ASSERT_NEAR(0.0, Int1D.integrate(f3, 1), 1e-5);
    ASSERT_NEAR(4.0/3.0, Int2D.integrate(f4, 2), 1e-5);
    ASSERT_NEAR(-32, Int3D.integrate(f5, 3), 1e-5);
    
    return 0;
}
