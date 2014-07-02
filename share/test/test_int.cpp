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
    auto f6 = [](const double x, const double y) -> double 
    { return pow(x,8) * pow(y,4) + 4*x*y; };
    auto f7 = [](const double x, const double y) -> double { return x/y + 7; };
    auto f8 = [](const double x, const double y, const double z) -> double
    {
        return x*y+(z*z);
    };
    auto f9 = [](const double x, const double y, const double z) -> double
    {
        return x+y*z;
    };
    
    GaussQuad1D<double> Int1D;
    GaussQuad2D<double> Int2D;
    GaussQuad3D<double> Int3D;
    
    cout << "running tests ... ";
    ASSERT_NEAR(2.0/3.0, Int1D.integrate(f1, 2), 1e-5);
    ASSERT_NEAR(10.0, Int1D.integrate(f2, 1), 1e-5);
    ASSERT_NEAR(0.0, Int1D.integrate(f3, 1), 1e-5);
    ASSERT_NEAR(4.0/3.0, Int2D.integrate(f4, 2), 1e-5);
    ASSERT_NEAR(-32, Int3D.integrate(f5, 3), 1e-5);
    ASSERT_NEAR(4.0/45.0, Int2D.integrate(f6, 5), 1e-5);
    ASSERT_NEAR(28.0, Int2D.integrate(f7, 2), 1e-5);
    ASSERT_NEAR(8.0/3.0, Int3D.integrate(f8, 2), 1e-5);
    ASSERT_NEAR(0, Int3D.integrate(f9, 1), 1e-5);
    cout << "... tests passed." << endl;
    
    return 0;
}
