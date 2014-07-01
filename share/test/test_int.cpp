#include "../integration.hpp"
#include "assertion_helper.h"
#include <iostream>
#include <cassert>

using namespace std;

int main()
{
    auto f1 = [](const double x) -> double { return x*x; };
    auto f2 = [](const double x) -> double { return 3*x+5; };
    auto f3 = [](const double x) -> double { return x; };
    auto f4 = [](const double x, const double y) -> double { return x*x + y; };
    
    GaussQuad1D<double> Int1D();
    GaussQuad2D<double> Int2D();
    
    ASSERT_NEAR(2.0/3.0, Int1D.integrate(f1, 2), 1e-5);
    ASSERT_NEAR(10.0, Int1D.integrate(f2, 1), 1e-5);
    ASSERT_NEAR(0.0, Int1D.integrate(f3, 1), 1e-5);
    ASSERT_NEAR(4.0/3.0, Int2D.integrate(f4, 2), 1e-5);
    
    return 0;
}
