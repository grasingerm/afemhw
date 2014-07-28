#include <cassert>
#include <cmath>
#include <iostream>
#include <iomanip>

#include "point.hpp"
#include "test/debug.hpp"

using namespace pt;
using namespace std;

int main()
{
    setprecision(9);

    Point2D p(1.1, 2.2);
    Point2D q(-1.1, -2.2);
    Point2D o(0.0, 0.0);
    
    cout << "Diagnostics printing..." << endl;
    cout << p << " " << q << " " << o << endl;
    cout << p[X] << " " << p[Y] << endl;
    cout << 1e-4 << endl;
    cout << "...finished." << endl;
    
    Point2D a = p+q;
    
    cout << "Start tests..." << endl;
    
    assert(a == a);
    cout << p+q << " == " << o << "?" << endl;
    cout << fabs(p[0]+q[0]) << "," << fabs(p[1]+q[1]) << endl;
    assert((p+q) == o);
    assert(q+p == o);
    assert(o+o == o);
    cout << p+o << " == " << p << "?" << endl;
    assert(p+o == p);
    assert(o+p == p);
    assert(q+o == q);
    assert((p-q) == (2*p));
    assert(q-p == 2*q);
    assert(q-p == -2*q);
    
    cout << "... tests passed." << endl;

    return 0;
}
