#include <cassert>
#include <armadillo>
#include "test_helpers.h"

using namespace arma;

int main(int argc, char* argv[])
{
    printf("Setting up data...\n");

    /* testing vector norms */
    printf("Testing vector norms...");
    vec a(2);
    a << 0 << endr << 0 << endr;

    assert_equal(0.0, norm(a, 2));

    a.ones();
    assert_equal(1.0, norm(a, 2));

    a(0) = 3.0;
    a(1) = 4.0;
    assert_near(5.0, norm(a, 2), 1.e-5);

    a(0) = -5.0;
    a(1) = -12.0;
    assert_near(13.0, norm(a, 2), 1.e-5);
    printf(" ...passed.\n");

    /* clean up and exit */
    printf("All tests passed.\n\n");
    return 0;
}

