#include <cassert>
#include <cmath>
#include "element.h"
#include "test_helpers.h"

int main(int argc, char* argv[])
{
    printf("Setting up data...\n");
    /*
     * Initialize 2 node bar element
     * Bar2(x1, x2, g_dof1, g_dof2)
     */
    Bar2 bar(0, 1, 0, 1);
    double xi;

    printf("Starting tests...\n\n");

    /* test partition of unity */
    printf("Testing partition of unity...");
    for (xi = -1.0; xi <= 1.0; xi+=0.1)
        assert_near(bar.nodes[0].N(xi) + bar.nodes[1].N(xi), 1, 0.01);
    assert_near(bar.nodes[0].N(-1), 1, 0.001);
    assert_near(bar.nodes[0].N(1), -1, 0.001);
    assert_near(bar.nodes[1].N(-1), -1, 0.001);
    assert_near(bar.nodes[1].N(1), 1, 0.001);
    printf(" ...passed.\n");

    /* clean up and exit */
    printf("All tests passed.\n\n");
    return 0;
}

