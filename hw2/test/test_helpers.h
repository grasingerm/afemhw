#ifndef _TEST_HELPERS
#define _TEST_HELPERS

#include <cmath>

#define assert_equal(a, b) assert((a)==(b))
#define assert_near(a, b, tol) assert(abs(((a)-(b))/(a)) <= (tol))
#define assert_nequal(a, b) assert((a)!=(b))
#define assert_lt(a, b) assert((a)<(b))
#define assert_le(a, b) assert((a)<=(b))
#define assert_gt(a, b) assert((a)>(b))
#define assert_ge(a, b) assert((a)>=(b))

#endif

