#include "integration.hpp"

template <class T> gquad::GaussQuad1D
    ::integrate(function<T(const double)> f, const double x)
