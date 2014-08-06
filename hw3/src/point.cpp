#include <cstdlib>
#include <cmath>
#include "point.hpp"

#define TOLERANCE 1e-5

namespace pt
{

Point2D& Point2D::operator=(const Point2D& p)
{
    for (int i = 0; i < 3; i++) coords[i] = p.coords[i];
    return *this;
}

Point2D::Point2D(Point2D&& p)
{
    for (int i = 0; i < 3; i++) coords[i] = p.coords[i];
}

Point2D& Point2D::operator=(Point2D&& p)
{
    for (int i = 0; i < 3; i++) coords[i] = p.coords[i];
    return *this;
}

Point2D operator+(const Point2D& p, const Point2D& q)
{
    return Point2D(p[X]+q[X], p[Y]+q[Y]);
}

Point2D operator-(const Point2D& p)
{
    return Point2D(-p[X], -p[Y]);
}

Point2D operator-(const Point2D& p, const Point2D& q)
{
    return Point2D(p[X]-q[X], p[Y]-q[Y]);
}

bool operator==(const Point2D& p, const Point2D &q)
{
    for (int i = 0; i < 3; i++)
        if (fabs(p[i]-q[i]) > 1e-5) return false;

    return true;
}

Point2D operator*(const Point2D& p, const double multiplier)
{
    return Point2D(multiplier*p[X], multiplier*p[Y]);
}

double euclidean_distance(const Point2D& p, const Point2D& q)
{
    return sqrt((p[X]-q[X])*(p[X]-q[X]) + (p[Y]-q[Y])*(p[Y]-q[Y]));
}

arma::vec::fixed<2> vec_from_pts(const Point2D& p, const Point2D& q)
{
    arma::vec::fixed<2> line;
    line << p[X]-q[X] << p[Y]-q[Y];
    return line;
}

arma::vec::fixed<2> unit_normal(const arma::vec& u)
{
    arma::vec::fixed<2> N;
    double slope = u(1)/u(0);
    N << 1 << -1/slope;
    N /= norm_2D(N);
    return N;
}

}
