#ifndef __POINT_H__
#define __POINT_H__

#include <iostream>
#include <armadillo>
#include <array>
#include <cmath>

namespace pt
{

enum euclidean_component { X, Y, Z };

class Point2D
{
public:
    Point2D() {}
    Point2D(const double x, const double y) { coords[X] = x; coords[Y] = y; }
    
    Point2D(const Point2D& p) 
        { coords[X] = p.coords[X]; coords[Y] = p.coords[Y]; }
    Point2D& operator=(const Point2D&);
    
    Point2D(Point2D&&);
    Point2D& operator=(Point2D&&);
   
    inline const double& operator[] (const int c) const { return coords[c]; }
    inline double& operator[] (const int c) { return coords[c]; }
    
    friend std::ostream& operator<< (std::ostream& stream, const Point2D& p)
    {
        std::cout << "(" << p.coords[X] << "," << p.coords[Y] << ")";
        return stream;
    }
    
private:
    double coords[2];
};

/* operator overloading */
Point2D operator+(const Point2D&, const Point2D&);
Point2D operator-(const Point2D&, const Point2D&);
Point2D operator-(const Point2D&);
Point2D operator*(const Point2D&, const double);
inline Point2D operator*(const double m, Point2D& p) { return p*m; }
bool operator==(const Point2D&, const Point2D&);

double euclidean_distance(const Point2D&, const Point2D&);
arma::vec::fixed<2> vec_from_pts(const Point2D&, const Point2D&);
arma::vec::fixed<2> unit_normal(const arma::vec&);
inline double norm_2D(const arma::vec& u) { return sqrt(u(0)*u(0)+u(1)*u(1)); }

class Node2D
{
public:
    Node2D(const Point2D& p, const unsigned int gdof1, const unsigned int gdof2)
        : pt(p) { gdofs[0] = gdof1; gdofs[1] = gdof2; }
    std::array<unsigned int,2> gdofs;
}

}

#endif /* __POINT_H__ */
