#ifndef __POINT_H__
#define __POINT_H__

#include <iostream>

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

double euclidean_distance(Point2D&, Point2D&);

}

#endif /* __POINT_H__ */
