#include "primitives.h"

#include <Eigen/Dense>

using namespace std;

namespace VlasovTucker
{
Point::Point(Vector3d coords) :
    coords(coords) {}

double Point::operator[](int i) const
{
    return coords[i];
}

double& Point::operator[](int i)
{
    return coords[i];
}

Point Point::operator+(const Point& p) const
{
    return Point({coords[0] + p.coords[0],
                  coords[1] + p.coords[1],
                  coords[2] + p.coords[2]});
}

Point Point::operator-(const Point& p) const
{
    return Point({coords[0] - p.coords[0],
                  coords[1] - p.coords[1],
                  coords[2] - p.coords[2]});
}

Point operator-(const Point& p)
{
    return (-1) * p;
}

Point Point::operator/(double d) const
{
    return Point({coords[0] / d,
                  coords[1] / d,
                  coords[2] / d});
}

Point Point::operator*(double d) const
{
    return Point({coords[0] * d,
                  coords[1] * d,
                  coords[2] * d});
}

Point operator*(double d, const Point& p)
{
    return p * d;
}

bool Point::operator==(const Point& p) const
{
    return coords[0] == p.coords[0] &&
           coords[1] == p.coords[1] &&
           coords[2] == p.coords[2];
}

double Point::Abs() const
{
    return sqrt(coords[0] * coords[0] + 
                coords[1] * coords[1] +
                coords[2] * coords[2]);
}

double Point::DotProduct(const Point& p) const
{
    return coords[0] * p.coords[0] +
           coords[1] * p.coords[1] +
           coords[2] * p.coords[2];
}

Point Point::CrossProduct(const Point& p) const
{
    return Point({coords[1] * p.coords[2] - coords[2] * p.coords[1],
                  coords[2] * p.coords[0] - coords[0] * p.coords[2],
                  coords[0] * p.coords[1] - coords[1] * p.coords[0]});
}

ostream& operator<<(ostream& os, const Point& p)
{
    os << "Point: {" << p[0] << ", " << p[1]  << ", " << p[2] << "}";
    return os;
}

Face::Face(Point* p0, Point* p1, Point* p2) :
    points({p0, p1, p2})
{
    centroid = (*p0 + *p1 + *p2) / 3.0;        

    Point a = *p1 - *p0;
    Point b = *p2 - *p0;

    normal = a.CrossProduct(b);

    double nLength = normal.Abs();
    if (nLength == 0)
        throw runtime_error("");

    normal = normal / normal.Abs();

    area = a.CrossProduct(b).Abs() / 2.0;
}

ostream& operator<<(ostream& os, const Face& f)
{
    os << "Face: {" << *f.points[0] << ",\n"
       << "       " << *f.points[1] << ",\n"
       << "       " << *f.points[2];

    return os;
}

Tet::Tet(Point* p0, Point* p1, Point* p2, Point* p3) :
    points({p0, p1, p2, p3})
{
    centroid = (*p0 + *p1 + *p2 + *p3) / 4.0;

    volume = abs(Orientation()) / 6.0;
}

double Tet::Orientation() const
{
    Eigen::Matrix4d m;
    m << points[0]->coords[0], points[0]->coords[1], points[0]->coords[2], 1,
         points[1]->coords[0], points[1]->coords[1], points[1]->coords[2], 1,
         points[2]->coords[0], points[2]->coords[1], points[2]->coords[2], 1,
         points[3]->coords[0], points[3]->coords[1], points[3]->coords[2], 1;

    return m.determinant();
}

ostream& operator<<(ostream& os, const Tet& t)
{
    os << "Tet: {" << *t.points[0] << ",\n"
       << "      " << *t.points[1] << ",\n"
       << "      " << *t.points[2] << ",\n"
       << "      " << *t.points[3] << "}";

    return os;
}

// All cyclic permutations are equivalent
bool Triple::operator==(const Triple& triple) const
{
    return (p0 == triple.p0 && p1 == triple.p1 && p2 == triple.p2) ||
           (p0 == triple.p1 && p1 == triple.p2 && p2 == triple.p0) ||
           (p0 == triple.p2 && p1 == triple.p0 && p2 == triple.p1);
}

size_t HashTriple::operator()(const Triple& triple) const
{
    return hash<Point*>()(triple.p0) & 
           hash<Point*>()(triple.p1) & 
           hash<Point*>()(triple.p2);
}
}