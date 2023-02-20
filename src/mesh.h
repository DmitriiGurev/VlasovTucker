#pragma once

#include <iostream>
#include <cmath>
#include <string>
#include <vector>
#include <array>
#include <unordered_map>

class Point
{
public:
    Point() {}

    Point(double x, double y, double z) :
        x(x), y(y), z(z) {}

public:
    double x;
    double y;
    double z;

    int index;

public:
    Point operator+(const Point& p) const
    {
        return Point(x + p.x, y + p.y, z + p.z);
    }

    Point operator-(const Point& p) const
    {
        return Point(x - p.x, y - p.y, z - p.z);
    }

    Point operator/(double d) const
    {
        return Point(x / d, y / d, z / d);
    }

    bool operator==(const Point& p) const
    {
        return x == p.x &&
               y == p.y &&
               z == p.z;
    }

    double Abs() const
    {
        return sqrt(x * x + y * y + z * z);
    }

    double DotProduct(const Point& p) const
    {
        return x * p.x + y * p.y + z * p.z;
    }

    Point CrossProduct(const Point& p) const
    {
        return {y * p.z - z * p.y,
                z * p.x - x * p.z,
                x * p.y - y * p.x};
    }

    friend std::ostream& operator<<(std::ostream& os, const Point& p);
};

// Forward declaration
class Tet;

enum FaceType
{
    Internal,
    Boundary
};

class Face
{
public:
    Face() {}

    Face(Point* p0, Point* p1, Point* p2) :
        points({p0, p1, p2})
    {
        centroid = (*p0 + *p1 + *p2) / 3.0;        

        Point a = *p1 - *p0;
        Point b = *p2 - *p0;

        normal = {a.y * b.z - a.z * b.y,
                  a.z * b.x - a.x * b.z,
                  a.x * b.y - a.y * b.x};

        normal = normal / normal.Abs();

        area = a.CrossProduct(b).Abs() / 2.0;
    }

public:
    std::array<Point*, 3> points = {nullptr, nullptr, nullptr};

    // The tet this face is adjacent to
    Tet* adjTet = nullptr;
    // The face-index in the adjacent tet
    int adjTetInd;

    double area;
    Point centroid;
    Point normal;

    FaceType type = Internal;

    std::vector<std::string> bcTypes;

    int index;

    int entity = 0;

public:
    friend std::ostream& operator<<(std::ostream& os, const Face& f);
};

class Tet
{
public:
    Tet() {}

    Tet(Point* p0, Point* p1, Point* p2, Point* p3) :
        points({p0, p1, p2, p3})
    {
        centroid = (*p0 + *p1 + *p2 + *p3) / 4.0;

        volume = std::abs(Orientation()) / 6.0;
    }

    double Orientation() const;

public:
    std::array<Point*, 4> points  = {nullptr, nullptr, nullptr, nullptr};
    std::array<Face*, 4>  faces   = {nullptr, nullptr, nullptr, nullptr};
    std::array<Tet*, 4>   adjTets = {nullptr, nullptr, nullptr, nullptr};

    Point centroid;
    double volume;

    int index;

public:
    friend std::ostream& operator<<(std::ostream& os, const Tet& t);
};

class Mesh
{
public:
    Mesh(std::string fileName);

    ~Mesh();

public:
    std::vector<Point*> points;
    std::vector<Face*>  faces;
    std::vector<Tet*>   tets;

private:
    struct KeyTriple
    {
        Point* p0;
        Point* p1;
        Point* p2;
    
        // All cyclic permutations are equivalent
        bool operator==(const KeyTriple& keyTriple) const
        {
            return (p0 == keyTriple.p0 &&
                    p1 == keyTriple.p1 &&
                    p2 == keyTriple.p2) ||
                    (p0 == keyTriple.p1 &&
                    p1 == keyTriple.p2 &&
                    p2 == keyTriple.p0) ||
                    (p0 == keyTriple.p2 &&
                    p1 == keyTriple.p0 &&
                    p2 == keyTriple.p1);
        }
    };
    
    struct HashTriple
    {
        std::size_t operator()(const KeyTriple& keyTriple) const
        {
            return std::hash<Point*>()(keyTriple.p0) & 
                   std::hash<Point*>()(keyTriple.p1) & 
                   std::hash<Point*>()(keyTriple.p2);
        }
    };
    
    std::unordered_map<KeyTriple, Face*, HashTriple> _pointsToFaces;
};