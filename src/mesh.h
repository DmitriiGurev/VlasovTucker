#pragma once

#include <iostream>
#include <cmath>
#include <string>
#include <vector>
#include <array>
#include <unordered_map>

namespace VlasovTucker
{
class Point
{
public:
    Point() {}

    Point(std::array<double, 3> coords);

public:
    std::array<double, 3> coords;

    int index; // insertion index

public:
    Point operator+(const Point& p) const;
    Point operator-(const Point& p) const;
    Point operator/(double d) const;
    Point operator*(double d) const;

    bool operator==(const Point& p) const;

    double Abs() const;
    double DotProduct(const Point& p) const;
    Point CrossProduct(const Point& p) const;

    friend std::ostream& operator<<(std::ostream& os, const Point& p);
};

class Tet; // forward declaration

enum FaceType
{
    Internal,
    Boundary
};

class Face
{
public:
    Face() {}

    Face(Point* p0, Point* p1, Point* p2);

public:
    std::array<Point*, 3> points = {nullptr, nullptr, nullptr};

    Tet* adjTet = nullptr; // the tetrahedron this face is adjacent to
    int adjTetInd; // the face-index in the adjacent tetrahedron

    double area;
    Point centroid;
    Point normal;

    FaceType type = Internal;
    std::vector<std::string> bcTypes; // boundary conditions

    int index; // insertion index

    int entity = 0; // entity tag

public:
    friend std::ostream& operator<<(std::ostream& os, const Face& f);
};

class Tet
{
public:
    Tet() {}

    Tet(Point* p0, Point* p1, Point* p2, Point* p3);

    double Orientation() const;

public:
    std::array<Point*, 4> points  = {nullptr, nullptr, nullptr, nullptr};
    std::array<Face*, 4>  faces   = {nullptr, nullptr, nullptr, nullptr};
    std::array<Tet*, 4>   adjTets = {nullptr, nullptr, nullptr, nullptr};

    Point centroid;
    double volume;

    int index; // insertion index

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
    struct KeyTriple // key for the lookup table
    {
        Point* p0;
        Point* p1;
        Point* p2;
    
        // All cyclic permutations are equivalent
        bool operator==(const KeyTriple& keyTriple) const;
    };
    
    struct HashTriple // hashing function for the lookup table
    {
        std::size_t operator()(const KeyTriple& keyTriple) const;
    };
    
    std::unordered_map<KeyTriple, Face*, HashTriple> _pointsToFaces;
};
}