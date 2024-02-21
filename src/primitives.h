#pragma once

#include <iostream>
#include <string>
#include <vector>
#include <array>
#include <cmath>

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
    // Coordinate access
    double operator[](int i) const;
    double& operator[](int i);

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

enum class FaceType
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

    FaceType type = FaceType::Internal;
    std::vector<std::string> bcTypes; // boundary conditions

    int index; // insertion index

    int entity = -1; // entity tag

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

struct Triple
{
    Point* p0;
    Point* p1;
    Point* p2;

    // All cyclic permutations are equivalent
    bool operator==(const Triple& triple) const;
};

struct HashTriple // hashing function for the lookup table
{
    std::size_t operator()(const Triple& triple) const;
};
}