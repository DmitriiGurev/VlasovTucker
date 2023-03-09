#include "mesh.h"

#include <mshio/mshio.h>
#include <Eigen/Dense>

#include <cassert>
#include <map>
#include <algorithm>

using namespace std;

Mesh::Mesh(string fileName)
{
    // Read a .msh file
    mshio::MshSpec spec = mshio::load_msh(fileName);
        
    // Fill in the vector of nodes
    for (int i = 0; i < spec.nodes.num_nodes; i++)
	{
        Point* point = new Point({spec.nodes.entity_blocks[0].data[3 * i],
                                  spec.nodes.entity_blocks[0].data[3 * i + 1],
                                  spec.nodes.entity_blocks[0].data[3 * i + 2]});

        point->index = points.size(); // insertion index

        points.push_back(point);
    }

    // Fill in the vector of tetrahedra
    for (auto entityBlock : spec.elements.entity_blocks)
    {
        if (entityBlock.element_type == 4)
        {
            for (int i = 0; i < entityBlock.num_elements_in_block; i++)
            {
                int i0 = entityBlock.data[5 * i + 1] - 1;
                int i1 = entityBlock.data[5 * i + 2] - 1;
                int i2 = entityBlock.data[5 * i + 3] - 1;
                int i3 = entityBlock.data[5 * i + 4] - 1;

                Point* p0 = points[i0];
                Point* p1 = points[i1];
                Point* p2 = points[i2];
                Point* p3 = points[i3];

                Tet* tet = new Tet({p0, p1, p2, p3});  

                tet->index = tets.size(); // insertion index

                // Ensure that the orientation is correct
                assert(tet->Orientation() <= 0);

                tets.push_back(tet);

                // Create positively oriented faces
                tet->faces[0] = new Face({p1, p2, p3});
                tet->faces[1] = new Face({p0, p3, p2});
                tet->faces[2] = new Face({p0, p1, p3});
                tet->faces[3] = new Face({p0, p2, p1});

                for (int j = 0; j < 4; j++)
                {
                    Face* face = tet->faces[j];
                    // Fill in the adjacency information for the new faces
                    face->adjTet = tet;
                    face->adjTetInd = j;

                    face->index = faces.size(); // insertion index

                    faces.push_back(face);

                    // Fill in the lookup table
                    _pointsToFaces[{face->points[0],
                                    face->points[1],
                                    face->points[2]}] = face;
                }
            }
        }
    }

    // Read names of physical groups
    map<int, string> physGroupTagToName;
    for (auto physGroup : spec.physical_groups)
        physGroupTagToName[physGroup.tag] = physGroup.name;

    // Map entity tags to corresponding vectors of physical group names
    map<int, vector<string>> entityToPhysGroups;
    for (auto surface : spec.entities.surfaces)
    {
        for (auto physGroupTag : surface.physical_group_tags)  // only for surfaces
            entityToPhysGroups[surface.tag].push_back(physGroupTagToName[physGroupTag]);
    }

    // Scan the boundary faces (.msh does not contain non-boundary faces since they are not 
    // assigned to any physical group)
    for (auto entityBlock : spec.elements.entity_blocks)
    {
        if (entityBlock.element_type == 2)
        {
            for (int i = 0; i < entityBlock.num_elements_in_block; i++)
	        {
                int i0 = entityBlock.data[4 * i + 1] - 1;
                int i1 = entityBlock.data[4 * i + 2] - 1;
                int i2 = entityBlock.data[4 * i + 3] - 1;
       
                Face* face = _pointsToFaces[{points[i0], points[i1], points[i2]}];
                
                face->type = Boundary;
                face->bcTypes = entityToPhysGroups[entityBlock.entity_tag];
                face->entity = entityBlock.entity_tag;
            }
        }
    }

    // Fill in the tet-tet adjacency information
    for(auto face : faces)
    {
        Point* p0 = face->points[0];
        Point* p1 = face->points[1];
        Point* p2 = face->points[2];

        if (_pointsToFaces.count({p0, p2, p1}) > 0)
        {
            Face* invFace = _pointsToFaces[{p0, p2, p1}];
            face->adjTet->adjTets[face->adjTetInd] = invFace->adjTet;
        }
    }       

    // Group boundaries assigned to one periodic-type physical group
    map<char, vector<vector<Face*>>> pairsOfSymPlanes;
    for (auto pair : entityToPhysGroups)
    {
        int entityTag = pair.first;
        for (auto name : pair.second)
        {
            if (name.substr(0, 19) == "Boundary: Periodic ")
            {
                char tag = name[19];

                vector<Face*> plane;
                for (auto face : faces)
                {
                    if (face->type == Boundary && face->entity == entityTag)
                        plane.push_back(face);
                }
                pairsOfSymPlanes[tag].push_back(plane);
            }
        }
    }
    
    // Lexicographically sort the faces so that parallel planes are sorted identically
    for (auto& pair : pairsOfSymPlanes)
    {
        for (auto& plane : pair.second)
        {
            sort(plane.begin(), plane.end(),
                [&](const auto& lhs, const auto& rhs)
                {
                    if (lhs->centroid.x != rhs->centroid.x)
                    {
                        return lhs->centroid.x < rhs->centroid.x;
                    }
                    if (lhs->centroid.y != rhs->centroid.y)
                    {
                        return lhs->centroid.y < rhs->centroid.y;
                    }
                    if (lhs->centroid.z != rhs->centroid.z)
                    {
                        return lhs->centroid.z < rhs->centroid.z;
                    }
                    return false;
                });
        }
    }

    // Connect the tetrahedra adjacent to the periodic boundaries
    for (auto pair : pairsOfSymPlanes)
    {
        for (int i = 0; i < pair.second[0].size(); i++)
        {
            Face* face = pair.second[0][i];
            Face* oppFace = pair.second[1][i];
            
            face->adjTet->adjTets[face->adjTetInd] = oppFace->adjTet;
            oppFace->adjTet->adjTets[oppFace->adjTetInd] = face->adjTet;
        }
    }
}

Mesh::~Mesh()
{
    for (auto point : points)
        delete point;
    
    for (auto face : faces)
        delete face;

    for (auto tet : tets)
        delete tet;
}

// All cyclic permutations are equivalent
bool Mesh::KeyTriple::operator==(const KeyTriple& keyTriple) const
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

std::size_t Mesh::HashTriple::operator()(const KeyTriple& keyTriple) const
{
    return std::hash<Point*>()(keyTriple.p0) & 
            std::hash<Point*>()(keyTriple.p1) & 
            std::hash<Point*>()(keyTriple.p2);
}

Point::Point(double x, double y, double z) :
    x(x), y(y), z(z) {}

Point Point::operator+(const Point& p) const
{
    return Point(x + p.x, y + p.y, z + p.z);
}

Point Point::operator-(const Point& p) const
{
    return Point(x - p.x, y - p.y, z - p.z);
}

Point Point::operator/(double d) const
{
    return Point(x / d, y / d, z / d);
}

bool Point::operator==(const Point& p) const
{
    return x == p.x &&
            y == p.y &&
            z == p.z;
}

double Point::Abs() const
{
    return sqrt(x * x + y * y + z * z);
}

double Point::DotProduct(const Point& p) const
{
    return x * p.x + y * p.y + z * p.z;
}

Point Point::CrossProduct(const Point& p) const
{
    return {y * p.z - z * p.y,
            z * p.x - x * p.z,
            x * p.y - y * p.x};
}

std::ostream& operator<<(std::ostream& os, const Point& p)
{
    os << "Point: {" << p.x << ", " << p.y  << ", " << p.z << "}";
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
        throw std::runtime_error("");

    normal = normal / normal.Abs();

    area = a.CrossProduct(b).Abs() / 2.0;
}

std::ostream& operator<<(std::ostream& os, const Face& f)
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

    volume = std::abs(Orientation()) / 6.0;
}

double Tet::Orientation() const
{
    Eigen::Matrix4d m;
    m << points[0]->x, points[0]->y, points[0]->z, 1,
         points[1]->x, points[1]->y, points[1]->z, 1,
         points[2]->x, points[2]->y, points[2]->z, 1,
         points[3]->x, points[3]->y, points[3]->z, 1;

    return m.determinant();
}

std::ostream& operator<<(std::ostream& os, const Tet& t)
{
    os << "Tet: {" << *t.points[0] << ",\n"
       << "      " << *t.points[1] << ",\n"
       << "      " << *t.points[2] << ",\n"
       << "      " << *t.points[3] << "}";

    return os;
}