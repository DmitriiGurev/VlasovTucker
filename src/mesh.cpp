#include "mesh.h"

#include <mshio/mshio.h>
#include <Eigen/Dense>

#include <cassert>
#include <map>
#include <tuple>

using namespace std;

Mesh::Mesh(string fileName)
{
    mshio::MshSpec spec = mshio::load_msh(fileName);

    map<int, string> physGroupTagToName;
    for (auto physGroup : spec.physical_groups)
        physGroupTagToName[physGroup.tag] = physGroup.name;

    // Only for surfaces
    map<int, vector<string>> entityToPhysGroups;
    for (auto surface : spec.entities.surfaces)
    {
        for (auto physGroupTag : surface.physical_group_tags)
            entityToPhysGroups[surface.tag].push_back(physGroupTagToName[physGroupTag]);
    }
        
    for (int i = 0; i < spec.nodes.num_nodes; i++)
	{
        Point* point = new Point({spec.nodes.entity_blocks[0].data[3 * i],
                                  spec.nodes.entity_blocks[0].data[3 * i + 1],
                                  spec.nodes.entity_blocks[0].data[3 * i + 2]});
        point->index = points.size();
        points.push_back(point);
    }

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
                tet->index = tets.size();

                assert(tet->Orientation() <= 0);

                tets.push_back(tet);

                tet->faces[0] = new Face({p1, p2, p3});
                tet->faces[1] = new Face({p0, p3, p2});
                tet->faces[2] = new Face({p0, p1, p3});
                tet->faces[3] = new Face({p0, p2, p1});

                for (int j = 0; j < 4; j++)
                {
                    Face* face = tet->faces[j];
                    face->adjTet = tet;
                    face->adjTetInd = j;
                    face->index = faces.size();
                    faces.push_back(face);  
                    _pointsToFaces[{*face->points[0],
                                    *face->points[1],
                                    *face->points[2]}] = face;
                }
            }
        }
    }

    // Boundary faces
    for (auto entityBlock : spec.elements.entity_blocks)
    {
        if (entityBlock.element_type == 2)
        {
            for (int i = 0; i < entityBlock.num_elements_in_block; i++)
	        {
                int i0 = entityBlock.data[4 * i + 1] - 1;
                int i1 = entityBlock.data[4 * i + 2] - 1;
                int i2 = entityBlock.data[4 * i + 3] - 1;
       
                Face* face = _pointsToFaces[{*points[i0], *points[i1], *points[i2]}];
                face->type = Boundary;
                face->bcTypes = entityToPhysGroups[entityBlock.entity_tag];
            }
        }
    }

    for(auto face : faces)
    {
        Point* p0 = face->points[0];
        Point* p1 = face->points[1];
        Point* p2 = face->points[2];

        if (_pointsToFaces.count({*p0, *p2, *p1}) > 0)
        {
            Face* invFace = _pointsToFaces[{*p0, *p2, *p1}];
            face->adjTet->adjTets[face->adjTetInd] = invFace->adjTet;
        }

        if (face->type == Boundary)
        {
            if (find(face->bcTypes.begin(),
                    face->bcTypes.end(),
                    "Boundary: Periodic") !=
                    face->bcTypes.end())
            {
                for (Point t : {Point(1, 0, 0),
                                Point(0, 1, 0),
                                Point(0, 0, 1),
                                Point(-1, 0, 0),
                                Point(0, -1, 0),
                                Point(0, 0, -1)})
                    if (_pointsToFaces.count({*p0 + t, *p2 + t, *p1 + t}) > 0)
                    {
                        Face* invFace = _pointsToFaces[{*p0 + t, *p2 + t, *p1 + t}];
                        face->adjTet->adjTets[face->adjTetInd] = invFace->adjTet;
                    }
            }
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

    cout << "~Mesh\n";
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

std::ostream& operator<<(std::ostream& os, const Point& p)
{
    os << "Point: {" << p.x << ", " << p.y  << ", " << p.z << "}";
    return os;
}

std::ostream& operator<<(std::ostream& os, const Face& f)
{
    os << "Face: {" << *f.points[0] << ",\n"
       << "       " << *f.points[1] << ",\n"
       << "       " << *f.points[2];

    return os;
}

std::ostream& operator<<(std::ostream& os, const Tet& t)
{
    os << "Tet: {" << *t.points[0] << ",\n"
       << "      " << *t.points[1] << ",\n"
       << "      " << *t.points[2] << ",\n"
       << "      " << *t.points[3] << "}";

    return os;
}