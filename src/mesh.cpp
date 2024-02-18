#include "mesh.h"

#include <cassert>
#include <map>
#include <algorithm>

#include <fstream>

using namespace std;

namespace VlasovTucker
{
Mesh::Mesh(string mshFile)
{
    // Read a .msh file
    _mshSpec = mshio::load_msh(mshFile);

    // Read names of physical groups
    unordered_map<int, string> physGroupTagToName;
    for (auto physGroup : _mshSpec.physical_groups)
        physGroupTagToName[physGroup.tag] = physGroup.name;

    // Map entity tags to corresponding vectors of physical group names
    for (auto surface : _mshSpec.entities.surfaces)
    {
        for (auto physGroupTag : surface.physical_group_tags)  // only for surfaces
            _entityToPhysGroups[surface.tag].push_back(physGroupTagToName[physGroupTag]);
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

unordered_map<int, vector<string>> Mesh::BoundaryLabels() const
{
    return _entityToPhysGroups;
}

void Mesh::PrintBoundaryLabels() const
{
    for (auto pair : _entityToPhysGroups)
    {
        int number = pair.first;
        int nLabels = pair.second.size();

        cout << number << ": {";
        if (nLabels == 0)
        {
            cout << "}\n";
        }
        else
        {
            for (int i = 0; i < nLabels - 1; i++) 
                cout << pair.second[i] << ", ";
            cout << pair.second[nLabels - 1] << "}\n";
        }
    }
}

void Mesh::SetPeriodicBounaries(const std::vector<std::array<int, 2>>& periodicPairs)
{
    _periodicPairs = periodicPairs;
}

std::vector<std::array<int, 2>> Mesh::PeriodicBoundaries() const
{
    return _periodicPairs;
}

void Mesh::Reconstruct(double scaleFactor)
{
    // Fill the vector of nodes
    _ExtractPoints(scaleFactor);

    // Fill the vectors of tetrahedra and faces
    _ExtractTetsAndFaces();

    // Scan the boundary faces (.msh does not contain non-boundary faces since they are not 
    // assigned to any physical group)
    _LabelBoundaryFaces();

    // Fill the tet-tet adjacency information
    _FillAdjacencyInfo();

    // Connect the tetrahedra adjacent to the periodic boundaries (if any)
    _ConfigurePeriodicity();
}

void Mesh::_ExtractPoints(double scaleFactor)
{
    for (int i = 0; i < _mshSpec.nodes.num_nodes; i++)
	{
        Point* point = new Point({_mshSpec.nodes.entity_blocks[0].data[3 * i],
                                  _mshSpec.nodes.entity_blocks[0].data[3 * i + 1],
                                  _mshSpec.nodes.entity_blocks[0].data[3 * i + 2]});

        *point = *point * scaleFactor; // scaling

        point->index = points.size(); // insertion index

        points.push_back(point);
    }
}

void Mesh::_ExtractTetsAndFaces()
{
    for (auto entityBlock : _mshSpec.elements.entity_blocks)
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

                // Make sure that the orientation is correct
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
                    // Fill the adjacency information for the new faces
                    face->adjTet = tet;
                    face->adjTetInd = j;

                    face->index = faces.size(); // insertion index

                    faces.push_back(face);

                    // Fill the lookup table
                    _pointsToFaces[{face->points[0],
                                    face->points[1],
                                    face->points[2]}] = face;
                }
            }
        }
    }
}

void Mesh::_LabelBoundaryFaces()
{
    for (auto entityBlock : _mshSpec.elements.entity_blocks)
    {
        if (entityBlock.element_type == 2)
        {
            for (int i = 0; i < entityBlock.num_elements_in_block; i++)
	        {
                int i0 = entityBlock.data[4 * i + 1] - 1;
                int i1 = entityBlock.data[4 * i + 2] - 1;
                int i2 = entityBlock.data[4 * i + 3] - 1;
       
                Face* face = _pointsToFaces[{points[i0], points[i1], points[i2]}];
                
                face->type = FaceType::Boundary;
                face->bcTypes = _entityToPhysGroups[entityBlock.entity_tag];
                face->entity = entityBlock.entity_tag;
            }
        }
    }
}

void Mesh::_FillAdjacencyInfo()
{
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
}

void SortFacesInPlane(vector<Face*>& plane)
{
    auto ApproxEqual = [](double a, double b)
    {
        return abs(a - b) < 1e-8;
    };

    sort(plane.begin(), plane.end(), [&](const auto& lhs, const auto& rhs)
    {
        if (!ApproxEqual(lhs->centroid.coords[0], rhs->centroid.coords[0]))
        {
            return lhs->centroid.coords[0] < rhs->centroid.coords[0];
        }
        else if (!ApproxEqual(lhs->centroid.coords[1], rhs->centroid.coords[1]))
        {
            return lhs->centroid.coords[1] < rhs->centroid.coords[1];
        }
        else if (!ApproxEqual(lhs->centroid.coords[2], rhs->centroid.coords[2]))
        {
            return lhs->centroid.coords[2] < rhs->centroid.coords[2];
        }
        else
        {
            return false;
        }
    });
}

void Mesh::_ConfigurePeriodicity()
{
    if (_periodicPairs.empty())
        return;

    // Collect faces into periodic planes
    vector<array<vector<Face*>, 2>> pairsOfPeriodicPlanes;
    for (auto pairInds : _periodicPairs)
    {
        array<vector<Face*>, 2> pairOfPlanes;
        for (int i : {0, 1})
        {
            vector<Face*> plane;
            for (auto face : faces)
            {
                if ((face->type == FaceType::Boundary) && (face->entity == pairInds[i]))
                    plane.push_back(face);
            }
            pairOfPlanes[i] = plane;
        }

        // TODO: Throw a runtime error with an explanation
        assert(pairOfPlanes[0].size() == pairOfPlanes[1].size());

        pairsOfPeriodicPlanes.push_back(pairOfPlanes);
    }

    // Lexicographically sort the faces so that parallel planes were ordered identically
    for (auto& pair : pairsOfPeriodicPlanes)
    {
        for (int i : {0, 1})
            SortFacesInPlane(pair[i]);
    }

    for (auto pairOfPlanes : pairsOfPeriodicPlanes)
    {
        for (int i = 0; i < pairOfPlanes[0].size(); i++)
        {
            Face* face = pairOfPlanes[0][i];
            Face* oppFace = pairOfPlanes[1][i];
            
            face->adjTet->adjTets[face->adjTetInd] = oppFace->adjTet;
            oppFace->adjTet->adjTets[oppFace->adjTetInd] = face->adjTet;
        }
    }
}
}