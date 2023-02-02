#pragma once

#include "header.h"
 
enum bc_type
{
    SYMMETRYZ = 0,
    INLET = 	1,
    OUTLET = 	2,
    WALL = 		3,
    SYMMETRYY = 4,
    SYMMETRYX = 5 
};

class Mesh {
public:
	Mesh(const std::string& path, double scale = 1.0);
	void read(const std::string& path, double scale = 1.0);

	double compute_tetra_volume(std::vector<std::vector<double>> tetra);
	std::vector<std::vector<int>> get_faces_for_cell(int ic);

	int nbf;
	int nbc;
	int nc;
	int nv;
	int nf;

	std::vector<std::vector<int>> 	 bcface_vert_lists;
	std::vector<std::set<int>> 		 bcface_vert_set_lists;
	std::vector<int> 				 bcface_bctype;
	std::vector<std::vector<int>> 	 vert_list_for_cell;
	std::vector<std::vector<int>> 	 bf_for_each_bc;
	std::vector<std::vector<double>> vert_coo;
	std::vector<std::vector<double>> cell_center_coo;
	std::vector<double>				 cell_volumes;
	std::vector<std::vector<int>> 	 cell_list_for_vertex;
	std::vector<int>				 cell_num_for_vertex;
	std::vector<std::vector<int>> 	 cell_neighbors_list;
	std::vector<std::vector<int>> 	 face_vert_list;
	std::vector<std::vector<int>> 	 cell_face_list;
	std::vector<double>				 face_areas;
	std::vector<std::vector<double>> face_normals;
	std::vector<std::vector<int>> 	 cell_face_normal_direction;
	std::vector<std::vector<int>> 	 bound_face_info;
	std::vector<double>				 cell_diam;
	std::vector<std::vector<double>> face_centers;

	void write_tecplot(std::vector<std::vector<double>>data, std::string filename,
		std::vector<std::string> var_names, double time = 0.0);
};