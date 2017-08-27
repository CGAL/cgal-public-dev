#ifndef HEXEXTR_H
#define HEXEXTR_H
#pragma once
#include"typedefs.h"

class HexExtr{
  public:
void save_parametrized_mesh();
//defintion in hexextr.cpp
    HexExtr(std::string);
    int set_dart_info(LCC_3&);
    void load_mesh(std::string);
    void extract(std::string);
    int calculate_cell_type(LCC_3&, Dart_handle);
    void set_parametrization(LCC_3&);
    void save_mesh(std::string);

// definition in preprocessing.cpp
    Aff_transformation get_parametrization_matrix(Point, Point, Point, Point, Point, Point, Point, Point);
    Aff_transformation extract_transition_function(Dart_handle);
    void preprocess();
    void set_chiral_symmetry_aff_transformations();

//definiton in sanitization.cpp
    void sanitize();
    void truncate_precision();
    void check_singularity();
    void fix_singularity(Dart_handle&);
    Point_3 get_projected_parameters(Dart_handle&, Dart_handle&);
    Dart_handle get_singular_edge(Dart_handle&);
    bool is_singular(Dart_handle&);
    Aff_transformation get_transition(Dart_handle&, Dart_handle&);
    void propagate_parameters(Dart_handle&);

//definition in hexahedron_extraction.cpp
    void extract_hexes();
    Aff_transformation* find_tet_parametrization(Point);
    bool does_intersect(Tetrahedron_3, Point_3);

//definition in connection_extraction.cpp
    void extract_connections();

//definition in post_processing.cpp
    void refine();
    bool post_processing_req();
    void post_processing();
    void merge_vertices();
    void annihilate_darts();
    int flipped_cells();
    bool are_hexes();
    bool are_quad_strips();
    bool are_quads();
    
// input_tet_mesh holds the initial parametrized tet mesh, and output_mesh holds the hexahedral mesh extracted after running extract() 
    LCC_3 input_tet_mesh, output_mesh;

// contains six orthogonal directions in 3D space given by (1,0,0), (0,1,0), (0,0,1), (-1,0,0), (0,-1,0), (0,0,-1) 
    std::vector<Direction> directions;

//identity tranformation
    Aff_transformation identity;

// chiral cubical symmetry group used in calculation of g
    std::vector<Aff_transformation> G;

// g is a 2D vector of transformations: it gives the transformation to be applied on one tet to convert into adjacent tet.
    std::vector<std::vector<Aff_transformation>> g;

// Indexed by the unique number assigned to each tet in input_tet_mesh, this gives the transformation to be applied to the vertices of a tet to tranform to parametrized space
    std::vector<Aff_transformation> parametrization_matrices;

// Keeps track of points in output_mesh from which make_hexahedron() was called (to prevent duplicate hexes)
    std::unordered_map<Point_3, Dart_handle> hex_handles;

// Keeps track of all points in output_mesh, by mapping the points of hexahedra found in parametrized space to the cartesian space.
    std::unordered_map<Point_3, Point_3> output_points;
       
};

#endif
