#ifndef HEXEXTR_H
#define HEXEXTR_H
#pragma once
#include"typedefs.h"

class HexExtr{
  public:
    HexExtr(std::string);
    Aff_transformation get_parametrization_matrix(Point, Point, Point, Point, Point, Point, Point, Point);
    Aff_transformation extract_transition_function(Dart_handle, const std::vector<Aff_transformation>&);
    int set_dart_info();
    void load_mesh(std::string);
    void preprocess();
    void set_chiral_symmetry_aff_transformations();
    void sanitize(std::vector<std::vector<Aff_transformation>>&);
    void truncate_precision(std::vector<std::vector<Aff_transformation>>&);
    void check_singularity(std::vector<std::vector<Aff_transformation>>&);
    void fix_singularity(Dart_handle&, std::vector<std::vector<Aff_transformation>>&);
    Point_3 get_projected_parameters(Dart_handle&, Dart_handle&, std::vector<std::vector<Aff_transformation>>&);
    Dart_handle get_singular_edge(Dart_handle&, std::vector<std::vector<Aff_transformation>>&);
    bool is_singular(Dart_handle&, std::vector<std::vector<Aff_transformation>>&);
    Aff_transformation get_transition(Dart_handle&, Dart_handle&, std::vector<std::vector<Aff_transformation>>&);
    void propagate_parameters(Dart_handle&, std::vector<std::vector<Aff_transformation>>&);
    void extract_hexes(std::vector<Aff_transformation>&);
    Aff_transformation* find_tet_parametrization(Point, std::vector<Aff_transformation>&);
    bool does_intersect(Tetrahedron_3, Point_3);
    int calculate_cell_type(LCC_3&, Dart_handle);
    void extract_connections();
    void extract();
    void refine();
    bool post_processing_req();
    void post_processing();
    void merge_vertices();
    void annihilate_darts();
    int flipped_cells();
    bool are_hexes();
    bool are_quad_strips();
    bool are_quads();
    
    std::unordered_map<Point_3, Dart_handle> hex_handles;
    std::unordered_map<Point_3, Point_3> output_points;
    std::vector<std::vector<Aff_transformation>> g;
    std::vector<Aff_transformation> parametrization_matrices;
    std::vector<Direction> directions;
    Aff_transformation identity;
    LCC_3 input_tet_mesh, output_mesh;
    std::vector<Aff_transformation> G; //chiral cubical symmetry group
       
};


#endif
