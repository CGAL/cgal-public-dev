#ifndef CGAL_LOD_PARAMETERS_H
#define CGAL_LOD_PARAMETERS_H

// STL includes.
#include <string>

namespace CGAL {
namespace Levels_of_detail {  

  template<typename FT>
  struct Parameters {

  public:

    // Input.

    // Path to the input data file.
    std::string data;

    // Label indices defined in the ply header: 
    // ground, building boundary, building interior, vegetation.
    std::string gi, bi, ii, vi;


    // Main parameters.

    // Scale in meters used to detect LOD.
    FT scale;

    // Noise level in meters used to detect LOD.
    FT noise_level;


    // Detecting building boundaries.

    // Alpha shape 2 size in meters.
    FT alpha_shape_size_2;

    // Grid 2 cell width in meters.
    FT grid_cell_width_2;

    // Region growing 2.
    FT region_growing_scale_2; // meters / number of points
    FT region_growing_noise_level_2; // meters
    FT region_growing_angle_2; // degrees
    FT region_growing_min_length_2; // meters


    // Computing building footprints.

    // Kinetic partitioning in 2D.
    FT kinetic_min_face_width_2; // meters
    std::size_t kinetic_max_intersections_2; // number

    // Tagging buildings.
    std::size_t min_faces_per_building_2; // number


    // Computing tree footprints.

    // Clustering.
    FT tree_grid_cell_width_2; // meters
    FT min_tree_height; // meters

    // Estimation.
    FT min_tree_radius; // meters

    // Creation.
    std::size_t min_faces_per_tree_2; // number


    // Extrusion.
    std::size_t extrusion_type; // 0 - min, 1 - avg, 2 - max - default


    // Detecting building roofs.

    // Region growing 3.
    FT region_growing_scale_3; // meters / number of points
    FT region_growing_noise_level_3; // meters
    FT region_growing_angle_3; // degrees
    FT region_growing_min_area_3; // meters

    // Roof cleaner.
    FT min_roof_size; // meters


    // Computing building roofs.

    // Kinetic partitioning in 3D.
    std::size_t kinetic_max_intersections_3; // number

    // Graph cut in 3D.
    FT graph_cut_beta_3; // floating in [0, 1]


    // Fitting tree models.
    FT tree_precision; // meters


    // Estimating smooth ground.
    FT ground_precision; // meters

    // Constructor.
    Parameters() : 
    data(""),
    gi("0"), bi("1"), ii("2"), vi("3"),
    scale(FT(4)),
    noise_level(FT(2)),
    alpha_shape_size_2(scale / FT(2)),
    grid_cell_width_2(scale / FT(4)),
    region_growing_scale_2(FT(12)),
    region_growing_noise_level_2(noise_level),
    region_growing_angle_2(FT(25)),
    region_growing_min_length_2(scale),
    kinetic_min_face_width_2(scale / FT(2)),
    kinetic_max_intersections_2(2),
    min_faces_per_building_2(2),
    tree_grid_cell_width_2(scale),
    min_tree_height(noise_level * FT(3) / FT(2)),
    min_tree_radius(noise_level),
    min_faces_per_tree_2(12),
    extrusion_type(2),
    region_growing_scale_3(region_growing_scale_2),
    region_growing_noise_level_3(region_growing_noise_level_2),
    region_growing_angle_3(region_growing_angle_2),
    region_growing_min_area_3(scale),
    min_roof_size(scale / FT(2)),
    kinetic_max_intersections_3(2),
    graph_cut_beta_3(FT(1) / FT(10)),
    tree_precision(scale),
    ground_precision(scale)
    { }

    // Update all parameters, which depend on scale and noise_level.
    void update_dependent() {

      alpha_shape_size_2 = scale / FT(2);
      grid_cell_width_2 = scale / FT(4);
        
      region_growing_noise_level_2 = noise_level;
      region_growing_min_length_2 = scale;
      
      kinetic_min_face_width_2 = scale / FT(2);

      tree_grid_cell_width_2 = scale;
      min_tree_height = noise_level * FT(3) / FT(2);
      min_tree_radius = noise_level;

      region_growing_scale_3 = region_growing_scale_2;
      region_growing_noise_level_3 = region_growing_noise_level_2;
      region_growing_angle_3 = region_growing_angle_2;
      region_growing_min_area_3 = scale;

      min_roof_size = scale / FT(2);
      tree_precision = scale;
      ground_precision = scale;
    }

    void save(const std::string path) const {

      const std::string file_path = path + "info.lod";
      std::ofstream file(file_path.c_str(), std::ios_base::out);

      if (!file) {
          
        std::cerr << std::endl << 
          "ERROR: Error saving file with parameters info!" 
        << std::endl << std::endl;

        exit(EXIT_FAILURE);
      }

      file << "Input: " << std::endl;
      file << "data : " << data << std::endl;
      file << std::endl;

      file << "Label indices: " << std::endl;
      file << "gi (ground) : " << gi << std::endl;
      file << "bi (building boundary) : " << bi << std::endl;
      file << "ii (building interior) : " << ii << std::endl;
      file << "vi (vegetation) : " << vi << std::endl;
      file << std::endl;

      file << "Main parameters: " << std::endl;
      file << "scale (meters) : " << scale << std::endl;
      file << "noise_level (meters) : " << noise_level << std::endl;
      file << std::endl;

      file << "Detecting building boudnaries: " << std::endl;
      file << 
        "alpha_shape_size_2 (meters) : " 
      << alpha_shape_size_2 << std::endl;
      file << 
        "grid_cell_width_2 (meters) : " 
      << grid_cell_width_2 << std::endl;
      file << 
        "region_growing_scale_2 (meters / number of points) : " 
      << region_growing_scale_2 << std::endl;
      file << 
        "region_growing_noise_level_2 (meters) : " 
      << region_growing_noise_level_2 << std::endl;
      file << 
        "region_growing_angle_2 (degrees) : " 
      << region_growing_angle_2 << std::endl;
      file << 
        "region_growing_min_length_2 (meters) : " 
      << region_growing_min_length_2 << std::endl;
      file << std::endl;

      file << "Computing building footprints: " << std::endl;
      file <<
        "kinetic_min_face_width_2 (meters) : "
      << kinetic_min_face_width_2 << std::endl;
      file <<
        "kinetic_max_intersections_2 (number) : "
      << kinetic_max_intersections_2 << std::endl;
      file <<
        "min_faces_per_building_2 (number) : "
      << min_faces_per_building_2 << std::endl;
      file << std::endl;

      file << "Computing tree footprints: " << std::endl;
      file <<
        "tree_grid_cell_width_2 (meters) : "
      << tree_grid_cell_width_2 << std::endl;
      file <<
        "min_tree_height (meters) : "
      << min_tree_height << std::endl;
      file <<
        "min_tree_radius (meters) : "
      << min_tree_radius << std::endl;
      file <<
        "min_faces_per_tree_2 (number) : "
      << min_faces_per_tree_2 << std::endl;
      file << std::endl;

      file << "Extrusion: " << std::endl;
      file <<
        "extrusion_type (0 - min, 1 - avg, 2 - max - default) : "
      << extrusion_type << std::endl;
      file << std::endl;

      file << "Detecting building roofs: " << std::endl;
      file << 
        "region_growing_scale_3 (meters / number of points) : " 
      << region_growing_scale_3 << std::endl;
      file << 
        "region_growing_noise_level_3 (meters) : " 
      << region_growing_noise_level_3 << std::endl;
      file << 
        "region_growing_angle_3 (degrees) : " 
      << region_growing_angle_3 << std::endl;
      file << 
        "region_growing_min_area_3 (meters) : " 
      << region_growing_min_area_3 << std::endl;
      file << 
        "min_roof_size (meters) : " 
      << min_roof_size << std::endl;
      file << std::endl;

      file << "Computing building roofs: " << std::endl;
      file <<
        "kinetic_max_intersections_3 (number) : "
      << kinetic_max_intersections_3 << std::endl;
      file << 
        "graph_cut_beta_3 (floating in [0, 1]) : " 
      << graph_cut_beta_3 << std::endl;
      file << std::endl;

      file << "Fitting tree models: " << std::endl;
      file <<
        "tree_precision (meters) : "
      << tree_precision << std::endl;
      file << std::endl;

      file << "Estimating smooth ground: " << std::endl;
      file <<
        "ground_precision (meters) : "
      << ground_precision << std::endl;

      file.close();
    }

  }; // Parameters

} // Levels_of_detail
} // CGAL

#endif // CGAL_LOD_PARAMETERS_H
