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
    FT region_growing_search_size_2; // meters / number of points
    FT region_growing_noise_level_2; // meters
    FT region_growing_angle_2; // degrees
    FT region_growing_minimum_length_2; // meters


    // Constructor.
    Parameters() : 
    data(""),
    gi("0"), bi("1"), ii("2"), vi("3"),
    scale(FT(4)),
    noise_level(FT(2)),
    alpha_shape_size_2(scale / FT(2)),
    grid_cell_width_2(scale / FT(4)),
    region_growing_search_size_2(FT(12)),
    region_growing_noise_level_2(noise_level),
    region_growing_angle_2(FT(25)),
    region_growing_minimum_length_2(scale)
    { }

    // Update all parameters, which depend on scale and noise_level.
    void update_dependent() {

      alpha_shape_size_2 = scale / FT(2);
      grid_cell_width_2 = scale / FT(4);
        
      region_growing_noise_level_2 = noise_level;
      region_growing_minimum_length_2 = scale;
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
        "region_growing_search_size_2 (meters / number of points) : " 
      << region_growing_search_size_2 << std::endl;
      file << 
        "region_growing_noise_level_2 (meters) : " 
      << region_growing_noise_level_2 << std::endl;
      file << 
        "region_growing_angle_2 (degrees) : " 
      << region_growing_angle_2 << std::endl;
      file << 
        "region_growing_minimum_length_2 (meters) : " 
      << region_growing_minimum_length_2 << std::endl;

      file.close();
    }

  }; // Parameters

} // Levels_of_detail
} // CGAL

#endif // CGAL_LOD_PARAMETERS_H
