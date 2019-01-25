#ifndef CGAL_LOD_PARAMETERS_H
#define CGAL_LOD_PARAMETERS_H

// STL includes.
#include <string>

namespace CGAL {
namespace Levels_of_detail {  

  template<typename FT>
  struct Parameters {

  public:

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

    // Alpha shape size in meters.
    FT alpha_shape_size;

    // Grid cell width in meters.
    FT grid_cell_width;

    // Region growing.
    FT region_growing_scale; // meters
    FT region_growing_noise_level; // meters
    FT region_growing_normal_threshold; // degrees
    FT region_growing_minimum_length; // meters


    // Constructor.
    Parameters() : 
    data(""),
    gi("0"), bi("1"), ii("2"), vi("3"),
    scale(FT(4)),
    noise_level(FT(2)),
    alpha_shape_size(scale / FT(2)),
    grid_cell_width(scale / FT(4)),
    region_growing_scale(scale),
    region_growing_noise_level(noise_level),
    region_growing_normal_threshold(FT(25)),
    region_growing_minimum_length(scale)
    { }

    // Update all parameters, which depend on scale and noise_level.
    void update_dependent() {

      alpha_shape_size = scale / FT(2);
      grid_cell_width = scale / FT(4);
        
      region_growing_scale = scale;
      region_growing_noise_level = noise_level;
      region_growing_minimum_length = scale;
    }

    // Set main parameters.
    void set_main(const FT scale_, const FT noise_level_) {

      scale = scale_;
      noise_level = noise_level_;
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
      file << "alpha_shape_size (meters) : " << alpha_shape_size << std::endl;
      file << "grid_cell_width (meters) : " << grid_cell_width << std::endl;
      file << "region_growing_scale (meters) : " << region_growing_scale << std::endl;
      file << "region_growing_noise_level (meters) : " << region_growing_noise_level << std::endl;
      file << "region_growing_normal_threshold (degrees) : " << region_growing_normal_threshold << std::endl;
      file << "region_growing_minimum_length (meters) : " << region_growing_minimum_length << std::endl;

      file.close();
    }

  }; // Parameters

} // Levels_of_detail
} // CGAL

#endif // CGAL_LOD_PARAMETERS_H
