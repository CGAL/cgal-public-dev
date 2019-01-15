#ifndef CGAL_LOD_PARAMETERS_H
#define CGAL_LOD_PARAMETERS_H

// STL includes.
#include <string>

namespace CGAL {

  namespace Levels_of_detail {  

    template<typename FT>
    struct Parameters {

    public:

      // Set to true if you want to print extra information.
      bool verbose;

      // Path to the input data file.
      std::string data;
      
      // Scale in meters used to detect LOD.
      FT scale;

      // Extent in meters to each found line or plane, in other words
      // this is a wall width.
      FT extent;

      // Label indices defined in the ply header: 
      // ground, building boundary, building interior, vegetation.
      std::string gi, bi, ii, vi;

      Parameters() : 
      verbose(true),
      data(""),
      scale(FT(4)),
      extent(FT(2)),
      gi("0"), bi("1"), ii("2"), vi("3") 
      { }

      // Update all parameters, which depend on scale and extent.
      void update_dependent() {

      }

      // Set main parameters.
      void set_main(const FT scale_, const FT extent_) {

        scale  = scale_;
        extent = extent_;
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
        file << "-data : " << data << std::endl;
        file << std::endl;

        file << "Label indices: " << std::endl;
        file << "-gi (ground) : " << gi << std::endl;
        file << "-bi (building boundary) : " << bi << std::endl;
        file << "-ii (building interior) : " << ii << std::endl;
        file << "-vi (vegetation) : " << vi << std::endl;
        file << std::endl;

        file << "Main parameters: " << std::endl;
        file << "-scale : " << scale << std::endl;
        file << "-extent : " << extent << std::endl;
        file << std::endl;

        file << "Info: " << std::endl;
        file << "-verbose : " << verbose << std::endl;

        file.close();
      }

    }; // Parameters

  } // Levels_of_detail

} // CGAL

#endif // CGAL_LOD_PARAMETERS_H
