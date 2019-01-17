#ifndef CGAL_LOD_SAVER_H
#define CGAL_LOD_SAVER_H

#if defined(WIN32) || defined(_WIN32)
#define _NL_ "\r\n"
#else
#define _NL_ "\n"
#endif

// STL includes.
#include <vector>
#include <string>
#include <sstream>
#include <fstream>

// CGAL includes.
#include <CGAL/IO/Color.h>
#include <CGAL/Point_set_3.h>
#include <CGAL/Point_set_3/IO.h>

namespace CGAL {

  namespace Levels_of_detail {

    template<typename GeometricTraits>
    class Saver {

    public:
      
      using Kernel = GeometricTraits;
      using Point_3 = typename Kernel::Point_3;
      using Point_set = Point_set_3<Point_3>;

      using Color_map = typename Point_set::template Property_map<unsigned char>;

      Saver() { 
        out.precision(20); 
      }

      void clear() { 
        out.str(std::string()); 
      }

      void export_planar_ground(
        const std::vector<Point_3> &vertices,
        const std::string file_path) {

        clear();

        const std::size_t num_vertices = vertices.size();
        const std::size_t num_faces = 1;

        add_ply_header(num_vertices, num_faces);

        for (std::size_t i = 0; i < num_vertices; ++i)
          out << vertices[i] << std::endl;
        
        out << num_vertices << " ";
        for (std::size_t i = 0; i < num_vertices; ++i)
          out << i << " ";
        out << Color(128, 64, 0) << std::endl;

        save(file_path + ".ply");
      }

      void export_point_set(
        const Point_set &point_set,
        const std::string file_path) {

        clear();

        const std::size_t num_points = point_set.size();

        add_ply_header(num_points);

        const bool use_color =
        point_set. template property_map<unsigned char>("r").second;
        Color_map r, g, b;

        if (use_color) {
          r = point_set. template property_map<unsigned char>("r").first;
          g = point_set. template property_map<unsigned char>("g").first;
          b = point_set. template property_map<unsigned char>("b").first;
        }

        for (auto it = point_set.begin(); it != point_set.end(); ++it) {

          std::string color = "";
          if (use_color) {

            color = 
            std::to_string(get(r, *it)) + " " +
            std::to_string(get(g, *it)) + " " +
            std::to_string(get(b, *it));

          } else color = "0 0 0";

          out << get(point_set.point_map(), *it) << " " 
          << color << std::endl;
          
        }
        save(file_path + ".ply");
      }

    private:
      std::stringstream out;

      inline std::string data() const { 
        return out.str(); 
      }

      void save(const std::string file_path) const {
        std::ofstream file(file_path.c_str(), std::ios_base::out);

        if (!file)
          std::cerr << std::endl
                    << "ERROR: Error saving file " << file_path << std::endl
                    << std::endl;

        file << data();
        file.close();
      }

      void add_ply_header(const std::size_t num_points) {

        out << 
				"ply" 				         +  std::string(_NL_) + ""               			<< 
				"format ascii 1.0"     +  std::string(_NL_) + ""     			          << 
				"element vertex "      << num_points       << "" + std::string(_NL_) + "" << 
				"property double x"    +  std::string(_NL_) + ""    			          << 
				"property double y"    +  std::string(_NL_) + ""    			          << 
				"property double z"    +  std::string(_NL_) + "" 				            <<
				"property uchar red"   +  std::string(_NL_) + "" 				            <<
				"property uchar green" +  std::string(_NL_) + "" 				            <<
				"property uchar blue"  +  std::string(_NL_) + "" 				            <<
				"end_header"           +  std::string(_NL_) + "";
      }

      void add_ply_header(
        const std::size_t num_vertices, 
        const std::size_t num_faces) {

        out << 
				"ply" 				         +  std::string(_NL_) + ""               			<< 
				"format ascii 1.0"     +  std::string(_NL_) + ""     			          << 
				"element vertex "      << num_vertices     << "" + std::string(_NL_) + "" << 
				"property double x"    +  std::string(_NL_) + ""    			          << 
				"property double y"    +  std::string(_NL_) + ""    			          << 
				"property double z"    +  std::string(_NL_) + "" 				            <<
				"element face "        << num_faces        << "" + std::string(_NL_) + "" << 
				"property list uchar int vertex_indices"         + std::string(_NL_) + "" <<
				"property uchar red"   +  std::string(_NL_) + "" 				            <<
				"property uchar green" +  std::string(_NL_) + "" 				            <<
				"property uchar blue"  +  std::string(_NL_) + "" 				            <<
				"end_header"           +  std::string(_NL_) + "";
      }
      
    }; // Saver

  } // namespace Levels_of_detail

} // namespace CGAL

#endif // CGAL_LOD_SAVER_H
