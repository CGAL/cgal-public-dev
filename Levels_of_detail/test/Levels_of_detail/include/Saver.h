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

  template<typename GeomTraits>
  class Saver {

  public:
    using Traits = GeomTraits;
    using Point_3 = typename Traits::Point_3;
    using Point_set = Point_set_3<Point_3>;
    using Points = std::vector<Point_3>;
    using Polylines = std::vector<Points>;

    using Color_map = typename Point_set::template Property_map<unsigned char>;

    Saver() { 
      out.precision(20); 
    }

    void clear() { 
      out.str(std::string()); 
    }

    void export_planar_ground(
      const Points& points,
      const std::string file_path) {

      clear();
      const std::size_t num_points = points.size();
      const std::size_t num_faces = 1;
      add_ply_header(num_points, num_faces);

      for (std::size_t i = 0; i < num_points; ++i)
        out << points[i] << std::endl;
        
      out << num_points << " ";
      for (std::size_t i = 0; i < num_points; ++i)
        out << i << " ";
      out << Color(128, 64, 0) << std::endl;

      save(file_path + ".ply");
    }

    void export_point_set(
      const Point_set& point_set,
      const std::string file_path) {

      clear();
      out << point_set;
      save(file_path + ".ply");
    }

    void export_polylines(
      const Polylines& polylines,
      const std::string file_path) {

      clear();

      for (std::size_t i = 0; i < polylines.size(); ++i) {
        const auto &polyline = polylines[i];

        out << polyline.size() << " ";
        for (std::size_t j = 0; j < polyline.size(); ++j)
          out << polyline[j] << " ";
        out << std::endl;
      }
      save(file_path + ".polylines");
    }

    template<
    typename Vertices, 
    typename Faces, 
    typename Colors>
    void export_polygon_soup(
      const Vertices& vertices, 
      const Faces& faces,
      const Colors& fcolors, 
      const std::string file_path) {

      clear();
      std::size_t num_vertices = vertices.size();
      std::size_t num_faces = faces.size();
      add_ply_header(num_vertices, num_faces);

      for (std::size_t i = 0; i < vertices.size(); ++i)
        out << vertices[i] << std::endl;

      for (std::size_t i = 0; i < faces.size(); ++i) {
        out << faces[i].size() << " ";

        for (std::size_t j = 0; j < faces[i].size(); ++j)
          out << faces[i][j] << " ";
        out << fcolors[i] << std::endl;
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
        std::cerr << std::endl << 
          "ERROR: Error saving file " << file_path 
        << std::endl << std::endl;

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
