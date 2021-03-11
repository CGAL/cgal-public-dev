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
#include <CGAL/Random.h>
#include <CGAL/IO/Color.h>
#include <CGAL/Point_set_3.h>
#include <CGAL/Point_set_3/IO.h>

namespace CGAL {
namespace Levels_of_detail {

  template<typename GeomTraits>
  class Saver {

  public:
    using Traits = GeomTraits;
    using FT = typename Traits::FT;
    using Point_2 = typename Traits::Point_2;
    using Point_3 = typename Traits::Point_3;
    using Vector_3 = typename Traits::Vector_3;
    using Point_set = Point_set_3<Point_3>;
    using Points = std::vector<Point_3>;
    using Polylines = std::vector<Points>;
    using Segment_2 = typename Traits::Segment_2;

    using Color_map = typename Point_set::template Property_map<unsigned char>;

    Saver() {
      out.precision(20);
    }

    void clear() {
      out.str(std::string());
    }

    void export_points(
      const std::vector< std::vector<Point_3> >& points,
      const std::string file_path) {

      if (points.size() == 0)
        return;

      clear();
      std::size_t num_points = 0;
      for (const auto& item : points)
        num_points += item.size();
      add_ply_header_points(num_points);

      for (std::size_t i = 0; i < points.size(); ++i) {
        CGAL::Random rnd(i);
        const auto r = 64 + rnd.get_int(0, 192);
        const auto g = 64 + rnd.get_int(0, 192);
        const auto b = 64 + rnd.get_int(0, 192);

        for (const auto& p : points[i])
          out << p << " " << r << " " << g << " " << b << std::endl;
      }
      save(file_path + ".ply");
    }

    void export_points(
      const std::vector<Point_3>& points,
      const std::vector< std::vector<std::size_t> >& regions,
      const std::string file_path) {

      if (points.size() == 0)
        return;

      clear();
      const std::size_t num_points = points.size();
      add_ply_header_points(num_points);

      for (std::size_t i = 0; i < regions.size(); ++i) {
        CGAL::Random rnd(i);
        const auto r = 64 + rnd.get_int(0, 192);
        const auto g = 64 + rnd.get_int(0, 192);
        const auto b = 64 + rnd.get_int(0, 192);

        for (const std::size_t idx : regions[i]) {
          const auto& p = points[idx];
          out << p << " " << r << " " << g << " " << b << std::endl;
        }
      }
      save(file_path + ".ply");
    }

    template<typename Color>
    void export_points(
      const std::vector<Point_3>& points,
      const Color color,
      const std::string file_path) {

      if (points.size() == 0)
        return;

      clear();
      const std::size_t num_points = points.size();
      add_ply_header_points(num_points);

      for (const auto& p : points)
        out << p << " " << color << std::endl;
      save(file_path + ".ply");
    }

    void export_points_with_normals(
      const std::vector<Point_3>& points,
      const std::vector<Vector_3>& normals,
      const std::string file_path) {

      if (points.size() == 0)
        return;

      clear();
      const std::size_t size = points.size();
      add_ply_header_normals(size);

      for (std::size_t i = 0; i < size; ++i)
        out << points[i] << " " << normals[i] << std::endl;
      save(file_path + ".ply");
    }

    void export_point_set(
      const Point_set& point_set,
      const std::string file_path) {

      if (point_set.size() == 0)
        return;

      clear();
      out << point_set;
      save(file_path + ".ply");
    }

    void export_polylines(
      const Polylines& polylines,
      const std::string file_path) {

      if (polylines.size() == 0)
        return;

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

    template<typename Color>
    void export_polygon_soup(
      const std::vector< std::vector<Point_3> >& polygons,
      const Color color,
      const std::string file_path) {

      clear();
      std::size_t num_vertices = 0;
      for (const auto& polygon : polygons)
        num_vertices += polygon.size();
      std::size_t num_faces = polygons.size();
      add_ply_header_mesh(num_vertices, num_faces);

      for (const auto& polygon : polygons)
        for (const auto& p : polygon)
          out << p << std::endl;

      std::size_t i = 0;
      for (const auto& polygon : polygons) {
        out << polygon.size() << " ";
        for (const auto& p : polygon)
          out << i++ << " ";
        out << color << std::endl;
      }

      save(file_path + ".ply");
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

      if (vertices.size() == 0 || faces.size() == 0)
        return;

      clear();
      std::size_t num_vertices = vertices.size();
      std::size_t num_faces = faces.size();
      add_ply_header_mesh(num_vertices, num_faces);

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

    void save_polylines(
      const std::vector<Segment_2>& segments,
      const std::string name) {

      std::vector< std::vector<Point_3> > polylines(segments.size());
      for (std::size_t i = 0; i < segments.size(); ++i) {
        const Point_2& s = segments[i].source();
        const Point_2& t = segments[i].target();

        polylines[i].push_back(Point_3(s.x(), s.y(), FT(0)));
        polylines[i].push_back(Point_3(t.x(), t.y(), FT(0)));
      }

      export_polylines(polylines, name);
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

    void add_ply_header_points(
      const std::size_t size) {

      out <<
			"ply" 				         +  std::string(_NL_) + ""               			 <<
			"format ascii 1.0"     +  std::string(_NL_) + ""     			           <<
			"element vertex "      << size        << "" + std::string(_NL_) + "" <<
			"property double x"    +  std::string(_NL_) + ""    			           <<
			"property double y"    +  std::string(_NL_) + ""    			           <<
			"property double z"    +  std::string(_NL_) + "" 				             <<
			"property uchar red"   +  std::string(_NL_) + "" 				             <<
			"property uchar green" +  std::string(_NL_) + "" 				             <<
			"property uchar blue"  +  std::string(_NL_) + "" 				             <<
			"end_header"           +  std::string(_NL_) + "";
    }

    void add_ply_header_normals(
      const std::size_t size) {

      out <<
			"ply" 				         +  std::string(_NL_) + ""               			 <<
			"format ascii 1.0"     +  std::string(_NL_) + ""     			           <<
			"element vertex "      << size        << "" + std::string(_NL_) + "" <<
			"property double x"    +  std::string(_NL_) + ""    			           <<
			"property double y"    +  std::string(_NL_) + ""    			           <<
			"property double z"    +  std::string(_NL_) + "" 				             <<
      "property double nx"   +  std::string(_NL_) + ""    			           <<
			"property double ny"   +  std::string(_NL_) + ""    			           <<
			"property double nz"   +  std::string(_NL_) + "" 				             <<
			"end_header"           +  std::string(_NL_) + "";
    }

    void add_ply_header_mesh(
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
