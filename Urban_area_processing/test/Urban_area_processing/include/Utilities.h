#ifndef CGAL_UAP_UTILITIES_H
#define CGAL_UAP_UTILITIES_H

// STL includes.
#include <tuple>
#include <string>
#include <sstream>
#include <utility>
#include <unordered_map>

// CGAL includes.
#include <CGAL/array.h>
#include <CGAL/Random.h>
#include <CGAL/IO/Color.h>
#include <CGAL/property_map.h>
#include <CGAL/Point_set_3.h>

// UAP includes.
#include <CGAL/Urban_area_processing/enum.h>

// Colors:
// 255 102 51  - roofs
// 0 143 0     - crowns
// 127 105 99  - ground
// 148 82 0    - trunks
// 214 214 214 - walls

namespace CGAL {
namespace Urban_area_processing {

  void merge_points(
    const std::vector<std::size_t>& a,
    const std::vector<std::size_t>& b,
    std::vector<std::size_t>& c) {

    c.clear();
    c.reserve(a.size() + b.size());
    for (const std::size_t idx : a)
      c.push_back(idx);
    for (const std::size_t idx : b)
      c.push_back(idx);
  }

  template<
  typename Point_set,
  typename Semantic_map>
  void get_labeled_points(
    const Point_set& point_set,
    const Semantic_map& semantic_map,
    const Semantic_label output_label,
    std::vector<std::size_t>& indices) {

    indices.clear();
    for (std::size_t i = 0; i < point_set.size(); ++i) {
      const Semantic_label label = 
      get(semantic_map, *(point_set.begin() + i));
      if (label == output_label)
        indices.push_back(i);
    }
  }

  template<typename Point_set>
  bool are_label_data_defined(
    const Point_set& point_set) {
    return point_set. template property_map<int>("label").second;
  }

  template<
  typename Point_set,
  typename Label_map>
  void load_input_data(
    const std::string path,
    Point_set& point_set,
    Label_map& label_map) {
    
    std::cout << std::endl << "Input data: " << std::endl;
    std::ifstream file(path.c_str(), std::ios_base::binary);
    file >> point_set; 
    file.close();
    std::cout << "File contains " << point_set.size() 
    << " points" << std::endl;

    if (are_label_data_defined(point_set)) {
      std::cout << 
        "Label data are defined!" 
      << std::endl << std::endl;
      label_map = point_set. template property_map<int>("label").first;
    } else {
      std::cerr << 
        "Label data are not defined!" 
      << std::endl << std::endl;
      exit(EXIT_FAILURE);
    }
  }

  template<typename LabelMap>
  struct Semantic_from_label_map {
    using Label_map = LabelMap;

    using key_type = typename boost::property_traits<Label_map>::key_type;
    using value_type = Semantic_label;
    using reference = value_type;
    using category = boost::readable_property_map_tag;

    using Label_to_semantic_map = std::unordered_map<int, Semantic_label>;

    Label_map m_label_map;
    Label_to_semantic_map m_label_to_semantic_map;

    Semantic_from_label_map() { }

    Semantic_from_label_map(
      Label_map label_map, 
      const std::string gi_str,
      const std::string bi_str,
      const std::string ii_str,
      const std::string vi_str) : 
    m_label_map(label_map) { 
      std::cout << "Setting semantic labels:" << std::endl;

      std::istringstream gi(gi_str);
      std::istringstream bi(bi_str);
      std::istringstream ii(ii_str);
      std::istringstream vi(vi_str);

      int idx;
      while (gi >> idx) {
        std::cout << idx << " is ground" << std::endl;
        m_label_to_semantic_map.insert(
          std::make_pair(idx, Semantic_label::GROUND));
      }
      while (bi >> idx) {
        std::cout << idx << " is building boundary" << std::endl;
        m_label_to_semantic_map.insert(
          std::make_pair(idx, Semantic_label::BUILDING_BOUNDARY));
      }
      while (ii >> idx) {
        std::cout << idx << " is building interior" << std::endl;
        m_label_to_semantic_map.insert(
          std::make_pair(idx, Semantic_label::BUILDING_INTERIOR));
      }
      while (vi >> idx) {
        std::cout << idx << " is vegetation" << std::endl;
        m_label_to_semantic_map.insert(
          std::make_pair(idx, Semantic_label::VEGETATION));
      }
      std::cout << std::endl;
    }

    friend value_type get(
      const Semantic_from_label_map& semantic_map, 
      const key_type& key) {

      const int label = get(semantic_map.m_label_map, key);
      const auto it = semantic_map.m_label_to_semantic_map.find(label);

      if (it == semantic_map.m_label_to_semantic_map.end())
        return Semantic_label::UNCLASSIFIED;
      return it->second;
    }
  }; // Semantic_from_label_map

  template<typename Point_set>
  void get_building_points(
    const std::string in_path,
    const std::string gi,
    const std::string bi,
    const std::string ii,
    const std::string vi,
    const std::string out_path,
    Point_set& point_set,
    std::vector<std::size_t>& building_points) {

    using Label_map = 
      typename Point_set:: template Property_map<int>;
    using Semantic_map = Semantic_from_label_map<Label_map>;

    Label_map label_map;
    load_input_data(in_path, point_set, label_map);
    Semantic_map semantic_map(label_map, gi, bi, ii, vi);

    std::vector<std::size_t> interior_points;
    std::vector<std::size_t> boundary_points;
    get_labeled_points(point_set, semantic_map, 
      Semantic_label::BUILDING_INTERIOR, interior_points);
    get_labeled_points(point_set, semantic_map, 
      Semantic_label::BUILDING_BOUNDARY, boundary_points);
    
    building_points.clear();
    merge_points(interior_points, boundary_points, building_points);
  }

  template<typename GeomTraits>
  struct Point_inserter {
    using Traits = GeomTraits;
    using Point_3 = typename Traits::Point_3;
      
    using argument_type0 = std::pair<Point_3, std::size_t>;
    using argument_type1 = std::pair<Point_3, Semantic_label>;
    using result_type = void;

    using Point_set = Point_set_3<Point_3>;
    using Color_map = typename Point_set:: template Property_map<unsigned char>;

    Point_set& m_point_set;
    Color_map m_red, m_green, m_blue;
    Point_inserter(Point_set& point_set) : 
    m_point_set(point_set) {
      m_red = 
      m_point_set.template add_property_map<unsigned char>("r", 0).first;
      m_green = 
      m_point_set.template add_property_map<unsigned char>("g", 0).first;
      m_blue = 
      m_point_set.template add_property_map<unsigned char>("b", 0).first;
    }

    void operator()(const argument_type0& arg) {
      const auto it = m_point_set.insert(arg.first);
      if (arg.second == std::size_t(-1)) 
        return;
      Random rand(arg.second);
      m_red[*it] = static_cast<unsigned char>(64 + rand.get_int(0, 192));
      m_green[*it] = static_cast<unsigned char>(64 + rand.get_int(0, 192));
      m_blue[*it] = static_cast<unsigned char>(64 + rand.get_int(0, 192));
    }

    void operator()(const argument_type1& arg) {
      const auto it = m_point_set.insert(arg.first);
      unsigned char r, g, b;
      switch (arg.second) {
        case Semantic_label::GROUND: {
          r = 127; g = 105; b = 99; break; }
        case Semantic_label::BUILDING_BOUNDARY: {
          r = 214; g = 214; b = 214; break; }
        case Semantic_label::BUILDING_INTERIOR: {
          r = 255; g = 102; b = 51; break; }
        case Semantic_label::VEGETATION: {
          r = 0; g = 143; b = 0; break; }
        default: {
          r = 0; g = 0; b = 0; break; }
      }
      m_red[*it] = r; m_green[*it] = g; m_blue[*it] = b;
    }
  }; // Point_inserter

  template<typename GeomTraits>
  struct Polyline_inserter {
    using Traits = GeomTraits;
    using Point_3 = typename Traits::Point_3;
    using Segment_3 = typename Traits::Segment_3;

    using argument_type0 = Segment_3;
    using argument_type1 = std::pair<Segment_3, std::size_t>;
    using result_type = void;

    using Polyline = std::vector<Point_3>;
    using Polylines = std::vector<Polyline>;

    Polylines &m_polylines;
    Polyline_inserter(Polylines& polylines) : 
    m_polylines(polylines) 
    { }

    void add_segment(const Segment_3& segment) {
      m_polylines.push_back(Polyline());
      m_polylines.back().push_back(segment.source());
      m_polylines.back().push_back(segment.target());
    }

    void operator()(const argument_type0& arg) {
      add_segment(arg);
    }

    void operator()(const argument_type1& arg) {
      add_segment(arg.first);
    }
  }; // Polyline_inserter

  template<typename GeomTraits>
  struct Polygon_inserter {
    using Traits = GeomTraits;
    using Color = CGAL::Color;
    using Indices = std::vector<std::size_t>;
    using Visibility_label = CGAL::Urban_area_processing::Visibility_label;

    using argument_type0 = std::pair<Indices, Visibility_label>;
    using argument_type1 = std::pair<Indices, Urban_object_type>;
    using argument_type2 = std::pair<Indices, std::size_t>;
    
    using result_type = void;

    std::vector<Indices>& m_polygons;
    std::vector<Color>& m_colors;
    Polygon_inserter(
      std::vector<Indices>& polygons,
      std::vector<Color>& colors) : 
    m_polygons(polygons), 
    m_colors(colors)
    { }

    result_type operator()(const argument_type0& arg) {
      m_polygons.push_back(arg.first);
      unsigned char r, g, b;
      switch (arg.second) {
        case Visibility_label::OUTSIDE: {
          r = 255; g = 77; b = 77; break; }
        case Visibility_label::INSIDE: {
          r = 77; g = 255; b = 77; break; }
        default: {
          r = 0; g = 0; b = 0; break; }
      }
      m_colors.push_back(Color(r, g, b));
    }

    result_type operator()(const argument_type1& arg) {
      m_polygons.push_back(arg.first);
      CGAL::Color color;
      switch (arg.second) {
        case Urban_object_type::GROUND:
        color = CGAL::Color(127, 105, 99);  
        break;
        case Urban_object_type::BUILDING_WALL:
        color = CGAL::Color(214, 214, 214);  
        break;
        case Urban_object_type::BUILDING_ROOF:
        color = CGAL::Color(255, 102, 51);  
        break;
        case Urban_object_type::TREE_TRUNK:
        color = CGAL::Color(148, 82, 0);
        break;
        case Urban_object_type::TREE_CROWN:
        color = CGAL::Color(0, 143, 0);
        break;
        default:
        color = CGAL::Color(0, 0, 0);  
        break;
      }
      m_colors.push_back(color);
    }

    result_type operator()(const argument_type2& arg) {
      m_polygons.push_back(arg.first);

      if (arg.second == std::size_t(-1)) {
        m_colors.push_back(Color(255, 77, 77));
        return;
      }

      Random rand(arg.second);
      const auto r = static_cast<unsigned char>(64 + rand.get_int(0, 192));
      const auto g = static_cast<unsigned char>(64 + rand.get_int(0, 192));
      const auto b = static_cast<unsigned char>(64 + rand.get_int(0, 192));
      m_colors.push_back(Color(r, g, b));
    }
  }; // Polygon_inserter

} // Urban_area_processing
} // CGAL

#endif // CGAL_UAP_UTILITIES_H
