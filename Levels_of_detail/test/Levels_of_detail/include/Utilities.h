#ifndef CGAL_LOD_UTILITIES_H
#define CGAL_LOD_UTILITIES_H

// STL includes.
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

// LOD includes.
#include <CGAL/Levels_of_detail/enumerations.h>

namespace CGAL {
namespace Levels_of_detail {

  template<typename LabelMap>
  struct Semantic_from_label_map {

  public:
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
        return Semantic_label::UNASSIGNED;

      return it->second;
    }

  }; // Semantic_from_label_map

  template<typename GeomTraits>
  struct Insert_point_colored_by_index {
      
  public:
    using Traits = GeomTraits;
    using Point_3 = typename Traits::Point_3;
      
    using argument_type = std::pair<Point_3, long>;
    using result_type = void;

    using Point_set = Point_set_3<Point_3>;
    using Color_map = typename Point_set:: template Property_map<unsigned char>;

    Point_set& m_point_set;
    Color_map m_red, m_green, m_blue;

    Insert_point_colored_by_index(Point_set& point_set) : 
    m_point_set(point_set) {
        
      m_red = 
      m_point_set.template add_property_map<unsigned char>("r", 0).first;
        
      m_green = 
      m_point_set.template add_property_map<unsigned char>("g", 0).first;
        
      m_blue = 
      m_point_set.template add_property_map<unsigned char>("b", 0).first;
    }

    void operator()(const argument_type& arg) {
        
      const auto it = m_point_set.insert(arg.first);
      if (arg.second == -1) 
        return;

      Random rand(arg.second);

      m_red[*it] = static_cast<unsigned char>(64 + rand.get_int(0, 192));
      m_green[*it] = static_cast<unsigned char>(64 + rand.get_int(0, 192));
      m_blue[*it] = static_cast<unsigned char>(64 + rand.get_int(0, 192));
    }
  }; // Insert_point_colored_by_index

  template<typename GeomTraits>
  struct Add_polyline_from_segment {

  public:
    using Traits = GeomTraits;
    using Point_3 = typename Traits::Point_3;

    using argument_type = typename Traits::Segment_3;
    using result_type = void;

    using Polyline = std::vector<Point_3>;
    using Polylines = std::vector<Polyline>;

    Polylines &m_polylines;

    Add_polyline_from_segment(Polylines& polylines) : 
    m_polylines(polylines) 
    { }

    void operator()(const argument_type& segment) {
      
      m_polylines.push_back(Polyline());
      m_polylines.back().push_back(segment.source());
      m_polylines.back().push_back(segment.target());
    }
  }; // Add_polyline_from_segment

  struct Add_polygon_with_color {
  
  public:
    using Color = CGAL::Color;
    using Indices = std::vector<std::size_t>;
    using Visibility_label = CGAL::Levels_of_detail::Visibility_label;

    using argument_type = std::pair<Indices, Visibility_label>;
    using result_type = void;

    std::vector<Indices>& m_polygons;
    std::vector<Color>& m_colors;
    const bool m_use_visibility;

    Add_polygon_with_color(
      std::vector<Indices>& polygons,
      std::vector<Color>& colors,
      const bool use_visibility) : 
    m_polygons(polygons), 
    m_colors(colors),
    m_use_visibility(use_visibility) 
    { }

    result_type operator()(const argument_type& arg) {
      
      m_polygons.push_back(arg.first);
      if (!m_use_visibility) {
        
        m_colors.push_back(Color(77, 77, 255));
        return;
      }

      unsigned char r, g, b;
      switch (arg.second) {

        case Visibility_label::OUTSIDE: {
          r = 255; g = 77; b = 77;
          break;
        }

        case Visibility_label::INSIDE: {
          r = 77; g = 255; b = 77;
          break;
        }

        default: {
          r = 0; g = 0; b = 0;
          break;
        }
      }
      m_colors.push_back(Color(r, g, b));
    }
  }; // Add_polygon_with_color

  struct Add_triangle_with_color {
    
  public:
    using Color = CGAL::Color;
    using Indices = std::vector<std::size_t>;

    using argument_type = 
    std::pair<CGAL::cpp11::array<std::size_t, 3>, std::size_t>;
    using result_type = void;

    std::vector<Indices>& m_triangles;
    std::vector<Color>& m_colors;

    Add_triangle_with_color(
      std::vector<Indices>& triangles,
      std::vector<Color>& colors) : 
    m_triangles(triangles), 
    m_colors(colors) 
    { }

    result_type operator()(const argument_type& arg) {
      
      m_triangles.push_back(std::vector<std::size_t>(3));
      for (std::size_t i = 0; i < 3; ++i)
        m_triangles.back()[i] = arg.first[i];

      unsigned char r, g, b;
      CGAL::Random rand(arg.second);

      r = static_cast<unsigned char>(64 + rand.get_int(0, 192));
      g = static_cast<unsigned char>(64 + rand.get_int(0, 192));
      b = static_cast<unsigned char>(64 + rand.get_int(0, 192));

      m_colors.push_back(Color(r, g, b));
    }
  }; // Add_triangle_with_color

} // Levels_of_detail
} // CGAL

#endif // CGAL_LOD_UTILITIES_H
