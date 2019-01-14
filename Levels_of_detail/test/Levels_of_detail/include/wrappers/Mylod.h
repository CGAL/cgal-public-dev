#ifndef CGAL_LEVELS_OF_DETAIL_MYLOD_H
#define CGAL_LEVELS_OF_DETAIL_MYLOD_H

#if defined(WIN32) || defined(_WIN32)
#define _SR_ "\\"
#else
#define _SR_ "/"
#endif

// STL includes.
#include <iostream>
#include <string>

// Boost includes.
#include <boost/function_output_iterator.hpp>

// CGAL includes.
#include <CGAL/Point_set_3.h>
#include <CGAL/Random.h>

// LOD includes.
#include <CGAL/Levels_of_detail.h>

// Local includes.
#include "../debugging/Mylog.h"
#include "../terminal/Myterminal_parser.h"

namespace CGAL {

namespace Levels_of_detail {

template <class GeomTraits> class Mylod {

public:
  using Kernel = GeomTraits;
  using FT = typename Kernel::FT;
  using Point_3 = typename Kernel::Point_3;

  using Input_parameters = char **;

  using Log = Mylog;
  using Terminal_parser = Myterminal_parser<FT>;

  using Point_set = CGAL::Point_set_3<Point_3>;
  using Point_map = typename Point_set::Point_map;

  typedef std::map<int, CGAL::Level_of_detail::Semantic_label> Map_l2sl;
  typedef boost::shared_ptr<Map_l2sl> Map_l2sl_ptr;

  struct Semantic_map_from_labels {
    typedef typename Point_set::Index key_type;
    typedef typename CGAL::Level_of_detail::Semantic_label value_type;
    typedef typename CGAL::Level_of_detail::Semantic_label reference;
    typedef boost::readable_property_map_tag category;

    Point_set *points;
    typename Point_set::template Property_map<int> label_map;
    Map_l2sl_ptr map_l2sl;

    Semantic_map_from_labels() {}
    Semantic_map_from_labels(Point_set *points)
        : points(points), map_l2sl(new Map_l2sl()) {
      label_map = points->template property_map<int>("label").first;
    }

    friend value_type get(const Semantic_map_from_labels &map,
                          const key_type &key) {
      int l = map.label_map[key];

      typename Map_l2sl::const_iterator found = map.map_l2sl->find(l);
      if (found == map.map_l2sl->end())
        return CGAL::Level_of_detail::Semantic_label::UNASSIGNED;

      return found->second;
    }
  };

  struct Insert_point_colored_by_index {
    typedef std::pair<Point_3, int> argument_type;
    typedef void result_type;

    Point_set &points;
    typename Point_set::template Property_map<unsigned char> red, green, blue;

    Insert_point_colored_by_index(Point_set &points) : points(points) {
      red = points.template add_property_map<unsigned char>("r", 0).first;
      green = points.template add_property_map<unsigned char>("g", 0).first;
      blue = points.template add_property_map<unsigned char>("b", 0).first;
      // points.check_colors();
    }

    void operator()(const argument_type &a) {
      typename Point_set::iterator pt = points.insert(a.first);
      if (a.second == -1)
        return;

      CGAL::Random rand(a.second);
      red[*pt] = static_cast<unsigned char>(64 + rand.get_int(0, 192));
      green[*pt] = static_cast<unsigned char>(64 + rand.get_int(0, 192));
      blue[*pt] = static_cast<unsigned char>(64 + rand.get_int(0, 192));
    }
  };

  typedef std::vector<Point_3> Polyline;
  typedef std::list<Polyline> Polylines_container;

  struct Add_polyline_from_segment {
    typedef typename Kernel::Segment_3 argument_type;
    typedef void result_type;

    Polylines_container &polylines;

    Add_polyline_from_segment(Polylines_container &polylines)
        : polylines(polylines) {}

    void operator()(const argument_type &a) {
      polylines.push_back(Polyline());
      polylines.back().push_back(a.source());
      polylines.back().push_back(a.target());
    }
  };

  struct Add_polygon_with_visibility_color {
    typedef std::pair<std::vector<std::size_t>,
                      CGAL::Level_of_detail::Visibility_label>
        argument_type;
    typedef void result_type;

    std::vector<std::vector<std::size_t>> &polygons;
    std::vector<CGAL::Color> &fcolors;

    Add_polygon_with_visibility_color(
        std::vector<std::vector<std::size_t>> &polygons,
        std::vector<CGAL::Color> &fcolors)
        : polygons(polygons), fcolors(fcolors) {}

    void operator()(const argument_type &a) {
      polygons.push_back(a.first);

      unsigned char r, g, b;
      if (a.second == CGAL::Level_of_detail::Visibility_label::OUTSIDE) {
        r = 186;
        g = 189;
        b = 182;
      } else if (a.second == CGAL::Level_of_detail::Visibility_label::INSIDE) {
        r = 245;
        g = 121;
        b = 0;
      } else {
        r = 77;
        g = 131;
        b = 186;
      }

      fcolors.push_back(CGAL::Color(r, g, b));
    }
  };

  struct Add_triangle_with_building_color {
    typedef std::pair<CGAL::cpp11::array<std::size_t, 3>, int> argument_type;
    typedef void result_type;

    std::vector<std::vector<std::size_t>> &polygons;
    std::vector<CGAL::Color> &fcolors;

    Add_triangle_with_building_color(
        std::vector<std::vector<std::size_t>> &polygons,
        std::vector<CGAL::Color> &fcolors)
        : polygons(polygons), fcolors(fcolors) {}

    void operator()(const argument_type &a) {
      polygons.push_back(std::vector<std::size_t>(3));
      for (std::size_t i = 0; i < 3; ++i)
        polygons.back()[i] = a.first[i];

      unsigned char r, g, b;

      if (a.second < 0) {
        r = 128;
        g = 128;
        b = 128;
      } else {
        CGAL::Random rand(a.second);
        r = static_cast<unsigned char>(64 + rand.get_int(0, 192));
        g = static_cast<unsigned char>(64 + rand.get_int(0, 192));
        b = static_cast<unsigned char>(64 + rand.get_int(0, 192));
      }

      fcolors.push_back(CGAL::Color(r, g, b));
    }
  };

  struct array_to_vector {
    std::vector<std::vector<std::size_t>> &vectors;

    array_to_vector(std::vector<std::vector<std::size_t>> &vectors)
        : vectors(vectors) {}

    void operator()(const CGAL::cpp11::array<std::size_t, 3> &ar) {
      vectors.push_back(std::vector<std::size_t>(3));
      vectors.back()[0] = ar[0];
      vectors.back()[1] = ar[1];
      vectors.back()[2] = ar[2];
    }
  };

  struct Add_triangle_with_metaface_id {
    typedef std::pair<CGAL::cpp11::array<std::size_t, 3>, int> argument_type;
    typedef void result_type;

    std::vector<typename Kernel::Point_3> &vertices;
    std::vector<std::vector<std::size_t>> &polygons;
    Polylines_container &polylines;

    typedef std::map<std::pair<std::size_t, std::size_t>, int> Map_edges;
    Map_edges map_edges;

    Add_triangle_with_metaface_id(
        std::vector<typename Kernel::Point_3> &vertices,
        std::vector<std::vector<std::size_t>> &polygons,
        Polylines_container &polylines)
        : vertices(vertices), polygons(polygons), polylines(polylines) {}

    void operator()(const argument_type &a) {
      polygons.push_back(std::vector<std::size_t>(3));
      for (std::size_t i = 0; i < 3; ++i) {
        polygons.back()[i] = a.first[i];

        // handle edge
        std::size_t idx_a = a.first[i];
        std::size_t idx_b = a.first[(i + 1) % 3];
        if (idx_a > idx_b)
          std::swap(idx_a, idx_b);

        typename Map_edges::iterator iter;
        bool inserted = false;
        boost::tie(iter, inserted) = map_edges.insert(
            std::make_pair(std::make_pair(idx_a, idx_b), a.second));
        if (!inserted && iter->second != a.second) // edge between two metafaces
        {
          polylines.push_back(Polyline());
          polylines.back().push_back(vertices[idx_a]);
          polylines.back().push_back(vertices[idx_b]);
        }
      }
    }
  };

  using Visibility_from_semantic_map =
      CGAL::Level_of_detail::Visibility_from_semantic_map<
          Semantic_map_from_labels>;

  using LOD =
      LOD::Level_of_detail<Kernel, Point_set, Point_map,
                           Semantic_map_from_labels,
                           Visibility_from_semantic_map, CGAL::Tag_true>;

  using Parameters = typename LOD::Parameters;

  Mylod(const int num_input_parameters, const Input_parameters input_parameters,
        const std::string &logs_path)
      : m_terminal_parser(num_input_parameters, input_parameters, logs_path),
        m_logs_path(logs_path),
        m_logs_path_0_1(logs_path + "lod_0_1" + std::string(_SR_)),
        m_logs_path_2(logs_path + "lod_2" + std::string(_SR_)) {}

  void run() {

    // First, we set LOD parameters parsing a terminal input or to the default
    // values, including the path to the input data.
    set_parameters();

    // Then, we load data from a file.
    load_input_data();

    // Create a semantic map.
    Semantic_map_from_labels semantic_map(&m_point_set);

    std::istringstream gi("0");
    std::istringstream bi("3");
    std::istringstream ii("2");
    std::istringstream vi("1");
    int idx;
    while (gi >> idx) {
      std::cerr << idx << " is ground" << std::endl;
      semantic_map.map_l2sl->insert(
          std::make_pair(idx, CGAL::Level_of_detail::Semantic_label::GROUND));
    }
    while (bi >> idx) {
      std::cerr << idx << " is building boundary" << std::endl;
      semantic_map.map_l2sl->insert(std::make_pair(
          idx, CGAL::Level_of_detail::Semantic_label::BUILDING_BOUNDARY));
    }
    while (ii >> idx) {
      std::cerr << idx << " is building interior" << std::endl;
      semantic_map.map_l2sl->insert(std::make_pair(
          idx, CGAL::Level_of_detail::Semantic_label::BUILDING_INTERIOR));
    }
    while (vi >> idx) {
      std::cerr << idx << " is vegetation" << std::endl;
      semantic_map.map_l2sl->insert(std::make_pair(
          idx, CGAL::Level_of_detail::Semantic_label::VEGETATION));
    }

    // Create LOD.
    LOD lod(m_point_set, m_point_set.point_map(), semantic_map);

    // Compute ground plane.
    lod.compute_planar_ground();

    // Compute boundary points.
    lod.detect_building_boundaries(
        m_parameters.alpha_shape_size(),         // alpha shape size
        m_parameters.region_growing_2_epsilon(), // region growing epsilon
        m_parameters.region_growing_2_cluster_epsilon(),  // region growing
                                                          // cluster epsilon
        m_parameters.region_growing_2_normal_threshold(), // region growing
                                                          // normal threshold
        0.5, // region growing min wall length
        m_parameters.segment_regularizer_2_max_angle_in_degrees(),
        m_parameters.segment_regularizer_2_max_difference_in_meters(),
        m_parameters.grid_cell_width()); // grid cell width

    Point_set filtered;
    lod.output_filtered_boundary_points(filtered.point_back_inserter());

    Log log;
    log.save_points(filtered, filtered.point_map(),
                    m_logs_path_0_1 + "1_building_boundary_points");

    // Compute segments.
    Point_set segment_points;
    Insert_point_colored_by_index inserter(segment_points);

    lod.output_segmented_boundary_points(
        boost::make_function_output_iterator(inserter));

    log.save_points(segment_points, segment_points.point_map(),
                    m_logs_path_0_1 + "2_regions");

    // Compute lines.
    Polylines_container polylines;
    Add_polyline_from_segment adder(polylines);
    lod.output_boundary_edges(boost::make_function_output_iterator(adder));

    log.save_polylines<Polylines_container, typename Kernel::Segment_3>(
        polylines, m_logs_path_0_1 + "3_segments");

    // Detect trees.
    lod.detect_trees(m_parameters.tree_grid_cell_width(),
                     m_parameters.min_tree_height(),
                     m_parameters.min_tree_radius());

    Point_set trees_point_set;
    Insert_point_colored_by_index tree_inserter(trees_point_set);

    lod.output_tree_points(boost::make_function_output_iterator(tree_inserter));
    log.save_points(trees_point_set, trees_point_set.point_map(),
                    m_logs_path_0_1 + "4_trees");

    // Partition.
    lod.partition(m_parameters.kinetic_partitioning_2_min_face_width(),
                  m_parameters.kinetic_partitioning_2_num_intersections(),
                  true);

    std::vector<typename Kernel::Point_3> vertices;
    std::vector<std::vector<std::size_t>> polygons;
    std::vector<CGAL::Color> fcolors;

    lod.output_partition_to_polygon_soup(std::back_inserter(vertices),
                                         std::back_inserter(polygons));

    CGAL::Random rand(time(0));
    for (std::size_t i = 0; i < polygons.size(); ++i) {
      unsigned char r, g, b;
      r = static_cast<unsigned char>(64 + rand.get_int(0, 192));
      g = static_cast<unsigned char>(64 + rand.get_int(0, 192));
      b = static_cast<unsigned char>(64 + rand.get_int(0, 192));
      fcolors.push_back(CGAL::Color(r, g, b));
    }

    log.save_polygon_soup(vertices, polygons, fcolors,
                          m_logs_path_0_1 + "5_partition");

    polygons.clear();
    fcolors.clear();
    Add_polygon_with_visibility_color vis_adder(polygons, fcolors);

    lod.output_partition_with_visibility_to_polygon_soup(
        std::back_inserter(vertices),
        boost::make_function_output_iterator(vis_adder));

    log.save_polygon_soup(vertices, polygons, fcolors,
                          m_logs_path_0_1 + "6_visibility");

    // LOD0.
    lod.compute_footprints(m_parameters.segment_constraints_threshold(),
                           m_parameters.min_num_building_floor_faces());

    Add_triangle_with_building_color bu_adder(polygons, fcolors);

    vertices.clear();
    polygons.clear();
    fcolors.clear();

    lod.output_building_footprints_to_triangle_soup(
        std::back_inserter(vertices),
        boost::make_function_output_iterator(bu_adder));

    log.save_polygon_soup(vertices, polygons, fcolors,
                          m_logs_path_0_1 + "7_buildings");

    Polylines_container walls;
    Add_polyline_from_segment wal_adder(walls);
    lod.output_building_footprints_to_segment_soup(
        boost::make_function_output_iterator(wal_adder));

    log.save_polylines<Polylines_container, typename Kernel::Segment_3>(
        walls, m_logs_path_0_1 + "8_walls");

    vertices.clear();
    polygons.clear();
    fcolors.clear();

    std::size_t first_building_facet;
    std::size_t first_vegetation_facet;
    std::size_t first_wall_facet;

    boost::tie(first_building_facet, first_vegetation_facet) =
        lod.output_lod0_to_triangle_soup(
            std::back_inserter(vertices),
            boost::make_function_output_iterator(array_to_vector(polygons)));

    // Fill colors according to facet type
    for (std::size_t i = 0; i < first_building_facet; ++i)
      fcolors.push_back(CGAL::Color(186, 189, 182));
    for (std::size_t i = first_building_facet; i < first_vegetation_facet; ++i)
      fcolors.push_back(CGAL::Color(245, 121, 0));
    for (std::size_t i = first_vegetation_facet; i < polygons.size(); ++i)
      fcolors.push_back(CGAL::Color(138, 226, 52));

    log.save_polygon_soup(vertices, polygons, fcolors, m_logs_path + "LOD0");

    // LOD1.
    lod.extrude_footprints();

    vertices.clear();
    polygons.clear();
    fcolors.clear();

    Polylines_container poly;
    Add_triangle_with_metaface_id lod1_adder(vertices, polygons, poly);

    std::tie(first_building_facet, first_wall_facet, first_vegetation_facet) =
        lod.output_lod1_to_triangle_soup(
            std::back_inserter(vertices),
            boost::make_function_output_iterator(lod1_adder));

    // Fill colors according to facet type
    for (std::size_t i = 0; i < first_building_facet; ++i)
      fcolors.push_back(CGAL::Color(186, 189, 182));
    for (std::size_t i = first_building_facet; i < first_wall_facet; ++i)
      fcolors.push_back(CGAL::Color(245, 121, 0));
    for (std::size_t i = first_wall_facet; i < first_vegetation_facet; ++i)
      fcolors.push_back(CGAL::Color(77, 131, 186));
    for (std::size_t i = first_vegetation_facet; i < polygons.size(); ++i)
      fcolors.push_back(CGAL::Color(138, 226, 52));

    log.save_polygon_soup(vertices, polygons, fcolors, m_logs_path + "LOD1");
  }

private:
  Terminal_parser m_terminal_parser;
  Parameters m_parameters;
  Point_set m_point_set;

  std::string m_logs_path;
  std::string m_logs_path_0_1;
  std::string m_logs_path_2;

  void set_parameters() {

    // Set all parameters that can be loaded from the terminal.
    // add_str_parameter  - adds a string-type parameter
    // add_val_parameter  - adds a scalar-type parameter
    // add_bool_parameter - adds a boolean parameter

    std::cout << "Input parameters: " << std::endl;

    m_terminal_parser.add_str_parameter("-data", m_parameters.path_to_input());

    m_terminal_parser.add_bool_parameter("-silent", m_parameters.silent());
    m_terminal_parser.add_bool_parameter("-verbose", m_parameters.verbose());
    m_terminal_parser.add_bool_parameter("-no_simplification",
                                         m_parameters.no_simplification());
    m_terminal_parser.add_bool_parameter("-no_regularization",
                                         m_parameters.no_regularization());
    m_terminal_parser.add_bool_parameter(
        "-no_consistent_visibility", m_parameters.no_consistent_visibility());

    m_terminal_parser.add_val_parameter("-scale", m_parameters.scale());
    m_terminal_parser.add_val_parameter("-eps", m_parameters.epsilon());

    m_parameters.update_dependent();

    m_terminal_parser.add_val_parameter("-alpha",
                                        m_parameters.alpha_shape_size());
    m_terminal_parser.add_val_parameter("-cell",
                                        m_parameters.grid_cell_width());

    m_terminal_parser.add_val_parameter(
        "-rg_nt_2d", m_parameters.region_growing_2_normal_threshold());
    m_terminal_parser.add_val_parameter(
        "-rg_min_2d", m_parameters.region_growing_2_min_points());
    m_terminal_parser.add_val_parameter(
        "-rg_eps_2d", m_parameters.region_growing_2_epsilon());
    m_terminal_parser.add_val_parameter(
        "-rg_ce_2d", m_parameters.region_growing_2_cluster_epsilon());

    m_terminal_parser.add_val_parameter(
        "-angle", m_parameters.segment_regularizer_2_max_angle_in_degrees());
  }

  void load_input_data() {

    std::cout << std::endl << "Input data: " << std::endl;

    std::ifstream loader(m_parameters.path_to_input().c_str(),
                         std::ios_base::in);
    loader >> m_point_set;

    loader.close();
    std::cout << "LOD input data are loaded: number of points: "
              << m_point_set.number_of_points() << std::endl;

    if (m_point_set.template property_map<int>("label").second) {
      std::cout << "LOD algorithm is applicable to these data!" << std::endl;
    }
  }
};

} // namespace Levels_of_detail

} // namespace CGAL

#endif // CGAL_LEVELS_OF_DETAIL_MYLOD_H
