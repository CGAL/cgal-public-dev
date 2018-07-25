#ifndef CGAL_LEVEL_OF_DETAIL_BUILDING_HEIGHT_ESTIMATOR_H
#define CGAL_LEVEL_OF_DETAIL_BUILDING_HEIGHT_ESTIMATOR_H

// STL includes.
#include <map>
#include <list>

// LOD includes.
#include <CGAL/Level_of_detail/Enumerations.h>
#include <CGAL/Level_of_detail/internal/Buildings/Building.h>
#include <CGAL/Level_of_detail/internal/Buildings/Average_height_fitter.h>
#include <CGAL/Level_of_detail/internal/Buildings/Fixed_10_meters_height_fitter.h>

namespace CGAL {

namespace Level_of_detail {

template<class GeomTraits, class Triangulation>
class Building_height_estimator {
			
public:
  using Kernel        = GeomTraits;

  typename Kernel::Compute_squared_distance_3 squared_distance_3;

  using FT      = typename Kernel::FT;
  using Point_3 = typename Kernel::Point_3;
  using Plane_3 = typename Kernel::Plane_3;
            
  using Triangulation_faces_iterator = typename Triangulation::Finite_faces_iterator;
  using Triangulation_face_handle   = typename Triangulation::Face_handle;

  using Building = Building<Kernel, Triangulation_face_handle>;
            
  using Roof_points          = std::vector<Point_3>;
  using Building_roof_points = std::map<int, Roof_points>;
  using Building_heights     = std::map<int, FT>;

  using Roof_points_iterator          = typename Roof_points::const_iterator;
  using Building_roof_points_iterator = typename Building_roof_points::const_iterator;
            
  using Average_height_fitter         = LOD::Average_height_fitter<Kernel>;
  using Fixed_10_meters_height_fitter = LOD::Fixed_10_meters_height_fitter<Kernel>;

  template<class PointMap>
  Building_height_estimator(
    const Triangulation  &triangulation,
    const PointMap       &point_map,
    const Plane_3        &ground_plane,
    const Flat_roof_type flat_roof_type) :

    m_triangulation(triangulation),
    m_ground_plane(ground_plane),
    m_flat_roof_type(flat_roof_type) { 

    collect_building_roof_points(point_map);
    fit_roofs();
  }

  inline FT get_local_ground_height() const {
    return m_ground_plane.point().z();
  }

  inline FT get_building_height(const int building_index) const {

    if (!does_building_height_exist(building_index)) return FT(0);
    return m_building_heights.at(building_index);
  }


  void estimate (Building& building) const {
    building.local_ground_height() = get_local_ground_height();

    CGAL_precondition(building.index() >= 0);
    building.height() = get_building_height(building.index());
    if (building.height() == FT(0)) building.is_valid() = false;
  }

private:
  const Triangulation &m_triangulation;
  const Plane_3       &m_ground_plane;
        
  const Flat_roof_type m_flat_roof_type;
            
  Building_roof_points m_building_roof_points;
  Building_heights     m_building_heights;

  Average_height_fitter         m_average_height_fitter;
  Fixed_10_meters_height_fitter m_fixed_10_meters_height_fitter;

  bool does_building_height_exist(const int building_index) const {
    return m_building_heights.find(building_index) != m_building_heights.end();
  }

  template<class PointMap>
  void collect_building_roof_points(const PointMap &point_map) {

    for (Triangulation_faces_iterator tf_it = m_triangulation.finite_faces_begin();
         tf_it != m_triangulation.finite_faces_end(); ++tf_it) {
      const int building_number = tf_it->info().group_number();
                    
      if (building_number < 0) continue;
      const auto &elements = tf_it->info().elements();

      for (auto element = elements.begin(); element != elements.end(); ++element)
        m_building_roof_points[building_number].push_back(get(point_map, **element));
    }
  }

  void fit_roofs() {

    for (Building_roof_points_iterator bp_it = m_building_roof_points.begin();
         bp_it != m_building_roof_points.end(); ++bp_it) {
      const Roof_points &roof_points = (*bp_it).second;

      clear_fitter();
      for (Roof_points_iterator rp_it = roof_points.begin(); rp_it != roof_points.end(); ++rp_it) {

        const Point_3 &p = *rp_it;
        const Point_3 &q = m_ground_plane.projection(p);

        const FT height = static_cast<FT>(CGAL::sqrt(CGAL::to_double(squared_distance_3(p, q))));
        if (is_valid_height(height)) add_height_to_fitter(height);
      }

      const int building_index           = (*bp_it).first;
      m_building_heights[building_index] = get_result_from_fitter();
    }
  }

  inline bool is_valid_height(const FT height) const {
    return std::isfinite(CGAL::to_double(height));
  }

  void clear_fitter() {
    switch (m_flat_roof_type) {

      case Flat_roof_type::AVERAGE:
        m_average_height_fitter.clear();
        break;

      case Flat_roof_type::FIXED_10_METERS:
        m_fixed_10_meters_height_fitter.clear();
        break;

      default:
        m_average_height_fitter.clear();
        break;
    }
  }

  void add_height_to_fitter(const FT height) {
    switch (m_flat_roof_type) {

      case Flat_roof_type::AVERAGE:
        m_average_height_fitter.add_height(height);
        break;

      case Flat_roof_type::FIXED_10_METERS:
        m_fixed_10_meters_height_fitter.add_height(height);
        break;

      default:
        m_average_height_fitter.add_height(height);
        break;
    }
  }

  FT get_result_from_fitter() const {
    switch (m_flat_roof_type) {

      case Flat_roof_type::AVERAGE:
        return m_average_height_fitter.get_result();
        break;

      case Flat_roof_type::FIXED_10_METERS:
        return m_fixed_10_meters_height_fitter.get_result();
        break;

      default:
        return m_average_height_fitter.get_result();
        break;
    }
  }
};

} // Level_of_detail

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_BUILDING_HEIGHT_ESTIMATOR_H
