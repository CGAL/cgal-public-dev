// Copyright (c) 2020 SARL GeometryFactory (France).
// All rights reserved.
//
// This file is a part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY, AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
//
// Author(s)     : Dmitry Anisimov, Simon Giraudot, Pierre Alliez, Florent Lafarge, and Andreas Fabri

#ifndef CGAL_URBAN_AREA_PROCESSING_INTERNAL_UTILS_H
#define CGAL_URBAN_AREA_PROCESSING_INTERNAL_UTILS_H

// #include <CGAL/license/Urban_area_processing.h>

namespace CGAL {
namespace Urban_area_processing {
namespace internal {

  template<
  typename GeomTraits,
  typename FaceHandle>
  class Triangulation_neighbor_query_2 {
  public:
    using Traits = GeomTraits;
    using Face_handle = FaceHandle;

    Triangulation_neighbor_query_2(
      const std::vector<Face_handle>& faces) :
    m_faces(faces)
    { }

    void operator()(
      const std::size_t face_index,
      std::vector<std::size_t>& neighbors) const {

      neighbors.clear();
      const auto face = *(m_faces.begin() + face_index);
      for (std::size_t k = 0; k < 3; ++k) {
        const auto neighbor = face->neighbor(k);
        neighbors.push_back(neighbor->info().index);
      }
    }

  private:
    const std::vector<Face_handle>& m_faces;
  };

  template<
  typename GeomTraits,
  typename FaceHandle>
  class Triangulation_labeled_region_2 {
  public:
    using Traits = GeomTraits;
    using Face_handle = FaceHandle;

    Triangulation_labeled_region_2(
      const std::vector<Face_handle>& faces,
      const std::size_t ref_label,
      const std::size_t min_faces) :
    m_faces(faces),
    m_ref_label(ref_label),
    m_min_faces(min_faces)
    { }

    bool is_part_of_region(
      const std::size_t,
      const std::size_t face_index, 
      const std::vector<std::size_t>& region) const {

      const bool is_valid = check_region(region);
      if (!is_valid) return false;

      const auto face = *(m_faces.begin() + face_index);
      return face->info().label == m_ref_label;
    }

    inline bool is_valid_region(
      const std::vector<std::size_t>& region) const {

      const bool is_valid = check_region(region);
      return is_valid && (region.size() >= m_min_faces);
    }

    void update(
      const std::vector<std::size_t>&) { }

  private:
    const std::vector<Face_handle>& m_faces;
    const std::size_t m_ref_label;
    const std::size_t m_min_faces;

    bool check_region(
      const std::vector<std::size_t>& region) const {
      
      if (region.size() == 0)
        return false;

      if (region.size() == 1) {
        const std::size_t face_index = region[0];
        const auto face = *(m_faces.begin() + face_index);
        return face->info().label == m_ref_label;
      }
      return true;
    }
  };

  template<typename Point>
  class Indexer {
  
    using Local_traits = Exact_predicates_inexact_constructions_kernel;
    using Local_point_3 = typename Local_traits::Point_3;

  public:
    std::size_t operator()(const Point& point) {
      const Local_point_3 p = Local_point_3(
        CGAL::to_double(point.x()),
        CGAL::to_double(point.y()),
        CGAL::to_double(point.z()));

      const auto pair = m_indices.insert(
        std::make_pair(
          p, 
          m_indices.size()));
      const auto& item = pair.first;
      const std::size_t idx = item->second;
      return idx;
    }
    void clear() { m_indices.clear(); }

  private:
    std::map<Local_point_3, std::size_t> m_indices;
  };

  template<typename GeomTraits> 
  class Default_sqrt {
    
  private:
    using Traits = GeomTraits;
    using FT = typename Traits::FT;

  public:
    FT operator()(const FT value) const { 
      
      CGAL_precondition(value >= FT(0));
      return static_cast<FT>(CGAL::sqrt(CGAL::to_double(value)));
    }
  };

  BOOST_MPL_HAS_XXX_TRAIT_NAMED_DEF(Has_nested_type_Sqrt, Sqrt, false)

  // Case: do_not_use_default = false.
  template<typename GeomTraits, 
  bool do_not_use_default = Has_nested_type_Sqrt<GeomTraits>::value>
  class Get_sqrt {
        
  public:
    using Traits = GeomTraits;
    using Sqrt = Default_sqrt<Traits>;

    static Sqrt sqrt_object(const Traits& ) { 
      return Sqrt();
    }
  };

  // Case: do_not_use_default = true.
  template<typename GeomTraits>
  class Get_sqrt<GeomTraits, true> {
        
  public:
    using Traits = GeomTraits;
    using Sqrt = typename Traits::Sqrt;

    static Sqrt sqrt_object(const Traits& traits) { 
      return traits.sqrt_object();
    }
  };

  template<typename FT>
  static FT max_value() {
    return FT(1000000000000);
  }

  template<typename FT>
  static FT tolerance() {
    return FT(1) / FT(100000);
  }

  template<
  typename Face_handle,
  typename Point_3>
  void point_3(
    const Face_handle& fh, 
    const std::size_t idx,
    Point_3& p) {
    p = Point_3(
      fh->vertex(idx)->point().x(), 
      fh->vertex(idx)->point().y(), 
      fh->info().z[idx]);
  }

  template<
  typename Face_handle,
  typename Triangle_3>
  void triangle_3(
    const Face_handle& fh, 
    Triangle_3& triangle) {
    
    using Traits = typename Kernel_traits<Triangle_3>::Kernel;
    using Point_3 = typename Traits::Point_3;
    
    Point_3 p0, p1, p2;
    point_3(fh, 0, p0);
    point_3(fh, 1, p1);
    point_3(fh, 2, p2);
    triangle = Triangle_3(p0, p1, p2);
  }

  template<typename Point_3>
  typename Kernel_traits<Point_3>::Kernel::Point_2
  point_2_from_point_3(const Point_3& point_3) {
    return typename Kernel_traits<Point_3>::Kernel::Point_2(
      point_3.x(), point_3.y());
  }

} // internal
} // Urban_area_processing
} // CGAL

#endif // CGAL_URBAN_AREA_PROCESSING_INTERNAL_UTILS_H
