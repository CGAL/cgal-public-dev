// Copyright (c) 2019 INRIA Sophia-Antipolis (France).
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
//
// Author(s)     : Dmitry Anisimov, Simon Giraudot, Pierre Alliez, Florent Lafarge, and Andreas Fabri
//

#ifndef CGAL_LEVELS_OF_DETAIL_INTERNAL_LOD2_IMAGE_RECONSTRUCTION_H
#define CGAL_LEVELS_OF_DETAIL_INTERNAL_LOD2_IMAGE_RECONSTRUCTION_H

// STL includes.
#include <map>
#include <set>
#include <list>
#include <queue>
#include <vector>
#include <utility>
#include <memory>

// CGAL includes.
#include <CGAL/property_map.h>

// Internal includes.
#include <CGAL/Levels_of_detail/internal/utils.h>
#include <CGAL/Levels_of_detail/internal/struct.h>

// Other includes.
#include <CGAL/Levels_of_detail/internal/Reconstruction/Image_creator.h>
#include <CGAL/Levels_of_detail/internal/Reconstruction/Image_data_structure.h>

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

template<
typename GeomTraits,
typename ImagePointer>
class LOD2_image_reconstruction {

public:
  using Traits = GeomTraits;
  using Image_ptr = ImagePointer;

  using FT = typename Traits::FT;
  using Segment_2 = typename Traits::Segment_2;
  using Plane_3 = typename Traits::Plane_3;

  using Triangulation = internal::Triangulation<Traits>;
  using Partition_2 = internal::Partition_2<Traits>;

  using Image_creator = internal::Image_creator<Traits, Image_ptr>;
  using Image_data_structure = internal::Image_data_structure<Traits>;

  LOD2_image_reconstruction(
    const std::vector<Segment_2>& boundary,
    const Triangulation& lod0,
    ImagePointer& image_ptr,
    Partition_2& partition_2,
    const FT noise_level_2,
    const FT min_length_2,
    const FT angle_bound_2,
    const FT ordinate_bound_2) :
  m_boundary(boundary),
  m_lod0(lod0),
  m_image_ptr(image_ptr),
  m_partition_2(partition_2),
  m_noise_level_2(noise_level_2),
  m_min_length_2(min_length_2),
  m_angle_bound_2(angle_bound_2),
  m_ordinate_bound_2(ordinate_bound_2),
  m_pi(static_cast<FT>(CGAL_PI)),
  m_image_creator( 
    m_image_ptr, m_boundary, m_noise_level_2) { 

    m_partition_2.clear();
  }

  void build() {
    m_image_creator.create_image();
    m_image_creator.clean_image();
    m_image_creator.create_label_pairs();
    m_image_creator.create_ridges();
    m_image_creator.create_contours();

    m_data_structure_ptr = std::make_shared<Image_data_structure>(
      m_boundary, 
      m_image_creator.get_ridges(),
      m_image_creator.get_image());
      
    m_data_structure_ptr->clear();
    m_data_structure_ptr->build();
  }

  void create_triangulation() {
  
  }

  void compute_visibility() {

  }

  void label_faces() {
  
  }

  void get_roof_planes(
    std::vector<Plane_3>& roof_planes) {

    const auto& plane_map = m_image_ptr->get_plane_map();
    roof_planes.clear();
    roof_planes.reserve(plane_map.size());
    
    for (const auto& pair : plane_map) {
      const auto& plane = pair.second;
      roof_planes.push_back(plane);
    }
  }

private:
  const std::vector<Segment_2>& m_boundary;
  const Triangulation& m_lod0;
  Image_ptr& m_image_ptr;
  Partition_2& m_partition_2;
  const FT m_noise_level_2;
  const FT m_min_length_2;
  const FT m_angle_bound_2;
  const FT m_ordinate_bound_2;
  const FT m_pi;

  Image_creator m_image_creator;
  std::shared_ptr<Image_data_structure> m_data_structure_ptr;
};

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_INTERNAL_LOD2_IMAGE_RECONSTRUCTION_H
