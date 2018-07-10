// Copyright (c) 2018  Liangliang Nan
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
//
// Author(s) : Liangliang Nan

#ifndef CGAL_POLYGONAL_SURFACE_RECONSTRUCTION_POINT_SET_WITH_PLANES_H
#define CGAL_POLYGONAL_SURFACE_RECONSTRUCTION_POINT_SET_WITH_PLANES_H

#include <CGAL/Point_set_3.h>
#include <CGAL/linear_least_squares_fitting_3.h>
#include <CGAL/license/Point_set_processing_3.h>
#include <CGAL/property_map.h>

#include <vector>


/*!
\file Point_set_with_planes.h
*/

namespace CGAL {


	namespace internal {

		// forward declaration
		template <typename Kernel>
		class Point_set_with_planes;


		/**
		*	A group of points (represented by their indices) belonging to a planar segment in a point set.
		*/
		template <typename Kernel>
		class Planar_segment : public std::vector<std::size_t>
		{
		public:
			typedef typename Kernel::Point_3            Point;
			typedef typename Kernel::Plane_3			Plane;
			typedef Point_set_with_planes<Kernel>		Point_set;

		public:

			// \param point_set the point set that owns this planar segment.
			Planar_segment(Point_set* point_set) : point_set_(point_set), supporting_plane_(nullptr) {}
			~Planar_segment() {}

			Point_set* point_set() { return point_set_; }
			void set_point_set(Point_set* point_set) { point_set_ = point_set; }

			// fits and returns the supporting plane of this planar segment
			Plane* fit_supporting_plane() {
				const typename Point_set::Point_map& points = point_set_->point_map();
				std::list<Point> pts;
				for (std::size_t i = 0; i < size(); ++i) {
					std::size_t idx = at(i);
					pts.push_back(points[idx]);
				}

				if (supporting_plane_)
					delete supporting_plane_;
				supporting_plane_ = new Plane;
				CGAL::linear_least_squares_fitting_3(pts.begin(), pts.end(), *supporting_plane_, CGAL::Dimension_tag<0>());
				return supporting_plane_;
			}

			// returns the supporting plane of this planar segment.
			// Note: returned plane is valid only if fit_supporting_plane() has been called.
			Plane* supporting_plane() const { return supporting_plane_; }

		private:
			Point_set * point_set_;
			Plane *		supporting_plane_; // the hypothesis generator owns this plane and manages the memory
		};


		/**
		*	An enriched point set that stores the extracted planar segments
		*/
		template <typename Kernel>
		class Point_set_with_planes : public Point_set_3<typename Kernel::Point_3>
		{
		public:

			typedef Point_set_3<typename Kernel::Point_3>	Parent_class;
			typedef Point_set_with_planes<Kernel>			This_class;

			typedef typename Kernel::FT						FT;
			typedef typename Kernel::Point_3				Point;
			typedef typename Kernel::Vector_3				Vector;
			typedef typename Kernel::Plane_3				Plane;
			typedef Planar_segment<Kernel>					Planar_segment;

		public:

			template <
				typename PointRange,
				typename NamedParameters
			>
			Point_set_with_planes(
				const PointRange& points,
				const NamedParameters& np) 
			{
				// basic geometric types
				typedef typename Point_set_processing_3::GetPointMap<PointRange, NamedParameters>::type PointMap;
				typedef typename Point_set_processing_3::GetNormalMap<PointRange, NamedParameters>::type NormalMap;
				typedef typename Point_set_processing_3::GetPlaneIndexMap<NamedParameters>::type PlaneIndexMap;

				CGAL_static_assertion_msg(!(boost::is_same<NormalMap,
					typename Point_set_processing_3::GetNormalMap<PointRange, NamedParameters>::NoMap>::value),
					"Error: no normal map");
				CGAL_static_assertion_msg(!(boost::is_same<PlaneIndexMap,
					typename Point_set_processing_3::GetPlaneIndexMap<NamedParameters>::NoMap>::value),
					"Error: no plane index map");

				PointMap point_map = boost::choose_param(get_param(np, internal_np::point_map), PointMap());
				NormalMap normal_map = boost::choose_param(get_param(np, internal_np::normal_map), NormalMap());
				PlaneIndexMap index_map = boost::choose_param(get_param(np, internal_np::plane_index_map), PlaneIndexMap());

				resize(points.size());
				add_normal_map();

				// get to know the number of plane from the plane indices
				int max_plane_index = 0;
				std::size_t idx = 0;
				for (typename PointRange::const_iterator it = points.begin(); it != points.end(); ++it) {
					int plane_index = get(index_map, idx);
					if (plane_index > max_plane_index)
						max_plane_index = plane_index;
					++idx;
				}
				std::size_t num_plane = max_plane_index + 1; // the first one has index 0

				for (std::size_t i = 0; i < num_plane; ++i)
					planar_segments_.push_back(new Planar_segment(this));

				idx = 0;
				for (typename PointRange::const_iterator it = points.begin(); it != points.end(); ++it) {
					m_points[idx] = get(point_map, *it);
					m_normals[idx] = get(normal_map, *it);

					int plane_index = get(index_map, idx);
					if (plane_index != -1) {
						Planar_segment* ps = planar_segments_[plane_index];
						ps->push_back(idx);
					}

					++idx;
				}
			}

			~Point_set_with_planes() {
				for (std::size_t i = 0; i < planar_segments_.size(); ++i)
					delete planar_segments_[i];
			}

			std::vector< Planar_segment* >& planar_segments() { return planar_segments_; }
			const std::vector< Planar_segment* >& planar_segments() const { return planar_segments_; }

		private:
			std::vector< Planar_segment* > planar_segments_;
		};

	} //namespace internal

} //namespace CGAL


#endif // CGAL_POLYGONAL_SURFACE_RECONSTRUCTION_POINT_SET_WITH_PLANES_H