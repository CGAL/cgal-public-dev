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

#ifndef CGAL_WRITE_POINT_SET_WITH_SEGMENTS_H
#define CGAL_WRITE_POINT_SET_WITH_SEGMENTS_H

#include <iostream>
#include <vector>

namespace CGAL {

	/**
		Writes a point set (coordinates + normals) to an ASCII vg format.
		\tparam Point_set_with_segments is a model of `Point_set_with_segments`.
		\param stream output stream.
		\param point_set the point set to be saved to the stream.
		\return true on success.
	*/
	template <typename Point_set_with_segments>
	bool write_point_set_with_segments(std::ostream& stream, const Point_set_with_segments& point_set);


	//////////////////////////////////////////////////////////////////////////


	// implementation


	/// \cond SKIP_IN_MANUAL
	namespace internal {

		template <typename Planar_segment>
		std::vector<float> get_segment_parameters(const Planar_segment* s) {
			const unsigned int num = 4;
			std::vector<float> para(num);
			const typename Planar_segment::Plane* plane = s->supporting_plane();
			para[0] = static_cast<float>(plane->a());
			para[1] = static_cast<float>(plane->b());
			para[2] = static_cast<float>(plane->c());
			para[3] = static_cast<float>(plane->d());
			return para;
		}

		template <typename Planar_segment>
		void write_segment(std::ostream& output, const Planar_segment* s) {
			const int type = 0;
			output << "group_type: " << type << std::endl;

			const std::vector<float>& para = get_segment_parameters(s);
			output << "num_group_parameters: " << para.size() << std::endl;
			output << "group_parameters: ";
			for (std::size_t i = 0; i < para.size(); ++i)
				output << para[i] << " ";
			output << std::endl;

			std::string label("no_label");
			output << "group_label: " << label << std::endl;

			output << "group_color: " << 0.5 << " " << 0.5 << " " << 0.8 << std::endl;

			std::size_t num_point = s->size();
			output << "group_num_point: " << num_point << std::endl;

			for (std::size_t i = 0; i < s->size(); ++i) {
				output << s->at(i) << " ";
			}
			output << std::endl;
		}
	}

	// the full description of the vg format is here:
	// https://github.com/LiangliangNan/PolyFit/blob/master/ReadMe-data.md

	template <typename Point_set_with_segments>
	bool write_point_set_with_segments(std::ostream& output, const Point_set_with_segments& point_set) {
		if (!output)
			return false;

		// write positions

		const std::size_t num = point_set.number_of_points();
		output << "num_points: " << num << std::endl;

		const typename Point_set_with_segments::Point_map& points = point_set.point_map();
		for (std::size_t i = 0; i < num; ++i)
			output << points[i] << " ";
		output << std::endl;

		// skip colors

		output << "num_colors: " << 0 << std::endl;	

		// write normals (if exist)

		std::size_t num_normals = 0;
		if (point_set.has_normal_map())
			num_normals = num;
		output << "num_normals: " << num_normals << std::endl;

		const typename Point_set_with_segments::Vector_map& normals = point_set.normal_map();
		for (std::size_t i = 0; i < num; ++i)
			output << normals[i] << " ";
		output << std::endl;

		// write planar segments (if exist)

		typedef typename Point_set_with_segments::Planar_segment	Planar_segment;
		const std::vector<Planar_segment*>& segments = point_set.planar_segments();
		output << "num_groups: " << segments.size() << std::endl;
		for (std::size_t i = 0; i < segments.size(); ++i) {
			const Planar_segment* s = segments[i];
			internal::write_segment(output, s);

			// children
			output << "num_children: " << 0 << std::endl; // skip
		}
		return true;
	}
	/// \endcond


} //namespace CGAL

#endif // CGAL_WRITE_POINT_SET_WITH_SEGMENTS_H
