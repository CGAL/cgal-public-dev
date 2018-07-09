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

#ifndef CGAL_READ_WRITE_POINT_SET_WITH_SEGMENTS_H
#define CGAL_READ_WRITE_POINT_SET_WITH_SEGMENTS_H

#include <iostream>
#include <vector>

namespace CGAL {

	/**
	   Reads a point set (coordinates + normals) from an ASCII vg format.
	   \tparam Point_set_with_planes is a model of `Point_set_with_planes`.
	   \param stream input stream.
	   \param point_set the point set to store the data.
	   \return true on success.
	*/
	template <typename Point_set_with_planes>
	bool read_point_set_with_segments(std::istream& stream, Point_set_with_planes& point_set);


	//////////////////////////////////////////////////////////////////////////


	// implementation


	/// \cond SKIP_IN_MANUAL
	namespace internal {

		template <typename Planar_segment>
		Planar_segment* read_segment(std::istream& input) {
			std::string dumy;
			int type;
			input >> dumy >> type;

			int num;
			input >> dumy >> num;
			CGAL_assertion(num == 4);
			std::vector<float> para(num);
			input >> dumy;
			for (int i = 0; i < num; ++i)
				input >> para[i];

			std::string label;
			input >> dumy >> label;

			float r, g, b;
			input >> dumy >> r >> g >> b;

			int num_points;
			input >> dumy >> num_points;

			Planar_segment* s = new Planar_segment;
			for (int i = 0; i < num_points; ++i) {
				int idx;
				input >> idx;
				s->push_back(idx);
			}

			// ignore label
			// ignore color

			return s;
		}

	}

	// the full description of the vg format is here:
	// https://github.com/LiangliangNan/PolyFit/blob/master/ReadMe-data.md

	template <typename Point_set_with_planes>
	bool read_point_set_with_segments(std::istream& input, Point_set_with_planes& point_set) {
		if (!input)
			return false;

		std::string dumy;
		std::size_t num;

		input >> dumy >> num;
		point_set.resize(num);

                typename Point_set_with_planes::Point_map points = point_set.point_map();
		for (std::size_t i = 0; i < num; ++i)
			input >> points[i];

		input >> dumy >> num;
		float rgb;
		for (std::size_t i = 0; i < num; ++i) {
			for (int j = 0; j < 3; ++j)
				input >> rgb;
		}

		input >> dumy >> num;

		point_set.add_normal_map();
                typename Point_set_with_planes::Vector_map normals = point_set.normal_map();

		for (std::size_t i = 0; i < num; ++i)
			input >> normals[i];

		std::size_t num_segments = 0;
		input >> dumy >> num_segments;

		typedef typename Point_set_with_planes::Planar_segment	Planar_segment;
		std::vector<Planar_segment*>& segments = point_set.planar_segments();
		for (std::size_t i = 0; i < num_segments; ++i) {
			Planar_segment* s(internal::read_segment<Planar_segment>(input));

			if (!s->empty()) {
				s->set_point_set(&point_set);
				s->fit_supporting_plane();
				segments.push_back(s);
			}

			int num_children = 0; // skip
			input >> dumy >> num_children;
		}
		return true;
	}
	/// \endcond


} //namespace CGAL

#endif // CGAL_READ_WRITE_POINT_SET_WITH_SEGMENTS_H
