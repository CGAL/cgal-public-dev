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
// Author(s)     : Liangliang Nan

#ifndef CGAL_POLYGONAL_SURFACE_RECONSTRUCTION_POINT_SET_WITH_SEGMENTS_H
#define CGAL_POLYGONAL_SURFACE_RECONSTRUCTION_POINT_SET_WITH_SEGMENTS_H


#include <fstream>
#include <vector>

#include <CGAL/Point_set_3.h>
#include <CGAL/linear_least_squares_fitting_3.h>


/*!
\file Point_set_with_segments.h
*/

namespace CGAL {


	// forward declaration
	template <typename Kernel>
	class Point_set_with_segments;


	/** \ingroup PkgPolygonalSurfaceReconstruction
	*
	*	A group of points (represented by their indices) belonging to a planar segment in a point set.
	*/
	template <typename Kernel>
	class Planar_segment : public std::vector<std::size_t>
	{
	public:
		typedef typename Kernel::Point_3            Point;
		typedef typename Kernel::Plane_3			Plane;
		typedef Point_set_with_segments<Kernel>     Point_set;

	public:

		// \param point_set the point set that owns this planar segment.
		Planar_segment(Point_set* point_set = 0) : point_set_(point_set), supporting_plane_(nullptr) {}
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
		Plane *		supporting_plane_; // the hypothesis generator owns this plane and will manages the memory
	};


	/** \ingroup PkgPolygonalSurfaceReconstruction
	*	An enriched point set that stores the extracted planar segments
	*/
	template <typename Kernel>
	class Point_set_with_segments : public Point_set_3<typename Kernel::Point_3>
	{
	public:

		typedef Point_set_3<typename Kernel::Point_3>	Parent_class;
		typedef Point_set_with_segments<Kernel>			This_class;
		typedef typename Kernel::FT						FT;
		typedef typename Kernel::Point_3				Point;
		typedef typename Kernel::Vector_3				Vector;
		typedef Planar_segment<Kernel>					Planar_segment;

	public:
		Point_set_with_segments() {}
		~Point_set_with_segments() {
			for (std::size_t i = 0; i < planar_segments_.size(); ++i)
				delete planar_segments_[i];
		}

		std::vector< Planar_segment* >& planar_segments() { return planar_segments_; }
		const std::vector< Planar_segment* >& planar_segments() const { return planar_segments_; }

		/*
		// ASCII vg file format.
		// These functions reads/writes a subset of PolyFit's vg format. The vg format stores a
		// point set (possibly with normals, colors) and the extracted primitives/segmentation.
		// PolyFit's original vg format supports first different primitives and it also supports
		// general point set segmentations. Since we are dealing planes only, these functions
		// read/write planar segments only. The original vg format is described here:
		// https://github.com/LiangliangNan/PolyFit/blob/master/ReadMe-data.md
		//
		num_points: num
		x  y  z
		...

		num_colors: num
		r g b
		...

		num_normals: num
		nx  ny  nz

		num_groups: num

		group_type: type (integer: 	VG_PLANE = 0, VG_CYLINDER = 1, VG_SPHERE = 2, VG_CONE = 3, VG_TORUS = 4, VG_GENERAL = 5)
		num_group_parameters: NUM_GROUP_PARAMETERS   // number of floating point values (integer)
		group_parameters: float[NUM_GROUP_PARAMETERS]
		group_label: label  // the first group info
		group_color: color (r, g, b)
		group_num_points: num
		idx ...

		num_children: num

		group_type: type (integer: 	VG_PLANE = 0, VG_CYLINDER = 1, VG_SPHERE = 2, VG_CONE = 3, VG_TORUS = 4, VG_GENERAL = 5)
		num_group_parameters: NUM_GROUP_PARAMETERS   // number of floating point values (integer)
		group_parameters: float[NUM_GROUP_PARAMETERS]
		group_label: label  // 0th child of group 0
		group_color: color (r, g, b)
		group_num_points: num
		idx ...

		group_type: type (integer: 	VG_PLANE = 0, VG_CYLINDER = 1, VG_SPHERE = 2, VG_CONE = 3, VG_TORUS = 4, VG_GENERAL = 5)
		num_group_parameters: NUM_GROUP_PARAMETERS   // number of floating point values (integer)
		group_parameters: float[NUM_GROUP_PARAMETERS]
		group_label: label  // 1st child of group 0
		group_color: color (r, g, b)
		group_num_points: num
		idx ...
		*/
		bool read(const std::string& file_name);
		bool save(const std::string& file_name);

	private:
		std::vector< Planar_segment* > planar_segments_;
	};


	//////////////////////////////////////////////////////////////////////////


	namespace {

		template <typename Planar_segment>
		std::vector<float> get_segment_parameters(const Planar_segment* s) {
			int num = 4;
			std::vector<float> para(num);

			const typename Planar_segment::Plane* plane = s->supporting_plane();
			para[0] = static_cast<float>(plane->a());
			para[1] = static_cast<float>(plane->b());
			para[2] = static_cast<float>(plane->c());
			para[3] = static_cast<float>(plane->d());

			return para;
		}


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

		template <typename Planar_segment>
		void write_segment(std::ostream& output, const Planar_segment* s) {
			int type = 0;
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


	template <typename Kernel>
	bool Point_set_with_segments<Kernel>::read(const std::string& file_name) {
		std::ifstream input(file_name.c_str());
		if (input.fail())
			return false;

		std::string dumy;
		std::size_t num;

		input >> dumy >> num;
		Parent_class::resize(num);

		Parent_class::add_normal_map();

		for (std::size_t i = 0; i < num; ++i)
			input >> Parent_class::m_points[i];

		input >> dumy >> num;
		float rgb;
		for (std::size_t i = 0; i < num; ++i) {
			for (int j = 0; j < 3; ++j)
				input >> rgb;
		}

		input >> dumy >> num;
		for (std::size_t i = 0; i < num; ++i)
			input >> Parent_class::m_normals[i];

		//////////////////////////////////////////////////////////////////////////

		std::size_t num_segments = 0;
		input >> dumy >> num_segments;
		for (std::size_t i = 0; i < num_segments; ++i) {
			Planar_segment* s(read_segment<Planar_segment>(input));

			if (!s->empty()) {
				s->set_point_set(this);
				s->fit_supporting_plane();
				planar_segments_.push_back(s);
			}

			int num_children = 0; // skip
			input >> dumy >> num_children;
		}
		return true;
	}


	template <typename Kernel>
	bool Point_set_with_segments<Kernel>::save(const std::string& file_name) {
		// open file
		std::ofstream output(file_name.c_str());
		if (output.fail())
			return false;

		output << "num_points: " << Parent_class::number_of_points() << std::endl;
		for (std::size_t i = 0; i < Parent_class::number_of_points(); ++i)
			output << Parent_class::m_points[i] << " ";
		output << std::endl;

		output << "num_colors: " << 0 << std::endl;		// skip colors
		output << "num_normals: " << 0 << std::endl;	// skip colors

		output << "num_groups: " << planar_segments_.size() << std::endl;
		for (std::size_t i = 0; i < planar_segments_.size(); ++i) {
			const Planar_segment* s = planar_segments_[i];
			write_segment(output, s);

			// children
			output << "num_children: " << 0 << std::endl; // skip
		}
		return true;
	}

} //namespace CGAL


#endif // CGAL_POLYGONAL_SURFACE_RECONSTRUCTION_POINT_SET_WITH_SEGMENTS_H
