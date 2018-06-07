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

#ifndef CGAL_POLYGONAL_SURFACE_RECONSTRUCTION_H
#define CGAL_POLYGONAL_SURFACE_RECONSTRUCTION_H

#include <CGAL/Point_set_with_segments.h>
#include <CGAL/Shape_detection_3.h>
#include <CGAL/algo/hypothesis.h>
#include <CGAL/algo/compute_confidences.h>
#include <CGAL/algo/face_selection.h>

/*!
\file Polygonal_surface_reconstruction.h
*/

namespace CGAL {

	/*!
	\ingroup PkgPolygonalSurfaceReconstruction

	\brief Implementation of the Polygonal Surface Reconstruction method.

	Given a set of 3D points with normals (either oriented or unoriented) sampled
	from the outer boundary of a piecewise planar object, the Polygonal Surface
	Reconstruction method \cgalCite{nan2017polyfit} outputs a simplified and
	watertight surface mesh interpolating the input point set.

	The reconstruction consists of the following steps:
	1) extracting planes from the input point set (can be skipped if planes are
	known or provided by other means);
	2) generating a set of face candidates by intersecting the extracted planar primitives;
	3) an optimal subset of the candidate faces is selected through optimization
	under hard constraints that enforce the final model to be manifold and watertight.

	*/
	template <class Kernel>
	class Polygonal_surface_reconstruction
	{
	public:

		/// \name Types 

		typedef typename Kernel::FT							FT;			///< number type.
		typedef typename Kernel::Point_3					Point;		///< point type.
		typedef typename Kernel::Vector_3					Vector;		///< vector type.
		typedef std::pair<Point, Vector>					Point_with_normal;
		typedef std::vector<Point_with_normal>				Point_list;

		typedef CGAL::Planar_segment<Kernel>				Planar_segment;

		typedef CGAL::Point_set_with_segments<Kernel>		Point_set_with_segments;
		///< enriched point set to access the extracted planar segments

		// Public methods
	public:

		/// \name Creation 
		/*!
		Creates a Polygonal Surface Reconstruction object
		*/
		Polygonal_surface_reconstruction() : hypothesis_(nullptr) {}
		/// \name Operations

		~Polygonal_surface_reconstruction() {
			if (hypothesis_)
				delete hypothesis_;
		}

		/*!
		Reconstruct from the input points represented by `input_range`.
		This is a one-shot reconstruction. Use the the step-by-step functions if you want to reuse the intermediate results.
		\tparam Surface_mesh is a model of `Surface_mesh`.
		\return `true` if plane extraction succeeded, `false` otherwise.
		*/
		template <typename Surface_mesh>
		bool reconstruct(
			const Point_list& points,			///< input point set with normals.
			Surface_mesh& output_mesh,			///< the final reconstruction results
			double wt_fitting = 0.43,			///< weight for the data fitting term.
			double wt_coverage = 0.27,			///< weight for the point coverage term.
			double wt_complexity = 0.30			///< weight for the model complexity term.
		);


		/** extracts planar segments from the input point set using RANSAC.
		\return `true` if plane extraction succeeded, `false` otherwise.
		*/
		bool extract_planes(
			const Point_list& points,			///< input point set with normals.
			Point_set_with_segments& segments	///< the extracted planar segments.
		);


		/** generates candidate faces.
		\tparam Surface_mesh is a model of `Surface_mesh`.
		\return `false` if error occurs, `true` otherwise.
		*/
		template <typename Surface_mesh>
		bool generate_candidate_faces(
			const Point_set_with_segments& segments,	///< point set with planar segments
			Surface_mesh& candidate_faces 				///< candidate faces by pairwise intersection of the planar segments
		);


		/** computes confidence values for each face.
		// - supporting point number:	stored as property 'f:num_supporting_points'
		// - face area:					stored as property 'f:face_area'
		// - covered area:				stored as property 'f:covered_area'
		\tparam Surface_mesh is a model of `Surface_mesh`.
		*/
		template <typename Surface_mesh>
		void compute_confidences(
			const Point_set_with_segments& segments,	///< point set with planar segments
			Surface_mesh& candidate_faces 				///< candidate faces by pairwise intersection of the planar segments
		);


		/** select the optimal subset of the faces to assemble a watertight polygonal mesh model.
		\tparam Surface_mesh is a model of `Surface_mesh`.
		\return `true` if optimization succeeded, `false` otherwise.
		*/
		template <typename Surface_mesh>
		bool select_faces(
			const Surface_mesh& candidate_faces,	///< candidate faces with face confidence values stored in property map "f:confidence"
			Surface_mesh& output_mesh,				///< the final reconstruction results
			double wt_fitting = 0.43,				///< weight for the data fitting term.
			double wt_coverage = 0.27,				///< weight for the point coverage term.
			double wt_complexity = 0.30				///< weight for the model complexity term.
		);

		// Data members.
	private:
		Hypothesis<Kernel> * hypothesis_;

	private: // disallow copying
		Polygonal_surface_reconstruction(const Polygonal_surface_reconstruction& psr);

	}; // end of Polygonal_surface_reconstruction


	   //////////////////////////////////////////////////////////////////////////s

	   // implementations

	template <class Kernel>
	template <typename Surface_mesh>
	bool Polygonal_surface_reconstruction<Kernel>::reconstruct(
		const Point_list& points,
		Surface_mesh& output_mesh,
		double wt_fitting /* = 0.43 */,
		double wt_coverage /* = 0.27 */,
		double wt_complexity /* = 0.30 */)
	{
		Point_set_with_segments segments;
		if (!extract_planes(points, segments))
			return false;

		Surface_mesh candidate_faces;
		if (!generate_candidate_faces(segments, candidate_faces))
			return false;

		compute_confidences(segments, candidate_faces);

		output_mesh.clear();
		return select_faces(candidate_faces, output_mesh, wt_fitting, wt_coverage, wt_complexity);
	}


	template <class Kernel>
	bool Polygonal_surface_reconstruction<Kernel>::extract_planes(
		const Point_list& points,
		Point_set_with_segments& segments)
	{
		typedef CGAL::First_of_pair_property_map<Point_with_normal>  Point_map;
		typedef CGAL::Second_of_pair_property_map<Point_with_normal> Normal_map;
		typedef CGAL::Shape_detection_3::Shape_detection_traits
			<Kernel, Point_list, Point_map, Normal_map>              Traits;
		typedef CGAL::Shape_detection_3::Efficient_RANSAC<Traits>    Efficient_ransac;
		typedef CGAL::Shape_detection_3::Plane<Traits>               Plane;

		// Instantiates the Efficient_ransac engine.
		Efficient_ransac ransac;

		// Provides the input data.
#if 1 // Why RANSAC requires an non-const input?
		Point_list* point_set = const_cast<Point_list*>(&points);
		ransac.set_input(*point_set);
#else		
		ransac.set_input(points);
#endif

		// Registers planar shapes via template method.
		ransac.template add_shape_factory<Plane>();

		// Detects registered shapes with default parameters.
		ransac.detect();

		const typename Efficient_ransac::Shape_range& shapes = ransac.shapes();

		// update the output.

		// clean first
		segments.clear();

		// upload points
		std::size_t num = points.size();
		for (std::size_t i = 0; i < num; ++i)
			segments.insert(points[i].first);

		// ignore the colors
		// ignore the normals

		// now the planar segments

		typename Efficient_ransac::Shape_range::const_iterator it = shapes.begin();
		for (; it != shapes.end(); ++it) {
			boost::shared_ptr<typename Efficient_ransac::Shape> shape = *it;
			const std::vector<std::size_t>& indices = (*it)->indices_of_assigned_points();
			Planar_segment* s = new Planar_segment;
			s->set_point_set(&segments);
			s->insert(s->end(), indices.begin(), indices.end());
			s->fit_supporting_plane();
			segments.planar_segments().push_back(s);
		}

		return !shapes.empty();
	}


	template <class Kernel>
	template <typename Surface_mesh>
	bool Polygonal_surface_reconstruction<Kernel>::generate_candidate_faces(
		const Point_set_with_segments& point_set,
		Surface_mesh& candidate_faces
	)
	{
		const std::vector< Planar_segment* >& planar_segments = point_set.planar_segments();
		if (planar_segments.size() < 4) {
			std::cerr << "not enough (" << planar_segments.size() << ") planar segments to"
				<< " reconstruct a closed polygonal mesh" << std::endl;
			return false;
		}

		if (!hypothesis_)
			hypothesis_ = new Hypothesis<Kernel>(&point_set);
		hypothesis_->generate(candidate_faces);

		return candidate_faces.num_faces() > 4;
	}


	template <class Kernel>
	template <typename Surface_mesh>
	void Polygonal_surface_reconstruction<Kernel>::compute_confidences(
		const Point_set_with_segments& segments,
		Surface_mesh& candidate_faces
	)
	{
		Candidate_confidences<Kernel> conf;
		conf.compute(segments, candidate_faces);
	}


	template <class Kernel>
	template <typename Surface_mesh>
	bool Polygonal_surface_reconstruction<Kernel>::select_faces(
		const Surface_mesh& candidate_faces,
		Surface_mesh& output_mesh,
		double wt_fitting /* = 0.43*/,
		double wt_coverage /* = 0.27*/,
		double wt_complexity /* = 0.30*/
	)
	{
		typedef typename Hypothesis<Kernel>::Adjacency Adjacency;
		const Adjacency& adjacency = hypothesis_->extract_adjacency(candidate_faces);

		Face_selection<Kernel> sel;
		return sel.optimize(candidate_faces, adjacency, output_mesh, wt_fitting, wt_coverage, wt_complexity);
	}

} //namespace CGAL

#endif // CGAL_POLYGONAL_SURFACE_RECONSTRUCTION_H
