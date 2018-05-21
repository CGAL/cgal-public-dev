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

#include <CGAL/disable_warnings.h>
#include <CGAL/Point_set_with_segments.h>
#include <CGAL/Shape_detection_3.h>

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

		/// @{
		/// \cond SKIP_IN_MANUAL
		typedef typename Kernel::FT							FT;			///< number type.
		typedef typename Kernel::Point_3					Point;		///< point type.
		typedef typename Kernel::Vector_3					Vector;		///< vector type.
		typedef std::pair<Point, Vector>					Point_with_normal;
		typedef std::vector<Point_with_normal>				Point_list;
		/// \endcondt

		typedef  CGAL::Point_set_with_segments<Kernel>		Point_set_with_segments;
		///< enriched point set to access the extracted planar segments

		// Public methods
	public:

		/// \name Creation 
		/*!
		Creates a Polygonal Surface Reconstruction object
		*/
		Polygonal_surface_reconstruction() {}

		/// \name Operations

		/*!
		/** reconstruct from the input points represented by `input_range`.
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
			const Point_list& points,					///< input point set with normals.
			Point_set_with_segments& planar_segments	///< the extracted planar segments.
		);


		/** generates candidate faces and computes confidence values for each face.
		The confidence values of the faces are stored as as property map with name "f:confidence".
		\tparam Surface_mesh is a model of `Surface_mesh`.
		\return `false` if error occurs, `true` otherwise.
		*/
		template <typename Surface_mesh>
		bool generate_candidate_faces(
			const Point_set_with_segments& planar_segments,	///< the input planar segments
			Surface_mesh& candidate_faces 					///< candidate faces by pairwise intersection of the planar segments
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


	private: // disallow copying
		Polygonal_surface_reconstruction(const Polygonal_surface_reconstruction& psr);

	}; // end of Polygonal_surface_reconstruction


	//////////////////////////////////////////////////////////////////////////s

	// implementations

	template <class Kernel>
	template <typename Surface_mesh>
	bool Polygonal_surface_reconstruction<Kernel>::reconstruct(
		const Point_list& point_list,
		Surface_mesh& output_mesh,
		double wt_fitting /* = 0.43 */,
		double wt_coverage /* = 0.27 */,
		double wt_complexity /* = 0.30 */)
	{
		Point_set_with_segments planar_segments;
		if (!extract_planes(point_list, planar_segments))
			return false;

		Surface_mesh candidate_faces;
		if (!generate_candidate_faces(planar_segments, candidate_faces))
			return false;

		output_mesh.clear();
		return select_faces(candidate_faces, output_mesh, wt_fitting, wt_coverage, wt_complexity);
	}


	template <class Kernel>
	bool Polygonal_surface_reconstruction<Kernel>::extract_planes(
		const Point_list& points,
		Point_set_with_segments& planar_segments)
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

		// Prints number of detected shapes.
		std::cout << ransac.shapes().end() - ransac.shapes().begin() << " shapes detected." << std::endl;

		return true;
	}


	template <class Kernel>
	template <typename Surface_mesh>
	bool Polygonal_surface_reconstruction<Kernel>::generate_candidate_faces(
		const Point_set_with_segments& planar_segments,
		Surface_mesh& candidate_faces
	)
	{
		return true;
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
		return true;
	}

} //namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_POLYGONAL_SURFACE_RECONSTRUCTION_H
