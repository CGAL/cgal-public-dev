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
#include <CGAL/Polygonal_surface_reconstruction_traits.h>
#include <CGAL/Point_set_with_segments.h>

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

	\tparam Traits a model of `Polygonal_surface_reconstruction_traits`
	*/
	template <class Traits>
	class Polygonal_surface_reconstruction
	{
	public:

		/// \name Types 

		/// @{
		/// \cond SKIP_IN_MANUAL
		typedef typename Traits::Input_range::iterator		Input_iterator;
		typedef typename Traits::FT							FT;			///< number type.
		typedef typename Traits::Point_3					Point;		///< point type.
		typedef typename Traits::Vector_3					Vector;		///< vector type.
		/// \endcond

		typedef typename Traits::Input_range				Input_range;
		///< model of the concept `Range` with random access iterators, providing input points and normals
		/// through the following two property maps.

		typedef typename Traits::Point_map					Point_map;
		///< property map to access the location of an input point.

		typedef typename Traits::Normal_map					Normal_map;
		///< property map to access the unoriented normal of an input point

		typedef typename Traits::Point_set_with_segments	Point_set_with_segments;
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
			const Input_range& input_range,			///< range of input data.
			Surface_mesh& output_mesh,				///< the final reconstruction results
			double wt_fitting = 0.43,				///< weight for the data fitting term.
			double wt_coverage = 0.27,				///< weight for the point coverage term.
			double wt_complexity = 0.30,			///< weight for the model complexity term.
			Point_map point_map = Point_map(),		///< property map to access the position of an input point.
			Normal_map normal_map = Normal_map()	///< property map to access the normal of an input point.
		);


		/** extracts planar segments from the input point set using RANSAC.
		\return `true` if plane extraction succeeded, `false` otherwise.
		*/
		bool extract_planes(
			const Input_range& input_range,				///< range of input data.
			Point_set_with_segments& planar_segments,	///< the extracted planar segments
			Point_map point_map = Point_map(),			///< property map to access the position of an input point.
			Normal_map normal_map = Normal_map()		///< property map to access the normal of an input point.
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

	template <class Gt>
	template <typename Surface_mesh>
	bool Polygonal_surface_reconstruction<Gt>::reconstruct(
		const Input_range& input_range,
		Surface_mesh& output_mesh,
		double wt_fitting /* = 0.43 */,
		double wt_coverage /* = 0.27 */,
		double wt_complexity /* = 0.30 */,
		Point_map point_map /* = Point_map() */,
		Normal_map normal_map /* = Normal_map() */)
	{
		Point_set_with_segments planar_segments;
		if (!extract_planes(input_range, planar_segments, point_map, normal_map))
			return false;

		Surface_mesh candidate_faces;
		if (!generate_candidate_faces(planar_segments, candidate_faces))
			return false;

		output_mesh.clear();
		return select_faces(candidate_faces, output_mesh, wt_fitting, wt_coverage, wt_complexity);
	}


	template <class Gt>
	bool Polygonal_surface_reconstruction<Gt>::extract_planes(
		const Input_range& input_range,
		Point_set_with_segments& planar_segments, 
		Point_map point_map /* = Point_map() */, 
		Normal_map normal_map /* = Normal_map() */ )
	{
		return true;
	}


	template <class Gt>
	template <typename Surface_mesh>
	bool Polygonal_surface_reconstruction<Gt>::generate_candidate_faces(
		const Point_set_with_segments& planar_segments,
		Surface_mesh& candidate_faces
	)
	{
		return true;
	}

	template <class Gt>
	template <typename Surface_mesh>
	bool Polygonal_surface_reconstruction<Gt>::select_faces(
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
