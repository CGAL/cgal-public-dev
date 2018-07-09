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

#ifndef CGAL_POLYGONAL_SURFACE_RECONSTRUCTION_H
#define CGAL_POLYGONAL_SURFACE_RECONSTRUCTION_H

#include <CGAL/Point_set_with_segments.h>
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

		typedef typename Kernel::FT				FT;			///< number type.
		typedef typename Kernel::Point_3		Point;		///< point type.
		typedef typename Kernel::Vector_3		Vector;		///< vector type.
		typedef std::pair<Point, Vector>		Point_with_normal;
		typedef std::vector<Point_with_normal>	Point_list;

	private:

		typedef CGAL::Planar_segment<Kernel>			Planar_segment;

		typedef CGAL::Point_set_with_segments<Kernel>	Point_set_with_segments;
		///< enriched point set to access the extracted planar segments

		// Public methods
	public:

		/// \name Creation 
		/*!
		Creates a Polygonal Surface Reconstruction object
		*/
		Polygonal_surface_reconstruction(
			const Point_set_with_segments& point_set_with_planes ///< input point set with planes.
		);

		/// \name Operations

		~Polygonal_surface_reconstruction() {
			delete hypothesis_;
		}

		/** select the optimal subset of the faces to assemble a watertight polygonal mesh model.
		\tparam PolygonMesh a model of `FaceGraph`
		\return `true` if optimization succeeded, `false` otherwise.
		*/
		template <typename PolygonMesh>
		bool reconstruct(
			PolygonMesh& output_mesh,		///< the final reconstruction result
			double wt_fitting = 0.43,		///< weight for the data fitting term.
			double wt_coverage = 0.27,		///< weight for the point coverage term.
			double wt_complexity = 0.30		///< weight for the model complexity term.
		);

		/*! Give the user the possibility to get the candidate faces.
		\tparam PolygonMesh a model of `FaceGraph`.
		*/
		template <typename PolygonMesh>
		void output_candidate_faces(PolygonMesh& candidate_faces) const { 
			candidate_faces = hypothesis_->candidate_faces(); 
		}

		// Data members.
	private:
		internal::Hypothesis<Kernel> * hypothesis_;

	private: // disallow copying
		Polygonal_surface_reconstruction(const Polygonal_surface_reconstruction& psr);

	}; // end of Polygonal_surface_reconstruction


	   //////////////////////////////////////////////////////////////////////////s

	   // implementations
	template <class Kernel>
	Polygonal_surface_reconstruction<Kernel>::Polygonal_surface_reconstruction(
		const Point_set_with_segments& point_set_with_planes
	) 
	{
		const std::vector< Planar_segment* >& planar_segments = point_set_with_planes.planar_segments();
		if (planar_segments.size() < 4) {
			std::cerr << "not enough (" << planar_segments.size() << ") planar segments to"
				<< " reconstruct a closed surface mesh" << std::endl;
			return;
		}

		hypothesis_ = new internal::Hypothesis<Kernel>(&point_set_with_planes);
		hypothesis_->generate();

		typedef typename internal::Hypothesis<Kernel>::Polygon_mesh		Polygon_mesh;
		Polygon_mesh& candidate_faces = hypothesis_->candidate_faces();

		typedef internal::Candidate_confidences<Kernel>		Candidate_confidences;
		Candidate_confidences conf;
		conf.compute(point_set_with_planes, candidate_faces);
	}


	template <class Kernel>
	template <typename PolygonMesh>
	bool Polygonal_surface_reconstruction<Kernel>::reconstruct(
		PolygonMesh& output_mesh,
		double wt_fitting /* = 0.43 */,
		double wt_coverage /* = 0.27 */,
		double wt_complexity /* = 0.30 */)
	{
		output_mesh.clear();

		typedef internal::Hypothesis<Kernel>		Hypothesis;
		typedef typename Hypothesis::Polygon_mesh	Polygon_mesh;
		const Polygon_mesh& candidate_faces = hypothesis_->candidate_faces();

		typedef typename Hypothesis::Adjacency Adjacency;
		const Adjacency& adjacency = hypothesis_->extract_adjacency();

		typedef internal::Face_selection<Kernel>	Face_selection;

		Face_selection sel;
		return sel.optimize(candidate_faces, adjacency, output_mesh, wt_fitting, wt_coverage, wt_complexity);
	}

} //namespace CGAL

#endif // CGAL_POLYGONAL_SURFACE_RECONSTRUCTION_H
