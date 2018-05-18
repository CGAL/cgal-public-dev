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
#include <CGAL/Point_set_3.h>
#include <CGAL/Surface_mesh.h>

/*!
  \file Polygonal_surface_reconstruction.h
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
		typedef typename Kernel::Plane_3					Plane;
		typedef typename Point_set_with_segments<Kernel>	Point_set;

	public:

		// \param point_set the point set that owns this planar segment.
		Planar_segment(Point_set* point_set) : point_set_(point_set) {}
		~Planar_segment() {}

		const Plane& supporting_plane() const { return plane_; }
		void set_supporting_plane(const Plane& plane) { plane_ = plane; }

	private:
		Point_set * point_set_;
		Plane		type_;
	};


	/** \ingroup PkgPolygonalSurfaceReconstruction
	*	An enriched point set that stores the extracted planar segments
	*/
	template <typename Kernel>
	class Point_set_with_segments : public Point_set_3<typename Kernel::Point_3>
	{
	public:
		typedef typename Planar_segment<Kernel>		PlanarSegment;

	public:
		Point_set_with_segments() {}
		~Point_set_with_segments() {}

		std::vector< std::unique_ptr<PlanarSegment> >& planar_segments() { return planar_segments_; }
		const std::vector< std::unique_ptr<PlanarSegment> >& planar_segments() const { return planar_segments_; }

	private:
		std::vector< std::unique_ptr<PlanarSegment> > planar_segments_;
	};



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

	\tparam Gt Geometric traits class.

	*/
	template <class Gt>
	class Polygonal_surface_reconstruction
	{
	public:

		typedef Gt Geom_traits; ///< Geometric traits class

		typedef typename Geom_traits::FT					FT;			
		typedef typename Geom_traits::Point_3				Point;		
		typedef typename Geom_traits::Vector_3				Vector;		
		typedef typename Geom_traits::Plane_3				Plane;		

		typedef CGAL::Point_set_3<Point>					Point_set;
		typedef CGAL::Point_set_with_segments<Geom_traits>	Point_set_with_segments;
		typedef CGAL::Surface_mesh<Point>					Surface_mesh;

		typedef typename Surface_mesh::Face_iterator		Face_iterator;
		typedef typename Surface_mesh::Edge_iterator		Edge_iterator;
		typedef typename Surface_mesh::Vertex_iterator		Vertex_iterator;

		typedef typename Surface_mesh::Face_index			Face_descriptor;
		typedef typename Surface_mesh::Edge_index			Edge_descriptor;
		typedef typename Surface_mesh::Vertex_index			Vertex_descriptor;

		// Public methods
	public:

		/// \name Creation 
		/*!
		Creates a Polygonal Surface Reconstruction object
		*/
		Polygonal_surface_reconstruction() {}


		/// \name Operations

		/** reconstruct from the input point set with normal property map.
		\remark
		This is a one-shot reconstruction. Please call the step-by-step functions if you want to reuse the intermediate results.
		\tparam Normal_map is a model of `ReadablePropertyMap` with value type `Vector_3<Kernel>`.
		\return `true` if plane extraction succeeded, `false` otherwise.
		*/
		template <typename Normal_map>
		bool reconstruct(
			const Point_set& point_set,  ///< input point set
			Normal_map normal_map, 		///< property map: `value_type of InputIterator` -> `Vector` (oriented or unoriented normal of an input point).
			Surface_mesh& output_mesh,			///< the final reconstruction results
			double wt_fitting = 0.43,			///< weight for the data fitting term.
			double wt_coverage = 0.27,			///< weight for the point coverage term.
			double wt_complexity = 0.30			///< weight for the model complexity term.
		);


		/** extracts planar segments from the input point set using RANSAC. 
		\tparam Normal_map is a model of `ReadablePropertyMap` with value type `Vector_3<Kernel>`.
		\return `true` if plane extraction succeeded, `false` otherwise.
		*/
		template <typename Normal_map>
		bool extract_planes(
			const Point_set& point_set,  ///< input point set
			Normal_map normal_map, 		///< property map: `value_type of InputIterator` -> `Vector` (oriented or unoriented normal of an input point).
			Point_set_with_segments& planar_segments	///< the extracted planar segments
		);


		/** generates candidate faces and computes confidence values for each face.
		The confidence values of the faces are stored as as property map with name "f:confidence". 
		\tparam Surface_mesh is a model of `Surface_mesh`.
		\return `false` if error occurs, `true` otherwise.
		*/
		bool generate_candidate_faces(
			const Point_set_with_segments& planar_segments,	///< the input planar segments
			Surface_mesh& candidate_faces 				///< candidate faces by pairwise intersection of the planar segments
		);


		/** select the optimal subset of the faces to assemble a watertight polygonal mesh model.

		\return `true` if optimization succeeded, `false` otherwise.
		*/
		bool select_faces(
			const Surface_mesh& candidate_faces, ///< candidate faces with face confidence values stored in property map "f:confidence"
			Surface_mesh& output_mesh,			///< the final reconstruction results
			double wt_fitting = 0.43,			///< weight for the data fitting term.
			double wt_coverage = 0.27,			///< weight for the point coverage term.
			double wt_complexity = 0.30			///< weight for the model complexity term.
		);

		// Data members.
	private:


	private: // disallow copying
		Polygonal_surface_reconstruction(const Polygonal_surface_reconstruction& psr) {}

	}; // end of Polygonal_surface_reconstruction


	//////////////////////////////////////////////////////////////////////////s

	// implementations

	template <class Gt>
	template <typename Normal_map>
	bool Polygonal_surface_reconstruction<Gt>::reconstruct(
		const Point_set& point_set,
		Normal_map normal_map,
		Surface_mesh& output_mesh,
		double wt_fitting /* = 0.43*/,
		double wt_coverage /* = 0.27*/,
		double wt_complexity /* = 0.30*/)
	{
		Point_set_with_segments planar_segments;
		if (!extract_planes(point_set, normal_map, planar_segments))
			return false;

		Surface_mesh candidate_faces;
		if (!generate_candidate_faces(planar_segments, candidate_faces))
			return false;

		output_mesh.clear();
		return select_faces(candidate_faces, output_mesh, wt_fitting, wt_coverage, wt_complexity);
	}


	template <class Gt>
	template <typename Normal_map>
	bool Polygonal_surface_reconstruction<Gt>::extract_planes(
		const Point_set& point_set,
		Normal_map normal_map,
		Point_set_with_segments& planar_segments
	)
	{
		return true;
	}


	template <class Gt>
	bool Polygonal_surface_reconstruction<Gt>::generate_candidate_faces(
		const Point_set_with_segments& planar_segments,
		Surface_mesh& candidate_faces
	)
	{
		return true;
	}

	template <class Gt>
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
