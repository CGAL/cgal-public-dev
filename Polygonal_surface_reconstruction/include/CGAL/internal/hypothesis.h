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

#ifndef CGAL_POLYGONAL_SURFACE_RECONSTRUCTION_HYPOTHESIS_H
#define CGAL_POLYGONAL_SURFACE_RECONSTRUCTION_HYPOTHESIS_H

#include <CGAL/Point_set_with_segments.h>
#include <CGAL/Surface_mesh.h>

#include <map>
#include <set>

namespace CGAL {

	template <typename Kernel>
	class Hypothesis
	{
	private:
		typedef typename Kernel::FT						FT;
		typedef typename Kernel::Point_3				Point;
		typedef typename Kernel::Vector_3				Vector;
		typedef typename Kernel::Plane_3				Plane;
		typedef CGAL::Planar_segment<Kernel>			Planar_segment;
		typedef CGAL::Point_set_with_segments<Kernel>	Point_set;

		typedef CGAL::Surface_mesh<Point>				Mesh;
		typedef typename Mesh::Face_index				Face_descriptor;
		typedef typename Mesh::Edge_index				Edge_descriptor;
		typedef typename Mesh::Vertex_index				Vertex_descriptor;

	public:
		Hypothesis(const Point_set* point_set) : point_set_(point_set) {}
		~Hypothesis();

		void generate(Mesh& mesh);

	private:
		// merge near co-planar segments
		void refine_planes();

		// construct a mesh representing the bounding box of the point set
		void construct_bbox_mesh(Mesh& mesh);

		// construct a mesh from the segments bounded by the bounding box mesh
		void construct_proxy_mesh(const Mesh& bbox_mesh, Mesh& mesh);

		// pairwise intersection
		void pairwise_intersection(Mesh& mesh);

	private:
		const Point_set * point_set_;

		// The intersection of the planes can be unreliable when the planes are near parallel. 
		// Here are the tricks we use in our implementation:
		//   - We first test if an intersection exists for every pair of planes. We then collect 
		//     plane triplets such that every pair in the plane triplet intersect. This is achieved 
		//     by testing each plane against the known intersecting pairs.
		//	 - The 3D vertices of the final faces are obtained by computing the intersections of
		//     the plane triplets. To cope with limited floating point precision, each vertex is 
		//     identified by the pointers of (in an increasing order) of the three planes from 
		//     which it is computed. By doing so, two vertices with almost identical positions can 
		//     be distinguished. This turned out to be quite robust in handling very close and near
		//     parallel planes.

		// the supporting planes of all planar segments and the bounding box faces
		std::vector<Plane*>			supporting_planes_;

		// precomputed intersecting point of all plane triplets
		std::map<Plane*, std::map<Plane*, std::map<Plane*, Point> > >	triplet_intersection_; 
	};


} //namespace CGAL

#include <CGAL/internal/hypothesis_impl.h>

#endif // CGAL_POLYGONAL_SURFACE_RECONSTRUCTION_HYPOTHESIS_H
