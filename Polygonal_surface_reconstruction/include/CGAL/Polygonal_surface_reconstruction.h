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

#include <CGAL/algo/hypothesis.h>
#include <CGAL/algo/compute_confidences.h>
#include <CGAL/algo/point_set_with_segments.h>
#include <CGAL/bounding_box.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/mip/mip_solver.h>

#include <map>

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
		typedef typename Kernel::Plane_3		Plane;		///< plane type.

	private:

		typedef CGAL::Planar_segment<Kernel>			Planar_segment;

		// enriched point set to access the extracted planar segments
		typedef CGAL::Point_set_with_segments<Kernel>	Point_set_with_segments;

		// polygon mesh storing the candidate faces
		typedef CGAL::Surface_mesh<Point>				Polygon_mesh; 

		typedef typename Polygon_mesh::Face_index		Face_descriptor;
		typedef typename Polygon_mesh::Vertex_index		Vertex_descriptor;
		typedef typename Polygon_mesh::Edge_index		Edge_descriptor;
		typedef typename Polygon_mesh::Halfedge_index	Halfedge_descriptor;

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

		/*! Give the user the possibility to access the intermediate candidate faces.
		\tparam PolygonMesh a model of `FaceGraph`.
		*/
		template <typename PolygonMesh>
		void output_candidate_faces(PolygonMesh& candidate_faces) const { 
			candidate_faces = hypothesis_->candidate_faces(); 
		}

		// Data members.
	private:
		internal::Hypothesis<Kernel> * hypothesis_;

		// the generated candidate faces stored as a polygon mesh
		Polygon_mesh candidate_faces_;

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
		hypothesis_->generate(candidate_faces_);

		typedef internal::Candidate_confidences<Kernel>		Candidate_confidences;
		Candidate_confidences conf;
		conf.compute(point_set_with_planes, candidate_faces_);
	}


	template <class Kernel>
	template <typename PolygonMesh>
	bool Polygonal_surface_reconstruction<Kernel>::reconstruct(
		PolygonMesh& output_mesh,
		double wt_fitting /* = 0.43 */,
		double wt_coverage /* = 0.27 */,
		double wt_complexity /* = 0.30 */)
	{
		if (candidate_faces_.num_faces() < 4) {
			std::cerr << "not enough (" << candidate_faces_.num_faces() << ") candidate faces to"
				<< " reconstruct a closed surface mesh" << std::endl;
			return false;
		}

		typedef internal::Hypothesis<Kernel>::Adjacency Adjacency;
		const Adjacency& adjacency = hypothesis_->extract_adjacency(candidate_faces_);

		output_mesh = candidate_faces_;

		// the number of supporting points of each face
		typename Polygon_mesh::template Property_map<Face_descriptor, std::size_t> face_num_supporting_points =
			output_mesh.template add_property_map<Face_descriptor, std::size_t>("f:num_supporting_points").first;

		// the area of each face
		typename Polygon_mesh::template Property_map<Face_descriptor, FT> face_areas =
			output_mesh.template add_property_map<Face_descriptor, FT>("f:face_area").first;

		// the point covered area of each face
		typename Polygon_mesh::template Property_map<Face_descriptor, FT> face_covered_areas =
			output_mesh.template add_property_map<Face_descriptor, FT>("f:covered_area").first;

		// the supporting plane of each face
		typename Polygon_mesh::template Property_map<Face_descriptor, const Plane*> face_supporting_planes =
			output_mesh.template add_property_map<Face_descriptor, const Plane*>("f:supp_plane").first;

		// give each face an index
		typename Polygon_mesh::template Property_map<Face_descriptor, std::size_t> face_indices =
			output_mesh.template add_property_map<Face_descriptor, std::size_t>("f:index").first;

		double total_points = 0.0;
		std::size_t idx = 0;
		BOOST_FOREACH(Face_descriptor f, output_mesh.faces()) {
			total_points += face_num_supporting_points[f];
			face_indices[f] = idx;
			++idx;
		}

		MIP_model model;

		// add variables

		// binary variables:
		// x[0] ... x[num_faces - 1] : binary labels of all the input faces
		// x[num_faces] ... x[num_faces + num_edges - 1] : binary labels of all the intersecting edges (remain or not)
		// x[num_faces + num_edges] ... x[num_faces + num_edges + num_edges] : binary labels of corner edges (sharp edge of not)

		std::size_t num_faces = output_mesh.number_of_faces();
		std::size_t num_edges(0);

		typedef typename internal::Hypothesis<Kernel>::Intersection	Intersection;

		std::map<const Intersection*, std::size_t> edge_usage_status;	// keep or remove an intersecting edges
		for (std::size_t i = 0; i < adjacency.size(); ++i) {
			const Intersection& fan = adjacency[i];
			if (fan.size() == 4) {
				std::size_t var_idx = num_faces + num_edges;
				edge_usage_status[&fan] = var_idx;
				++num_edges;
			}
		}

		std::size_t total_variables = num_faces + num_edges + num_edges;

		const std::vector<Variable*>& variables = model.create_n_variables(total_variables);
		for (std::size_t i = 0; i < total_variables; ++i) {
			Variable* v = variables[i];
			v->set_variable_type(Variable::BINARY);
		}

		// add objective

		const typename Polygon_mesh::template Property_map<Vertex_descriptor, Point>& coords = output_mesh.points();
		std::vector<Point> vertices(output_mesh.number_of_vertices());
		idx = 0;
		BOOST_FOREACH(Vertex_descriptor v, output_mesh.vertices()) {
			vertices[idx] = coords[v];
			++idx;
		}

		typedef typename CGAL::Iso_cuboid_3<Kernel>		Box;

		const Box& box = CGAL::bounding_box(vertices.begin(), vertices.end());
		FT dx = box.xmax() - box.xmin();
		FT dy = box.ymax() - box.ymin();
		FT dz = box.zmax() - box.zmin();
		FT box_area = FT(2.0) * (dx * dy + dy * dz + dz * dx);

		// choose a better scale: all actual values multiplied by total number of points
		double coeff_data_fitting = wt_fitting;
		double coeff_coverage = total_points * wt_coverage / box_area;
		double coeff_complexity = total_points * wt_complexity / double(adjacency.size());

		Linear_objective * objective = model.create_objective(Linear_objective::MINIMIZE);

		std::map<const Intersection*, std::size_t> edge_sharp_status;	// the edge is sharp or not
		std::size_t num_sharp_edges = 0;
		for (std::size_t i = 0; i < adjacency.size(); ++i) {
			const Intersection& fan = adjacency[i];
			if (fan.size() == 4) {
				std::size_t var_idx = num_faces + num_edges + num_sharp_edges;
				edge_sharp_status[&fan] = var_idx;

				// accumulate model complexity term
				objective->add_coefficient(variables[var_idx], coeff_complexity);
				++num_sharp_edges;
			}
		}
		CGAL_assertion(num_edges == num_sharp_edges);

		BOOST_FOREACH(Face_descriptor f, output_mesh.faces()) {
			std::size_t var_idx = face_indices[f];

			// accumulate data fitting term
			std::size_t num = face_num_supporting_points[f];
			objective->add_coefficient(variables[var_idx], -coeff_data_fitting * num);

			// accumulate model coverage term
			double uncovered_area = (face_areas[f] - face_covered_areas[f]);
			objective->add_coefficient(variables[var_idx], coeff_coverage * uncovered_area);
		}

		// Add constraints: the number of faces associated with an edge must be either 2 or 0
		std::size_t var_edge_used_idx = 0;
		for (std::size_t i = 0; i < adjacency.size(); ++i) {
			Linear_constraint* c = model.create_constraint(0.0, 0.0);
			const Intersection& fan = adjacency[i];
			for (std::size_t j = 0; j < fan.size(); ++j) {
				Face_descriptor f = output_mesh.face(fan[j]);
				std::size_t var_idx = face_indices[f];
				c->add_coefficient(variables[var_idx], 1.0);
			}

			if (fan.size() == 4) {
				std::size_t var_idx = num_faces + var_edge_used_idx;
				c->add_coefficient(variables[var_idx], -2.0);  // 
				++var_edge_used_idx;
			}
			else { // boundary edge
				   // will be set to 0 (i.e., we don't allow open surface)
			}
		}

		// Add constraints: for the sharp edges. The explanation of posing this constraint can be found here:
		// https://user-images.githubusercontent.com/15526536/30185644-12085a9c-942b-11e7-831d-290dd2a4d50c.png
		double M = 1.0;
		for (std::size_t i = 0; i < adjacency.size(); ++i) {
			const Intersection& fan = adjacency[i];
			if (fan.size() != 4)
				continue;

			// if an edge is sharp, the edge must be selected first:
			// X[var_edge_usage_idx] >= X[var_edge_sharp_idx]	
			Linear_constraint* c = model.create_constraint(0.0);
			std::size_t var_edge_usage_idx = edge_usage_status[&fan];
			c->add_coefficient(variables[var_edge_usage_idx], 1.0);
			std::size_t var_edge_sharp_idx = edge_sharp_status[&fan];
			c->add_coefficient(variables[var_edge_sharp_idx], -1.0);

			for (std::size_t j = 0; j < fan.size(); ++j) {
				Face_descriptor f1 = output_mesh.face(fan[j]);
				const Plane* plane1 = face_supporting_planes[f1];
				std::size_t fid1 = face_indices[f1];
				for (std::size_t k = j + 1; k < fan.size(); ++k) {
					Face_descriptor f2 = output_mesh.face(fan[k]);
					const Plane* plane2 = face_supporting_planes[f2];
					std::size_t fid2 = face_indices[f2];

					if (plane1 != plane2) {
						// the constraint is:
						//X[var_edge_sharp_idx] + M * (3 - (X[fid1] + X[fid2] + X[var_edge_usage_idx])) >= 1
						// which equals to  
						//X[var_edge_sharp_idx] - M * X[fid1] - M * X[fid2] - M * X[var_edge_usage_idx] >= 1 - 3M
						c = model.create_constraint(1.0 - 3.0 * M);
						c->add_coefficient(variables[var_edge_sharp_idx], 1.0);
						c->add_coefficient(variables[fid1], -M);
						c->add_coefficient(variables[fid2], -M);
						c->add_coefficient(variables[var_edge_usage_idx], -M);
					}
				}
			}
		}

		// Optimize 

		MIP_solver solver;
		if (solver.solve(&model)) {

			// mark results
			const std::vector<double>& X = solver.solution();

			std::vector<Face_descriptor> to_delete;
			std::size_t f_idx(0);
			BOOST_FOREACH(Face_descriptor f, output_mesh.faces()) {
				if (static_cast<int>(std::round(X[f_idx])) == 0)
					to_delete.push_back(f);
				++f_idx;
			}

			for (std::size_t i = 0; i < to_delete.size(); ++i) {
				Face_descriptor f = to_delete[i];
				Halfedge_descriptor h = output_mesh.halfedge(f);
				Euler::remove_face(h, output_mesh);
			}

			CGAL_assertion(model.is_valid(true));

			// mark the sharp edges
			typename Polygon_mesh::template Property_map<Edge_descriptor, bool> edge_is_sharp =
				output_mesh.template add_property_map<Edge_descriptor, bool>("e:sharp_edges").first;
			BOOST_FOREACH(Edge_descriptor e, output_mesh.edges())
				edge_is_sharp[e] = false;

			for (std::size_t i = 0; i < adjacency.size(); ++i) {
				const Intersection& fan = adjacency[i];
				if (fan.size() != 4)
					continue;

				std::size_t idx_sharp_var = edge_sharp_status[&fan];
				if (static_cast<int>(X[idx_sharp_var]) == 1) {
					for (std::size_t j = 0; j < fan.size(); ++j) {
						Halfedge_descriptor h = fan[j];
						Face_descriptor f = output_mesh.face(h);
						if (f != Polygon_mesh::null_face()) { // some faces may be deleted
							std::size_t fid = face_indices[f];
							if (static_cast<int>(std::round(X[fid])) == 1) {
								Edge_descriptor e = output_mesh.edge(h);
								edge_is_sharp[e] = true;
								break;
							}
						}
					}
				}
			}
		}
		else {
			std::cerr << "solving the binary program failed." << std::endl;
			return false;
		}

		return true;
	}

} //namespace CGAL

#endif // CGAL_POLYGONAL_SURFACE_RECONSTRUCTION_H
