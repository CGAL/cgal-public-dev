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

#ifndef CGAL_POLYGONAL_SURFACE_RECONSTRUCTION_FACE_SELECTION_H
#define CGAL_POLYGONAL_SURFACE_RECONSTRUCTION_FACE_SELECTION_H

#include <CGAL/Surface_mesh.h>
#include <CGAL/algo/hypothesis.h>
#include <CGAL/mip/mip_solver.h>

#include <set>
#include <map>

namespace CGAL {

	/** \ingroup PkgPolygonalSurfaceReconstruction
	*
	* 	Determines which faces will appear to assemble the final model.
	*	This is achieved by solving the linear integer program formulated
	*   in the original paper.
	*/

	template <typename Kernel>
	class Face_selection
	{
	private:
		typedef typename Kernel::FT						FT;
		typedef typename Kernel::Point_3				Point;
		typedef typename Kernel::Plane_3				Plane;
		typedef typename CGAL::Iso_cuboid_3<Kernel>		Box;

		typedef CGAL::Surface_mesh<Point>				Mesh;
		typedef typename Mesh::Face_index				Face_descriptor;
		typedef typename Mesh::Edge_index				Edge_descriptor;
		typedef typename Mesh::Vertex_index				Vertex_descriptor;
		typedef typename Mesh::Halfedge_index			Halfedge_descriptor;

		typedef typename Hypothesis<Kernel>::Intersection	Intersection;
		typedef typename std::vector<Intersection>			Adjacency;

	public:
		Face_selection() {}
		~Face_selection() {}

		bool optimize(
			const Mesh& candidate_faces,	///< candidate faces with face confidence values stored in property map "f:confidence"
			const Adjacency& adjacency,		///< the adjacency information captured during pairwise intersection
			Mesh& output_mesh,				///< the final reconstruction results
			double wt_fitting,				///< weight for the data fitting term.
			double wt_coverage,				///< weight for the point coverage term.
			double wt_complexity			///< weight for the model complexity term.
		);
	};


	template <class Kernel>
	bool Face_selection<Kernel>::optimize(
		const Mesh& candidate_faces,
		const Adjacency& adjacency,
		Mesh& mesh,
		double wt_fitting,
		double wt_coverage,
		double wt_complexity)
	{
		if (candidate_faces.num_faces() < 4)
			return false;

		mesh = candidate_faces;

		// the number of supporting points of each face
		typename Mesh::template Property_map<Face_descriptor, std::size_t> face_num_supporting_points =
			mesh.template add_property_map<Face_descriptor, std::size_t>("f:num_supporting_points").first;

		// the area of each face
		typename Mesh::template Property_map<Face_descriptor, FT> face_areas =
			mesh.template add_property_map<Face_descriptor, FT>("f:face_area").first;

		// the point covered area of each face
		typename Mesh::template Property_map<Face_descriptor, FT> face_covered_areas =
			mesh.template add_property_map<Face_descriptor, FT>("f:covered_area").first;

		// the supporting plane of each face
		typename Mesh::template Property_map<Face_descriptor, const Plane*> face_supporting_planes =
			mesh.template add_property_map<Face_descriptor, const Plane*>("f:supp_plane").first;

		// give each face an index
		typename Mesh::template Property_map<Face_descriptor, std::size_t> face_indices =
			mesh.template add_property_map<Face_descriptor, std::size_t>("f:index").first;

		double total_points = 0.0;
		std::size_t idx = 0;
		BOOST_FOREACH(Face_descriptor f, mesh.faces()) {
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

		std::size_t num_faces = mesh.number_of_faces();
		std::size_t num_edges(0);
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

		const typename Mesh::template Property_map<Vertex_descriptor, Point>& coords = mesh.points();
		std::vector<Point> vertices(mesh.number_of_vertices());
		idx = 0;
		BOOST_FOREACH(Vertex_descriptor v, mesh.vertices()) {
			vertices[idx] = coords[v];
			++idx;
		}
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

		BOOST_FOREACH(Face_descriptor f, mesh.faces()) {
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
				Face_descriptor f = mesh.face(fan[j]);
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
				Face_descriptor f1 = mesh.face(fan[j]);
				const Plane* plane1 = face_supporting_planes[f1];
				std::size_t fid1 = face_indices[f1];
				for (std::size_t k = j + 1; k < fan.size(); ++k) {
					Face_descriptor f2 = mesh.face(fan[k]);
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
			BOOST_FOREACH(Face_descriptor f, mesh.faces()) {
				if (static_cast<int>(std::round(X[f_idx])) == 0)
					to_delete.push_back(f);
				++f_idx;
			}

			for (std::size_t i = 0; i < to_delete.size(); ++i) {
				Face_descriptor f = to_delete[i];
				Halfedge_descriptor h = mesh.halfedge(f);
				Euler::remove_face(h, mesh);
			}

			CGAL_assertion(model.is_valid(true));

			// mark the sharp edges
			typename Mesh::template Property_map<Edge_descriptor, bool> edge_is_sharp =
				mesh.template add_property_map<Edge_descriptor, bool>("e:sharp_edges").first;
			BOOST_FOREACH(Edge_descriptor e, mesh.edges())
				edge_is_sharp[e] = false;

			for (std::size_t i = 0; i < adjacency.size(); ++i) {
				const Intersection& fan = adjacency[i];
				if (fan.size() != 4)
					continue;

				std::size_t idx_sharp_var = edge_sharp_status[&fan];
				if (static_cast<int>(X[idx_sharp_var]) == 1) {
					for (std::size_t j = 0; j < fan.size(); ++j) {
						Halfedge_descriptor h = fan[j];
						Face_descriptor f = mesh.face(h);
						if (f != Mesh::null_face()) { // some faces may be deleted
							std::size_t fid = face_indices[f];
							if (static_cast<int>(std::round(X[fid])) == 1) {
								Edge_descriptor e = mesh.edge(h);
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


#endif // CGAL_POLYGONAL_SURFACE_RECONSTRUCTION_FACE_SELECTION_H
