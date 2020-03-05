#pragma once
#include <list>
#include <vector>
#include "partition_objects.h"

namespace Skippy {
	class Partition_Vertex_Octree
	{
	public:
		Partition_Vertex_Octree(const std::vector<double> & dims);

	private:
		Partition_Vertex_Octree(const double _x_min, const double _x_max, const double _y_min, const double _y_max, const double _z_min, const double _z_max,
			Partition_Vertex_Octree* _parent);

	public:
		~Partition_Vertex_Octree();

		void add(Partition_Vertex* v);

		void remove(Partition_Vertex* v, bool destroy);

		bool exists(CGAL_Point_3 & M, const std::list<int> & P) const;

		Partition_Vertex* partial_match(const CGAL_Point_3 & M, const std::list<int> & P) const;

		Partition_Vertex* match(const CGAL_Point_3 & M, const std::list<int> & P) const;

		//void get_vertices_for_bounding_box_facet(const int i, std::list<Partition_Vertex *> & R) const;

		void get_vertices_for_boundary_plane(const int H_id, const int H_or, const CGAL_Plane & H, std::list<Partition_Vertex*> & R) const;

		void get_vertices_for_bounding_box_edge(const int i, const int j, 
			const int r_i, const int r_j, const CGAL_Plane & P_i, const CGAL_Plane & P_j,
			std::list<Partition_Vertex *> & R) const;

		void get_all_vertices(std::list<Partition_Vertex*> & V) const;

		void get_all_vertices_sorted_by_identifier(std::vector<Partition_Vertex*> & V) const;

	private:
		Partition_Vertex_Octree* assign(Partition_Vertex* v) const;

		void search(const CGAL_Point_3 & M, const CGAL_Inexact_Point_3 & hint_M, std::list<Partition_Vertex *> & R) const;

		void search(const std::vector<double> & dims, std::list<Partition_Vertex *> & R) const;

		void set_query(std::vector<double> & query, const int P_or, const CGAL_Plane & P) const;

		void set_query(std::vector<double> & query, const int i, const CGAL_Plane & P_i, const int j, const CGAL_Plane & P_j) const;

		static Partition_Vertex_Octree* remove(Partition_Vertex_Octree* T);

		bool is_empty_leaf() const;

		bool children_are_empty_leaves() const;

	public:
		bool is_node;
		const double eps;
		const double x_min;
		const double x_max;
		const double y_min;
		const double y_max;
		const double z_min;
		const double z_max;

		Partition_Vertex_Octree* parent;
		std::vector<Partition_Vertex_Octree*> children;

		std::list<Partition_Vertex*> vertices;
	};
}