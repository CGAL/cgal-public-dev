#pragma once
#include <list>
#include <vector>
#include <cassert>
#include <iterator>
#include "defs.h"
#include "octree_base_vertex.h"


namespace Skippy 
{
	class Octree_Base
	{
	public:
		Octree_Base(const std::vector<double> & dims);

	protected:
		Octree_Base(const double _x_min, const double _x_max, const double _y_min, const double _y_max, const double _z_min, const double _z_max, Octree_Base* _parent);

	public:
		virtual ~Octree_Base();

		void add(Octree_Base_Vertex* v);

		void remove(Octree_Base_Vertex* v, bool destroy);

		void get_all_vertices(std::list<Octree_Base_Vertex*> & V) const;

	private:
		Octree_Base* assign(Octree_Base_Vertex* v) const;

	public:
		void search(const CGAL_Point_3 & M, const CGAL_Inexact_Point_3 & hint_M, std::list<Octree_Base_Vertex*> & R) const;

		void search(const std::vector<double> & dims, std::list<Octree_Base_Vertex*> & R) const;

		static Octree_Base* remove(Octree_Base* T);

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

		Octree_Base* parent;
		std::vector<Octree_Base*> children;

		std::list<Octree_Base_Vertex*> vertices;
	};
}