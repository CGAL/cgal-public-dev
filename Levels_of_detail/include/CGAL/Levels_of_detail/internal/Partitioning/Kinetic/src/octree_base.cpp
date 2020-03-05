#include "../include/octree_base.h"

namespace Skippy
{
	Octree_Base::Octree_Base(const std::vector<double> & dims)
		: Octree_Base(dims[0], dims[1], dims[2], dims[3], dims[4], dims[5], nullptr)
	{

	}


	Octree_Base::Octree_Base(const double _x_min, const double _x_max, const double _y_min, const double _y_max,
		const double _z_min, const double _z_max, Octree_Base* _parent)
		:
		is_node(false),
		eps(0.0001),
		x_min(_x_min),
		x_max(_x_max),
		y_min(_y_min),
		y_max(_y_max),
		z_min(_z_min),
		z_max(_z_max),
		parent(_parent)
	{
		children = std::vector<Octree_Base*>(8, nullptr);
		vertices = std::list<Octree_Base_Vertex*>();
	}


	Octree_Base::~Octree_Base()
	{
		if (is_node) {
			for (size_t i = 0; i < children.size(); i++) {
				delete children[i];
			}
			children.clear();
		}

		else {
			for (std::list<Octree_Base_Vertex*>::iterator it_v = vertices.begin(); it_v != vertices.end(); it_v++) {
				delete (*it_v);
			}
		}
	}


	void Octree_Base::add(Octree_Base_Vertex* v)
	{
		double d_eps = 0.0004;

		if (is_node) {
			// If the current object is a node, we insert the vertex into the right subtree
			this->assign(v)->add(v);
		}

		else {
			const double dx = x_max - x_min, dy = y_max - y_min, dz = z_max - z_min;
			bool capacity_exceeded = (vertices.size() > 20);

			// If the current node hasn't exceed its maximal capacity, or if it can no longer be divided,
			// we insert the vertex
			if (!capacity_exceeded || (capacity_exceeded && (dz < d_eps || dy < d_eps || dx < d_eps))) {
				vertices.push_back(v);

			} else {
				is_node = true;

				// Otherwise, we divide the leaf into eight subtrees
				children[0] = new Octree_Base(x_min, x_min + dx / 2, y_min, y_min + dy / 2, z_min, z_min + dz / 2, this);
				children[1] = new Octree_Base(x_min + dx / 2, x_max, y_min, y_min + dy / 2, z_min, z_min + dz / 2, this);
				children[2] = new Octree_Base(x_min, x_min + dx / 2, y_min + dy / 2, y_max, z_min, z_min + dz / 2, this);
				children[3] = new Octree_Base(x_min + dx / 2, x_max, y_min + dy / 2, y_max, z_min, z_min + dz / 2, this);
				children[4] = new Octree_Base(x_min, x_min + dx / 2, y_min, y_min + dy / 2, z_min + dz / 2, z_max, this);
				children[5] = new Octree_Base(x_min + dx / 2, x_max, y_min, y_min + dy / 2, z_min + dz / 2, z_max, this);
				children[6] = new Octree_Base(x_min, x_min + dx / 2, y_min + dy / 2, y_max, z_min + dz / 2, z_max, this);
				children[7] = new Octree_Base(x_min + dx / 2, x_max, y_min + dy / 2, y_max, z_min + dz / 2, z_max, this);

				// We divide the list of points and assign them to one of the subtrees according to their coordinates
				for (std::list<Octree_Base_Vertex*>::iterator it_v = vertices.begin(); it_v != vertices.end(); it_v++) {
					Octree_Base_Vertex* u = (*it_v);
					this->assign(u)->add(u);
				}
				vertices.clear();

				// Finally, we add v to the corresponding node
				this->assign(v)->add(v);
			}
		}
	}


	void Octree_Base::remove(Octree_Base_Vertex* v, bool destroy)
	{
		// Assumes this function is called from the octree's root
		assert(parent == nullptr);

		// First we search for the subtree where v is located.
		// Once it's done, we remove v from the list of vertices.

		Octree_Base* T = this;
		while (T->is_node) {
			T = T->assign(v);
		}

		for (std::list<Octree_Base_Vertex*>::iterator it_v = T->vertices.begin(); it_v != T->vertices.end(); it_v++) {
			if ((*it_v) == v) {
				if (destroy) delete v;
				it_v = T->vertices.erase(it_v);
				break;
			}
		}

		// Second, we check if T is an empty leaf
		if (T->is_empty_leaf()) {

			// Loops on the successive parent nodes of T
			// At the end of this iterative process, either T_p == nullptr, because the root is reached,
			// or T == T_p because at least one of the children of T is not an empty leaf

			Octree_Base* T_p = T->parent;
			while (T_p != nullptr && T_p != T) {
				T = remove(T_p);
				if (T != T_p) {
					if (T != nullptr) {
						T_p = T->parent;
					} else {
						T_p = nullptr;
					}
				}
			}
		}
	}


	Octree_Base* Octree_Base::remove(Octree_Base* T)
	{
		assert(T != nullptr);

		// Given an octree, we check if it is a node with eight empty leaves.
		// If so, this octree is turned into a leaf, and we return a pointer to the octree's parent,
		// so that the same check can be performed on it.

		if (T->is_node && T->children_are_empty_leaves()) {

			// Obtains a pointer to the parent tree
			Octree_Base* T_par = T->parent;

			// Deletes the subtrees and returns a pointer to the parent tree
			for (int i = 0; i < 8; i++) {
				delete T->children[i];
				T->children[i] = nullptr;
			}
			T->is_node = false;

			// Returns a pointer to the parent node
			return T_par;

		} else {
			// Well, we do nothing
			return T;
		}
	}


	void Octree_Base::get_all_vertices(std::list<Octree_Base_Vertex*> & V) const
	{
		if (is_node) {
			for (int i = 0; i < 8; i++) children[i]->get_all_vertices(V);
		} else {
			std::copy(vertices.begin(), vertices.end(), std::back_inserter(V));
		}
	}


	bool Octree_Base::is_empty_leaf() const
	{
		return (!is_node && vertices.empty());
	}


	bool Octree_Base::children_are_empty_leaves() const
	{
		for (int i = 0; i < children.size(); i++) {
			if (children[i]->is_empty_leaf()) return false;
		}
		return true;
	}


	Octree_Base* Octree_Base::assign(Octree_Base_Vertex* v) const
	{
		const CGAL_Inexact_Point_3 & M = v->hint_M;
		const double dx = x_max - x_min;
		const double dy = y_max - y_min;
		const double dz = z_max - z_min;

		int i = 0;
		i += (M.x() < x_min + dx / 2 ? 0 : 1);
		i += (M.y() < y_min + dy / 2 ? 0 : 2);
		i += (M.z() < z_min + dz / 2 ? 0 : 4);

		return children[i];
	}


	void Octree_Base::search(const CGAL_Point_3 & M, const CGAL_Inexact_Point_3 & hint_M, std::list<Octree_Base_Vertex*> & R) const
	{
		const double a_min = hint_M.x() - eps, a_max = hint_M.x() + eps;
		const double b_min = hint_M.y() - eps, b_max = hint_M.y() + eps;
		const double c_min = hint_M.z() - eps, c_max = hint_M.z() + eps;

		if (((jin(x_min, a_min, x_max) || jin(x_min, a_max, x_max)) || (jin(a_min, x_min, a_max) || jin(a_min, x_max, a_max)))
			&& ((jin(y_min, b_min, y_max) || jin(y_min, b_max, y_max)) || (jin(b_min, y_min, b_max) || jin(b_min, y_max, b_max)))
			&& ((jin(z_min, c_min, z_max) || jin(z_min, c_max, z_max)) || (jin(c_min, z_min, c_max) || jin(c_min, z_max, c_max)))) {

			if (is_node) {
				// If the quadtree intercepts the area of search, and if this tree is a node,
				// we call the function recursively on each of the four subtrees
				for (int i = 0; i < 8; i++) children[i]->search(M, hint_M, R);

			} else {
				// Otherwise, if the quadtree is a leaf, we loop on the list of points
				// If a vertex is contained in the area of search, then we append it to the list of vertices
				for (std::list<Octree_Base_Vertex *>::const_iterator it_v = vertices.begin(); it_v != vertices.end(); it_v++) {
					Octree_Base_Vertex* v = (*it_v);
					double dx = v->hint_M.x() - hint_M.x(), dy = v->hint_M.y() - hint_M.y(), dz = v->hint_M.z() - hint_M.z();
					double squared_length = dx * dx + dy * dy + dz * dz;
					if (squared_length < eps) {
						R.push_back(v);
					}
				}
			}
		}
	}


	void Octree_Base::search(const std::vector<double> & dims, std::list<Octree_Base_Vertex *> & R) const
	{
		if (((jin(x_min, dims[0], x_max) || jin(x_min, dims[1], x_max)) || (jin(dims[0], x_min, dims[1]) || jin(dims[0], x_max, dims[1])))
			&& ((jin(y_min, dims[2], y_max) || jin(y_min, dims[3], y_max)) || (jin(dims[2], y_min, dims[3]) || jin(dims[2], y_max, dims[3])))
			&& ((jin(z_min, dims[4], z_max) || jin(z_min, dims[5], z_max)) || (jin(dims[4], z_min, dims[5]) || jin(dims[4], z_max, dims[5])))) {

			if (is_node) {
				// If the quadtree intercepts the area of search, and if this tree is a node,
				// we call the function recursively on each of the four subtrees
				for (int i = 0; i < 8; i++) children[i]->search(dims, R);

			} else {
				// Otherwise, if the quadtree is a leaf, we loop on the list of points
				// If a vertex is contained in the area of search, then we append it to the list of vertices
				for (std::list<Octree_Base_Vertex *>::const_iterator it_v = vertices.begin(); it_v != vertices.end(); it_v++) {
					Octree_Base_Vertex* v = (*it_v);
					const double x = v->hint_M.x(), y = v->hint_M.y(), z = v->hint_M.z();
					if (jin(dims[0], x, dims[1]) && jin(dims[2], y, dims[3]) && jin(dims[4], z, dims[5])) {
						R.push_back(v);
					}
				}
			}
		}
	}
}