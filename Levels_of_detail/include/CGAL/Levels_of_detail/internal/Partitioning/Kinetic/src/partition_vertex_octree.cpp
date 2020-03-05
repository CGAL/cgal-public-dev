#include "../include/partition_vertex_octree.h"

namespace Skippy {
	Partition_Vertex_Octree::Partition_Vertex_Octree(const double _x_min, const double _x_max, const double _y_min, const double _y_max,
		const double _z_min, const double _z_max, Partition_Vertex_Octree* _parent)
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
		children = std::vector<Partition_Vertex_Octree*>(8, nullptr);
		vertices = std::list<Partition_Vertex*>();
	}


	Partition_Vertex_Octree::Partition_Vertex_Octree(const std::vector<double> & dims)
		: Partition_Vertex_Octree(dims[0], dims[1], dims[2], dims[3], dims[4], dims[5], nullptr)
	{

	}


	Partition_Vertex_Octree::~Partition_Vertex_Octree()
	{
		if (is_node) {
			for (size_t i = 0; i < children.size(); i++) {
				delete children[i];
			}
			children.clear();
		}

		else {
			for (std::list<Partition_Vertex*>::iterator it_v = vertices.begin(); it_v != vertices.end(); it_v++) {
				delete (*it_v);
			}
		}
	}



	void Partition_Vertex_Octree::add(Partition_Vertex* v)
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
				children[0] = new Partition_Vertex_Octree(x_min, x_min + dx / 2, y_min, y_min + dy / 2, z_min, z_min + dz / 2, this);
				children[1] = new Partition_Vertex_Octree(x_min + dx / 2, x_max, y_min, y_min + dy / 2, z_min, z_min + dz / 2, this);
				children[2] = new Partition_Vertex_Octree(x_min, x_min + dx / 2, y_min + dy / 2, y_max, z_min, z_min + dz / 2, this);
				children[3] = new Partition_Vertex_Octree(x_min + dx / 2, x_max, y_min + dy / 2, y_max, z_min, z_min + dz / 2, this);
				children[4] = new Partition_Vertex_Octree(x_min, x_min + dx / 2, y_min, y_min + dy / 2, z_min + dz / 2, z_max, this);
				children[5] = new Partition_Vertex_Octree(x_min + dx / 2, x_max, y_min, y_min + dy / 2, z_min + dz / 2, z_max, this);
				children[6] = new Partition_Vertex_Octree(x_min, x_min + dx / 2, y_min + dy / 2, y_max, z_min + dz / 2, z_max, this);
				children[7] = new Partition_Vertex_Octree(x_min + dx / 2, x_max, y_min + dy / 2, y_max, z_min + dz / 2, z_max, this);

				// We divide the list of points and assign them to one of the subtrees according to their coordinates
				for (std::list<Partition_Vertex *>::iterator it_v = vertices.begin(); it_v != vertices.end(); it_v++) {
					Partition_Vertex* u = (*it_v);
					this->assign(u)->add(u);
				}
				vertices.clear();

				// Finally, we add v to the corresponding node
				this->assign(v)->add(v);
			}
		}
	}



	void Partition_Vertex_Octree::remove(Partition_Vertex* v, bool destroy)
	{
		// Assumes this function is called from the octree's root
		assert(parent == nullptr);

		// First we search for the subtree where v is located.
		// Once it's done, we remove v from the list of vertices.

		Partition_Vertex_Octree* T = this;
		while (T->is_node) {
			T = T->assign(v);
		}

		for (std::list<Partition_Vertex*>::iterator it_v = T->vertices.begin(); it_v != T->vertices.end(); it_v++) {
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

			Partition_Vertex_Octree* T_p = T->parent;
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



	Partition_Vertex_Octree* Partition_Vertex_Octree::remove(Partition_Vertex_Octree* T)
	{
		assert(T != nullptr);

		// Given an octree, we check if it is a node with eight empty leaves.
		// If so, this octree is turned into a leaf, and we return a pointer to the octree's parent,
		// so that the same check can be performed on it.

		if (T->is_node && T->children_are_empty_leaves()) {

			// Obtains a pointer to the parent tree
			Partition_Vertex_Octree* T_par = T->parent;

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



	bool Partition_Vertex_Octree::is_empty_leaf() const
	{
		return (!is_node && vertices.empty());
	}


	bool Partition_Vertex_Octree::children_are_empty_leaves() const
	{
		for (int i = 0; i < children.size(); i++) {
			if (children[i]->is_empty_leaf()) return false;
		}
		return true;
	}



	Partition_Vertex_Octree* Partition_Vertex_Octree::assign(Partition_Vertex* v) const
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


	void Partition_Vertex_Octree::search(const CGAL_Point_3 & M, const CGAL_Inexact_Point_3 & hint_M, std::list<Partition_Vertex *> & R) const
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
				for (std::list<Partition_Vertex *>::const_iterator it_v = vertices.begin(); it_v != vertices.end(); it_v++) {
					Partition_Vertex* v = (*it_v);
					if (CGAL::squared_distance(v->hint_M, hint_M) < eps) {
						R.push_back(v);
					}
				}
			}
		}
	}


	void Partition_Vertex_Octree::search(const std::vector<double> & dims, std::list<Partition_Vertex *> & R) const
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
				for (std::list<Partition_Vertex *>::const_iterator it_v = vertices.begin(); it_v != vertices.end(); it_v++) {
					Partition_Vertex* v = (*it_v);
					const double x = v->hint_M.x(), y = v->hint_M.y(), z = v->hint_M.z();
					if (jin(dims[0], x, dims[1]) && jin(dims[2], y, dims[3]) && jin(dims[4], z, dims[5])) {
						R.push_back(v);
					}
				}
			}
		}
	}


	bool Partition_Vertex_Octree::exists(CGAL_Point_3 & M, const std::list<int> & P) const
	{
		return (match(M, P) != nullptr);
	}



	Partition_Vertex* Partition_Vertex_Octree::match(const CGAL_Point_3 & M, const std::list<int> & P) const
	{
		// We test the existence of a Partition_Vertex with the same coordinates as M
		// and located at the intersection of the points listed in P
		std::list<Partition_Vertex*> V;

		double prec = 1.0 / (1 << 30) / (1 << 10);
		FT::set_relative_precision_of_to_double(prec);

		double x = to_double(M.x()), y = to_double(M.y()), z = to_double(M.z());
		CGAL_Inexact_Point_3 hint_M(x, y, z);
		search(M, hint_M, V);

		// We got a list of candidates, then we check if one of them has the same list P
		for (std::list<Partition_Vertex*>::iterator it_v = V.begin(); it_v != V.end(); it_v++) {
			if ((*it_v)->definitions_match(P)) {
				return (*it_v);
			}
		}

		return nullptr;
	}



	Partition_Vertex* Partition_Vertex_Octree::partial_match(const CGAL_Point_3 & M, const std::list<int> & P) const
	{
		// We test the existence of a Partition_Vertex with the same coordinates as M
		// and located at the intersection of the points listed in P
		std::list<Partition_Vertex*> V;

		double prec = 1.0 / (1 << 30) / (1 << 10);
		FT::set_relative_precision_of_to_double(prec);

		double x = to_double(M.x()), y = to_double(M.y()), z = to_double(M.z());
		CGAL_Inexact_Point_3 hint_M(x, y, z);
		search(M, hint_M, V);

		// We got a list of candidates, then we check if one of them has the same list P
		for (std::list<Partition_Vertex*>::iterator it_v = V.begin(); it_v != V.end(); it_v++) {
			if ((*it_v)->definitions_partially_match(P)) {
				return (*it_v);
			}
		}

		return nullptr;
	}



	void Partition_Vertex_Octree::set_query(std::vector<double> & query, const int i, const CGAL_Plane & P_i, const int j, const CGAL_Plane & P_j) const
	{
		// Initializes a query that corresponds to the entire bounding box
		double _query[6] = { x_min - eps, x_max + eps, y_min - eps, y_max + eps, z_min - eps, z_max + eps };
		query = std::vector<double>(_query, _query + sizeof(_query) / sizeof(double));

		// Refines the query by restricting query[.] to the dimensions of the planes i and j
		for (int k = 0; k <= 1; k++) {
			int l = (k == 0 ? i : j);
			const CGAL_Plane & P_l = (k == 0 ? P_i : P_j);
			
			double d = -CGAL::to_double(P_l.d());

			if (l <= 1) {
				query[0] = d - eps, query[1] = d + eps;
			} else if (l <= 3) {
				query[2] = d - eps, query[3] = d + eps;
			} else {
				query[4] = d - eps, query[5] = d + eps;
			}
		}
	}



	void Partition_Vertex_Octree::set_query(std::vector<double> & query, const int P_or, const CGAL_Plane & P) const
	{
		double prec = 1.0 / (1 << 30) / (1 << 10);
		FT::set_relative_precision_of_to_double(prec);

		// Initializes a query that corresponds to the entire bounding box
		double _query[6] = { x_min - eps, x_max + eps, y_min - eps, y_max + eps, z_min - eps, z_max + eps };
		query = std::vector<double>(_query, _query + sizeof(_query) / sizeof(double));

		double D = -CGAL::to_double(P.d());

		// Refines the query by restricting query[.] to 2D dimensions
		// The idea is that, depending on the orientation of P, 
		// two coordinates of the bounding box defining the query are modified.

		// For example, P_or = 0 -> restriction along the x-axis -> query[0 1] are reset.

		query[(P_or << 1) + 0] = D - eps;
		query[(P_or << 1) + 1] = D + eps;
	}



	void Partition_Vertex_Octree::get_vertices_for_boundary_plane(const int H_id, const int H_or, const CGAL_Plane & H, std::list<Partition_Vertex *> & R) const
	{
		R.clear();

		// Constructs a query to get all the points of the slicing plane i
		std::vector<double> query;
		set_query(query, H_or, H);

		// Searches points
		search(query, R);

		// Loops on the list R, and returns the result that correspond to the plane i only
		std::list<Partition_Vertex *>::iterator it_v = R.begin();
		while (it_v != R.end()) {
			if ((*it_v)->belongs_to_plane(H_id)) {
				++it_v;
			} else {
				it_v = R.erase(it_v);
			}
		}
	}



	/*void Partition_Vertex_Octree::get_vertices_for_bounding_box_facet(const int i, std::list<Partition_Vertex *> & R) const
	{
		R.clear();

		// Constructs a query to get all the points of facet i
		std::vector<double> query;
		set_query(query, i);

		// Searches points
		search(query, R);

		// Loops on the list R, and returns the result that correspond to the plane i only
		std::list<Partition_Vertex *>::iterator it_v = R.begin();
		while (it_v != R.end()) {
			if ((*it_v)->belongs_to_plane(i)) {
				++it_v;
			} else {
				it_v = R.erase(it_v);
			}
		}
	}*/



	void Partition_Vertex_Octree::get_vertices_for_bounding_box_edge(const int i, const int j, 
		const int r_i, const int r_j, const CGAL_Plane & P_i, const CGAL_Plane & P_j,
		std::list<Partition_Vertex *> & R) const
	{
		R.clear();

		// Constructs a query to get all the points at the intersection of bounding box facets i and j
		std::vector<double> query;
		set_query(query, i, P_i, j, P_j);

		// Searches points
		search(query, R);

		// Loops on the list R, and returns the results that correspond to the planes i and j only
		std::list<Partition_Vertex *>::iterator it_v = R.begin();
		while (it_v != R.end()) {
			if ((*it_v)->belongs_to_planes(r_i, r_j)) {
				++it_v;
			} else {
				it_v = R.erase(it_v);
			}
		}
	}


	void Partition_Vertex_Octree::get_all_vertices(std::list<Partition_Vertex*> & V) const
	{
		if (is_node) {
			for (int i = 0; i < 8; i++) children[i]->get_all_vertices(V);
		} else {
			std::copy(vertices.begin(), vertices.end(), std::back_inserter(V));
		}
	}


	void Partition_Vertex_Octree::get_all_vertices_sorted_by_identifier(std::vector<Partition_Vertex*> & V) const
	{
		std::list<Partition_Vertex*> l_V;
		get_all_vertices(l_V);

		V.clear();
		std::copy(l_V.begin(), l_V.end(), std::back_inserter(V));
		std::sort(V.begin(), V.end(), sort_by_vertex_id);
	}
}