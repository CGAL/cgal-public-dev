#include "../include/support_plane.h"
#include "../include/polygon_vertex.h"
#include "../include/polygon_edge.h"
#include "../include/intersection_line.h"
#include "../include/segment.h"
#include "../include/polygon_group.h"
#include "../include/universe.h"
#include "../include/event_queue.h"
#include "../include/parameters.h"
#include "../include/vars.h"
#include "../include/stats.h"
#include <iostream>
#include <CGAL/convex_hull_2.h>


namespace Skippy 
{
	using CGAL::to_double;


	Support_Plane::Support_Plane(const int plane_id, const CGAL_Plane & plane_equation, const CGAL_Color & plane_color)
		: id (++Counters::id_planes),
		real_id (plane_id),
		plane (plane_equation),
		color (plane_color)
	{
		Universe::map_of_planes.push_back(this);

		running_vertices = 0;
		polygon_set = nullptr;

		x_min = FLT_MAX, x_max = -FLT_MAX;
		y_min = x_min, y_max = x_max;
	}



	Support_Plane::~Support_Plane()
	{
		// When the destructor of Support_Planes [6 .. N] is called, 
		// most data has been previously cleared after the construction of the final partition
		// and the call to function clear_polygons().
		
		concurrences.clear();
		std::vector<Intersection_Line*>::iterator it_l;
		for (it_l = lines.begin(); it_l != lines.end(); it_l++) {
			delete (*it_l);
		}

		std::map<int, Planar_Segment*>::iterator it_pl_s;
		for (it_pl_s = borders.begin(); it_pl_s != borders.end(); it_pl_s++) {
			delete it_pl_s->second;
		}

		if (Universe::params->use_landmarks && id < 6) {
			delete_landmarks();
		}
	}



	void Support_Plane::clear_polygons()
	{
		// We suppose that id >= 6.
		// This function is called as Partition_Facets have been built.

		// Destroys directions
		for (std::vector<Polygon_Directions*>::iterator it_d = polygon_directions.begin(); it_d != polygon_directions.end(); it_d++) {
			delete (*it_d);
		}

		// Destroys groups of polygons
		for (std::list<Polygon_Group*>::iterator it_g = groups.begin(); it_g != groups.end(); it_g++) {
			delete (*it_g);
		}
		groups.clear();

		// Destroys polygons, vertices and edges
		delete polygon_set;
		vertices_r.clear();

		// Destroys segments on intersection lines.
		for (size_t i = 0 ; i < lines.size() ; ++i) {
			if (!lines[i]->is_border) {
				lines[i]->clear_segments();
			}
		}

		// Destroys
		if (Universe::params->use_landmarks) {
			delete_landmarks();
		}
	}



	CGAL_Point_2 Support_Plane::project(const CGAL_Point_3 & P) const
	{
		// Assuming that P is a point on this supporting plane,
		// returns its 2D coordinates in the local 2D frame.

		return plane.to_2d(P);
	}


	CGAL_Point_3 Support_Plane::backproject(const CGAL_Point_2 & P) const
	{
		// Assuming that P is a point on this supporting plane,
		// returns its 3D coordinates in the global 3D frame.

		return plane.to_3d(P);
	}



	void Support_Plane::init_intersection_lines(const std::map<int, CGAL_Line_3> & lines_3d)
	{
		for (std::map<int, CGAL_Line_3>::const_iterator it_l = lines_3d.begin(); it_l != lines_3d.end(); it_l++) {

			// Gets a definition of the 3D line
			const CGAL_Point_3 A = it_l->second.point();
			const CGAL_Point_3 B = A + it_l->second.to_vector();

			// Gets a definition of the projected 2D line
			const CGAL_Point_2 A_proj = plane.to_2d(A);
			const CGAL_Point_2 B_proj = plane.to_2d(B);
			CGAL_Line_2 L(A_proj, B_proj);

			// Checks if this line is equal to another line
			// i.e. if their directions are collinear, and a point of one line belongs to the other
			std::vector<Intersection_Line*>::iterator it_m;
			for (it_m = lines.begin(); it_m != lines.end(); it_m++) {
				const CGAL_Line_2 & M = (*it_m)->line;
				if (CGAL::determinant(L.to_vector(), M.to_vector()) == 0 && (M.a() * A_proj.x() + M.b() * A_proj.y() + M.c() == 0)) {
					break;
				}
			}

			if (it_m != lines.end()) {
				// This line is therefore shared by another plane
				(*it_m)->mark_as_intersected(it_l->first);
				std::cerr << "Warning : P[" << id << "] : the intersection of three planes is a line" << std::endl;

			} else {
				// We are going to insert a new line, but before that we consider K_proj = A_proj + v_ab.
				// If K = backproject(K_proj) is above L and below the plane represented by L, or vice-versa,
				// then we revert the coefficients of L.

				const CGAL_Point_2 K_proj = A_proj + CGAL_Vector_2(L.a(), L.b());
				const CGAL_Point_3 K = plane.to_3d(K_proj);
				
				const CGAL_Plane & P = Universe::map_of_planes[it_l->first]->plane;

				const FT k_l = L.a() * K_proj.x() + L.b() * K_proj.y() + L.c();
				const FT k_p = P.a() * K.x() + P.b() * K.y() + P.c() * K.z() + P.d();
				if (k_l * k_p < 0) L = CGAL_Line_2(-L.a(), -L.b(), -L.c());

				Intersection_Line* I = new Intersection_Line(id, L, it_l->first);
				lines.push_back(I);
			}
		}

		if (Universe::params->use_landmarks) {
			init_landmarks();
		}
	}



	Intersection_Line* Support_Plane::get_line_by_identifier(const int id_object) const
	{
		int j = id_object - lines[0]->id_object;
		Intersection_Line* I = lines[j];

		assert(I->id_object == id_object);
		return I;
	}


	Intersection_Line* Support_Plane::get_line_for_plane(const int id_plane) const
	{
		for (std::vector<Intersection_Line*>::const_iterator it_l = lines.begin(); it_l != lines.end(); it_l++) {
			if ((*it_l)->intersects(id_plane)) {
				return (*it_l);
			}
		}
		return nullptr;
	}



	void Support_Plane::get_concurrent_lines(const Intersection_Line* I_1, const Intersection_Line* I_2, std::list<Intersection_Line*> & L)
	{
		int p = jmin(I_1->id_object, I_2->id_object), q = jmax(I_1->id_object, I_2->id_object);
		L.clear();

		// All triplets (t_0, t_1, t_2) are such that t_0 < t_1 < t_2, and we know that p < q.
		// We make use of this representation to select the interesting triplets and stop the search
		// once t_0 > p.

		for (std::set<Triplet, Triplet_Comparator>::const_iterator it_t = concurrences.begin(); it_t != concurrences.end(); it_t++) {
			const Triplet & T = (*it_t);

			const int t_0 = std::get<0>(T), t_1 = std::get<1>(T), t_2 = std::get<2>(T);

			if (t_0 < p) {
				if (t_1 == p && t_2 == q) {
					L.push_back(get_line_by_identifier(t_0));
				}
			} else if (t_0 == p) {
				if (t_1 < q) {
					if (t_2 == q) L.push_back(get_line_by_identifier(t_1));
				} else if (t_1 == q) {
					L.push_back(get_line_by_identifier(t_2));
				} else {
					// t_1 > q : Nothing else to search for
					return;
				}
			} else {
				return;
			}
		}
	}




	void Support_Plane::get_indices_of_intersecting_planes(const Intersection_Line* I_1, const Intersection_Line* I_2, std::list<int> & P)
	{
		// We consider the point located at the intersection of this plane, I_1 and I_2,
		// and we determine the indices of all the planes that intersect in that location.

		P.clear();
		P.push_back(real_id);

		// If I_1 and I_2 represent more than one plane, then indices of such planes are added to the result.
		for (int p_id : I_1->planes) P.push_back(Universe::map_of_planes[p_id]->real_id);
		for (int p_id : I_2->planes) P.push_back(Universe::map_of_planes[p_id]->real_id);

		// However, if there exists another line I_3 which intersects I_1 at the same location as I_2,
		// then planes represented by I_3 are also added to the list P.

		std::list<Intersection_Line*> I_L;
		get_concurrent_lines(I_1, I_2, I_L);

		for (std::list<Intersection_Line*>::iterator it_l = I_L.begin(); it_l != I_L.end(); it_l++) {
			Intersection_Line* L = (*it_l);
			for (int p_id : L->planes) P.push_back(Universe::map_of_planes[p_id]->real_id);
		}

		// Sorts indices
		std::vector<int> V_P;
		V_P.reserve(P.size());

		std::copy(P.begin(), P.end(), std::back_inserter(V_P));
		std::sort(V_P.begin(), V_P.end());

		// Moves back sorted indices to P
		P.clear();
		std::copy(V_P.begin(), V_P.end(), std::back_inserter(P));
	}


#if 0
	void Support_Plane::get_indices_of_intersecting_planes(const Intersection_Line* I_1, const Intersection_Line* I_2, std::list<int> & P)
	{
		// We consider the point located at the intersection of this plane, I_1 and I_2,
		// and we determine the indices of all the planes that intersect in that location.

		P.clear();
		P.push_back(id);

		// If I_1 and I_2 represent more than one plane, then indices of such planes are added to the result.
		std::copy(I_1->planes.begin(), I_1->planes.end(), std::back_inserter(P));
		std::copy(I_2->planes.begin(), I_2->planes.end(), std::back_inserter(P));

		// However, if there exists another line I_3 which intersects I_1 at the same location as I_2,
		// then planes represented by I_3 are also added to the list P.

		std::list<Intersection_Line*> I_L;
		get_concurrent_lines(I_1, I_2, I_L);

		for (std::list<Intersection_Line*>::iterator it_l = I_L.begin(); it_l != I_L.end(); it_l++) {
			Intersection_Line* L = (*it_l);
			std::copy(L->planes.begin(), L->planes.end(), std::back_inserter(P));
		}

		// Sorts indices
		std::vector<int> V_P;
		V_P.reserve(P.size());

		std::copy(P.begin(), P.end(), std::back_inserter(V_P));
		std::sort(V_P.begin(), V_P.end());

		// Moves back sorted indices to P
		P.clear();
		std::copy(V_P.begin(), V_P.end(), std::back_inserter(P));
	}
#endif


	void Support_Plane::insert_triplets_of_concurrent_lines(const std::vector<Intersection_Line*> & I)
	{
		size_t n = I.size();

		std::vector<int> J;
		J.reserve(n);

		// Gets a sorted vector of indices of concurrent lines

		for (int i = 0; i < n; i++) J.push_back(I[i]->id_object);
		std::sort(J.begin(), J.end());

		// From J, we extract a list of possible triplets
		// We test if the first triplet is already present in the list of triplets of concurrent lines
		// If so, then it's not necessary to insert any of the elements, since this intersection of n lines has already been processed before

		if (concurrences.find(std::make_tuple(J[0], J[1], J[2])) != concurrences.end()) return;

		for (int i = 0; i < n - 2; i++) {
			for (int j = i + 1; j < n - 1; j++) {
				for (int k = j + 1; k < n; k++) {
					concurrences.insert(std::make_tuple(J[i], J[j], J[k]));
				}
			}
		}
	}


	void Support_Plane::init_landmarks()
	{
		size_t n = lines.size();
		landmarks = new CGAL_Point_2**[n - 1];
		for (size_t i = 0 ; i < n - 1 ; ++i) {
			landmarks[i] = new CGAL_Point_2*[n - 1 - i];
			memset(landmarks[i], 0, (n - 1 - i) * sizeof(CGAL_Point_2*));
		}
	}



	void Support_Plane::delete_landmarks()
	{
		size_t n = lines.size();

		for (size_t i = 0 ; i < n - 1 ; ++i) {
			for (size_t j = 0 ; j < n - 1 - i ; ++j) {
				if (landmarks[i][j] != nullptr) delete landmarks[i][j];
			}
			delete[] landmarks[i];
		}
		delete[] landmarks;
	}


	const CGAL_Point_2 & Support_Plane::get_landmark(Intersection_Line* I, Intersection_Line* J)
	{
		int i = I->id_object - lines[0]->id_object;
		int j = J->id_object - lines[0]->id_object;

		int min_ij = jmin(i, j), max_ij = jmax(i, j);
		int p = min_ij, q = max_ij - min_ij - 1;

		if (landmarks[p][q] == nullptr) {
			CGAL_Point_2 M;
			const CGAL_Line_2 & I_L = I->line;
			const CGAL_Line_2 & J_L = J->line;
			CGAL::cpp11::result_of<K::Intersect_2(CGAL_Line_2, CGAL_Line_2)>::type object = CGAL::intersection(I_L, J_L);
			if (object) {
				if (const CGAL_Point_2* ptr = boost::get<CGAL_Point_2>(&*object)) {
					M = *ptr;
				}
			}
			landmarks[p][q] = new CGAL_Point_2(M);
		}

		return (*landmarks[p][q]);
	}



	CGAL_Point_2 Support_Plane::get_intersection_point(Intersection_Line* I, Intersection_Line* J) const
	{
		CGAL_Point_2 M;
		const CGAL_Line_2 & I_L = I->line;
		const CGAL_Line_2 & J_L = J->line;
		CGAL::cpp11::result_of<K::Intersect_2(CGAL_Line_2, CGAL_Line_2)>::type object = CGAL::intersection(I_L, J_L);
		if (object) {
			if (const CGAL_Point_2* ptr = boost::get<CGAL_Point_2>(&*object)) {
				M = *ptr;
			}
		}

		return M;
	}


	void Support_Plane::init_bounding_polygon(const std::list<CGAL_Point_3> & bounding_polygon_3d, const std::vector<std::list<int> > & bounding_facets)
	{
		// Converts the coordinates of the 3D polygon to the 2D local frame
		// and finds the coordinates of a squared bounding box for the support plane

		std::vector<CGAL_Point_2> points;

		FT ft_x_min = FT(FLT_MAX), ft_x_max = FT(-FLT_MAX);
		FT ft_y_min = ft_x_min, ft_y_max = ft_x_max;

		for (std::list<CGAL_Point_3>::const_iterator it_p = bounding_polygon_3d.begin(); it_p != bounding_polygon_3d.end(); it_p++) {
			CGAL_Point_2 pt = project(*it_p);

			FT x = pt.x(), y = pt.y();
			if (x < ft_x_min) ft_x_min = x;
			if (x > ft_x_max) ft_x_max = x;
			if (y < ft_y_min) ft_y_min = y;
			if (y > ft_y_max) ft_y_max = y;

			points.push_back(pt);
		}

		x_min = to_double(ft_x_min), x_max = to_double(ft_x_max);
		y_min = to_double(ft_y_min), y_max = to_double(ft_y_max);

		// Builds the polygon
		int n = int(points.size());

		for (int i = 0; i < n; i++) {

			// We are going to define the Planar_Segment [P_i, P_{i + 1}]
			// It is attached to the facets whose indices are listed at BF[i].
			// In case BF[i] is a list of size greater than 1, segments are duplicated.
			for (std::list<int>::const_iterator it_f = bounding_facets[i].begin(); it_f != bounding_facets[i].end(); it_f++) {
				int f = (*it_f);

				Intersection_Line* L = get_line_for_plane(f);
				if (L == nullptr) {
					L = get_line_for_plane(f);
				}
				Constraint CL = std::make_pair(L, (f % 2 == 0 ? PLUS : MINUS));
				Planar_Segment* s_f = new Planar_Segment(id, CL, points[i], points[(i + 1) % n]);
			}
		}

		// Discards the lines that don't intersect the bounding polygon
		for (std::vector<Intersection_Line*>::iterator it_l = lines.begin(); it_l != lines.end(); it_l++) {
			Intersection_Line* I = (*it_l);
			const CGAL_Line_2 & L = I->line;

			bool positive_vertices_exist = false;
			bool negative_vertices_exist = false;

			// We are going to test if the points that define the bounding polygon are located on both sides of L.
			// If not, it can be discarded, because no Polygon_Vertex will ever reach it

			for (int i = 0; i < n; i++) {
				const CGAL_Point_2 & M = points[i];
				FT epsilon = L.a() * M.x() + L.b() * M.y() + L.c();

				if (epsilon == 0) {
					// L intersects the extrema of the bouding polygon so we already know it shouldn't be discarded
					positive_vertices_exist = negative_vertices_exist = true;
					break;
				}

				if (epsilon > 0) {
					positive_vertices_exist = true;
				} else {
					negative_vertices_exist = true;
				}

				if (positive_vertices_exist && negative_vertices_exist) break;
			}

			I->set_inside(positive_vertices_exist && negative_vertices_exist);
		}

		// Sets a vector with lines that are inside the bounding polygon
		lines_inside.clear();
		lines_inside.reserve(lines.size());
		for (std::vector<Intersection_Line*>::const_iterator it_l = lines.begin(); it_l != lines.end(); ++it_l) {
			if ((*it_l)->is_inside) {
				lines_inside.push_back(*it_l);
			}
		}
		lines_inside.shrink_to_fit();
	}



	bool Support_Plane::assert_inside(Polygon_Vertex_R* v, const FT & t)
	{
		CGAL_Point_2 V_t = v->pt(t);

		for (std::map<int, Planar_Segment *>::iterator it_s = borders.begin() ; it_s != borders.end() ; ++it_s) {
			Intersection_Line* I = it_s->second->support;
			Sign eps_0 = I->sign(v->pt(v->t_init));
			Sign eps = I->sign(V_t);

			if (eps_0 != ZERO && eps != ZERO && eps_0 != eps) {
				return false;
			}
		}

		return true;
	}


	void Support_Plane::search_lines_in_neighborhood(Polygon_Vertex_R* v, const FT & t_1, const FT & t_2,
		std::list<std::tuple<Intersection_Line*, bool, FT> > & L)
	{
		// As Polygon_Vertex v propagates, we want to answer the following question :
		// which Intersection_Lines are susceptible to get intersected by it between t_1 and t_2 ?

		// Step 1
		// Finds the coordinates of the bounding box for v between t_1 and t_2

		CGAL_Point_2 M_1 = v->pt(t_1);
		CGAL_Point_2 M_2 = v->pt(t_2);

		double prec = 1.0 / (1 << 30) / (1 << 10);
		FT::set_relative_precision_of_to_double(prec);

		double x_1 = 0, x_2 = 0, y_1 = 0, y_2 = 0;
		if (Universe::params->fast_schedule) {
			x_1 = to_double(M_1.x()), x_2 = to_double(M_2.x()), y_1 = to_double(M_1.y()), y_2 = to_double(M_2.y());
		}

		for (int j = 0; j < lines_inside.size(); j++) {
			Intersection_Line* I_j = lines_inside[j];

			if (Universe::params->fast_schedule) {
				const double & h_a = I_j->hint_a, &h_b = I_j->hint_b, &h_c = I_j->hint_c;
				const double dh_1 = h_a * x_1 + h_b * y_1 + h_c;
				const double dh_2 = h_a * x_2 + h_b * y_2 + h_c;

				if (dh_1 * dh_2 > 0.1) {
					continue;
				} else if (dh_1 * dh_2 < -0.1) {
					L.push_back(std::make_tuple(I_j, false, 0));
				} else {
					std::pair<FT, bool> R = v->get_intersection_time(I_j);
					if (R.second) {
						if (R.first >= t_1 && R.first < t_2) {
							L.push_back(std::make_tuple(I_j, true, R.first));
						}
					}
				}

			} else {
				std::pair<FT, bool> R = v->get_intersection_time(I_j);
				if (R.second) {
					if (R.first >= t_1 && R.first < t_2) {
						L.push_back(std::make_tuple(I_j, true, R.first));
					}
				}
			}
		}
	}



	void Support_Plane::init_polygon_set()
	{
		polygon_set = new Polygon_Set(id, lines);
	}



	void Support_Plane::init_polygon(const std::vector<CGAL_Point_3> & polygon)
	{
		std::vector<CGAL_Point_2> P;
		project_polygon(polygon, P);

		init_polygon(P);
	}



	void Support_Plane::init_polygon(const std::vector<CGAL_Point_2> & polygon)
	{
		// Gets a polygon and splits it hierarchically
		Polygon_Tree* T = initialize_tree_of_polygons(polygon);
		if (T == nullptr) return;

		// Loops on all polygons of T and insert them into cells
		initialize_cells(T);
		delete T;
	}



	void Support_Plane::init_polygon(const CGAL_Point_2 & initial_barycenter,
		const std::vector<std::pair<CGAL_Point_2, CGAL_Vector_2> > & initial_directions)
	{
		Polygon_Tree* T = initialize_tree_of_polygons(initial_barycenter, initial_directions);
		if (T == nullptr) return;

		initialize_cells(T);
		delete T;
	}



	void Support_Plane::initialize_cells(Polygon_Tree* T)
	{
		std::list<Polygon*> LP;
		T->get_polygons(LP);

		// Loops on all the initial polygons
		// Each polygon is used to initialize a Polygon_Node structure
		for (std::list<Polygon*>::iterator it_p = LP.begin(); it_p != LP.end(); it_p++) {
			Polygon* P = (*it_p);

			// Gets a signature for the polygon
			const CGAL_Point_2 M = P->get_barycenter(0);
			std::vector<bool> S = Polygon_Node::make_signature(M, lines);

			// Inserts a couple (S, P) in the set of polygon cells
			polygon_set->insert(S, P);
		}

		// If we built polygons without running vertices,
		// then we only keep trace of the constraints that delimit them.

		for (auto it_n = polygon_set->cells_begin(); it_n != polygon_set->cells_end(); ++it_n) {
			Polygon_Node* N = it_n->second;
			for (std::list<Polygon*>::const_iterator it_p = N->polygons_begin(); it_p != N->polygons_end(); ++it_p) {
				if ((*it_p)->running_vertices == 0) {
					(*it_p)->forget();
				}
			}
		}

		T->remove_reference_to_polygons();
	}



	void Support_Plane::set_initial_propagation_directions(const std::vector<CGAL_Point_2> & R, CGAL_Point_2 & O, std::vector<std::pair<CGAL_Point_2, CGAL_Vector_2> > & D) const
	{
		// Computes the barycenter of all points

		const int n = int(R.size());
		FT x_bar = 0, y_bar = 0;

		for (int i = 0; i < n; i++) {
			const CGAL_Point_2 & pt = R[i];
			x_bar += pt.x(), y_bar += pt.y();
		}
		O = CGAL_Point_2(x_bar / n, y_bar / n);

		// The directions of all the points M[i] are given by the vector OM[i]
		// We finally get pairs of (M[i], OM[i]) that we store in the vector D

		D.reserve(n);
		for (int i = 0; i < n; i++) {
			CGAL_Vector_2 dR = R[i] - O;
			assert(dR.x() != 0 || dR.y() != 0);
			D.push_back(std::make_pair(R[i], R[i] - O));
		}
	}



	void Support_Plane::project_polygon(const std::vector<CGAL_Point_3> & P_0, std::vector<CGAL_Point_2> & P) const
	{
		std::vector<CGAL_Point_2> Q;
		Q.reserve(P_0.size());

		for (int i = 0; i < int(P_0.size()); i++) {
			Q.push_back(project(P_0[i]));
		}

		CGAL::convex_hull_2(Q.begin(), Q.end(), std::back_inserter(P));
	}



	std::tuple<int, int, int> Support_Plane::get_subvolume(const CGAL_Point_3 & M,
		const int gx, const int gy, const int gz,
		const std::vector<FT> & grad_x, const std::vector<FT> & grad_y, const std::vector<FT> & grad_z)
	{
		// We have a grid of subvolumes, whose graduations along all axes
		// are represented by vectors grad_x, grad_y, grad_z. 

		int i = 0, j = 0, k = 0;
		FT x = M.x(), y = M.y(), z = M.z();

		// We simply loop on each vector to get i s.t. grad_x[i] <= M.x() < grad_x[i + 1].
		// Same for variables j and k.
		
		while (i < gx - 1) {
			if (grad_x[i + 1] < x) ++i;
			else break;
		}

		while (j < gy - 1) {
			if (grad_y[j + 1] < y) ++j;
			else break;
		}

		while (k < gz - 1) {
			if (grad_z[k + 1] < z) ++k;
			else break;
		}

		if (i >= gx || j >= gy || k >= gz) {
			throw std::logic_error("Error : inconsistent value returned from get_subvolume().");
		}

		return std::make_tuple(i, j, k);
	}



	void Support_Plane::project_and_decompose_with_respect_to_slices(const std::vector<CGAL_Point_3> & polygon,
		const int gx, const int gy, const int gz,
		CGAL_Point_2 & initial_barycenter,
		std::vector<std::pair<CGAL_Point_2, CGAL_Vector_2> > & initial_directions,
		std::set<std::tuple<int, int, int> > & subvolumes)
	{
		// This function is called by Kinetic_Propagation_Multiple,
		// as we exhibit the subpolygons that propagate in each subvolume.

		int slices = gx + gy + gz - 3;
		int ns = int(Universe::map_of_planes.size()) - slices;

		// Step 1.
		// First of all the polygon gets projected onto the plane.

		std::vector<CGAL_Point_2> P;
		project_polygon(polygon, P);

		// Then, a tree of polygons is obtained.
		// We only take into account the lines associated to the slicing planes.

		set_initial_propagation_directions(P, initial_barycenter, initial_directions);

		// Step 2.
		// Tries to avoid the precise initialization and the splitting of a polygon tree
		// and the slow computations it requires, by computing the roughest bounding box
		// and determining if it already falls into a subvolume (i, j, k).

		FT x_min = FT(FLT_MAX), x_max = -x_min;
		FT y_min = x_min, y_max = -y_min;
		FT z_min = x_min, z_max = -z_min;

		for (size_t i = 0 ; i < polygon.size() ; ++i) {
			const CGAL_Point_3 & pt = polygon[i];
			FT x = pt.x(), y = pt.y(), z = pt.z();
			if (x < x_min) x_min = x;
			if (x > x_max) x_max = x;
			if (y < y_min) y_min = y;
			if (y > y_max) y_max = y;
			if (z < z_min) z_min = z;
			if (z > z_max) z_max = z;
		}

		// Gets graduations along each axis (slicing planes + boundaries)

		std::vector<FT> grad_x(1 + gx), grad_y(1 + gy), grad_z(1 + gz);
		
		grad_x[0] = -Universe::map_of_planes[0]->plane.d() / Universe::map_of_planes[0]->plane.a();
		grad_y[0] = -Universe::map_of_planes[2]->plane.d() / Universe::map_of_planes[2]->plane.b();
		grad_z[0] = -Universe::map_of_planes[4]->plane.d() / Universe::map_of_planes[4]->plane.c();

		grad_x[gx] = -Universe::map_of_planes[1]->plane.d() / Universe::map_of_planes[1]->plane.a();
		grad_y[gy] = -Universe::map_of_planes[3]->plane.d() / Universe::map_of_planes[3]->plane.b();
		grad_z[gz] = -Universe::map_of_planes[5]->plane.d() / Universe::map_of_planes[5]->plane.c();
		
		for (size_t i = 0 ; i < gx - 1 ; ++i) {
			const CGAL_Plane & H = Universe::map_of_planes[ns + i]->plane;
			grad_x[i + 1] = - H.d() / H.a();
		}
		
		for (size_t j = 0 ; j < gy - 1 ; ++j) {
			const CGAL_Plane & H = Universe::map_of_planes[ns + gx - 1 + j]->plane;
			grad_y[j + 1] = - H.d() / H.b();
		}
		
		for (size_t k = 0 ; k < gz - 1 ; ++k) {
			const CGAL_Plane & H = Universe::map_of_planes[ns + gx - 1 + gy - 1 + k]->plane;
			grad_z[k + 1] = - H.d() / H.c();
		}

		// Gets the id of the subvolumes to which the lowest and highest point of the polygon belong.
		
		std::tuple<int, int, int> tuple_pt_min = get_subvolume(CGAL_Point_3(x_min, y_min, z_min), gx, gy, gz, grad_x, grad_y, grad_z);
		std::tuple<int, int, int> tuple_pt_max = get_subvolume(CGAL_Point_3(x_max, y_max, z_max), gx, gy, gz, grad_x, grad_y, grad_z);
		bool full_decomposition = (tuple_pt_min != tuple_pt_max);

		if (tuple_pt_min == tuple_pt_max) {
			// The primitive belongs to only one subvolume.
			// We insert an identifier tuple into 'subvolumes' and return
			subvolumes.insert(tuple_pt_min);

		} else {
			// Step 3.
			// Performs a full decomposition of each primitive relatively
			// to all the lines that correspond to the slicing planes. 

			Polygon_Directions* D = new Polygon_Directions(initial_barycenter, initial_directions);
			Polygon_Tree* polygon_tree = new Polygon_Tree(id, 0, D);

			for (std::vector<Intersection_Line *>::const_iterator it_l = lines.begin(); it_l != lines.end(); it_l++) {
				Intersection_Line* I = (*it_l);
				bool is_slice = false;

				for (int p : I->planes) {
					if (p >= ns) {
						is_slice = true;
						break;
					}
				}

				std::list<Polygon_Vertex *> intersection_pts;
				polygon_tree->split(I, intersection_pts, 0, 0, false);
			}

			// Gets all subpolygons and assigns them to a subvolume.

			std::list<Polygon*> LP;
			polygon_tree->get_polygons(LP);

			for (std::list<Polygon*>::const_iterator it_p = LP.begin(); it_p != LP.end(); ++it_p) {

				Polygon* subpolygon = (*it_p);
				CGAL_Point_2 m = subpolygon->get_barycenter(0);
				CGAL_Point_3 M = backproject(m);
				std::tuple<int, int, int> T = get_subvolume(M, gx, gy, gz, grad_x, grad_y, grad_z);
				
				/*while (i < gx - 1) {
					const CGAL_Plane & H = Universe::map_of_planes[ns + i]->plane;
					if (H.a() * M.x() + H.b() * M.y() + H.c() * M.z() + H.d() > 0) ++i;
					else break;
				}*/

				subvolumes.insert(T);
			}

			delete polygon_tree;
			delete D;
		}
	}



	Polygon_Tree* Support_Plane::initialize_tree_of_polygons(const std::vector<CGAL_Point_2> & polygon)
	{
		// Projects and regularizes polygon
		const std::vector<CGAL_Point_2> & P = polygon;
		std::vector<Intersection_Line*> reference_lines;

		// Gets the locations and directions of all initial vertices
		CGAL_Point_2 initial_barycenter;
		std::vector< std::pair<CGAL_Point_2, CGAL_Vector_2> > initial_directions;
		set_initial_propagation_directions(P, initial_barycenter, initial_directions);

		Polygon_Tree* T = initialize_tree_of_polygons(initial_barycenter, initial_directions);
		return T;
	}



	Polygon_Tree* Support_Plane::initialize_tree_of_polygons(const CGAL_Point_2 & initial_barycenter,
		const std::vector<std::pair<CGAL_Point_2, CGAL_Vector_2> > & initial_directions)
	{
		// Part 1.
		// Builds an initial polygon, itself used to initialize a structure of polygonal tree
		// This tree will be later hierarchically decomposed, depending on possible intersections with Intersection_Lines

		int seed = int(polygon_directions.size());
		Polygon_Directions* D = new Polygon_Directions(initial_barycenter, initial_directions);
		polygon_directions.push_back(D);

		Polygon_Tree* polygon_tree = new Polygon_Tree(id, seed, D);


		// Part 2.
		// Divides the polygon into subpolygons it is intersects with any of the intersection lines.
		// In case of intersection, we keep trace of the constrained vertices that delimit the intersection.
		// Later, we will initialize segments using them.

		typedef std::tuple<Polygon_Vertex*, Polygon_Vertex*, Polygon_Vertex*, Polygon_Vertex*> Quadruplet;

		std::map<int, Quadruplet> polyline_intersections;
		std::map<int, std::map<int, CGAL_Point_2> > mutual_intersections;

		// Part 2.1
		// If we make a polygon propagate in multiple subvolumes, then a part of it may be outside the current bounding box.
		// If so we only keep the inner portion and initialize a 2-uplet of segments.

		if (Universe::params->use_grid) {
			for (std::vector<Intersection_Line*>::const_iterator it_l = lines.begin() ; it_l != lines.end() ; ++it_l) {
				Intersection_Line* I = (*it_l);
				int i = I->id_object;
				if (!I->is_border) break;

				// Determines what the sign of the constraint must be
				// on the side of the subpolygon that we keep (if any)

				Sign I_eps = ZERO;
				for (int p : I->planes) {
					if (p <= 5) { 
						I_eps = (p % 2 == 0 ? PLUS : MINUS); 
						break; 
					}
				}

				// If the polygon is intersected by the line I,
				// then we split the tree and iterate by selecting one subtree
				
				// For now, we don't try to collect the intersection points.
				// Indeed, the polygon may be subdivided again and therefore its intersection with I may change.

				if (polygon_tree->split(I, 0, Universe::params->K, false)) {
					Polygon_Tree* selected_subtree = (I_eps == PLUS ? polygon_tree->subtree_plus : polygon_tree->subtree_minus);
					Polygon* subpolygon = selected_subtree->polygon;

					selected_subtree->polygon = nullptr;
					delete polygon_tree->subtree_plus;
					delete polygon_tree->subtree_minus;
					polygon_tree->polygon = subpolygon;
					polygon_tree->is_node = false;
				}
			}

			// Builds 2-uplets.

			Polygon* polygon = polygon_tree->polygon;
			for (std::vector<Intersection_Line*>::const_iterator it_l = lines.begin() ; it_l != lines.end() ; ++it_l) {
				Intersection_Line* I = (*it_l);
				int i = I->id_object;
				if (!I->is_border) break;

				Polygon_Vertex *v1 = nullptr, *v2 = nullptr;
				for (auto it_e = polygon->edges.begin() ; it_e != polygon->edges.end() ; ++it_e) {
					if ((*it_e)->is_constrained_by(I)) {
						v1 = (*it_e)->v1, v2 = (*it_e)->v2;
						break;
					}
				}

				if (v1 == nullptr) continue;

				if (v1->sign_of_constraint(I) == PLUS) {
					polyline_intersections[I->id_object] = std::make_tuple(v1, nullptr, v2, nullptr);
				} else {
					polyline_intersections[I->id_object] = std::make_tuple(nullptr, v1, nullptr, v2);
				}

				update_mutual_intersections(v1, i, mutual_intersections);
				update_mutual_intersections(v2, i, mutual_intersections);
			}
		}

		// Part 2.2
		// Processes inner lines.

		for (std::vector<Intersection_Line *>::const_iterator it_l = lines.begin(); it_l != lines.end(); it_l++) {
			Intersection_Line* I = (*it_l);
			if (I->is_border) continue;

			int i = I->id_object;
			std::list<Polygon_Vertex *> intersection_pts;
			std::vector<Polygon_Vertex*> ranked_pts;

			if (polygon_tree->split(I, intersection_pts, 0, Universe::params->K, true)) {
				Polygon_Tree::sort_intersection_points(I, intersection_pts, 0, ranked_pts);

				// If the tree had not been divided by another line during the initialization process,
				// or if only one leaf of the tree is split by the line, then 'intersection' only contains at least 4 vertices.
				// (2 intersections points, that are duplicated since there are 2 sides to take into account)

				// We return the 4 furthest points along the line I.

				size_t n = ranked_pts.size();
				Polygon_Vertex *v1_p = nullptr, *v1_n = nullptr, *v2_p = nullptr, *v2_n = nullptr;

				if (ranked_pts[0]->sign_of_constraint(I) == PLUS) {
					v1_p = ranked_pts[0]; v1_n = ranked_pts[1];
				} else {
					v1_p = ranked_pts[1]; v1_n = ranked_pts[0];
				}
				if (ranked_pts[n - 1]->sign_of_constraint(I) == PLUS) {
					v2_p = ranked_pts[n - 1]; v2_n = ranked_pts[n - 2];
				} else {
					v2_p = ranked_pts[n - 2]; v2_n = ranked_pts[n - 1];
				}

				polyline_intersections[i] = std::make_tuple(v1_p, v1_n, v2_p, v2_n);
				update_mutual_intersections(ranked_pts, i, mutual_intersections);
			}
		}


		// Part 3.
		// Initializes segments

		for (std::map<int, Quadruplet>::const_iterator it_pli = polyline_intersections.cbegin(); it_pli != polyline_intersections.cend(); it_pli++) {
			int i = it_pli->first;
			Intersection_Line* I = get_line_by_identifier(i);

			const Quadruplet & Q = it_pli->second;
			Polygon_Vertex *v1_p = std::get<0>(Q), *v1_n = std::get<1>(Q), *v2_p = std::get<2>(Q), *v2_n = std::get<3>(Q);

			const std::map<int, CGAL_Point_2> & L = mutual_intersections[i];

			for (int s = 0 ; s < 2 ; ++s) {
				Polygon_Vertex *v1, *v2;
				Constraint C;

				if (s == 0) {
					v1 = v1_p, v2 = v2_p;
					C = Constraint(I, PLUS);
				} else {
					v1 = v1_n, v2 = v2_n;
					C = Constraint(I, MINUS);
				}

				if (v1 == nullptr || v2 == nullptr) continue;

				// There are now 4 possible cases.
				// Depending on the type of vertices v1 and v2 actually are, we initialize a segment.

				Polygon_Vertex_R* v1_r = v1->to_r(), *v2_r = v2->to_r();
				Polygon_Vertex_S* v1_s = v1->to_s(), *v2_s = v2->to_s();

				const CGAL_Point_2 V1 = v1->get_M(), &V2 = v2->get_M();
				const CGAL_Point_3 W1 = backproject(V1), W2 = backproject(V2);

				if (v1_r == nullptr && v2_r == nullptr) {
					// Initiliazes a constant segment
					Constraint C_1 = v1_s->get_other_constraint(C);
					Constraint C_2 = v2_s->get_other_constraint(C);
					std::list<Intersection_Line*> C_L;
					for (auto it_l = L.begin() ; it_l != L.end() ; ++it_l) {
						int l = it_l->first;
						if (l != C_1.first->id_object && l != C_2.first->id_object) {
							C_L.push_back(get_line_by_identifier(l));
						}
					}
					init_1_uplets_of_constant_segments(I, v1_s, v2_s, C_1, C_2, C_L, 0);

				} else if (v1_r == nullptr && v2_r != nullptr) {
					// Initializes a unidirectional segment
					Constraint C_1 = v1_s->get_other_constraint(C);
					init_1_uplets_of_unidirectional_segments(I, V1, V2, W2, v2_r, C_1, 0);
					for (auto it_l = L.begin() ; it_l != L.end() ; ++it_l) {
						int l = it_l->first;
						if (l != C_1.first->id_object) {
							v2_r->indicate_line_initially_crossed_by_segments(get_line_by_identifier(l));
						}
					}

				} else if (v1_r != nullptr && v2_r == nullptr) {
					// Same
					Constraint C_2 = v2_s->get_other_constraint(C);
					init_1_uplets_of_unidirectional_segments(I, V2, V1, W1, v1_r, C_2, 0);
					for (auto it_l = L.begin() ; it_l != L.end() ; ++it_l) {
						int l = it_l->first;
						if (l != C_2.first->id_object) {
							v1_r->indicate_line_initially_crossed_by_segments(get_line_by_identifier(l));
						}
					}

				} else {
					if (L.empty()) {
						// Initializes a bidirectional segment
						init_2_uplets_of_bidirectional_segments(I, V1, W1, V2, W2, v1_r, v2_r, 0);

					} else {
						// Or a 2-uplet of unidirectional segments
						std::map<int, CGAL_Point_2>::const_iterator it_l = L.cbegin();
						Intersection_Line* J = get_line_by_identifier(it_l->first);
						const CGAL_Point_2 & M_ij = it_l->second;

						Sign eps_1 = J->sign(v1_r->get_M()), eps_2 = (eps_1 == PLUS ? MINUS : PLUS);
						assert(eps_1 != ZERO && eps_2 != ZERO);
						Constraint C_1(J, eps_1), C_2(J, eps_2);

						init_1_uplets_of_unidirectional_segments(I, M_ij, V1, W1, v1_r, C_1, 0);
						init_1_uplets_of_unidirectional_segments(I, M_ij, V2, W2, v2_r, C_2, 0);
						++it_l;

						while (it_l != L.end()) {
							Intersection_Line* K = get_line_by_identifier(it_l->first);
							const CGAL_Point_2 & M_ik = it_l->second;
							int r = locate_intersection(M_ij, V1, V2, M_ik);
							if (r == 0) {
								// Intersects both segments [M_ij v1] and [M_ij v2]
								std::cerr << "Warning : not implemented case in the initialization process -- no guarantee of result" << std::endl;
							} else if (r == 1) {
								// Intersects segment [M_ij v1]
								v1_r->indicate_line_initially_crossed_by_segments(K);
							} else {
								// Intersects segment [M_ij v2]
								v2_r->indicate_line_initially_crossed_by_segments(K);
							}
							++it_l;
						}
					}					
				}
			}
		}

		return polygon_tree;
	}



	void Support_Plane::update_mutual_intersections(const std::vector<Polygon_Vertex*> & ranked_points,
		const int current_line,
		std::map<int, std::map<int, CGAL_Point_2> > & mutual_intersections) const
	{
		size_t n = ranked_points.size();

		for (size_t q = 0; q < n; ++q) {
			update_mutual_intersections(ranked_points[q], current_line, mutual_intersections);
		}
	}



	void Support_Plane::update_mutual_intersections(Polygon_Vertex* v, const int current_line, 
		std::map<int, std::map<int, CGAL_Point_2> > & mutual_intersections) const
	{
		if (Polygon_Vertex_S* v_s = v->to_s()) {
			// Vertex v_s has two constraints : I, of index current_line,
			// and another line H processed before.
			// If we haven't met this intersection of lines yet, then we mark down
			// its location in the map mutual_intersections.
			int h = v_s->get_constraint().first->id_object;
			if (h == current_line) h = v_s->get_second_constraint().first->id_object;

			if (mutual_intersections[current_line].find(h) == mutual_intersections[current_line].end()) {
				CGAL_Point_2 M = v_s->get_M();
				mutual_intersections[current_line][h] = M;
				mutual_intersections[h][current_line] = M;
			}
		}
	}


#if 0
	Polygon_Tree* Support_Plane::initialize_tree_of_polygons(const std::vector<CGAL_Point_2> & polygon)
	{
		// Part 1.
		// Initializes a 2D polygon

		// Projects and regularizes polygon
		const std::vector<CGAL_Point_2> & P = polygon;
		std::vector<Intersection_Line*> reference_lines;

		// Gets the locations and directions of all initial vertices
		CGAL_Point_2 initial_barycenter;
		std::vector< std::pair<CGAL_Point_2, CGAL_Vector_2> > initial_directions;
		set_initial_propagation_directions(P, initial_barycenter, initial_directions);

		// Builds an initial polygon, itself used to initialize a structure of polygonal tree
		// This tree will be later hierarchically decomposed, depending on possible intersections with Intersection_Lines

		int seed = int(polygon_directions.size());
		Polygon_Directions* D = new Polygon_Directions(initial_barycenter, initial_directions);
		polygon_directions.push_back(D);

		Polygon_Tree* polygon_tree = new Polygon_Tree(id, seed, D);


		// Part 2.
		// Divides the polygon into subpolygons it is intersects with any of the intersection lines.
		// In case of intersection, we keep trace of the constrained vertices that delimit the intersection.
		// Later, we will initialize segments using them.

		typedef std::tuple<Polygon_Vertex_R*, Polygon_Vertex_R*, Polygon_Vertex_R*, Polygon_Vertex_R*> Quadruplet;

		std::map<int, Quadruplet> polyline_intersections;
		std::map<int, std::map<int, CGAL_Point_2> > mutual_intersections;

		for (std::vector<Intersection_Line *>::const_iterator it_l = lines.begin(); it_l != lines.end(); it_l++) {
			Intersection_Line* I = (*it_l);

			int i = I->id_object;
			std::list<Polygon_Vertex *> intersection_pts;

			if (polygon_tree->split(I, intersection_pts, 0, Universe::params->K, true)) {

				// If the tree had not been divided by another line during the initialization process,
				// or if only one leaf of the tree is split by the line, then 'intersection' only contains 4 active vertices.
				// (2 intersections points, that are duplicated since there are 2 sides to take into account.)

				// Otherwise it contains several 4-uplets of vertices that correspond to the intersections of the line and each of the intersected subpolygons of the tree. 
				// In the latter case, only 4 vertices are active (others correspond to the intersections of two lines, that's why they are disabled).

				Polygon_Vertex_R *v1_p = nullptr, *v1_n = nullptr, *v2_p = nullptr, *v2_n = nullptr;

				for (std::list<Polygon_Vertex *>::iterator it_v = intersection_pts.begin(); it_v != intersection_pts.end(); it_v++) {
					
					if (Polygon_Vertex_S* v_s = (*it_v)->to_s()) {
						// If the vertex is not active, then it is at the intersection of I and another line H processed before
						// If we haven't met it yet, then we mark down it location in map crossed_lines_I
						// and we notify the line H that an intersection with I exists
						int h = v_s->get_constraint().first->id_object;
						if (h == i) h = v_s->get_second_constraint().first->id_object;
						if (mutual_intersections[i].find(h) == mutual_intersections[i].end()) {
							mutual_intersections[i][h] = v_s->get_M();
							mutual_intersections[h][i] = v_s->get_M();
						}
						continue;

					} else if (Polygon_Vertex_R* v_r = (*it_v)->to_r()) {
						// We find two pairs of opposite vertices
						if (v_r->get_constraint().second == PLUS) {
							set_element_in_quadruplet(v1_n, v2_n, v1_p, v2_p, v_r);
						} else {
							set_element_in_quadruplet(v1_p, v2_p, v1_n, v2_n, v_r);
						}
					}
				}

				polyline_intersections[i] = std::make_tuple(v1_p, v1_n, v2_p, v2_n);
			}
		}

		// Part 3.
		// Initializes segments

		for (std::map<int, Quadruplet>::const_iterator it_pli = polyline_intersections.cbegin(); it_pli != polyline_intersections.cend(); it_pli++) {
			int i = it_pli->first;

			// Exhibits the 4-uplet of vertices

			const Quadruplet & Q = it_pli->second;
			Polygon_Vertex_R *v1_p = std::get<0>(Q), *v1_n = std::get<1>(Q), *v2_p = std::get<2>(Q), *v2_n = std::get<3>(Q);
			Intersection_Line* I = v1_p->get_constraint().first;

			bool v1_paired = Polygon_Vertex_R::are_paired_vertices(v1_p, v1_n);
			bool v2_paired = Polygon_Vertex_R::are_paired_vertices(v2_p, v2_n);

			const CGAL_Point_2 & V1 = v1_p->get_M(), &V2 = v2_p->get_M();
			const CGAL_Point_3 W1 = backproject(V1), W2 = backproject(V2);

			const std::map<int, CGAL_Point_2> & L = mutual_intersections[i];

			if (L.empty()) {
				// If line I doesn't intersect any other line,
				// then we initialize a 4-uplet of opposite bidirectional segments without dilimitation
				// w.r.t. the other lines of the support plane.
				if (v1_paired && v2_paired) {
					init_4_uplets_of_bidirectional_segments(I, V1, W1, V2, W2, v1_p, v1_n, v2_p, v2_n, 0);
				} else {
					init_2_uplets_of_bidirectional_segments(I, V1, W1, V2, W2, v1_p, v2_p, 0);
					init_2_uplets_of_bidirectional_segments(I, V1, W1, V2, W2, v1_n, v2_n, 0);
				}

			} else {
				// I intersects (I_1 .. I_n).
				// We use I_1 to initialize unidirectional segments with an initial constraint,
				// and we assign the other lines (I_2 .. I_n) to one of the 2 2-uplets of segments thus created,
				// by marking them as intersected by the segments.

				// Sets J, which plays the role of I_1.
				// It intersects the line I in M_ij.
				std::map<int, CGAL_Point_2>::const_iterator it_l = L.cbegin();
				Intersection_Line* J = get_line_by_identifier(it_l->first);
				const CGAL_Point_2 & M_ij = it_l->second;
				//const CGAL_Point_3 W_ij = backproject(M_ij);

				Sign eps_1 = J->sign(v1_p->get_M());
				Sign eps_2 = (eps_1 == PLUS ? MINUS : PLUS);
				assert(eps_1 != ZERO && eps_2 != ZERO);
				Constraint C_1(J, eps_1), C_2(J, eps_2);

				// We initialize two 2-uplets of segments [M_ij v1], [M_ij v2]
				// J is used to set the initial constraint of these segments.

				if (v1_paired) {
					init_2_uplets_of_unidirectional_segments(I, M_ij, V1, W1, v1_p, v1_n, C_1, 0);
				} else {
					init_1_uplets_of_unidirectional_segments(I, M_ij, V1, W1, v1_p, C_1, 0);
					init_1_uplets_of_unidirectional_segments(I, M_ij, V1, W1, v1_n, C_1, 0);
				}

				if (v2_paired) {
					init_2_uplets_of_unidirectional_segments(I, M_ij, V2, W2, v2_p, v2_n, C_2, 0);
				} else {
					init_1_uplets_of_unidirectional_segments(I, M_ij, V2, W2, v2_p, C_2, 0);
					init_1_uplets_of_unidirectional_segments(I, M_ij, V2, W2, v2_n, C_2, 0);
				}

				while (++it_l != L.cend()) {
					// Sets sucessive lines K, which play the roles of I_2 .. I_n.
					// They also intersect I in a set of points M_ik, 
					// we need to find if its intersects [M_ij v1], [M_ij v2], or both.
					Intersection_Line* K = get_line_by_identifier(it_l->first);
					const CGAL_Point_2 & M_ik = it_l->second;

					int r = locate_intersection(M_ij, V1, V2, M_ik);
					if (r == 0) {
						// Intersects both segments [M_ij v1] and [M_ij v2]
						// Nothing to do ?
						std::cerr << "Warning : not implemented case in the initialization process -- no guarantee of result" << std::endl;
					} else if (r == 1) {
						// Intersects segment [M_ij v1]
						v1_p->indicate_line_initially_crossed_by_segments(K);
					} else {
						// Intersects segment [M_ij v2]
						v2_p->indicate_line_initially_crossed_by_segments(K);

					}
				}

				v1_p->copy_crossed_lines(v1_n);
				v2_p->copy_crossed_lines(v2_n);
			}
		}

		return polygon_tree;
	}
#endif


	int Support_Plane::locate_intersection(const CGAL_Point_2 & M_0, const CGAL_Point_2 & M_1, const CGAL_Point_2 & M_2, const CGAL_Point_2 & P) const
	{
		// By calling this function, we have two segments s_1 = [M_1 M_0], s_2 = [M_0 M_2]
		// and we want to determine if the point P belongs to s_1, s_2 or both.
		// We suppose that P belongs to [M_1 M_2].

		const FT x_0 = M_0.x(), x_1 = M_1.x(), x_2 = M_2.x();
		if (x_1 != x_2) {
			// [M_1 M_2] is not vertical : the x coordinate can be used.
			// Three possibilities :
			// - P.x() = x_0 => P = M_0 => return 0
			// - P.x() is in ]x_0 x_1]  => return 1
			// - P.x() is in ]x_0 x_2]  => return 2
			const FT x = P.x();
			if (x_0 == x) {
				return 0;
			} else if ((x_0 <= x && x <= x_1) || (x_1 <= x && x <= x_0)) {
				return 1;
			} else {
				return 2;
			}
		} else {
			// [M_1 M_2] is vertical.
			// We should use the y coordinate, but with a similar reasoning.
			const FT y_0 = M_0.y(), y_1 = M_1.y(), y = P.y();
			if (y_0 == y) {
				return 0;
			} else if ((y_0 <= y && y <= y_1) || (y_1 <= y && y <= y_0)) {
				return 1;
			} else {
				return 2;
			}
		}
	}



	void Support_Plane::set_element_in_quadruplet(Polygon_Vertex_R* & v1_os, Polygon_Vertex_R* & v2_os, Polygon_Vertex_R* & v1_ts, Polygon_Vertex_R* & v2_ts, Polygon_Vertex_R* v) const
	{
		if (v1_os != nullptr && v2_os == nullptr) {
			if ((v1_os->get_M() - v->get_M()) == CGAL::NULL_VECTOR) {
				v1_ts = v;
				return;
			}
		} else if (v1_os == nullptr && v2_os != nullptr) {
			if ((v2_os->get_M() - v->get_M()) == CGAL::NULL_VECTOR) {
				v2_ts = v;
				return;
			}
		} else if (v1_os != nullptr && v2_os != nullptr) {
			if ((v1_os->get_M() - v->get_M()).squared_length() < (v2_os->get_M() - v->get_M()).squared_length()) {
				v1_ts = v;
				return;
			} else {
				v2_ts = v;
				return;
			}
		}
		if (v1_ts == nullptr) {
			v1_ts = v;
			return;
		} else if (v2_ts == nullptr && (v1_ts->dM * v->dM < 0)) {
			v2_ts = v;
			return;
		}
	}



	void Support_Plane::init_schedule()
	{
		// Loops on the list of all active vertices and schedule events

		for (std::map<int, Polygon_Vertex_R*>::iterator it_v = vertices_r.begin(); it_v != vertices_r.end(); it_v++) {
			Polygon_Vertex_R* v = it_v->second;

			// Two possible cases :
			// - either the vertex is independant, and we compute its events normally,
			// - or it is paired to another vertex, and in this case, we duplicate the events of its 'master' vertex

			if (v->is_independent()) {
				if (v->unconstrained()) {
					v->schedule_events();
				} else {
					const Constraint & C = v->get_constraint();
					Intersection_Line* I = C.first;
					v->schedule_events(I, nullptr);
				}
			}

			else {
				Polygon_Vertex_R* v_ts = v->get_paired_vertex();
				v->schedule_events(v_ts);
			}
		}
	}



	void Support_Plane::print_events(const std::list<Event_Vertex_Line*> & E, const std::string & type) const
	{
		for (std::list<Event_Vertex_Line*>::const_iterator it_e = E.begin() ; it_e != E.end() ; ++it_e) {
			std::string prefix = (it_e == E.begin() ? "* " : "  ");
			Event_Vertex_Line* e = (*it_e);
			std::cout << prefix << " [" << e->plane << "] " << e->intersectant << " " << e->intersected << " " << e->t_intersectant << " " << type << std::endl;
		}
	}



	void Support_Plane::process_event(Event_Vertex_Line* e)
	{
		int intersectant = e->intersectant;
		int intersected = e->intersected;
		const FT t_intersectant = e->t_intersectant;

		Polygon_Vertex_R* v = vertices_r[intersectant];

		// Step 1.
		// Gets simulateneous events that are related to the same subpolygon as v's,
		// and counts the number of vertices intersecting lines at time t_intersectant.

		std::list<Event_Vertex_Line*> E;
		std::list<Polygon_Vertex_R*> V;
	
		get_references_to_simulatenous_events_for_this_polygon(v, t_intersectant, E);
		get_number_of_vertices_simultaneously_intersecting_lines(v, e, E, V);
		pop_fitered_simultaneous_events_from_queue(E, V);
		E.push_front(e);

		// Step 2.
		// Associates the appropriate handler to the list of simultaneous events E.

		try {
			switch (V.size()) {
			case 1: process_event_1(E, V); break;
			case 2: process_event_2(E, V); break;
			case 3: process_event_3(E, V); break;
			default: throw std::logic_error("Error : Couldn't handle a set of simultaneous events");
			}

		} catch (std::exception & except) {
			// Prints the support plane
			std::cout << "Error raised while processing [" << id << "] " << intersectant << " " << intersected << " " << t_intersectant << std::endl;
			draw(to_double(t_intersectant), 0.5, 20000, 0.5, 0.5, 10);
			throw except;
		}

		if (Universe::params->rt_check) {
			real_time_check(t_intersectant);
		}
	}



	void Support_Plane::process_event_1(const std::list<Event_Vertex_Line*> & E, const std::list<Polygon_Vertex_R*> & V)
	{
		// We associate a list of simultaneous events that involve only one vertex to an appropriate handler.
		
		Event_Vertex_Line* e_vl = E.front();

		int intersectant = e_vl->intersectant;
		int intersected = e_vl->intersected;
		FT t_intersectant = e_vl->t_intersectant;

		Polygon_Vertex_R* v = V.front();

		// Determines if v is constrained by a line,
		// and determines the location of the event on the plane.

		CGAL_Point_2 V_t;
		Constraint C_0;
		get_constraint_and_location_at_time_t(v, intersected, t_intersectant, C_0, V_t);

		std::list<Intersection_Line*> I;
		for (Event_Vertex_Line* e : E) I.push_back(get_line_by_identifier(e->intersected));

		try {
			Intersection_Line* I_0 = C_0.first;

			if (I_0 != nullptr) {
				// Case C2.
				// At time t, a constrained vertex v intersects a set of lines I.
				if (Universe::params->print_schedule) print_events(E, "C2");
				constrained_vertex_intersects_line(E, v, V_t);

			} else if (Polygon_Vertex_R* v_n = v->get_constrained_neighbor(I_0, I)) {
				// Case B.
				// At time t, an unconstrained vertex intersects a set of lines I,
				// and a vertex constrained by one of the lines of I.
				if (Universe::params->print_schedule) print_events(E, "B");
				unconstrained_vertex_intersects_line_kth_time(E, v, V_t);

			} else {
				// Case A.
				// At time t, an unconstrained vertex intersects a set of lines I.
				if (Universe::params->print_schedule) print_events(E, "A");
				unconstrained_vertex_intersects_line(E, v, V_t);
			}
		} catch (std::exception & except) {
			throw except;
		}
	}



	void Support_Plane::process_event_2(const std::list<Event_Vertex_Line*> & E, const std::list<Polygon_Vertex_R*> & V)
	{
		// We associate a list of simultaneous events, involving two neighbor vertices, to an appropriate handler.

		Polygon_Vertex_R *v1 = V.front(), *v2 = V.back();

		std::list<Event_Vertex_Line*> E_1, E_2;
		split_simultaneous_events(E, v1, v2, E_1, E_2);

		try {
			if (Polygon_Vertex_R::constrained_vertices_meet(E_1, E_2, v1, v2)) {
				// Case C1.
				// At time t, two constrained vertices intersect in a region corner,
				// defined as the intersection of at least two intersection lines.
				if (Universe::params->print_schedule) print_events(E, "C1");
				constrained_vertices_intersect(E_1, E_2, v1, v2);
			} else {
				// Case D1.
				// Two vertices, forming an edge of non-null length, intersect one common line.
				if (Universe::params->print_schedule) print_events(E, "D1");
				edge_intersects_line(E_1, E_2, v1, v2);
			}

		} catch (std::exception & except) {
			throw except;
		}
	}


	void Support_Plane::process_event_3(const std::list<Event_Vertex_Line*> & E_ref, const std::list<Polygon_Vertex_R*> & V_ref)
	{
		// We associate a list of simultaneous events, involving three neighbor vertices, to an appropriate handler.

		Polygon_Vertex_R* v1, *v, *v2;
		std::list<Event_Vertex_Line*> E_1, E, E_2;

		try {
			if (Polygon_Edge::two_edges_intersect_two_lines(E_ref, V_ref, E_1, E, E_2, v1, v, v2)) {
				// Case D2.
				// Two adjacent edges of a subpolygon intersect two adjacent edges.
				if (Universe::params->print_schedule) print_events(E_ref, "D2");
				two_edges_intersect_two_lines(E_1, E, E_2, v1, v, v2);
			} else {
				throw std::logic_error("Couldn't handle a set of simultaneous events");
			}

		} catch (std::exception & except) {
			throw except;
		}
	}



	void Support_Plane::get_constraint_and_location_at_time_t(Polygon_Vertex_R* v, const int intersected, const FT & t_intersectant, Constraint & C_res, CGAL_Point_2 & V_res)
	{
		if (v->unconstrained()) {
			C_res = Constraint();
			V_res = v->pt(t_intersectant);
		} else {
			C_res = v->get_constraint();
			Intersection_Line* J = get_line_by_identifier(intersected);
			if (Universe::params->use_landmarks) {
				V_res = get_landmark(C_res.first, J);
			} else {
				V_res = get_intersection_point(C_res.first, J);
			}
		}
	}



	void Support_Plane::get_constraint_and_location_at_time_t(Polygon_Vertex_R* v, Intersection_Line* I_intersect, const FT & t_intersectant, Constraint & C_res, CGAL_Point_2 & V_res)
	{
		if (v->unconstrained()) {
			C_res = Constraint();
			V_res = v->pt(t_intersectant);
		} else {
			C_res = v->get_constraint();
			if (Universe::params->use_landmarks) {
				V_res = get_landmark(C_res.first, I_intersect);
			} else {
				V_res = get_intersection_point(C_res.first, I_intersect);
			}
		}
	}



	void Support_Plane::get_references_to_simulatenous_events_for_this_polygon(Polygon_Vertex_R* v, const FT & t, std::list<Event_Vertex_Line*> & E)
	{
		// We are looking for references to events which occur at time t,
		// and are related to the same polygon as v's.
		Polygon* P_ref = v->get_polygon();

		// Searches events related to this plane, occuring at time t.
		Event_Queue* Q = Universe::event_queue;
		std::list<Event_Vertex_Line*> E_0;
		Q->get_references_to_simultaneous_events_for_given_plane(id, t, E_0);
		
		// Filters events.
		for (std::list<Event_Vertex_Line*>::iterator it_e = E_0.begin() ; it_e != E_0.end() ; ++it_e) {
			Event_Vertex_Line* e_0 = (*it_e);
			Polygon_Vertex_R* v_0 = vertices_r[e_0->intersectant];
			if (v == v_0 || v_0->get_polygon() == P_ref) {
				E.push_back(e_0);
			}
		}
	}



	int Support_Plane::get_number_of_vertices_simultaneously_intersecting_lines(Polygon_Vertex_R* v, Event_Vertex_Line* e,
		const std::list<Event_Vertex_Line*> & E, std::list<Polygon_Vertex_R*> & V) const
	{
		// Given a list of events E, happening simultaneously to another event e,
		// we count the number of distinct vertices intersecting lines.

		V.push_back(v);

		if (E.empty()) return 1;

		// We loop on v's neighbors, on v->e1's side.
		// We iterate every time we find an event associated to v1.
		// Otherwise we break the loop and do the same operation on v->e2's side.

		int n = int(v->get_polygon()->vertices.size());
		int n1 = loop_on_neighbor_vertices_with_simultaneous_events(v, v->e1, e, E, V);
		if (n1 == n - 1) {
			return n;
		} else {
			int n2 = loop_on_neighbor_vertices_with_simultaneous_events(v, v->e2, e, E, V);
			return 1 + n1 + n2;
		}
	}



	int Support_Plane::loop_on_neighbor_vertices_with_simultaneous_events(Polygon_Vertex_R* v_init, Polygon_Edge* e_init, Event_Vertex_Line* e, const std::list<Event_Vertex_Line*> & E, std::list<Polygon_Vertex_R*> & V) const
	{
		int n = 0;

		// Given a vertex, we loop on the neighbor vertices along one of its edges,
		// as long as the following conditions are satisfied :
		// - v1 is a running vertex
		// - v1 is not v_init (which means that we've looped on all the vertices)
		// - v1 intersects at least one of the lines intersected by v1_prev.

		Polygon_Vertex* v1_prev = v_init;
		Polygon_Vertex* v1 = e_init->other_vertex(v_init);

		std::list<int> I_prev;
		I_prev.push_back(e->intersected);
		if (v_init->constrained()) {
			I_prev.push_back(v_init->get_constraint().first->id_object);
		}
		for (Event_Vertex_Line* e_vl : E) {
			if (e_vl->intersectant == v_init->id_object) {
				I_prev.push_back(e_vl->intersected);
			}
		}

		while (v1->to_r() != nullptr && v1 != v_init) {

			int id_1 = v1->id_object;
			std::list<int> I;
			for (std::list<Event_Vertex_Line*>::const_iterator it_e = E.begin() ; it_e != E.end() ; ++it_e) {
				if ((*it_e)->intersectant == id_1) {
					I.push_back((*it_e)->intersected);
				}
			}

			// Does v1 intersect at least one line that is also intersected by v1_prev ?

			bool lists_intersect = false;
			for (std::list<int>::iterator it = I.begin() ; it != I.end() ; ++it) {
				if (std::find(I_prev.begin(), I_prev.end(), *it) != I_prev.end()) {
					lists_intersect = true;
				}
			}

			if (lists_intersect) {
				V.push_back(v1->to_r());
				++n;
			} else {
				break;
			}

			Polygon_Vertex* v1_next = v1->e1->other_vertex(v1);
			if (v1_next == v1_prev) v1_next = v1->e2->other_vertex(v1);
			
			I_prev = I;
			v1_prev = v1; 
			v1 = v1_next;
		}

		return n;
	}


	
	void Support_Plane::pop_fitered_simultaneous_events_from_queue(std::list<Event_Vertex_Line*> & E, const std::list<Polygon_Vertex_R*> & V) const
	{
		// We have identified a list of events E simultaneous to another event e,
		// and a list V of vertices whose collisions are going to be handled.

		// Events of E that involve elements of V are going to be removed from the queue.
		// Others are removed from the list E.

		Event_Queue* Q = Universe::event_queue;

		std::list<Event_Vertex_Line*>::iterator it_e = E.begin();
		while (it_e != E.end()) {
			Event_Vertex_Line* e = (*it_e);

			std::map<int, Polygon_Vertex_R*>::const_iterator it_v = vertices_r.find(e->intersectant);
			assert(it_v != vertices_r.end());

			Polygon_Vertex_R* v = it_v->second;
			if (std::find(V.begin(), V.end(), v) != V.end()) {
				Q->erase(e);
				++it_e;
			} else {
				it_e = E.erase(it_e);
			}
		}
	}



	void Support_Plane::split_simultaneous_events(const std::list<Event_Vertex_Line*> & E, Polygon_Vertex_R* v1, Polygon_Vertex_R* v2, 
		std::list<Event_Vertex_Line*> & E_1, std::list<Event_Vertex_Line*> & E_2) const
	{
		// We suppose that E contains a set of simultaneous events related to v1 and v2.
		// We decompose E into subsets E_1 and E_2 depending on e->intersectant.

		int id_1 = v1->id_object, id_2 = v2->id_object;

		for (std::list<Event_Vertex_Line*>::const_iterator it_e = E.begin() ; it_e != E.end() ; ++it_e) {
			Event_Vertex_Line* e = (*it_e);
			if (e->intersectant == id_1) {
				E_1.push_back(e);
			} else {
				E_2.push_back(e);
			}
		}
	}


#if 0
	void Support_Plane::process_event(Event_Vertex_Line* e)
	{
		int intersectant = e->intersectant;
		int intersected = e->intersected;
		const FT t_intersectant = e->t_intersectant;

		Polygon_Vertex_R* v = vertices_r[intersectant];

		Intersection_Line* I_0;
		CGAL_Point_2 V_t;

		if (v->unconstrained()) {
			I_0 = nullptr;
			V_t = v->pt(t_intersectant);
		} else {
			I_0 = v->get_constraint().first;
			Intersection_Line* I_1 = get_line_by_identifier(intersected);
			if (Universe::params->use_landmarks) {
				V_t = get_landmark(I_0, I_1);
			} else {
				V_t = get_intersection_point(I_0, I_1);
			}
		}
		
		// Intersection_Line* I_0 = (v->unconstrained() ? nullptr : v->get_constraint().first);
		// CGAL_Point_2 V_t = v->pt(t_intersectant);

		bool verbose = false;

		// Step 1.
		// In this function, we would like to associate the correct handle to the event e that has been popped from the queue.
		// We first factorize events, by obtaining the list of lines I intersected by vertex v at time t.

		std::list<Intersection_Line*> I;
		std::list<Event_Vertex_Line*> E;
		get_simultaneous_events(v, t_intersectant, V_t, E, e, I);
		v->decrement_queued_events(int(E.size()));

		// Step 2.
		// Performs different tests, to associate the current handle to E_VL.

		try {
			if (Polygon_Vertex_R* v_n = v->get_constrained_neighbor(I_0, I)) {
				if (I_0 == nullptr) {
					// Case B.
					// An unconstrained vertex v intersects a set of lines I_k,
					// and meets at the same time one of its neighbors v_n, constrained by one of these lines.
					// The vertex is going to be passed to the adjacent polygon, if it exists.
					if (verbose) std::cout << "[" << id << "] " << e->intersectant << " " << e->intersected << " " << e->t_intersectant << " B" << std::endl;
					unconstrained_vertex_intersects_line_kth_time(E, v, V_t);

				} else {
					// Case C1.
					// A constrained vertex intersects a set of lines I_k,
					// and meets at the same time one of its neighbors, constrained by one of these lines.
					// The two constrained vertices should be merged into a double-constrained and still vertex.
					if (verbose) std::cout << "[" << id << "] " << e->intersectant << " " << e->intersected << " " << e->t_intersectant << " C1" << std::endl;
					constrained_vertices_intersect(E, v, v_n, V_t);
				}

			} else {
				if (Polygon_Vertex_R* v_p = v->get_neighbor_intersecting_identical_line(I_0, I, t_intersectant)) {
					// Case D.
					// There exists a line intersected, at the same time, by two vertices v and v_p.
					// We should actually process a collision between an edge and a line.
					if (verbose) {
						std::cout << "[" << id << "] " << e->intersectant << " " << e->intersected << " " << e->t_intersectant << " D" << std::endl;
						std::cout << "[" << id << "] " << v_p->id_object << " " << e->intersected << " " << e->t_intersectant << " D *" << std::endl;
					}
					edge_intersects_line(E, v, v_p, V_t);

				} else if (I_0 == nullptr) {
					// Case A.
					// There exists no constrained neighbor, 
					// and no vertex intersecting one of the lines intersected by v at time t.
					// There is only a unconstrained vertex, intersecting a set of lines.
					if (verbose) std::cout << "[" << id << "] " << e->intersectant << " " << e->intersected << " " << e->t_intersectant << " A" << std::endl;
					unconstrained_vertex_intersects_line(E, v, V_t);

				} else {
					// Case C2.
					// There exists no neighbor constrained by one of the lines of I,
					// and there exists no vertex intersecting of the lines intersected by v at time t.
					if (verbose) std::cout << "[" << id << "] " << e->intersectant << " " << e->intersected << " " << e->t_intersectant << " C2" << std::endl;
					constrained_vertex_intersects_line(E, v, V_t);
				}
			}

		} catch (const std::exception & except) {
			// Prints the support plane
			std::cout << "Error raised while processing [" << id << "] " << intersectant << " " << intersected << " " << t_intersectant << std::endl;
			draw(to_double(t_intersectant), 0.5, 20000, 0.5, 0.5, 10);
			throw except;
		}

		if (Universe::params->rt_check) {
			real_time_check(t_intersectant);
		}
	}
#endif



	void Support_Plane::get_simultaneous_events(Polygon_Vertex_R* v, const FT & t, const CGAL_Point_2 & V_t,
		std::list<Event_Vertex_Line*> & E, Event_Vertex_Line* e, std::list<Intersection_Line*> & I)
	{
		Event_Queue* Q = Universe::event_queue;

		I.clear();
		E.clear();

		// We have popped an event e from the queue
		// Now, we determine if some other events, involving the same vertex, occur simultaneously

		std::list<Event_Vertex_Line*> E_0;
		Q->get_simultaneous_events_for_this_vertex(id, e->intersectant, e->t_intersectant, E_0);

		for (std::list<Event_Vertex_Line*>::iterator it_e = E_0.begin(); it_e != E_0.end(); it_e++) {
			Event_Vertex_Line* e_p = (*it_e);
			Intersection_Line* I_p = get_line_by_identifier(e_p->intersected);

			if (I_p->includes(V_t)) {
				// If the intersection of I_ref and v is the same as the intersection of I and v,
				// then e_ref and e occur simultaneously. We remove e from the queue and add it to E.
				// Q->erase(e_p);
				E.push_back(e_p);
				I.push_back(I_p);
			}
		}

		// Finally adds e_ref to E and I
		if (e != nullptr) {
			E.push_back(e);
			I.push_back(get_line_by_identifier(e->intersected));
		}
	}



	void Support_Plane::get_simultaneous_events_as_edge_intersects_line(Polygon_Vertex_R* v, const FT & t, const CGAL_Point_2 & V_t,
		Intersection_Line* I_L, std::list<Event_Vertex_Line*> & E, Event_Vertex_Line* & e_vl)
	{
		e_vl = nullptr;
		E.clear();

		Event_Queue* Q = Universe::event_queue;

		// We support that an edge e = (v0 v) intersects the line I, at time t.
		// We've already popped the event (v0 intersects I_L), we should now pop (v intersects I_L)

		Q->get_simultaneous_events_for_this_vertex(id, v->id_object, t, E);
		for (std::list<Event_Vertex_Line *>::iterator it_e = E.begin(); it_e != E.end(); it_e++) {
			if ((*it_e)->intersected == I_L->id_object) {
				e_vl = (*it_e);
				break;
			}
		}
	}




	void Support_Plane::append(Event_Vertex_Line* e_vl, std::string type)
	{
		FILE* file = fopen("events.txt", "a+");
		if (file != NULL) {
			fprintf(file, "[%i] %i %i %lf %s\n", id, e_vl->intersectant, e_vl->intersected, to_double(e_vl->t_intersectant), type.c_str());
			fclose(file);
		}
	}



	void Support_Plane::get_polygon_description(std::list<std::list<CGAL_Point_3> > & P, std::list<CGAL_Color> & C, const double t)
	{
		polygon_set->get_polygon_description(P, C, t);
	}


	void Support_Plane::get_polygon_description(Polygon_Vertex_Octree* V, std::list<std::list<int> > & P, const FT & t)
	{
		polygon_set->get_polygon_description(V, P, t);
	}


	void Support_Plane::group_final_polygons()
	{
		// We suppose that the propagation phase is over.
		// It is time to group adjacent polygons into subsets, 
		// which depends on the existence of a segment between two adjacent polygons.

		// Step 1.
		// We initialize a graph that connects every pair of adjacent polygons
		// that are not separated by a graph.
		
		polygon_set->build_graph();

		// Step 2.
		// We explore the edges of the graph, to get the connected components of the graph.
		// These connected components represent the desired Polygon_Groups.

		polygon_set->get_groups(groups);

		// Step 3.
		// Clears the graph.

		polygon_set->clear_graph();
	}



	void Support_Plane::check()
	{
		// We assume that this function is called as all vertices have stopped propagating.
		
		for (auto it_c = polygon_set->cells_begin() ; it_c != polygon_set->cells_end() ; ++it_c) {
			const Signature & S = it_c->first;
			Polygon_Node* C = it_c->second;

			// Test 1 : for each polygon contained in each cell,
			// we assert that a vertex hasn't silently crossed a line

			for (std::list<Polygon*>::const_iterator it_p = C->polygons_begin(); it_p != C->polygons_end(); ++it_p) {
				Polygon* P = (*it_p);
				for (std::list<Polygon_Vertex*>::const_iterator it_v = P->vertices.begin(); it_v != P->vertices.end(); ++it_v) {
					Polygon_Vertex* v = (*it_v);

					// Checks that a vertex hasn't silently crossed a line
					const CGAL_Point_2 & M = v->get_M();
					for (std::map<Intersection_Line*, int>::iterator it_l = polygon_set->dictionary.begin(); it_l != polygon_set->dictionary.end(); it_l++) {
						Intersection_Line* I = it_l->first;
						bool SL = S[it_l->second];
						FT d = I->a() * M.x() + I->b() * M.y() + I->c();
						if (SL && (d < 0) || (!SL && (d > 0))) {
							std::cout << "Plane " << id << " : error with vertex : " << (*it_v)->id_object << " : missed line detected." << std::endl;
							assert(0);
						}
					}
				}
			}

			// Test 2 : when a cell is on the outer ring of the polygon set,
			// we assert that there exists a couple of segments of opposite signs that border it.

			for (std::map<Signature, Polygon_Node*, Vector_Bool_Comparator>::const_iterator it_n = polygon_set->cells_begin() ; it_n != polygon_set->cells_end() ; ++it_n) {
				Polygon_Node* N = it_n->second;
				size_t cs = N->contours_size();
				for (size_t i = 0 ; i < cs ; ++i) {
					Constraint C_i = N->get_contour(i);

					Intersection_Line* I = C_i.first;
					if (I->is_border) continue;

					Signature S_adj = polygon_set->get_adjacent_polygons_signature(N->get_one_polygon(), I);
					if (polygon_set->exists(S_adj)) continue;

					Constraint C_i_prev = N->get_contour(i == 0 ? cs - 1 : i - 1);
					Constraint C_i_next = N->get_contour(i == cs - 1 ? 0 : i + 1);
					CGAL_Point_2 A = get_intersection_point(C_i_prev.first, C_i.first);
					CGAL_Point_2 B = get_intersection_point(C_i.first, C_i_next.first);
					CGAL_Point_2 M = CGAL::midpoint(A, B);

					if (!I->exist_segments_including_point_outside_intersections(M, FLT_MAX)) {
						std::cout << "Plane " << id << " : hole detected near points (" << A.x() << ", " << A.y() << ") and (" << B.x() << ", " << B.y() << ") " << std::endl;
						std::cout << "       i.e. at intersection of lines " << C_i_prev.first->id_object << ", " << C_i.first->id_object << ", " << C_i_next.first->id_object << std::endl;
						assert(0);
					}
				}
			}
		}
	}



	void Support_Plane::real_time_check(const FT & t)
	{
		// For each cell
		for (std::map<Signature, Polygon_Node*, Vector_Bool_Comparator>::const_iterator it_c = polygon_set->cells_begin(); it_c != polygon_set->cells_end(); it_c++) {
			const Signature & S = it_c->first;
			Polygon_Node* C = it_c->second;

			// For each polygon contained in the cell
			for (std::list<Polygon*>::const_iterator it_p = C->polygons_begin(); it_p != C->polygons_end(); it_p++) {
				Polygon* P = (*it_p);
				for (std::list<Polygon_Vertex*>::iterator it_v = P->vertices.begin(); it_v != P->vertices.end(); it_v++) {
					if (Polygon_Vertex_R* v = (*it_v)->to_r()) {

						// Checks that a vertex hasn't silently crossed a line
						CGAL_Point_2 V_t = v->pt(t);
						for (std::map<Intersection_Line*, int>::iterator it_l = polygon_set->dictionary.begin(); it_l != polygon_set->dictionary.end(); it_l++) {
							Intersection_Line* I = it_l->first;
							if (S[it_l->second]) {
								if (I->a() * V_t.x() + I->b() * V_t.y() + I->c() < 0) {
									std::cout << "Error with vertex : " << (*it_v)->id_object << " : missed line at t = " << t << std::endl;
									draw(to_double(t), 0.5, 20000, 0.5, 0.5, 10);
									exit(0);
								}
							} else {
								if (I->a() * V_t.x() + I->b() * V_t.y() + I->c() > 0) {
									std::cout << "Error with vertex : " << (*it_v)->id_object << " : missed line at t = " << t << std::endl;
									draw(to_double(t), 0.5, 20000, 0.5, 0.5, 10);
									exit(0);
								}
							}
						}
					}
				}
			}
		}
	}
}