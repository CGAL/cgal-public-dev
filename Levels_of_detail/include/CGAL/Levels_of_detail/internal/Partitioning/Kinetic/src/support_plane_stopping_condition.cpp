#include "../include/support_plane.h"
#include "../include/intersection_line.h"
#include "../include/polygon_vertex.h"
#include "../include/polygon_edge.h"

#include "../include/universe.h"
#include "../include/parameters.h"
#include "../include/octree_base_vertex.h"
#include "../include/point_cloud_vertex.h"


namespace Skippy {

	bool Support_Plane::K_based_stopping_condition(Polygon_Vertex_R* v, bool lookup) const
	{
		// In this function, we determine if vertex v should keep propagating.
		// We suppose that it is called as v intersects a couple of segments of opposite signs.

		// The boolean lookup is specific to events of type D.
		// Indeed some functions have been factorized and in the current algorithm,
		// we want to prevent K from being decremented twice.

		if (v->K == 0) {
			return false;
		} else {
			if (!lookup) {
				--v->K;
			}
			return true;
		}
	}



	bool Support_Plane::density_based_stopping_condition(Polygon_Vertex_R* v, const CGAL_Point_2 & V_t, const CGAL_Point_3 & W_t) const
	{
		if (Universe::params->stopping_condition == 1) {
			return density_based_stopping_condition_box(v, V_t, W_t);
		} else {
			return density_based_stopping_condition_cone(v, V_t, W_t);
		}
	}



	bool Support_Plane::density_based_stopping_condition_box(Polygon_Vertex_R* v, const CGAL_Point_2 & V_t, const CGAL_Point_3 & W_t) const
	{
		// At time t, a vertex v intersects a line I in V_t with a propagation speed dM.
		// We would like to define a rectangular cuboid and count the number of input points
		// that belong to this cuboid. Their faces should be expressed as planes.

		// Part 1.
		// We are interested in defining a local 3D frame in W_t whose two vectors
		// correspond to dM and the orthogonal vector to this support plane.

		const CGAL_Point_3 & O = W_t;
		CGAL_Vector_3 j = backproject(V_t + v->dM) - W_t;
		CGAL_Vector_3 k = plane.orthogonal_vector();
		CGAL_Vector_3 i = CGAL::cross_product(j, k);

		const FT & a = Universe::params->density_box_width;
		const FT & b = Universe::params->density_box_length;
		const FT & c = Universe::params->density_box_height;

		double k_x = CGAL::to_double(k.x());
		double k_y = CGAL::to_double(k.y());
		double k_z = CGAL::to_double(k.z());
		double k_l = sqrt(k_x * k_x + k_y * k_y + k_z * k_z);

		// Part 2.
		// In the local 3D frame (O ni nj nk) where ni, nj and nk are the normalized versions of i, j and k,
		// we consider the 8 points with the following coordinates :
		// P[0 .. 3] : A- = (-a/2, 0, -c/2), D- = (a/2, 0, -c/2), B- = (-a/2, b, -c/2), C- = (a/2, b, -c/2),
		// P[4 .. 7] : A+ = (-a/2, 0,  c/2), D+ = (a/2, 0,  c/2), B+ = (-a/2, b,  c/2), C+ = (a/2, b,  c/2).

		// As we work with i, j and k which are not normalized,
		// we should find approximate coordinates of these 8 points with the exact kernel.

		// We want |t_a| such that O + t_a * i is at distance 'a' from O and belongs to a plane parallel to (O j k),
		//         |t_b| such that O + t_b * j is at distance 'b' from O and belongs to a plane parallel to (O i k),
		//         |t_c| such that O + t_c * k is at distance 'c' from O and belongs to a plane parallel to (O i j).

		FT half_t_a = project_to_parallel_plane(a, i) / FT(2);
		FT t_b = project_to_parallel_plane(b, j);
		FT half_t_c = project_to_parallel_plane(c, k) / FT(2);

		std::vector<CGAL_Point_3> P;
		for (int r = 0; r < 8; r++) {
			FT x = (r % 2 == 0 ? -half_t_a : half_t_a);
			FT y = ((r / 2) % 2 == 0 ? 0 : t_b);
			FT z = (r > 3 ? -half_t_c : half_t_c);
			P.push_back(O + x * i + y * j + z * k);
		}

		// Step 3.
		// Reasons on the minimal/maximal coordinates of the 8 points along all axes
		// in order to formulate a query for the octree

		double _query[6] = { FLT_MAX, -FLT_MAX, FLT_MAX, -FLT_MAX, FLT_MAX, -FLT_MAX };
		std::vector<double> query(_query, _query + sizeof(_query) / sizeof(double));

		for (int r = 0; r < 8; r++) {
			double x = CGAL::to_double(P[r].x()), y = CGAL::to_double(P[r].y()), z = CGAL::to_double(P[r].z());
			if (x < query[0]) query[0] = x;
			if (x > query[1]) query[1] = x;
			if (y < query[2]) query[2] = y;
			if (y > query[3]) query[3] = y;
			if (z < query[4]) query[4] = z;
			if (z > query[5]) query[5] = z;
		}

		std::list<Octree_Base_Vertex*> C;
		assert(Universe::point_cloud_octree != nullptr);
		Universe::point_cloud_octree->search(query, C);

		if (int(C.size()) < Universe::params->density_pts) {
			return false;
		}

		// Step 4.
		// Filters the results returned by the previous query,
		// by discarding the points that are not inside the six planes that delimit the previous cuboid.
		// First we set the plane equations. For each plane the three first coefficients are provided by [ijk].xyz(),
		// the last one can be deduced from a point included in the plane.

		FT I_inf = -(i.x() * P[0].x() + i.y() * P[0].y() + i.z() * P[0].z()), I_sup = -(i.x() * P[7].x() + i.y() * P[7].y() + i.z() * P[7].z());
		FT J_inf = -(j.x() * P[0].x() + j.y() * P[0].y() + j.z() * P[0].z()), J_sup = -(j.x() * P[7].x() + j.y() * P[7].y() + j.z() * P[7].z());
		FT K_inf = -(k.x() * P[0].x() + k.y() * P[0].y() + k.z() * P[0].z()), K_sup = -(k.x() * P[7].x() + k.y() * P[7].y() + k.z() * P[7].z());

		// Now loops on the points.
		// Initially, we consider that all points are between each couple of plane (I_inf, I_sup), (J_inf, J_sup), (K_inf, K_sup).
		// We decrease the number of points if it appears that one of these relations is not satisfied.

		int inside = int(C.size());

		for (auto it_p = C.begin(); it_p != C.end(); ++it_p) {
			if (inside < Universe::params->density_pts) return false;

			Point_Cloud_Vertex* PCV = dynamic_cast<Point_Cloud_Vertex*>(*it_p);
			assert(PCV != nullptr);
			const CGAL_Point_3 & M = PCV->M;
			const CGAL_Inexact_Vector_3 & M_nor = PCV->hint_normal;

			FT d, d_inf, d_sup;

			// Fourth test : angle deviation
			double deviation = fabs(M_nor.x() * k_x + M_nor.y() * k_y + M_nor.z() * k_z) / k_l;
			if (deviation < Universe::params->density_deviation) {
				--inside;
				continue;
			}

			// First test : I_inf, I_sup
			d = -(i.x() * M.x() + i.y() * M.y() + i.z() * M.z());
			if (!jin(I_inf, d, I_sup) && !jin(I_sup, d, I_inf)) {
				--inside; 
				continue;
			}

			// Second test : J_inf, J_sup
			d = -(j.x() * M.x() + j.y() * M.y() + j.z() * M.z());
			if (!jin(J_inf, d, J_sup) && !jin(J_sup, d, J_inf)) {
				--inside; 
				continue;
			}

			// Third test : K_inf, K_sup
			d = -(k.x() * M.x() + k.y() * M.y() + k.z() * M.z());
			if (!jin(K_inf, d, K_sup) && !jin(K_sup, d, K_inf)) {
				--inside; 
				continue;
			}
		}

		return (inside >= Universe::params->density_pts);
	}



	bool Support_Plane::density_based_stopping_condition_cone(Polygon_Vertex_R* v, const CGAL_Point_2 & V_t, const CGAL_Point_3 & W_t) const
	{
		// At time t, a vertex v intersects a line I in V_t with a propagation speed dM.
		// We would like to define a cone and count the number of input points it includes.

		// Part 1.
		// We are interested in defining a local 3D frame in W_t whose two vectors
		// correspond to dM and the orthogonal vector to this support plane.

		const CGAL_Point_3 & O = W_t;
		CGAL_Vector_3 j = backproject(V_t + v->dM) - W_t;
		CGAL_Vector_3 k = plane.orthogonal_vector();
		CGAL_Vector_3 i = CGAL::cross_product(j, k);

		const FT & a = Universe::params->density_cone_base;
		const FT & b = Universe::params->density_cone_height;

		// Part 2.
		// In the local 3D frame (O ni nj nk) where ni, nj and nk are the normalized versions of i, j and k,
		// we consider the 8 points with the following coordinates :
		// P[0 .. 3] : A- = (-a/2, 0, -a/2), D- = (a/2, 0, -a/2), B- = (-a/2, b, -a/2), C- = (a/2, b, -a/2),
		// P[4 .. 7] : A+ = (-a/2, 0,  a/2), D+ = (a/2, 0,  a/2), B+ = (-a/2, b,  a/2), C+ = (a/2, b,  a/2).

		// As we work with i, j and k which are not normalized,
		// we should find approximate coordinates of these 8 points with the exact kernel.

		// We want |t_a| such that O + t_a * i is at distance 'a' from O and belongs to a plane parallel to (O j k),
		//         |t_b| such that O + t_b * j is at distance 'b' from O and belongs to a plane parallel to (O i k)
		
		FT half_t_a = project_to_parallel_plane(a, i) / FT(2);
		FT t_b = project_to_parallel_plane(b, j);

		std::vector<CGAL_Point_3> P;
		for (int r = 0; r < 8; r++) {
			FT x = (r % 2 == 0 ? -half_t_a : half_t_a);
			FT y = ((r / 2) % 2 == 0 ? 0 : t_b);
			FT z = (r > 3 ? -half_t_a : half_t_a);
			P.push_back(O + x * i + y * j + z * k);
		}

		// Step 3.
		// Reasons on the minimal/maximal coordinates of the 8 points along all axes
		// in order to formulate a query for the octree

		double _query[6] = { FLT_MAX, -FLT_MAX, FLT_MAX, -FLT_MAX, FLT_MAX, -FLT_MAX };
		std::vector<double> query(_query, _query + sizeof(_query) / sizeof(double));

		for (int r = 0; r < 8; r++) {
			double x = CGAL::to_double(P[r].x()), y = CGAL::to_double(P[r].y()), z = CGAL::to_double(P[r].z());
			if (x < query[0]) query[0] = x;
			if (x > query[1]) query[1] = x;
			if (y < query[2]) query[2] = y;
			if (y > query[3]) query[3] = y;
			if (z < query[4]) query[4] = z;
			if (z > query[5]) query[5] = z;
		}

		std::list<Octree_Base_Vertex*> C;
		assert(Universe::point_cloud_octree != nullptr);
		Universe::point_cloud_octree->search(query, C);

		if (int(C.size()) < Universe::params->density_pts) {
			return false;
		}

		// Step 4.
		// Filters the results returned by the previous query,
		// by discarding the points that are not inside the cone defined by W_t, a and b.

		int inside = int(C.size());

		for (auto it_p = C.begin() ; it_p != C.end() ; ++it_p) {
			if (inside < Universe::params->density_pts) return false;

			const CGAL_Point_3 & M = (*it_p)->M;

			// Finds the intersection of the line (W_t, j) and the plane Q defined by (M, j) :
			// it is defined by l s.t. N = W_t + lj belongs to Q.

			CGAL_Plane Q(M, j);
			FT l = -(Q.a() * W_t.x() + Q.b() * W_t.y() + Q.c() * W_t.z() + Q.d()) / (Q.a() * j.x() + Q.b() * j.y() + Q.c() * j.z());

			if (l < 0 || l > t_b) {
				--inside;
				continue;
			}

			CGAL_Point_3 N = W_t + l * j;
			CGAL_Vector_3 MN = N - M;

			// Computes the maximal distance between M and N

			FT r_mn = a * l / b;
			if (MN.squared_length() > r_mn * r_mn) {
				--inside;
				continue;
			}
		}

		return (inside >= Universe::params->density_pts);
	}



	FT Support_Plane::project_to_parallel_plane(const FT & D, const CGAL_Vector_3 & n) const
	{
		// For a point M (x, y, z) in a plane (P) ax + by + cz + d = 0,
		// We want N in a plane (P') ax + by + cz + d' = 0 at a distance D from M.
		// Such point is obtained by computing N = M + t * n (= vec(a, b, c)) where t = D / sqrt(a2 + b2 + c2).

		FT s = n.x() * n.x() + n.y() * n.y() + n.z() * n.z();

		double s_d = CGAL::to_double(s);
		double sqrt_s_d = sqrt(s_d);

		std::stringstream stream;
		stream << std::setprecision(15);
		stream << sqrt_s_d;

		// std::string str_sqrt_s = stream.str();
		FT sqrt_s;
		stream >> sqrt_s;

		return D / sqrt_s;
	}



	bool Support_Plane::propagation_continues_outside_intersections(Intersection_Line* I, Polygon_Vertex_R* v, bool lookup, const CGAL_Point_2 & V_t, const CGAL_Point_3 & W_t, const FT & t, const std::vector<bool> & S, const int seed) const
	{
		// A vertex v intersects a line I at time t in V_t.

		// First of all, we evaluate the three conditions when the propagation must, or must not continue.
		// If no definitive answer is returned, then we evaluate the stopping condition for v and return its value.

		if (I->is_border || polygon_set->exists(S, seed)) {
			// If I is a border, or if there's already a polygon beyond I, then v must not propagate.
			return false;

		} else {
			if (!I->exist_segments_including_point_outside_intersections(V_t, t)) {
				// If there is no segment on I, then v must propagate.
				return true;

			} else {
				// The consistency of the data structure is preserved.
				// At this stage we know that v intersects a segment, so we simply evaluate the stopping condition of v in this point.
				if (Universe::params->stopping_condition == 0) {
					return K_based_stopping_condition(v, lookup);
				} else {
					return density_based_stopping_condition(v, V_t, W_t);
				}
			}
		}
	}



	bool Support_Plane::propagation_continues_at_intersection(Intersection_Line* I, Polygon_Vertex_R* v, bool lookup, const CGAL_Point_2 & V_t,
		const CGAL_Point_3 & W_t, const Constraint & C, const FT & t, const std::vector<bool> & S, const int seed) const
	{
		// Similarly as before,
		// some conditions determine whether the vertex v must continue or must stop.
		// If none of these conditions returns a definite answer, then we evaluate the stopping condition of v.

		if (I->is_border || polygon_set->exists(S, seed)) {
			return false;

		} else {
			if (!I->exist_segments_including_point_at_intersection(V_t, C, t)) {
				return true;

			} else {
				if (v != nullptr) {
					// The consistency of the data structure is preserved in all cases.
					// It is now up to the stopping condition of v.
					if (Universe::params->stopping_condition == 0) {
						return K_based_stopping_condition(v, lookup);
					} else {
						return density_based_stopping_condition(v, V_t, W_t);
					}

				} else {
					// Obviously, the function is called in a context of hole prevention.
					// I is not a border, there doesn't exist a polygon in S but there exists segments in V_t
					// So we don't have to propagate
					return false;
				}	
			}
		}
	}



	bool Support_Plane::propagation_continues_at_intersection(Intersection_Line* I, Polygon_Vertex_R* v, bool lookup, const CGAL_Point_2 & V_t,
		const CGAL_Point_3 & W_t, const std::list<Constraint> & C_limits, const FT & t, const std::vector<bool> & S, const int seed) const
	{
		// Same function as above, in the case of a multi-line intersection

		if (I->is_border || polygon_set->exists(S, seed)) {
			return false;

		} else {
			if (!I->exist_segments_including_point_at_intersection(V_t, C_limits, t)) {
				return true;

			} else {
				if (v != nullptr) {
					if (Universe::params->stopping_condition == 0) {
						return K_based_stopping_condition(v, lookup);
					} else {
						return density_based_stopping_condition(v, V_t, W_t);
					}

				} else {
					return false;
				}	
			}
		}
	}
}