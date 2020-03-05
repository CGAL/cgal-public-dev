#include "../include/support_plane_objects.h"
#include "../include/segment.h"
#include "../include/intersection_line.h"
#include "../include/polygon_vertex.h"
#include "../include/polygon_edge.h"
#include "../include/support_plane.h"
#include "../include/universe.h"
#include "../include/event_queue.h"
#include "../include/parameters.h"


namespace Skippy 
{
	Polygon_Segment::Polygon_Segment(const int _id_plane, const int _seed, const Constraint & C_support)
		: Segment(_id_plane, C_support)
	{
		seed = _seed;

		C_init.first = nullptr;
		C_stop.first = nullptr;
		C_crossed = std::list<Intersection_Line *>();

		// This element has just been created and inserted in a list of segments by Segment()
		// So we access to the last element of the list to set the direct access iterator

		if (C_support.second == PLUS) {
			iterator = --C_support.first->segments_plus.end();
		} else {
			iterator = --C_support.first->segments_minus.end();
		}
	}


	Polygon_Segment::~Polygon_Segment() 
	{ 
	
	}


	CGAL_Point_2 Polygon_Segment::end() const
	{
		Support_Plane* SP = Universe::map_of_planes[id_plane];
		if (Universe::params->use_landmarks) {
			return SP->get_landmark(support, C_stop.first);
		} else {
			return SP->get_intersection_point(support, C_stop.first);
		}
	}


	/*const CGAL_Point_2 & Polygon_Segment::end() const
	{
		assert(C_stop.first != nullptr);

		// return Universe::map_of_planes[id_plane]->get_landmark(support, C_stop.first);
		return Universe::map_of_planes[id_plane]->get_intersection_point(support, C_stop.first);
	}*/



	void Polygon_Segment::insert_as_crossed_line(Intersection_Line* I)
	{
		C_crossed.push_back(I);
	}



	void Polygon_Segment::insert_as_crossed_lines(const std::list<Intersection_Line*> & L)
	{
		std::copy(L.begin(), L.end(), std::back_inserter(C_crossed));
	}



	bool Polygon_Segment::includes_point_at_intersection(const CGAL_Point_2 & M, const std::list<Constraint> & C_limits, const FT & t) const
	{
		// Same function as before, 
		// for the case when several lines intersect in M

		for (std::list<Constraint>::const_iterator it_c = C_limits.begin(); it_c != C_limits.end(); it_c++) {
			Constraint C = (*it_c);

			Intersection_Line* J = C.first;
			Sign J_eps = C.second;

			// If J is crossed by the segment, then M belongs to it
			std::list<Intersection_Line*>::const_iterator it_l = std::find(C_crossed.cbegin(), C_crossed.cend(), J);
			if (it_l != C_crossed.cend()) return true;

			// Compares J to the lines that delimit the segment.
			// If they are not nullptr, then the segment includes M 
			// iff the initial or stop constraint have the same sign as C

			Intersection_Line *I_init = C_init.first, *I_stop = C_stop.first;
			if (I_init == J) {
				return (C_init.second == J_eps);
			} else if (I_stop == J) {
				return (C_stop.second == J_eps);
			}
		}

		// Numerical computation
		return includes_point_on_support_line(M, t);
	}



	bool Polygon_Segment::includes_edge(const CGAL_Point_2 & V_1, const CGAL_Point_2 & V_2, const Constraint & C_1, const Constraint & C_2,
		const std::list<Intersection_Line*> & CL) const
	{
		// This function is called when all segments no longer propagate,
		// as we group adjacent polygons of the set to define the facets of the partition.

		// We assume that e = (v1 v2) is delimited by two constraints C_1 and C_2.

		Intersection_Line *I_1 = C_1.first, *I_2 = C_2.first;
		Sign eps_1 = C_1.second, eps_2 = C_2.second;

		// First of all we determine if C_1 and C_2 belong to the list of lines crossed by the segment.
		// If so, we return true.

		if (std::find(C_crossed.cbegin(), C_crossed.cend(), I_1) != C_crossed.cend()) return true;
		if (std::find(C_crossed.cbegin(), C_crossed.cend(), I_2) != C_crossed.cend()) return true;

		// After that we check if I_1 or I_2 are used to define C_init.
		// If so, we compare the sign of the constraint to C_1 and C_2.

		Intersection_Line* I_init = C_init.first;
		if (I_init == I_1) {
			return (C_init.second == eps_1);
		} else if (I_init == I_2) {
			return (C_init.second == eps_2);
		}

		// Same for C_stop.

		Intersection_Line* I_stop = C_stop.first;
		if (I_stop == I_1) {
			return (C_stop.second == eps_1);
		} else if (I_stop == I_2) {
			return (C_stop.second == eps_2);
		}

		// Finally, we focus on lines which are concurrent with this->support, I_1 and I_2.
		// If such lines don't exist, then we return false because we have already listed all the possible cases.

		if (CL.empty() && Universe::params->cl_flush == -1) {
			return false;
		} else {
			// Arithmetic computation
			const CGAL_Point_2 &A = origin(), &B = end();
			CGAL_Point_2 M = CGAL::midpoint(V_1, V_2);

			FT x_a = A.x(), x_b = B.x(), x = M.x();
			if (x_a != x_b) {
				return ((x_a <= x && x <= x_b) || (x_b <= x && x <= x_a));
			} else {
				FT y_a = A.y(), y_b = B.y(), y = M.y();
				return ((y_a <= y && y <= y_b) || (y_b <= y && y <= y_a));
			}

			/*if (std::find(CL.begin(), CL.end(), I_init) != CL.end() && std::find(CL.begin(), CL.end(), I_stop) != CL.end()) {
				// Arithmetic computation
				const CGAL_Point_2 &A = origin(), &B = end();
				CGAL_Point_2 M = CGAL::midpoint(V_1, V_2);

				FT x_a = A.x(), x_b = B.x(), x = M.x();
				if (x_a != x_b) {
					return ((x_a <= x && x <= x_b) || (x_b <= x && x <= x_a));
				} else {
					FT y_a = A.y(), y_b = B.y(), y = M.y();
					return ((y_a <= y && y <= y_b) || (y_b <= y && y <= y_a));
				}

			} else {
				return false;
			}*/
		}
	}


	std::list<Intersection_Line*>::const_iterator Polygon_Segment::crossed_lines_begin() const
	{
		return C_crossed.cbegin();
	}


	std::list<Intersection_Line*>::const_iterator Polygon_Segment::crossed_lines_end() const
	{
		return C_crossed.cend();
	}


	Polygon_Segment_R* Polygon_Segment::to_r()
	{
		return dynamic_cast<Polygon_Segment_R*>(this);
	}


	Polygon_Segment_S* Polygon_Segment::to_s()
	{
		return dynamic_cast<Polygon_Segment_S*>(this);
	}


	void Polygon_Segment::check() const
	{
		//std::cout << "TODO : check()" << std::endl;
	}



	void Polygon_Segment::who_is() const
	{
		std::cout << "*** Segment at fault : " << std::endl;
		std::cout << "id : " << id_object << " [" << id_plane << "] on line " << support->id_object << " which has " << support->planes.size() << " planes : " << std::endl;
		for (int i : support->planes) {
			std::cout << "# " << i << std::endl;
		}
		if (C_init.first == nullptr) {
			std::cout << "C_init = nullptr" << std::endl;
		} else {
			std::cout << "C_init = " << C_init.first->id_object << ", " << (C_init.second == 1 ? "PLUS" : "MINUS") << std::endl;
		}
		if (C_stop.first == nullptr) {
			std::cout << "C_stop = nullptr" << std::endl;
		} else {
			std::cout << "C_stop = " << C_stop.first->id_object << ", " << (C_stop.second == 1 ? "PLUS" : "MINUS") << std::endl;
		}
		std::cout << "|C_crossed| = " << C_crossed.size() << std::endl;
	}
}