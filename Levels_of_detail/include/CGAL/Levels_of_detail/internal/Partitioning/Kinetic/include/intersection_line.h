#pragma once
#include "support_plane_object.h"
#include "support_plane_objects.h"
#include "defs_cgal.h"


namespace Skippy
{
	class Intersection_Line : public Support_Plane_Object 
	{
	public:
		Intersection_Line(const int _id_plane, const CGAL_Line_2 & _line, int intersected);

		~Intersection_Line();

		void clear_segments();

		void combine_signs(std::list<Segment*> & S) const;

		void mark_as_intersected(int intersected);

		bool intersects(const int id_plane) const;

		void set_inside(const bool _is_inside);

		Sign sign(const CGAL_Point_2 & pt) const;

		Sign sign(Polygon_Vertex* v, const FT & t) const;


		bool includes(const CGAL_Point_2 & M) const;

		bool is_parallel(Intersection_Line* I) const;



		bool exist_segments_including_point_outside_intersections(const CGAL_Point_2 & V_t, const FT & t) const;

		bool exist_segments_including_point_at_intersection(const CGAL_Point_2 & V_t, const Constraint & C, const FT & t) const;

		bool exist_segments_including_point_at_intersection(const CGAL_Point_2 & V_t, const std::list<Constraint> & C_limits, const FT & t) const;

		// bool exists_segment_adjacent_to_edge(Polygon_Edge* e) const;

		bool exists_segment_adjacent_to_edge(const Constraint & C_prev, const Constraint & C, const Constraint & C_next) const;

		static bool represents_same_intersection(const std::pair<Constraint, Constraint> & v1, const std::pair<Constraint, Constraint> & v2);


		const FT & a() const;

		const FT & b() const;

		const FT & c() const;

		void reject_any_further_segment();

	private:
		void reject_any_further_segment(std::list<Segment*> & S);

		bool exists_segment_including_point_outside_intersections(const CGAL_Point_2 & V_t, const FT & t, const std::list<Segment*> & S) const;

		bool exists_segment_including_point_at_intersection(const CGAL_Point_2 & V_t, const Constraint & C, const FT & t, const std::list<Segment*> & S) const;

		bool exists_segment_including_point_at_intersection(const CGAL_Point_2 & V_t, const std::list<Constraint> & C_limits, const FT & t, const std::list<Segment*> & S) const;

	public:
		bool is_border;
		bool is_inside;
		bool reject_segments;

		const CGAL_Line_2 line;
		const FT _a;
		const FT _b;
		const FT _c;

		double hint_a;
		double hint_b;
		double hint_c;

		std::list<int> planes;
		std::list<Segment*> segments_minus;
		std::list<Segment*> segments_plus;
	};
}