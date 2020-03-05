#pragma once
#include "support_plane_object.h"
#include "support_plane_objects.h"
#include "defs_cgal.h"


namespace Skippy
{
	class Segment : public Support_Plane_Object {
	protected:
		Segment(const int _id_plane, const Constraint & C_support);

	public:
		virtual ~Segment();

		static bool closed_segment_includes(const CGAL_Point_2 & M, const CGAL_Point_2 & A, const CGAL_Point_2 & B);

		static bool half_closed_segment_includes(const CGAL_Point_2 & M, const CGAL_Point_2 & A, const CGAL_Point_2 & B);

	public:
		Intersection_Line* support;
	};


	class Planar_Segment : public Segment {
	public:
		Planar_Segment(const int _id_plane, const Constraint & C_support, CGAL_Point_2 & _A, CGAL_Point_2 & _B);

		~Planar_Segment();

		bool checks_if_belongs(const CGAL_Point_2 & V) const;

	public:
		CGAL_Point_2 A;
		CGAL_Point_2 B;
	};

	
	class Polygon_Segment : public Segment
	{
	protected:
		Polygon_Segment(const int _id_plane, const int _seed, const Constraint & C_support);

	public:
		virtual ~Polygon_Segment();

		//virtual const CGAL_Point_2 & origin() const = 0;

		//const CGAL_Point_2 & end() const;

		virtual CGAL_Point_2 origin() const = 0;

		CGAL_Point_2 end() const;

		void insert_as_crossed_line(Intersection_Line* I);

		void insert_as_crossed_lines(const std::list<Intersection_Line*> & L);

		virtual bool includes_point_on_support_line(const CGAL_Point_2 & M, const FT & t) const = 0;

		virtual bool includes_point_at_intersection(const CGAL_Point_2 & M, const Constraint & C, const FT & t) const = 0;

		bool includes_point_at_intersection(const CGAL_Point_2 & M, const std::list<Constraint> & C_limits, const FT & t) const;

		bool includes_edge(const CGAL_Point_2 & V_1, const CGAL_Point_2 & V_2, const Constraint & C_1, const Constraint & C_2, const std::list<Intersection_Line*> & CL) const;

		std::list<Intersection_Line*>::const_iterator crossed_lines_begin() const;

		std::list<Intersection_Line*>::const_iterator crossed_lines_end() const;

		virtual void check() const;

		virtual void who_is() const;

		Polygon_Segment_R* to_r();

		Polygon_Segment_S* to_s();

	public:
		int seed;
		Constraint C_init;
		Constraint C_stop;
		std::list<Intersection_Line*> C_crossed;
		std::list<Segment*>::iterator iterator;
	};



	class Polygon_Segment_R : public Polygon_Segment
	{
	protected:
		Polygon_Segment_R(const int _id_plane, const int _seed, const Constraint & C_support);
	public:
		Polygon_Segment_R(const int _id_plane, const int _seed, const FT & _t_init, const Constraint & C_support, Polygon_Vertex_R* & v, const CGAL_Point_2 & _O, const CGAL_Point_2 & _A, const CGAL_Vector_2 & _dA);

		Polygon_Segment_R(const int _id_plane, const int _seed, const FT & _t_init, const Constraint & C_support, Polygon_Vertex_R* & v, const CGAL_Point_2 & _A, const CGAL_Vector_2 & _dA);

		Polygon_Segment_R(const int _id_plane, const int _seed, const FT & _t_init, const Constraint & _C_init, const Constraint & C_support, Polygon_Vertex_R* & v, const CGAL_Point_2 & _A, const CGAL_Vector_2 & _dA);

		Polygon_Segment_R(const int _id_plane, const int _seed, const FT & _t_init, const Constraint & _C_init, const Constraint & C_support, Polygon_Vertex_R* & v, const CGAL_Point_2 & _O, const CGAL_Point_2 & _A, const CGAL_Vector_2 & _dA);

		~Polygon_Segment_R();

	public:
		Sign get_sign() const;

		// const CGAL_Point_2 & origin() const;
		
		CGAL_Point_2 origin() const;

		bool stopped() const;

		CGAL_Point_2 pt(const FT & t) const;

		void update_translation(const FT & t, const CGAL_Point_2 & A_t, const CGAL_Vector_2 & dA_t);

		void stop(const Constraint C, const FT & t);

		static void set_as_opposite_bidirectional_segments(Polygon_Segment_R* s_1, Polygon_Segment_R* s_2);

		bool exists_opposite_segment_without_initinal_constraint() const;

		Polygon_Segment_R* get_opposite() const;

		void set_pseudo_init_constraint(const CGAL_Point_2 & A, std::list<Intersection_Line*> & L);

		void set_pseudo_init_constraint(const Constraint C_pseudo_init);

		Constraint get_pseudo_init_constraint() const;

		bool includes_point_on_support_line(const CGAL_Point_2 & M, const FT & t) const;

		bool includes_point_at_intersection(const CGAL_Point_2 & M, const Constraint & C, const FT & t) const;

		void who_is() const;

		const Segment_Translation* get_current_translation() const;

	protected:
		Sign sign;
		CGAL_Point_2* A;
		Segment_Translation* Tr;
		std::list<Segment_Translation*> Tr_previous;
		Polygon_Segment_R* opposite;
	};



	class Polygon_Segment_S : public Polygon_Segment
	{
	public:
		Polygon_Segment_S(const int _id_plane, const int _seed, const Constraint & _C_init, const Constraint & _C_support, const Constraint & _C_stop, const CGAL_Point_2 & _A, const CGAL_Point_2 & _B);

		Polygon_Segment_S(const int _id_plane, const int _seed, const Constraint & _C_init, const Constraint & _C_support, const std::list<Intersection_Line*> & C_crossed, const Constraint & _C_stop, const bool merge_at_boundaries);

		Polygon_Segment_S(const int _id_plane, const int _seed, const Constraint & _C_init, const Constraint & _C_support, const std::list<Intersection_Line*> & C_crossed_1, const std::list<Intersection_Line*> & C_crossed_2, const Constraint & _C_stop);

		~Polygon_Segment_S();

	public:
		bool includes_point_on_support_line(const CGAL_Point_2 & M, const FT & t) const;

		bool includes_point_at_intersection(const CGAL_Point_2 & M, const Constraint & C, const FT & t) const;

		Constraint get_other_constraint(const Constraint & C) const;

	protected:
		//const CGAL_Point_2 & origin() const;

		CGAL_Point_2 origin() const;

		void init_by_merging(const Constraint & C_init, const Constraint & C_support, std::list<Intersection_Line*> C_crossed, const Constraint & C_stop);

		Polygon_Segment_S* get_adjacent_segment(const std::list<Segment*> & L, const Constraint & C) const;
	};
}