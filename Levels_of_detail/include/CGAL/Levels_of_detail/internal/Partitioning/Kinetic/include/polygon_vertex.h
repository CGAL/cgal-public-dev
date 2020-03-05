#pragma once
#include "support_plane_object.h"
#include "support_plane_objects.h"
#include "defs_cgal.h"


namespace Skippy 
{
	typedef enum { 
		RUNNING_VERTEX, 
		STILL_VERTEX 
	} Polygon_Vertex_Type;


	class Polygon_Vertex : public Support_Plane_Object
	{
	protected:
		Polygon_Vertex(const int _id_plane);

	public:
		virtual ~Polygon_Vertex();

	public:
		void add(Polygon_Edge* e);

		void remove(Polygon_Edge* e);

		bool has_neighbor(Polygon_Vertex* v) const;

		Polygon_Edge* other_edge(Polygon_Edge* e) const;

		//virtual const CGAL_Point_2 & get_M() const = 0;
		
		virtual CGAL_Point_2 get_M() const = 0;

		virtual CGAL_Point_2 pt(const FT & t) const = 0;

		bool unconstrained() const;

		bool constrained() const;

		bool is_constrained_by(Intersection_Line* I) const;

		Constraint get_constraint() const;

		Sign sign_of_constraint(Intersection_Line* I) const;

		Polygon_Vertex_R* to_r();

		Polygon_Vertex_S* to_s();

	protected:
		bool is_constrained_neighbor(const std::list<Intersection_Line*> & I) const;

	public:
		virtual void set_polygon(Polygon* _polygon) = 0;

		Polygon* get_polygon();
		
	public:
		Polygon_Vertex_Type type;
		Polygon_Edge* e1;
		Polygon_Edge* e2;

	protected:
		std::list<Constraint> constraints;
		Polygon* polygon;
	};



	class Polygon_Vertex_R : public Polygon_Vertex
	{
	public:
		Polygon_Vertex_R(const int _id_plane, const FT & _t_init, const CGAL_Point_2 & _M, const CGAL_Vector_2 & _dM, const int _K, Event_Flags flags);

		Polygon_Vertex_R(const int _id_plane, const FT & _t_init, const CGAL_Point_2 & _M, const CGAL_Vector_2 & _dM, const Constraint & C, Intersection_Line* I_discarded, Polygon_Vertex_R* v_n, const int _K, Event_Flags flags);

		Polygon_Vertex_R(const int _id_plane, const FT & _t_init, const CGAL_Point_2 & _M, const CGAL_Vector_2 & _dM, const Constraint & C, const std::list<Intersection_Line*> & I_discarded, Polygon_Vertex_R* v_n, const int _K, Event_Flags flags);

		Polygon_Vertex_R(Polygon_Vertex_R* v_ts, Event_Flags flags);

		~Polygon_Vertex_R();

	public:
		void set_polygon(Polygon* _polygon);

		void set_M(const CGAL_Point_2 & _M);
		
		CGAL_Point_2 get_M() const;

		// const CGAL_Point_2 & get_M() const;
		
		CGAL_Point_2 pt(const FT & t) const;

		void stop(const FT & t_stop);

		void transfer_segments(Polygon_Vertex_R* v_dest, const FT & t, const CGAL_Point_2 & V_t);

		void stop_segments(Intersection_Line* I, const int seed, const FT & t);

		void extend_segments(Intersection_Line* I) const;

		Polygon_Vertex_R* get_constrained_neighbor(Intersection_Line* I_0, const std::list<Intersection_Line*> & I) const;

		Polygon_Vertex_R* get_neighbor_intersecting_identical_line(Intersection_Line* I_0, const std::list<Intersection_Line*> & I, const FT & t) const;

		void indicate_line_initially_crossed_by_segments(Intersection_Line* I) const;

		void set_as_guided_segment(Polygon_Segment_R* s);

		void copy_crossed_lines(Polygon_Vertex_R* v_os) const;

		bool has_events() const;

		void set_tau();

		void schedule_events(Polygon_Vertex_R* v_ts);

		void schedule_events();

		void schedule_events(Intersection_Line* I_redirect, Polygon_Vertex_R* v_n);

		void schedule_events(const std::list<Intersection_Line*> & discarded, Polygon_Vertex_R* v_n);

		void reschedule_events();

		void reschedule_events(Polygon_Vertex_R* v_n);




		bool meets_constrained_neighbor_at_next_event(Polygon_Vertex_R* v_n, Intersection_Line* & I_n, std::list<Event_Vertex_Line*> & E) const;

		void get_upcoming_event(std::list<Event_Vertex_Line*> & E) const;

		void reschedule_events_after_neighbors(Polygon_Vertex_R* v_n, std::list<Event_Vertex_Line*> & E_VL);

		void reschedule_events_by_finding_lines(std::list<Event_Vertex_Line*> & E_VL);



		Event_Vertex_Line* get_event_for_line(const int I_object) const;

		void delete_event(Event_Vertex_Line* e_vl);

		void delete_events(const std::list<Event_Vertex_Line*> & E_VL);

		void decrement_queued_events(const int n);



		std::pair<FT, bool> get_intersection_time(Intersection_Line* I) const;

		static void set_paired_vertices(Polygon_Vertex_R* v_ts, Polygon_Vertex_R* v_os);

		static bool are_paired_vertices(Polygon_Vertex_R* v_ts, Polygon_Vertex_R* v_os);

		void set_paired_vertex(const Companion & _paired_vertex);

		void reset_paired_vertex();

		bool has_paired_vertex() const;

		bool is_independent() const;

		Polygon_Vertex_R* get_paired_vertex() const;

	protected:
		Event_Vertex_Line* make_event(Intersection_Line* I) const;

		Event_Vertex_Line* make_event(Intersection_Line* I, const FT & t) const;

		void copy_events(Polygon_Vertex_R* v_ts);

	public:
		void deschedule_events();

	public:
		Polygon_Vertex* to_base();

	public:
		static bool constrained_vertices_meet(const std::list<Event_Vertex_Line*> & E_1, const std::list<Event_Vertex_Line*> & E_2,
			Polygon_Vertex_R* v1, Polygon_Vertex_R* v2);

	protected:
		CGAL_Point_2 M;
		
	public:
		CGAL_Vector_2 dM;
		FT t_init;
		int K;

	protected:
		std::map<int, Event_Vertex_Line*> lines_to_events;
		int queued_events;

	protected:
		int k_th_interval;
		FT tau_inf;
		FT tau_sup;

	public:
		int crossed_lines;

	protected:
		std::list<Polygon_Segment_R*> guided_segments;
		Companion paired_vertex;
	};


	class Polygon_Vertex_S : public Polygon_Vertex
	{
	public:
		Polygon_Vertex_S(const int _id_plane, const Constraint & C_1, const Constraint & C_2);

		~Polygon_Vertex_S();

	public:
		void set_polygon(Polygon* _polygon);

		CGAL_Point_2 get_M() const;

		// const CGAL_Point_2 & get_M() const;
	
		CGAL_Point_2 pt(const FT & t) const;

		Constraint get_second_constraint() const;

		Constraint get_other_constraint(const Constraint & C) const;

		bool represents_same_intersection(Polygon_Vertex_S* v1, Polygon_Vertex_S* v2) const;

		bool represents_same_intersection(Polygon_Vertex_S* v) const;

		Polygon_Vertex* to_base();
	};
}