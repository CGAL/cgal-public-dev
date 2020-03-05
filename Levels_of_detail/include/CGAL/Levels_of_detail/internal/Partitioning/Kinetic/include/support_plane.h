#pragma once
#include "defs.h"
#include "event.h"
#include "support_plane_objects.h"
#include "polygon.h"
#include "polygon_set.h"



namespace Skippy 
{
	typedef std::tuple<int, int, int> Triplet;

	class Support_Plane
	{
	public:
		Support_Plane(const int plane_id, const CGAL_Plane & plane_equation, const CGAL_Color & plane_color);

		~Support_Plane();

	public:
		static void construct_bounding_polygon_of_support_plane(const CGAL_Point_3 & pt_min,
			const CGAL_Point_3 & pt_max,
			const std::vector<CGAL_Point_3> & box_corners,
			const std::vector<std::pair<int, int> > & box_edges,
			const CGAL_Plane & P,
			std::list<CGAL_Point_3> & bounding_polygon,
			std::vector<std::list<int> > & bounding_facets);

		static void construct_bounding_polygon_of_support_plane(const CGAL_Point_3 & pt_min,
			const CGAL_Point_3 & pt_max,
			const std::vector<CGAL_Point_3> & box_corners,
			const std::vector<std::pair<int, int> > & box_edges,
			const CGAL_Plane & P,
			std::list<CGAL_Point_3> & bounding_polygon);

	protected:
		static void find_next_object_colliding_plane(const CGAL_Point_3 & pt_min,
			const CGAL_Point_3 & pt_max,
			const std::vector<CGAL_Point_3> & box_corners,
			const std::vector<std::pair<int, int> > & box_edges,
			const CGAL_Plane & P,
			std::pair<bool, int> & next_object,
			CGAL_Point_3 & M);

		static bool find_next_object_colliding_plane(const CGAL_Point_3 & pt_min,
			const CGAL_Point_3 & pt_max,
			const std::vector<CGAL_Point_3> & box_corners,
			const std::vector<std::pair<int, int> > & box_edges,
			const CGAL_Plane & P,
			const std::vector<int> & V_ind,
			const std::vector<int> & E_ind,
			const std::pair<bool, int> & previous_object,
			std::pair<bool, int> & next_object,
			CGAL_Point_3 & M);

		static std::vector<int> box_edges_to_facets(const int i);

		static std::vector<int> box_vertices_to_facets(const int i);

		static std::vector<int> box_facets_to_vertices(const int i);

		static std::vector<int> box_facets_to_edges(const int i);

		static CGAL_Point_3 intersection_plane_and_edge_of_bounding_box(const CGAL_Point_3 & pt_min, const CGAL_Point_3 & pt_max, const CGAL_Plane & P, const int i);

	public:
		CGAL_Point_2 project(const CGAL_Point_3 & P) const;

		CGAL_Point_3 backproject(const CGAL_Point_2 & P) const;

		void init_intersection_lines(const std::map<int, CGAL_Line_3> & lines_3d);

		Intersection_Line* get_line_by_identifier(const int id_object) const;

		Intersection_Line* get_line_for_plane(const int id_plane) const;

		void get_concurrent_lines(const Intersection_Line* I_1, const Intersection_Line* I_2, std::list<Intersection_Line*> & L);

		void get_indices_of_intersecting_planes(const Intersection_Line* I_1, const Intersection_Line* I_2, std::list<int> & P);

		void insert_triplets_of_concurrent_lines(const std::vector<Intersection_Line*> & I);


		void init_landmarks();

		void delete_landmarks();

		const CGAL_Point_2 & get_landmark(Intersection_Line* I_1, Intersection_Line* I_2);

		CGAL_Point_2 get_intersection_point(Intersection_Line* I_1, Intersection_Line* I_2) const;


	public:
		bool assert_inside(Polygon_Vertex_R* v, const FT & t);

		void search_lines_in_neighborhood(Polygon_Vertex_R* v, const FT & t_1, const FT & t_2, std::list<std::tuple<Intersection_Line*, bool, FT> > & L);

	public:
		void init_bounding_polygon(const std::list<CGAL_Point_3> & bounding_polygon_3d, const std::vector<std::list<int> > & bounding_facets);

		void init_polygon_set();

		void init_polygon(const std::vector<CGAL_Point_3> & polygon);

		void init_polygon(const std::vector<CGAL_Point_2> & polygon);

		void init_polygon(const CGAL_Point_2 & initial_barycenter,
			const std::vector<std::pair<CGAL_Point_2, CGAL_Vector_2> > & initial_directions);

		void project_and_decompose_with_respect_to_slices(const std::vector<CGAL_Point_3> & polygon,
			const int gx, const int gy, const int gz,
			CGAL_Point_2 & initial_barycenter,
			std::vector<std::pair<CGAL_Point_2, CGAL_Vector_2> > & initial_directions,
			std::set<std::tuple<int, int, int> > & subdivs);

		std::tuple<int, int, int> get_subvolume(const CGAL_Point_3 & M,
			const int gx, const int gy, const int gz,
			const std::vector<FT> & grad_x, const std::vector<FT> & grad_y, const std::vector<FT> & grad_z);

		void clear_polygons();

	protected:
		Polygon_Tree* initialize_tree_of_polygons(const std::vector<CGAL_Point_2> & polygon);

		Polygon_Tree* initialize_tree_of_polygons(const CGAL_Point_2 & initial_barycenter,
			const std::vector<std::pair<CGAL_Point_2, CGAL_Vector_2> > & initial_directions);

		void initialize_cells(Polygon_Tree* T);

		void update_mutual_intersections(Polygon_Vertex* v, const int current_line, 
			std::map<int, std::map<int, CGAL_Point_2> > & mutual_intersections) const;

		void update_mutual_intersections(const std::vector<Polygon_Vertex*> & ranked_points,
			const int current_line,
			std::map<int, std::map<int, CGAL_Point_2> > & mutual_intersections) const;

		void project_polygon(const std::vector<CGAL_Point_3> & P_0, std::vector<CGAL_Point_2> & P) const;

		void set_initial_propagation_directions(const std::vector<CGAL_Point_2> & R, CGAL_Point_2 & O, std::vector<std::pair<CGAL_Point_2, CGAL_Vector_2> > & D) const;

		void set_element_in_quadruplet(Polygon_Vertex_R* & v1_os, Polygon_Vertex_R* & v2_os, Polygon_Vertex_R* & v1_ts, Polygon_Vertex_R* & v2_ts, Polygon_Vertex_R* v) const;

		int locate_intersection(const CGAL_Point_2 & M, const CGAL_Point_2 & M_0, const CGAL_Point_2 & M_1, const CGAL_Point_2 & P) const;

	protected:
		// Segment initializers

		// Case 1A. Typicalled called for an event of type C : a constrained vertex reaches the intersection of two lines.

		void init_1_uplets_of_unidirectional_segments_of_null_length(Intersection_Line* I, const CGAL_Point_2 & V_t,
			Polygon_Vertex_R* v, const Constraint C_init, const FT & t);

		void init_2_uplets_of_unidirectional_segments_of_null_length(Intersection_Line* I, const CGAL_Point_2 & V_t,
			Polygon_Vertex_R* v_ts, Polygon_Vertex_R* v_os, const Constraint C_init, const FT & t);

		// Case 1B. Typically called at the initialization or when an edge reaches a line.

		void init_1_uplets_of_unidirectional_segments(Intersection_Line* I, const CGAL_Point_2 & V_c, const CGAL_Point_2 & V_t, const CGAL_Point_3 & W_t,
			Polygon_Vertex_R* v, const Constraint C_init, const FT & t);

		void init_2_uplets_of_unidirectional_segments(Intersection_Line* I, const CGAL_Point_2 & V_c, const CGAL_Point_2 & V_t, const CGAL_Point_3 & W_t,
			Polygon_Vertex_R* v_ts, Polygon_Vertex_R* v_os, const Constraint C_init, const FT & t);

		// Case 2A. Typically called for an event of type A : a unconstrained vertex reaches a line.

		void init_2_uplets_of_bidirectional_segments_of_null_length(Intersection_Line* I, const CGAL_Point_2 & V_t, const CGAL_Point_3 & W_t,
			Polygon_Vertex_R* v1, Polygon_Vertex_R* v2, const FT & t);

		void init_4_uplets_of_bidirectional_segments_of_null_length(Intersection_Line* I, const CGAL_Point_2 & V_t, const CGAL_Point_3 & W_t,
			Polygon_Vertex_R* v1_1, Polygon_Vertex_R* v1_2, Polygon_Vertex_R* v2_1, Polygon_Vertex_R* v2_2, const FT & t);

		// Case 2B. Typically called at the initialization or when an edge reaches a line.

		void init_2_uplets_of_bidirectional_segments(Intersection_Line* I, const CGAL_Point_2 & V1_t, const CGAL_Point_3 & W1_t, const CGAL_Point_2 & V2_t, const CGAL_Point_3 & W2_t,
			Polygon_Vertex_R* v1, Polygon_Vertex_R* v2, const FT & t);

		void init_4_uplets_of_bidirectional_segments(Intersection_Line* I, const CGAL_Point_2 & V1_t, const CGAL_Point_3 & W1_t, const CGAL_Point_2 & V2_t, const CGAL_Point_3 & W2_t,
			Polygon_Vertex_R* v1_1, Polygon_Vertex_R* v1_2, Polygon_Vertex_R* v2_1, Polygon_Vertex_R* v2_2, const FT & t);

		// Case 3. Typically called when an edge reaches a line.

		void init_1_uplets_of_constant_segments(Intersection_Line* I, 
			Polygon_Vertex* v1_ts, Polygon_Vertex* v2_ts, const Constraint C_1, const Constraint C_2, const FT & t);

		void init_1_uplets_of_constant_segments(Intersection_Line* I, 
			Polygon_Vertex* v1_ts, Polygon_Vertex* v2_ts, const Constraint C_1, const Constraint C_2, const std::list<Intersection_Line*> & C_crossed, const FT & t);

		void init_2_uplets_of_constant_segments(Intersection_Line* I,
			Polygon_Vertex* v1_ts, Polygon_Vertex* v1_os, Polygon_Vertex* v2_ts, Polygon_Vertex* v2_os,
			const Constraint C_1, const Constraint C_2, const FT & t);

		void init_2_uplets_of_constant_segments(Intersection_Line* I,
			Polygon_Vertex* v1_ts, Polygon_Vertex* v1_os, Polygon_Vertex* v2_ts, Polygon_Vertex* v2_os,
			const Constraint C_1, const Constraint C_2, const std::list<Intersection_Line*> & C_crossed, const FT & t);

		void build_n_uplets_of_segments_as_edge_collides_line(Intersection_Line* I, const Constraint & C_1, const Constraint & C_2,
			const CGAL_Point_2 & V1_t, const CGAL_Point_3 & W1_t, const CGAL_Point_2 & V2_t, const CGAL_Point_3 & W2_t,
			Polygon_Vertex* v1_ts, Polygon_Vertex* v1_os, Polygon_Vertex* v2_ts, Polygon_Vertex* v2_os, const FT & t);

	public:
		void init_schedule();

		void process_event(Event_Vertex_Line* e);

		void append(Event_Vertex_Line* e_vl, std::string type);

		void real_time_check(const FT & t);

		void check();

	protected:
		void process_event_1(const std::list<Event_Vertex_Line*> & E, const std::list<Polygon_Vertex_R*> & V);

		void process_event_2(const std::list<Event_Vertex_Line*> & E, const std::list<Polygon_Vertex_R*> & V);
		
		void process_event_3(const std::list<Event_Vertex_Line*> & E, const std::list<Polygon_Vertex_R*> & V);

		void print_events(const std::list<Event_Vertex_Line*> & E, const std::string & type) const;

		void get_constraint_and_location_at_time_t(Polygon_Vertex_R* v, const int intersected, const FT & t_intersectant, Constraint & C_res, CGAL_Point_2 & V_res);

		void get_constraint_and_location_at_time_t(Polygon_Vertex_R* v, Intersection_Line* I_intersect, const FT & t_intersectant, Constraint & C_res, CGAL_Point_2 & V_res);

		void get_references_to_simulatenous_events_for_this_polygon(Polygon_Vertex_R* v, const FT & t, std::list<Event_Vertex_Line*> & E);

		int get_number_of_vertices_simultaneously_intersecting_lines(Polygon_Vertex_R* v, Event_Vertex_Line* e, const std::list<Event_Vertex_Line*> & E, std::list<Polygon_Vertex_R*> & V) const;

		int loop_on_neighbor_vertices_with_simultaneous_events(Polygon_Vertex_R* v_init, Polygon_Edge* e_init, Event_Vertex_Line* e, const std::list<Event_Vertex_Line*> & E, std::list<Polygon_Vertex_R*> & V) const;

		void pop_fitered_simultaneous_events_from_queue(std::list<Event_Vertex_Line*> & E, const std::list<Polygon_Vertex_R*> & V) const;

		void split_simultaneous_events(const std::list<Event_Vertex_Line*> & E, Polygon_Vertex_R* v1, Polygon_Vertex_R* v2, std::list<Event_Vertex_Line*> & E_1, std::list<Event_Vertex_Line*> & E_2) const;

	protected:
		void get_simultaneous_events(Polygon_Vertex_R* v, const FT & t, const CGAL_Point_2 & V_t, std::list<Event_Vertex_Line*> & E, Event_Vertex_Line* e, std::list<Intersection_Line*> & I);

		void get_simultaneous_events_as_edge_intersects_line(Polygon_Vertex_R* v, const FT & t, const CGAL_Point_2 & V_t, Intersection_Line* I_L,
			std::list<Event_Vertex_Line*> & E, Event_Vertex_Line* & e_vl);

		// Handle case A
		void unconstrained_vertex_intersects_line(const std::list<Event_Vertex_Line*> & E, Polygon_Vertex_R* v, const CGAL_Point_2 & V_t);

		// Handle case B
		void unconstrained_vertex_intersects_line_kth_time(const std::list<Event_Vertex_Line*> & E, Polygon_Vertex_R* v, const CGAL_Point_2 & V_t);

		// Handle case C1
		void constrained_vertices_intersect(const std::list<Event_Vertex_Line*> & E_1, const std::list<Event_Vertex_Line*> & E_2, Polygon_Vertex_R* v1, Polygon_Vertex_R* v2);

		// Handle case C2
		void constrained_vertex_intersects_line(const std::list<Event_Vertex_Line*> & E, Polygon_Vertex_R* v, const CGAL_Point_2 & V_t);

		void constrained_vertex_intersects_single_line(Event_Vertex_Line* e, Polygon_Vertex_R* v, const CGAL_Point_2 & V_t);

		void constrained_vertex_intersects_multiple_lines(const std::list<Event_Vertex_Line*> & E, Polygon_Vertex_R* v, const CGAL_Point_2 & V_t);

		// Helpers case C2
		void get_sorted_lines(Polygon_Vertex_R* v, const CGAL_Point_2 & V_t, const std::list<Event_Vertex_Line*> & E, std::vector<Intersection_Line*> & I);

		// void get_sorted_lines(Polygon_Vertex_R* v, const std::list<Event_Vertex_Line*> & E, std::vector<Intersection_Line*> & I);

		void get_sorted_constraints(Polygon_Vertex_R* v, const std::vector<Intersection_Line*> & I, std::vector<std::pair<Constraint, Constraint> > & C);

		void get_definitions_of_constrained_vertices(Polygon_Edge* e, const FT t, const std::vector<Intersection_Line*> & I, std::vector<std::pair<CGAL_Point_2, CGAL_Vector_2> > & D);

		bool constrained_vertex_intersects_multiple_lines_ts(Polygon_Vertex_R* v, const FT t, const CGAL_Point_2 & V_t, const CGAL_Point_3 & W_t,
			const std::vector<Intersection_Line*> & I, const std::vector<std::pair<Constraint, Constraint> > & C);

		void constrained_vertex_intersects_multiple_lines_os(Polygon_Vertex_R* v, const FT t, const CGAL_Point_2 & V_t, const CGAL_Point_3 & W_t,
			const std::vector<Intersection_Line*> & I, const std::vector<std::pair<Constraint, Constraint> > & C);

		void get_stop_constraints_at_multiple_intersections(const std::vector<std::pair<Constraint, Constraint> > & C, const int k, const bool H_ts, std::list<Constraint> & C_k);

		// Handle case D1
		void edge_intersects_line(const std::list<Event_Vertex_Line*> & E_1, const std::list<Event_Vertex_Line*> & E_2, Polygon_Vertex_R* v1, Polygon_Vertex_R* v2);

		// Helpers case D1
		void edge_propagates_frontally(Intersection_Line* I, Polygon_Vertex_R* v1, Polygon_Vertex_R* v2, const FT & t,
			const CGAL_Point_2 & V_1, const CGAL_Point_3 & W_1, const CGAL_Point_2 & V_2, const CGAL_Point_3 & W_2, const bool propagate_frontally);

		void edge_propagates_laterally(Intersection_Line* I, Polygon_Vertex_R* v, const FT & t, const CGAL_Point_2 & V_t);

		bool edge_is_propagating_frontally(Intersection_Line* I, Polygon_Vertex_R* v1, Polygon_Vertex_R* v2, const FT & t, bool lookup, 
			const Constraint & C_1, const Constraint & C_2,
			const CGAL_Point_2 & V_1, const CGAL_Point_3 & W_1, const CGAL_Point_2 & V_2, const CGAL_Point_3 & W_2);

		bool edge_is_propagating_laterally(Intersection_Line* I, Intersection_Line* I_l, Polygon_Vertex_R* v, const FT & t, const CGAL_Point_2 & V_t, const CGAL_Point_3 & W_t);

		// Handle case D2
		void two_edges_intersect_two_lines(const std::list<Event_Vertex_Line*> & E_1, const std::list<Event_Vertex_Line*> & E, const std::list<Event_Vertex_Line*> & E_2,
			Polygon_Vertex_R* v1, Polygon_Vertex_R* v, Polygon_Vertex_R* v2);

		// Evaluate stop conditions

		bool propagation_continues_outside_intersections(Intersection_Line* I, Polygon_Vertex_R* v, bool lookup, const CGAL_Point_2 & V_t,
			const CGAL_Point_3 & W_t, const FT & t, const std::vector<bool> & S, const int seed) const;

		bool propagation_continues_at_intersection(Intersection_Line* I, Polygon_Vertex_R* v, bool lookup, const CGAL_Point_2 & V_t,
			const CGAL_Point_3 & W_t, const Constraint & C, const FT & t, const std::vector<bool> & S, const int seed) const;

		bool propagation_continues_at_intersection(Intersection_Line* I, Polygon_Vertex_R* v, bool lookup, const CGAL_Point_2 & V_t,
			const CGAL_Point_3 & W_t, const std::list<Constraint> & C_limits, const FT & t, const std::vector<bool> & S, const int seed) const;

		bool K_based_stopping_condition(Polygon_Vertex_R* v, bool lookup) const;

		bool density_based_stopping_condition(Polygon_Vertex_R* v, const CGAL_Point_2 & V_t, const CGAL_Point_3 & W_t) const;

		bool density_based_stopping_condition_box(Polygon_Vertex_R* v, const CGAL_Point_2 & V_t, const CGAL_Point_3 & W_t) const;

		bool density_based_stopping_condition_cone(Polygon_Vertex_R* v, const CGAL_Point_2 & V_t, const CGAL_Point_3 & W_t) const;

		FT project_to_parallel_plane(const FT & D, const CGAL_Vector_3 & n) const;

	public:
		void get_polygon_description(std::list<std::list<CGAL_Point_3> > & P, std::list<CGAL_Color> & C, const double t);

		void get_polygon_description(Polygon_Vertex_Octree* V, std::list<std::list<int> > & P, const FT & t);

		void group_final_polygons();

	public:
		void draw(FT t, double dt, int size, double x_step, double y_step, int font_size);

	protected:
		void init_colormap(std::map<Intersection_Line *, CGAL_Color> & colors);

		void init_bounding_square(double & x_min, double & y_min, double & x_max, double & y_max, double & dx, double & dy);

		void init_markers(std::ostream & os, std::map<Intersection_Line *, CGAL_Color> & colors);

		void draw_grid(std::ostream & os, int rows, int cols, double x_min, double y_min, double dx, double dy, int font_size, double x_s, double y_s);

		void draw_bounding_polygon(std::ostream & os, int rows, int cols, double x_min, double y_min, double dx, double dy, int font_size, std::map<Intersection_Line*, CGAL_Color> & colors);

		void draw_intersection_lines(std::ostream & os, int rows, int cols, double x_min, double y_min, double dx, double dy, int font_size, std::map<Intersection_Line*, CGAL_Color> & colors);

		void draw_polygon_segments(std::ostream & os, int rows, int cols, double x_min, double y_min, double dx, double dy, int font_size, std::map<Intersection_Line*, CGAL_Color> & colors, FT t, double dt);

		void draw_polygon_segments(std::ostream & os, int rows, int cols, double x_min, double y_min, double dx, double dy, int font_size, std::map<Intersection_Line*, CGAL_Color> & colors, FT t, double dt, const std::list<Segment*> & segments);

		void draw_polygon(std::ostream & os, int rows, int cols, double x_min, double y_min, double dx, double dy,
			int font_size, std::map<Intersection_Line*, CGAL_Color> & colors, FT t, Polygon_Set* T, double dt);

	public:
		const int id;
		const int real_id;
		const CGAL_Plane plane;
		CGAL_Color color;

		double x_min, x_max, y_min, y_max;

		std::vector<Intersection_Line*> lines;
		std::map<int, Planar_Segment *> borders;
		std::map<int, Polygon_Vertex_R*> vertices_r;

		int running_vertices;
		CGAL_Point_2*** landmarks;
		
		std::vector<Polygon_Directions*> polygon_directions;
		Polygon_Set* polygon_set;

		std::list<Polygon_Group*> groups;
		std::vector<Intersection_Line*> lines_inside;

	protected:
		typedef struct _Triplet_Comparator {
			bool operator() (const Triplet & TL, const Triplet & TR) const {
				int tl_0 = std::get<0>(TL), tl_1 = std::get<1>(TL), tl_2 = std::get<2>(TL);
				int tr_0 = std::get<0>(TR), tr_1 = std::get<1>(TR), tr_2 = std::get<2>(TR);
				return (tl_0 < tr_0) || (tl_0 == tr_0 && tl_1 < tr_1) || (tl_0 == tr_0 && tl_1 == tr_1 && tl_2 < tr_2);
			}
		} Triplet_Comparator;

		std::set<Triplet, Triplet_Comparator> concurrences;
	};
}