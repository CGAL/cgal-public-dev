#pragma once
#include "support_plane_objects.h"
#include "polygon_vertex.h"

namespace Skippy 
{
	class Polygon_Node;

	class Polygon_Directions
	{
	public:
		Polygon_Directions(const CGAL_Point_2 & _O,
			const std::vector<std::pair<CGAL_Point_2, CGAL_Vector_2> > & _vertices);
			// const std::vector<Intersection_Line*> & _reference_lines);

		~Polygon_Directions();

		size_t size() const;

		std::pair<CGAL_Point_2, CGAL_Vector_2> & get_vertex(int i);

		// Intersection_Line* get_reference_line(int i);

		const CGAL_Point_2 & get_barycenter() const;

	protected:
		const CGAL_Point_2 O;
		std::vector<std::pair<CGAL_Point_2, CGAL_Vector_2> > vertices;
		// std::vector<Intersection_Line*> reference_lines;
	};


	class Polygon_Tree {
	public:
		Polygon_Tree();

		Polygon_Tree(const int id_plane, const int seed, Polygon_Directions* D);

		Polygon_Tree(const int id_plane, const int seed, Polygon_Directions* D, const CGAL_Vector_2 & OM, const FT & t);

		Polygon_Tree(Polygon_Tree* _parent, Polygon* _polygon);

		~Polygon_Tree();

		bool split(Intersection_Line* I, const FT & t, const int K, bool pair_constrained_vertices);

		bool split(Intersection_Line* I, std::list<Polygon_Vertex *> & V, const FT & t, const int K, bool pair_constrained_vertices);

		void turn_into_node(Intersection_Line* I, Polygon* subpolygon_plus, Polygon* subpolygon_minus, bool destroy_polygon);

		//void get_polygon_description(std::list<std::list<CGAL_Point_3> > & P, std::list<CGAL_Color> & C, const double t);

		void get_polygons(std::list<Polygon*> & P);

		Polygon* remove_reference_to_polygon();

		void remove_reference_to_polygons();

		static void sort_intersection_points(Intersection_Line* I, const std::list<Polygon_Vertex*> & V, const FT & t, std::vector<Polygon_Vertex*> & V_res);

	public:
		bool is_node;
		Intersection_Line* line;
		Polygon_Tree* parent;
		Polygon_Tree* subtree_plus;
		Polygon_Tree* subtree_minus;
		Polygon* polygon;
	};


	class Polygon {
	public:
		Polygon(const int _seed, const std::list<Polygon_Vertex*> & _vertices, const std::list<Polygon_Edge*> & _edges);

		Polygon(const int id_plane, const int _seed, Polygon_Directions* D);

		Polygon(const int id_plane, const int _seed, Polygon_Directions* D, const CGAL_Vector_2 & OM, const FT & t);

		~Polygon();

		inline void set_cell(Polygon_Node* _cell) { cell = _cell; }

		inline Polygon_Node* get_cell() { return cell; }

		void clear();

		bool intersection_exists(Intersection_Line* I, const FT & t, std::map<int, Sign> & vertices_signs,
			std::list<Polygon_Vertex *> & vertices_plus, std::list<Polygon_Edge *> & edges_plus, Polygon_Edge* & edge_plus_front, Polygon_Edge* & edge_plus_back,
			std::list<Polygon_Vertex *> & vertices_minus, std::list<Polygon_Edge *> & edges_minus, Polygon_Edge* & edge_minus_front, Polygon_Edge* & edge_minus_back,
			std::list<Polygon_Vertex *> & vertices_zero);

		void exhibit_sequence(Intersection_Line* I, Polygon_Vertex* v, Sign epsilon, std::map<int, Sign> & signs,
			std::list<Polygon_Vertex *> & sequence_vertices, std::list<Polygon_Edge *> & sequence_edges,
			Polygon_Edge* & excluded_edge_front, Polygon_Edge* & excluded_edge_back);

		void split_unconstrained_vertex(Intersection_Line* I, Polygon_Vertex_R* v, const FT & t, const CGAL_Point_2 & V_t, Polygon_Vertex_R* & v1, Polygon_Vertex_R* & v2, Constraint & C);

		void remove_vertices_colliding_on_intersection_line(Intersection_Line* I, Polygon_Vertex_R* v, const FT & t, const CGAL_Point_2 & V_t, Polygon_Vertex_R* & vc_ts, CGAL_Vector_2 & u_ref);

		Polygon_Vertex_R* get_tied_vertex(Intersection_Line* I, Polygon_Vertex_R* v, const FT & t, const CGAL_Point_2 & V_t);

		Polygon_Vertex_R* get_tied_vertex(Intersection_Line* I, Polygon_Vertex_R* v, const FT & t, const CGAL_Point_2 & V_t, const CGAL_Vector_2 & u);

		void insert_unconstrained_vertex(Intersection_Line* I, Polygon_Vertex_R* v, const FT & t, const CGAL_Point_2 & V_t, Polygon_Vertex_R* vc_ts, Polygon_Vertex_R* v_c);

		void remove_vertices_colliding_in_dead_end(Polygon_Vertex_R* v1, Polygon_Vertex_R* v2, const FT & t);

		void redirect_constrained_vertex(Intersection_Line* I, Polygon_Vertex_R* v, const FT & t, const CGAL_Point_2 & V_t,
			Polygon_Vertex_S* & v_de, Polygon_Vertex_R* & v_ts, Constraint & C_ts);

		void redirect_constrained_vertex(Intersection_Line* I, Polygon_Vertex_R* v, Polygon_Edge* e, const FT & t,
			const std::vector<std::pair<CGAL_Point_2, CGAL_Vector_2> > & D, const std::vector<std::pair<Constraint, Constraint> > & C,
			const std::list<Intersection_Line*> & I_discarded, Polygon_Vertex_R* & v_ts);

		void remove_intersectant_edge(Intersection_Line* I, Polygon_Vertex_R* v1, Polygon_Vertex_R* v2, const FT & t,
			const CGAL_Point_2 & V1_t, const CGAL_Point_2 & V2_t,
			Polygon_Edge* & e, Polygon_Vertex* & v1_ts, Polygon_Vertex* & v2_ts,
			Constraint & C_1, Constraint & C_2, Constraint & C_ts);

		void remove_intersectant_edges(Intersection_Line* I_1, Intersection_Line* I_2,
			Polygon_Vertex_R* v1, Polygon_Vertex_R* v, Polygon_Vertex_R* v2, const FT & t,
			const CGAL_Point_2 & V1_t, const CGAL_Point_2 & V_t, const CGAL_Point_2 & V2_t,
			Polygon_Vertex* & v1_ts, Polygon_Vertex_S* & v_de, Polygon_Vertex* & v2_ts);



		// bool is_adjacent_to(Polygon* P, Polygon_Edge* & e, Polygon_Edge* & p_e);

		CGAL_Point_2 get_barycenter(const FT & t) const;


		static Polygon* build_as_unconstrained_vertex_propagates(Polygon_Vertex_R* v, const FT & t,
			const int polygon_seed, Polygon_Vertex_R* v1_ts, Polygon_Vertex_R* v2_ts, const Constraint & C_ts,
			Polygon_Vertex_R* & v1_os, Polygon_Vertex_R* & v2_os);

		static Polygon* build_as_constrained_vertex_propagates(Polygon_Vertex_R* v, const FT & t,
			const int polygon_seed, Polygon_Vertex_S* vd_ts, Polygon_Vertex_R* v_ts,
			const Constraint & C_os, const Constraint & C_v_ts, Polygon_Vertex_R* & v_os);

		static Polygon* build_as_constrained_vertex_propagates(Polygon_Vertex_R* v, const FT & t,
			const int seed, const int k,
			const std::vector<std::pair<CGAL_Point_2, CGAL_Vector_2> > & D, const std::vector<std::pair<Constraint, Constraint> > & C,
			const std::list<Intersection_Line*> & I_discarded, Polygon_Vertex_R* & v_k, Polygon_Vertex_R* & v_l);

		static Polygon* build_as_edge_propagates(Polygon_Vertex_R* v1, Polygon_Vertex_R* v2, 
			Polygon_Vertex* v1_ts, Polygon_Vertex* v2_ts, Polygon_Edge* e, const FT & t,
			const int seed, const Constraint & C_ts, Polygon_Vertex* & v1_os, Polygon_Vertex* & v2_os);

		void get_active_constrained_vertices(Intersection_Line* I_1, Intersection_Line* I_2, Polygon_Vertex_R* & v_1, Polygon_Vertex_R* & v_2);

		void shift_vertices(Polygon_Vertex_R* v1, Polygon_Vertex_R* v2, const CGAL_Point_2 & M, const FT & t);

		static void shift_remaining_vertices_and_schedule_events(Polygon* P_1, Polygon* P_2, Intersection_Line* I_1, Intersection_Line* I_2, const CGAL_Point_2 & M, const FT & t);

		static void shift_remaining_vertices_and_schedule_events(const std::list<Polygon*> P, std::list<Intersection_Line*> I, const CGAL_Point_2 & M, const FT & t);

		void forget();

	private:
		template <typename T>
		inline void insert_object(T* object, std::list<T*> & container) {
			container.push_back(object);
			if (std::is_same<T*, Polygon_Vertex*>::value) {
				Polygon_Vertex* v = reinterpret_cast<Polygon_Vertex*>(object);
				v->set_polygon(this);
			}
		}

		template <typename T>
		inline void insert_object(T* object_1, T* object_2, std::list<T*> & container) {
			insert_object<T>(object_1, container);
			insert_object<T>(object_2, container);
		}

		template <typename T>
		inline void insert_object(T* object_1, T* object_2, T* object_3, std::list<T*> & container) {
			insert_object<T>(object_1, container);
			insert_object<T>(object_2, container);
			insert_object<T>(object_3, container);
		}

		template <typename T>
		inline void insert_object(T* object_1, T* object_2, T* object_3, T* object_4, std::list<T*> & container) {
			insert_object<T>(object_1, container);
			insert_object<T>(object_2, container);
			insert_object<T>(object_3, container);
			insert_object<T>(object_4, container);
		}

		template <typename T>
		inline void delete_object(T* object, std::list<T*> & container, bool calls_destructor) {
			typename std::list<T*>::iterator it_obj = std::find(container.begin(), container.end(), object);
			assert(it_obj != container.end());
			container.erase(it_obj);
			if (std::is_same<T*, Polygon_Vertex*>::value) {
				Polygon_Vertex* v = reinterpret_cast<Polygon_Vertex*>(object);
				v->set_polygon(nullptr);
			}
			if (calls_destructor) {
				delete object;
			}
		}

		template <typename T>
		inline void delete_object(T* object_1, T* object_2, std::list<T*> & container, bool calls_destructor) {
			delete_object<T>(object_1, container, calls_destructor);
			delete_object<T>(object_2, container, calls_destructor);
		}

		template <typename T>
		inline void delete_object(T* object_1, T* object_2, T* object_3, std::list<T*> & container, bool calls_destructor) {
			delete_object<T>(object_1, container, calls_destructor);
			delete_object<T>(object_2, container, calls_destructor);
			delete_object<T>(object_3, container, calls_destructor);
		}

		template <typename T>
		inline void delete_object(T* object_1, T* object_2, T* object_3, T* object_4, std::list<T*> & container, bool calls_destructor) {
			delete_object<T>(object_1, container, calls_destructor);
			delete_object<T>(object_2, container, calls_destructor);
			delete_object<T>(object_3, container, calls_destructor);
			delete_object<T>(object_4, container, calls_destructor);
		}
	public:
		std::list<Polygon_Vertex*> vertices;
		std::list<Polygon_Edge*> edges;
		int seed;
		int running_vertices;

	protected:
		Polygon_Node* cell;
	};
}
