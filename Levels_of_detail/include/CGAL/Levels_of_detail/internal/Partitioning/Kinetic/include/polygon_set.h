#pragma once
#include <set>
#include "polygon_node.h"
#include "polygon_group.h"
#include "polygon_vertex_octree.h"
#include <CGAL/Partition_traits_2.h>


namespace Skippy {

	class Polygon_Set
	{
	public:
		Polygon_Set(int _id_plane, const std::vector<Intersection_Line*> & L);

		~Polygon_Set();

		void insert(const Signature & S, Polygon* P);

		bool exists(const Signature & S);

		bool exists(const Signature & S, const int seed);

		std::map<Signature, Polygon_Node*, Vector_Bool_Comparator>::const_iterator cells_begin() const;
		
		std::map<Signature, Polygon_Node*, Vector_Bool_Comparator>::const_iterator cells_end() const;

		Polygon* get_adjacent_polygon(Polygon* P_ts, Intersection_Line* I);

		std::vector<bool> get_adjacent_polygons_signature(Polygon* P, Intersection_Line* I);

		std::vector<bool> get_adjacent_polygons_signature(Polygon* P, Intersection_Line* I_1, Intersection_Line* I_2);

		void get_signature_of_adjacent_cell(std::vector<bool> & S, Intersection_Line* I);

		void get_polygons(std::list<Polygon*> & P);

		/*void get_polygon_description(std::list<std::list<CGAL_Point_3> > & polygons, std::list<CGAL_Color> & colors, const FT & t);

		void get_sequence_of_3d_vertices(const CGAL::Partition_traits_2<K>::Polygon_2 & polygon, std::list<CGAL_Point_3> & P);

		bool exists_similar_edge_in_adjacent_cell(int d, const Signature & S, const Constraint & C, const Constraint & C_prev, const Constraint & C_next, const CGAL_Point_2 & A, const CGAL_Point_2 & B, const FT & t);

		bool exists_similar_edge_in_adjacent_cell(int d, const Signature & S, const Constraint & C, Polygon_Vertex_R* v, Polygon_Vertex_S* v_s, const FT & t);

		bool exists_similar_edge_in_adjacent_cell(int d, const Signature & S, Intersection_Line* I, const CGAL_Point_2 & A, const CGAL_Point_2 & B, const FT & t);*/

		void get_polygon_description(std::list<std::list<CGAL_Point_3> > & polygons, std::list<CGAL_Color> & colors, const double t);

		void get_sequence_of_3d_vertices(Polygon* polygon, const double t, std::list<CGAL_Point_3> & P);



		void get_polygon_description(Polygon_Vertex_Octree* V, std::list<std::list<int> > & P, const FT & t);

		void get_sequence_of_3d_vertices(Polygon* polygon, const FT & t, std::list<CGAL_Point_3> & P);

	public:
		void build_graph();

		void get_groups(std::list<Polygon_Group*> & G);

		void clear_graph();

	public:
		int id_plane;
		std::map<Intersection_Line*, int> dictionary;
	
	protected:
		std::map<Signature, Polygon_Node*, Vector_Bool_Comparator> nodes;
		std::list<Polygon_Link*> edges;
	};

}