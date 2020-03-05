#pragma once
#include "polygon.h"

namespace Skippy {

	typedef std::vector<bool> Signature;

	class Polygon_Link;

	class Polygon_Node
	{
	public:
		Polygon_Node(const Signature & S, Polygon* P);

		~Polygon_Node();

		static Signature make_signature(const CGAL_Point_2 & O, const std::vector<Intersection_Line*> & L);

		Polygon* get_one_polygon() const;

		void push(Polygon* P);

		Signature get_signature() const;

		size_t size() const;

		std::list<Polygon*>::const_iterator polygons_begin() const;
		
		std::list<Polygon*>::const_iterator polygons_end() const;

		bool contours_are_empty() const;

		void set_contours(const std::list<Constraint> & _contours);

		size_t contours_size() const;

		Constraint get_contour(size_t i) const;

		void set_contour_as_not_on_border(size_t i);

		void set_contour_as_not_on_border(Intersection_Line* I);

		void append_contours(std::list<std::tuple<Constraint, Constraint, Constraint> > & T) const;

		void insert_link(Polygon_Link* l);

		void clear_links();

		bool is_visited() const;

		void set_visited();

		bool is_queued() const;

		void set_queued();

		std::list<Polygon_Link*>::const_iterator links_begin() const;

		std::list<Polygon_Link*>::const_iterator links_end() const;

	protected:
		Signature signature;
		std::list<Polygon*> polygons;
		std::vector<Constraint> contours;
		bool* borders;

		std::list<Polygon_Link*> links;
		bool visited_node;
		bool queued_node;
	};


	class Polygon_Link
	{
	public:	
		Polygon_Link(Polygon_Node* v1, Polygon_Node* v2);

		~Polygon_Link();

		Polygon_Node* other_node(Polygon_Node* v) const;

	protected:
		Polygon_Node* v1;
		Polygon_Node* v2;
	};


	struct Vector_Bool_Comparator
	{
		bool operator() (const std::vector<bool> & SL, const std::vector<bool> & SR) const
		{
			return SL < SR;
		}
	};
}