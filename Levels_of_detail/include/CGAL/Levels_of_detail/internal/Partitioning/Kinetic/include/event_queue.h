#pragma once
#include <map>
#include <vector>

#include "event.h"


namespace Skippy {
	class Event_Queue
	{
	public:
		Event_Queue();

		~Event_Queue();

	public:
		bool has_events() const;

		Event_Vertex_Line* pop();

		void push(Event_Vertex_Line* e_vl);

		void push(const std::list<Event_Vertex_Line*> & E_VL);

		void push(Event_Vertex_Line* e_vl, const std::map<std::pair<FT, int>, std::list<Event_Vertex_Line*> >::iterator it_push);

		void erase(Event_Vertex_Line* e_vl);

		void erase(const std::list<Event_Vertex_Line*> & E_VL);

		void get_simultaneous_events_for_this_vertex(const int plane, const int intersectant, const FT & t_intersectant, std::list<Event_Vertex_Line*> & E);

		void get_references_to_simultaneous_events_for_given_plane(const int plane, const FT & t, std::list<Event_Vertex_Line*> & E) const;

	protected:
		void print() const;

	protected:
		std::map<std::pair<FT, int>, std::list<Event_Vertex_Line*> > queue;
	};
}