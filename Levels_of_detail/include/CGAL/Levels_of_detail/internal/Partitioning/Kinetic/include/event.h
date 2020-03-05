#pragma once
#include "defs_cgal.h"
#include <map>

namespace Skippy {

	typedef enum {
		NO_SCHEDULE = 0,
		SCHEDULE = 1,
		SORT = 2,
	} Event_Flags;


	class Event_Vertex_Line
	{
	public:
		Event_Vertex_Line(const int _plane, const int _intersectant, const int _intersected, const FT & _t_intersectant);

		~Event_Vertex_Line();

	public:
		int plane;
		int intersectant;
		int intersected;
		FT t_intersectant;

		bool is_queued;
		std::map<std::pair<FT, int>, std::list<Event_Vertex_Line*> >::iterator queue_iterator;
	};



	inline bool compare_events_by_intersectant_object(Event_Vertex_Line* e_i, Event_Vertex_Line* e_j) {
		return (e_i->intersectant < e_j->intersectant);
	}



	inline bool compare_events_by_intersected_object(Event_Vertex_Line* e_i, Event_Vertex_Line* e_j) {
		return (e_i->intersected < e_j->intersected);
	}



	inline bool compare_events_by_time(Event_Vertex_Line* e_i, Event_Vertex_Line* e_j) {
		return (e_i->t_intersectant < e_j->t_intersectant);
	}
}