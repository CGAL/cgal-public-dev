#include "../include/event.h"
#include <iterator>
#include <iostream>


namespace Skippy 
{

	Event_Vertex_Line::Event_Vertex_Line(const int _plane, const int _intersectant, const int _intersected, const FT & _t_intersectant)
		: plane (_plane),
		intersectant(_intersectant),
		intersected(_intersected),
		t_intersectant(_t_intersectant),
		is_queued(false)
	{
		// The iterator is left uninitialized.
		// Indeed, the event is not assigned yet to an entry of the queue.
	}


	Event_Vertex_Line::~Event_Vertex_Line()
	{
	}
}