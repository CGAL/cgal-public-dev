#include "../include/event_queue.h"
#include "../include/parameters.h"
#include <vector>
#include <iterator>
#include <algorithm>
#include <iostream>
#include <cassert>
#include "../include/universe.h"
#include "../include/stats.h"


namespace Skippy {

	using CGAL::to_double;


	Event_Queue::Event_Queue()
	{
		queue = std::map<std::pair<FT, int>, std::list<Event_Vertex_Line*> >();
	}



	Event_Queue::~Event_Queue()
	{
		// There is nothing to do when deleting the queue :
		// all the events that is contains have already been destroyed.
		// In other words, 'queue' is empty. 
	}



	void Event_Queue::print() const
	{
		for (auto it_le = queue.begin(); it_le != queue.end(); ++it_le) {
			std::cout << "** t = " << it_le->first.first << " :: ";
			const std::list<Event_Vertex_Line*> E = it_le->second;
			auto it_e = E.begin();
			while (it_e != E.end()) {
				Event_Vertex_Line* e = (*it_e);
				assert(e->queue_iterator == it_le);
				std::cout << e->intersectant << " " << e->intersected;
				if (++it_e != E.end()) {
					std::cout << ", ";
				} else {
					std::cout << std::endl;
				}
			}
		}
	}



	bool Event_Queue::has_events() const
	{
		// Tests if the queue is empty by calling the corresponding method.
		return queue.empty();
	}



	Event_Vertex_Line* Event_Queue::pop()
	{
		// print();

		if (queue.empty()) return nullptr;

		// Considers the first entry of the map,
		// and returns the first element of the list that is assigned to it.

		std::list<Event_Vertex_Line*> & E = queue.begin()->second;
		Event_Vertex_Line* e_vl = E.front();

		// The element is removed from that list,
		// and if it becomes empty, then the first entry of the queue is removed.

		E.erase(E.begin());
		if (E.empty()) queue.erase(queue.begin());

		// Finally returns the element

		e_vl->is_queued = false;
		return e_vl;
	}



	void Event_Queue::push(Event_Vertex_Line* e_vl)
	{
		// Appends an event to the queue.
		std::pair<std::map<std::pair<FT, int>, std::list<Event_Vertex_Line*> >::iterator, bool> it_pushed;
		it_pushed = queue.insert(std::make_pair(std::make_pair(e_vl->t_intersectant, e_vl->plane), std::list<Event_Vertex_Line*>(1, e_vl)));

		if (!it_pushed.second) {
			// If we reach this block, it means that the intersection couldn't immediately take place :
			// there exists events that are simultaneous to e_vl.
			// No problem, we just append e_vl to a list of simultaneous events.
			it_pushed.first->second.push_back(e_vl);
		}

		// Sets the event as queued,
		// and sets a pointer to the map entry where it is saved.
		e_vl->is_queued = true;
		e_vl->queue_iterator = it_pushed.first;
	}



	void Event_Queue::push(const std::list<Event_Vertex_Line*> & E_VL)
	{
		for (std::list<Event_Vertex_Line*>::const_iterator it_e_vl = E_VL.cbegin(); it_e_vl != E_VL.end(); ++it_e_vl) {
			push(*it_e_vl);
		}
	}



	void Event_Queue::push(Event_Vertex_Line* e_vl, const std::map<std::pair<FT, int>, std::list<Event_Vertex_Line*> >::iterator it_push)
	{
		// Appends an event to the queue, with prior knowledge of the entry of the map where e_vl should be inserted.
		// We suppose that it_loc points to a valid entry of the map, and there exists a list at it_loc->second.
		// All we have to do is to insert e_vl at the end of that list.

		assert(it_push != queue.end());
		it_push->second.push_back(e_vl);

		// Sets the event as queued,
		// and sets a pointer to the map entry where it is saved.
		e_vl->is_queued = true;
		e_vl->queue_iterator = it_push;
	}



	void Event_Queue::erase(Event_Vertex_Line* e_vl)
	{
		// To erase a given event, we first access to the entry of the queue where it is stored.
		std::map<std::pair<FT, int>, std::list<Event_Vertex_Line*> >::iterator it_loc = e_vl->queue_iterator;
		assert(it_loc != queue.end());

		// This entry corresponds to a list of elements.
		// We remove e_vl from that list.
		std::list<Event_Vertex_Line*> & E = it_loc->second;
		std::list<Event_Vertex_Line*>::iterator it_e = std::find(E.begin(), E.end(), e_vl);
		assert(it_e != E.end());
		E.erase(it_e);

		// If the list of simultaneous events gets null,
		// then we remove the entry of the map that corresponds to e_vl.
		if (E.empty()) queue.erase(it_loc);

		// The event is no longer queued.
		e_vl->is_queued = false;
		e_vl->queue_iterator = queue.end();
	}



	void Event_Queue::erase(const std::list<Event_Vertex_Line*> & E_VL)
	{
		for (std::list<Event_Vertex_Line*>::const_iterator it_e_vl = E_VL.cbegin(); it_e_vl != E_VL.end(); ++it_e_vl) {
			erase(*it_e_vl);
		}
	}



	void Event_Queue::get_simultaneous_events_for_this_vertex(const int plane, const int intersectant, const FT & t_intersectant, std::list<Event_Vertex_Line*> & E_VL)
	{
		// Gets an iterator to the first entry of the map,
		// which contains the events that come just after, or simultaneously to e_vl.
		// If there are no more events or not occuring at the same time, returns.

		std::map<std::pair<FT, int>, std::list<Event_Vertex_Line*> >::iterator it_queue_begin = queue.begin();
		if (it_queue_begin == queue.cend() || it_queue_begin->first.first != t_intersectant) {
			return;
		}

		// If there are simultaneous events, then checks the indices of the colliding vertices.
		// If they coincide, then pops the event from the queue.

		std::list<Event_Vertex_Line*> & E = it_queue_begin->second;
		std::list<Event_Vertex_Line*>::iterator it_e = E.begin();
		while (it_e != E.end()) {
			Event_Vertex_Line* e = (*it_e);
			if (e->plane == plane && e->intersectant == intersectant) {
				it_e = E.erase(it_e);
				e->is_queued = false;
				e->queue_iterator = queue.end();
				E_VL.push_back(e);
			} else {
				++it_e;
			}
		}

		if (E.empty()) {
			queue.erase(it_queue_begin);
		}
	}



	void Event_Queue::get_references_to_simultaneous_events_for_given_plane(const int plane, const FT & t, std::list<Event_Vertex_Line*> & E) const
	{
		// As this function is called, an event has been popped from the queue, from the entry (t, plane).
		// We are looking for simultaenous events that apply to the same polygon as that event.

		// We get an iterator to the first entry of the map.
		// If there exists such events, there must still exist an entry (t, plane).
		// If it is not the case, returns. Otherwise, copies events, which are not popped from the queue.

		std::map<std::pair<FT, int>, std::list<Event_Vertex_Line*> >::const_iterator it_queue_begin = queue.cbegin();
		if (it_queue_begin == queue.cend() || it_queue_begin->first.second != plane || it_queue_begin->first.first != t) {
			return;
		} else {
			const std::list<Event_Vertex_Line*> & elements = it_queue_begin->second;
			std::copy(elements.begin(), elements.end(), std::back_inserter(E));
		}
	}
}