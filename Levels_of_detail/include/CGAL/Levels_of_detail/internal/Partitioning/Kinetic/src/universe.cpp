#include "../include/universe.h"
#include <chrono>

namespace Skippy {

#pragma warning(suppress:4244)
	std::default_random_engine Universe::generator(0);

	Kinetic_Parameters* Universe::params = nullptr;

	std::vector<Support_Plane*> Universe::map_of_planes;

	int Universe::moving_objects = 0;

	Event_Queue* Universe::event_queue = nullptr;

	Octree_Base* Universe::point_cloud_octree = nullptr;

	Oriented_Bounding_Box* Universe::bounding_box = nullptr;
}