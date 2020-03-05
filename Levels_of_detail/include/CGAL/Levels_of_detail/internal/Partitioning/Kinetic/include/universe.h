#pragma once

#include "defs.h"
#include <map>
#include <list>
#include <vector>
#include <random>

#include "defs_cgal.h"
#include "octree_base.h"

namespace Skippy {

	class Kinetic_Parameters;

	class Polygon_Vertex_R;

	class Support_Plane;

	class Event_Queue;

	class Partition;

	class Oriented_Bounding_Box;

	namespace Universe
	{
		extern std::default_random_engine generator;

		extern Kinetic_Parameters* params;

		extern std::vector<Support_Plane*> map_of_planes;

		extern int moving_objects;

		extern Event_Queue* event_queue;

		extern Octree_Base* point_cloud_octree;

		extern Oriented_Bounding_Box* bounding_box;
	};
}