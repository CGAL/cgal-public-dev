#pragma once
#include <vector>

namespace Skippy {
	namespace Counters {
		extern int id_planes;
		// extern int id_objects;
		extern std::vector<int> id_objects;

		extern int id_partition_vertex;
		extern int id_partition_edge;
		extern int id_partition_facet;
		extern int id_partition_polyhedron;

		extern std::vector<int> par_v_local_ids;
		extern std::vector<int> par_e_local_ids;
	}
}