#include "../include/partition_simple.h"
#include "../include/partition.h"
#include "../include/octree_base.h"
#include "../include/vars.h"



namespace Skippy 
{
	Partition_Simple::Partition_Simple(const std::vector<CGAL_Plane> & plane_defs)
		: Partition(plane_defs)
	{
	}


	Partition_Simple::~Partition_Simple()
	{
	}


	void Partition_Simple::build()
	{
		std::list<Partition_Vertex*> V;
		std::list<Partition_Edge*> E;

		// Part 1 : support planes of polygons
		// Groups of adjacent polygons inside trees directly define planar facets
		build_inner_facets(V, E);

		// Part 2 : facets of the bounding box
		// We first need to create the missing elements of the bounding box (corner vertices and edges)
		// After that we loop on the edges to get facets
		build_missing_vertices_for_bounding_box(V);
		build_missing_edges_for_bounding_box(E);
		
		edges = E; 
		E.clear();

		build_facets_for_bounding_box();

		// Part 3 : cleans the partition
		Partition::remove_bivalent_vertices(V, edges);
		reset_indices(V, edges);

		// Part 4 : get final polyhedrons
		build_polyhedrons();
	}
}