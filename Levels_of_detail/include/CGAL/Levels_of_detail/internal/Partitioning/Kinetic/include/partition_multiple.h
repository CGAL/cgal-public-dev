#pragma once
#include "partition.h"


namespace Skippy
{
	class Partition_Multiple : public Partition
	{
	public:
		Partition_Multiple(const std::vector<CGAL_Plane> & plane_defs);

		~Partition_Multiple();

		void incremental_build();

		void finalize(const int slices);

	protected:
		void rebuild_edges_at_boundaries(const int slices);

		void check(const int id, const CGAL_Plane & P);

		void rebuild_edges(const int id, const CGAL_Plane & P, bool edges_may_intersect);

		void exhibit_and_merge_duplicate_paths(const int id, Partition_Vertex* v, Partition_Edge* e1, Partition_Edge* e2,
			std::map<int, std::list<Partition_Edge*>::iterator> & E_map);

		void exhibit_duplicate_path(Partition_Vertex* v_init, Partition_Edge* e_init,
			std::vector<Partition_Edge*> & E, std::vector<Partition_Vertex*> & V);

		Partition_Vertex* shorten_duplicate_paths(std::vector<Partition_Edge*> & E_1, std::vector<Partition_Vertex*> & V_1,
			std::vector<Partition_Edge*> & E_2, std::vector<Partition_Vertex*> & V_2);

		void sort_vertices_by_distance(Partition_Vertex* v_source, Partition_Vertex* v_dest,
			std::vector<Partition_Vertex*> & V_1, std::vector<Partition_Vertex*> & V_2,
			std::vector<std::tuple<Partition_Vertex*, int, FT> > & V);

		void merge_duplicate_paths(Partition_Vertex* v_source,
			std::vector<Partition_Edge*> & E_1, std::vector<Partition_Edge*> & E_2,
			std::vector<std::tuple<Partition_Vertex*, int, FT> > & V,
			std::map<int, std::list<Partition_Edge*>::iterator> & E_map);

		void remove_overlapping_edges(const int id, const CGAL_Plane & P,
			std::map<int, std::list<Partition_Edge*>::iterator> & E_map);

		void remove_self_intersecting_edges(const CGAL_Plane & P, std::list<Partition_Vertex*> & V,
			std::map<int, std::list<Partition_Edge*>::iterator> & E_map);

		void get_bounding_box_for_rtree(Partition_Edge* e, const std::map<int, CGAL_Inexact_Point_2> & hints_points_2d,
			const double & eps, Boost_Box & B);

		void make_contiguous_indices(const int slices);

		void build_facets_on_slicing_planes(const int slices);

		void reset_indices_of_facets();
	};
}