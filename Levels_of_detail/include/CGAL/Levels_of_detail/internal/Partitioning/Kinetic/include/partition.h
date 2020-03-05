#pragma once
#include "support_plane_objects.h"
#include "partition_objects.h"
#include "partition_vertex_octree.h"


namespace Skippy {

	class Partition
	{
	public:
		Partition(const std::vector<CGAL_Plane> & plane_defs);
		
		virtual ~Partition();

	protected:
		typedef std::pair<CGAL_Direction_2, Partition_HalfEdge> Direction_H;

		typedef std::pair<CGAL_Direction_2, Partition_Side> Direction_S;

		void build_inner_facets(std::list<Partition_Vertex*> & unique_vertices,
			std::list<Partition_Edge*> & unique_edges);

		Partition_Vertex* get_partition_vertex(const std::pair<Constraint, Constraint> & v, 
			std::map<std::pair<int, int>, Partition_Vertex*> & conversions,
			std::map<int, std::pair<Partition_Edge*, Partition_Edge*> > & adjacent_edges,
			std::list<Partition_Vertex*> & unique_vertices);

		Partition_Edge* get_partition_edge(Partition_Vertex* v1, Partition_Vertex* v2,
			std::map<int, std::pair<Partition_Edge*, Partition_Edge*> > & adjacent_edges,
			std::list<Partition_Edge*> & unique_edges);

		void get_sequence_of_partition_edges(std::map<int, std::pair<Partition_Edge*, Partition_Edge*> > & adjacent_edges, std::list<Partition_Edge*> & sorted_edges);


		void build_missing_vertices_for_bounding_box(std::list<Partition_Vertex*> & unique_vertices);

		void build_missing_edges_for_bounding_box(std::list<Partition_Edge*> & unique_edges);

		void build_facets_for_bounding_box();

		int get_orientation_of_normal_vector(const CGAL_Plane & P);

		void build_facets_on_boundaries(const int P_id, const CGAL_Plane & P);

		void build_halfedges(const int H_index, const int H_or, const CGAL_Plane & H, Partition_Edge* e, std::vector<std::vector<Direction_H> > & directions, bool** & queued);

		void add_halfedge(std::vector<std::vector<Direction_H> > & D, int id_v, const CGAL_Vector_2 & u_theta, Partition_Edge* e, bool v1_v2);


		void loop_and_build_facets(const int H_id, const int H_or, std::vector<std::vector<Direction_H> > & directions, bool** & queued);


		void remove_bivalent_vertices();

		void remove_bivalent_vertices(std::list<Partition_Vertex*> & V, std::list<Partition_Edge*> & E);



		void reset_indices(std::list<Partition_Vertex*> & V, std::list<Partition_Edge*> & E);

		void build_polyhedrons();

		void index_facets_sides(std::vector<std::vector<Direction_S> > & D);

		void get_local_frame(Partition_Edge* e, Partition_Facet* f, Partition_Vertex* & O, CGAL_Vector_3 & v_i, CGAL_Vector_3 & v_j, CGAL_Vector_3 & v_k);

		void add_side(std::vector<Direction_S> & D, const CGAL_Direction_2 & u_theta, Partition_Facet* F, bool positive_side_comes_first);

		CGAL_Direction_2 get_angle_of_projected_vertex(const CGAL_Point_3 & A, Partition_Vertex* O, const CGAL_Vector_3 & v_i, const CGAL_Vector_3 & v_j, const CGAL_Vector_3 & v_k);

		void init_table_of_queued_sides(bool** & queued);

		void delete_table_of_queued_sides(bool** & queued);

		void loop_and_build_polyhedrons(std::vector<std::vector<Direction_S> > & directions, bool** & queued);

		void debug();

		void set_built_while_processing_negative_side_of_slicing_plane(Partition_Edge* e);

	public:
		KINETIC_PARTITION_API void ply_facets(const std::string & filename, std::vector<CGAL_Color> & colors) const;

		KINETIC_PARTITION_API void ply_individual_polyhedron(const std::string & filename, const int P) const;

		KINETIC_PARTITION_API std::string export_partition(const std::vector<int> & polygons_to_planes);

	protected:
		std::string export_partition_without_reorientation(const std::vector<int> & polygons_to_planes) const;

		std::string export_partition_with_reorientation(const std::vector<int> & polygons_to_planes);

		std::string export_partition_elements(const std::vector<Partition_Vertex*> & vertices) const;

	public:
		KINETIC_PARTITION_API void export_partition(const std::vector<int> & polygons_to_planes, const std::string & filename);


		KINETIC_PARTITION_API void get_all_vertices(std::list<Partition_Vertex*> & V) const;

		KINETIC_PARTITION_API void get_all_vertices_sorted_by_identifier(std::vector<Partition_Vertex*> & V) const;


		KINETIC_PARTITION_API std::list<Partition_Facet*>::const_iterator planar_facets_begin(const int id) const;

		KINETIC_PARTITION_API std::list<Partition_Facet*>::const_iterator planar_facets_end(const int id) const;

		KINETIC_PARTITION_API std::list<Partition_Facet*>::iterator planar_facets_begin(const int id);

		KINETIC_PARTITION_API std::list<Partition_Facet*>::iterator planar_facets_end(const int id);


		KINETIC_PARTITION_API std::list<Partition_Polyhedron*>::iterator polyhedrons_begin();

		KINETIC_PARTITION_API std::list<Partition_Polyhedron*>::const_iterator polyhedrons_begin() const;

		KINETIC_PARTITION_API std::list<Partition_Polyhedron*>::iterator polyhedrons_end();

		KINETIC_PARTITION_API std::list<Partition_Polyhedron*>::const_iterator polyhedrons_end() const;

		KINETIC_PARTITION_API size_t polyhedrons_size() const;

		std::vector<std::list<Partition_Facet*> > facets;

		const std::vector<CGAL_Plane> & planes;

	protected:
		std::vector<double> dims;
		Partition_Vertex_Octree* octree_vertices;
		std::list<Partition_Edge*> edges;
		std::list<Partition_Polyhedron*> polyhedrons;
	};

}