#pragma once
#include "defs.h"
#include "propagation.h"



namespace Skippy 
{
	class Kinetic_Propagation_Multiple : public Kinetic_Propagation
	{
	protected:
		struct Support_Plane_Basic_Info 
		{
			Support_Plane_Basic_Info(int _id, 
				const std::vector<std::list<int> > & _BF, 
				const std::list<CGAL_Point_3> & _BP_3, 
				const std::list<int> & _PR) 
			{
				id = _id;
				BF = _BF;
				BP_3 = _BP_3;
				PR = _PR;
			}
			
			int id;
			std::vector<std::list<int> > BF;
			std::list<CGAL_Point_3> BP_3;
			std::list<int> PR;
		};

	public:
		KINETIC_PARTITION_API Kinetic_Propagation_Multiple();

		KINETIC_PARTITION_API Kinetic_Propagation_Multiple(const int N);

		KINETIC_PARTITION_API Kinetic_Propagation_Multiple(int argc, char *argv[]);

		KINETIC_PARTITION_API Kinetic_Propagation_Multiple(const std::string & filename, Preprocess process = NONE);

		KINETIC_PARTITION_API Kinetic_Propagation_Multiple(const std::vector<std::vector<CGAL_Inexact_Point_3> > & primitives, Preprocess process = NONE);

		KINETIC_PARTITION_API Kinetic_Propagation_Multiple(const std::vector<std::vector<CGAL_Point_3> > & primitives, Preprocess process = NONE);

		KINETIC_PARTITION_API virtual ~Kinetic_Propagation_Multiple();

		KINETIC_PARTITION_API void run();

		void delete_multiple_kinetic_data_structure();

	protected:
		void init_multiple_kinetic_data_structure(const int i, const int j, const int k, const int N,
			std::map<int, std::map<int, CGAL_Line_3> > & lines);

		int init_support_planes_slices_and_lines(CGAL_Point_3 & pt_min,
			CGAL_Point_3 & pt_max,
			std::vector<CGAL_Point_3> & box_corners,
			std::vector<std::pair<int, int> > & box_edges,
			std::map<int, std::map<int, CGAL_Line_3> > & lines);

		void init_slicing_planes();

		void init_intersection_lines_for_slices_and_facets(const int slices,
			std::map<int, std::map<int, CGAL_Line_3> > & lines);

		void init_intersection_lines_partially(const std::vector<int> & indices,
			std::map<int, std::map<int, CGAL_Line_3> > & lines);

		//void project_primitives_to_their_respective_planes();

		void decompose_primitives_in_subvolumes(const CGAL_Point_3 & pt_min,
			const CGAL_Point_3 & pt_max,
			const std::vector<CGAL_Point_3> & box_corners,
			const std::vector<std::pair<int, int> > & box_edges,
			const std::map<int, std::map<int, CGAL_Line_3> > & lines,
			const int gx, const int gy, const int gz);

		void init_current_bounding_box(const int i, const int j, const int k,
			const int N,
			std::vector<int> & bounding_planes_indices,
			CGAL_Point_3 & pt_min, 
			CGAL_Point_3 & pt_max,
			std::vector<CGAL_Point_3> & box_corners,
			std::vector<std::pair<int, int> > & box_edges);

		void init_current_support_planes(const int i, const int j, const int k,
			const int N,
			std::vector<int> & bounding_planes_indices,
			const CGAL_Point_3 & pt_min,
			const CGAL_Point_3 & pt_max,
			const std::vector<CGAL_Point_3> & box_corners,
			const std::vector<std::pair<int, int> > & box_edges,
			const std::map<int, std::map<int, CGAL_Line_3> > & lines,
			const std::vector<Support_Plane_Basic_Info> & considered_inner_planes);

		void exhibit_considered_inner_planes(const int i, const int j, const int k,
			const int N,
			std::vector<int> & bounding_planes_indices,
			const CGAL_Point_3 & pt_min,
			const CGAL_Point_3 & pt_max,
			const std::vector<CGAL_Point_3> & box_corners,
			const std::vector<std::pair<int, int> > & box_edges,
			std::vector<Support_Plane_Basic_Info> & considered_planes,
			std::vector<int> & all_planes_indices);

		void get_subset_of_3d_lines(const std::map<int, std::map<int, CGAL_Line_3> > & lines,
			std::map<int, std::map<int, CGAL_Line_3> > & subset);

		void build_polygons(const std::vector<Support_Plane_Basic_Info> & considered_planes) const;

	protected:
		std::vector<CGAL_Point_2> primitives_to_barycenters;
		std::vector<std::vector<std::pair<CGAL_Point_2, CGAL_Vector_2> > > primitives_to_directions;
		std::vector<std::set<std::tuple<int, int, int> > > primitives_to_subvolumes;
	};
}