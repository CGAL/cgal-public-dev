#include "../include/propagation_multiple.h"
#include "../include/propagation.h"
#include "../include/universe.h"
#include "../include/parameters.h"
#include "../include/support_plane.h"
#include "../include/vars.h"
#include "../include/event_queue.h"
#include "../include/partition_multiple.h"
#include "../include/stats.h"
#include <CGAL/Polygon_2.h>
#include <CGAL/intersection_2.h>
#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/convex_hull_2.h>



namespace Skippy
{
	Kinetic_Propagation_Multiple::Kinetic_Propagation_Multiple()
		: Kinetic_Propagation()
	{
	
	}



	Kinetic_Propagation_Multiple::Kinetic_Propagation_Multiple(const int N)
		: Kinetic_Propagation(N)
	{
	
	}



	Kinetic_Propagation_Multiple::Kinetic_Propagation_Multiple(int argc, char *argv[])
		: Kinetic_Propagation(argc, argv)
	{
	
	}



	Kinetic_Propagation_Multiple::Kinetic_Propagation_Multiple(const std::string & filename, Preprocess process)
		: Kinetic_Propagation(filename, process)
	{
	
	}



	Kinetic_Propagation_Multiple::Kinetic_Propagation_Multiple(const std::vector<std::vector<CGAL_Inexact_Point_3> > & primitives, Preprocess process)
		: Kinetic_Propagation(primitives, process)
	{
	
	}



	Kinetic_Propagation_Multiple::Kinetic_Propagation_Multiple(const std::vector<std::vector<CGAL_Point_3> > & primitives, Preprocess process)
		: Kinetic_Propagation(primitives, process)
	{
	
	}



	Kinetic_Propagation_Multiple::~Kinetic_Propagation_Multiple()
	{
	
	}



	void Kinetic_Propagation_Multiple::run()
	{
		bool perfs = Universe::params->perfs;

		if (polygons.empty()) return;

		std::cout << "Basename : " << Universe::params->basename << std::endl;

		clock_t t_0 = clock();

		Kinetic_Parameters* params = Universe::params;
		int gx = params->grid_x, gy = params->grid_y, gz = params->grid_z;
		int slices = gx + gy + gz - 3;

		try {

			// Part 1.
			// Initialization of the global kinetic data structure.

			// We set the different planes of the kinetic system :
			//   planes[0 .. 5]     : main bounding box facets
			//   planes[6 .. N - 1] : support planes of planar primitives
			//   planes[N .. N + M] : slicing planes, with M = (g_x - 1) + (g_y - 1) + (g_z - 1)
			// but also the 3D intersection lines between these planes.

			CGAL_Point_3 pt_min, pt_max;
			std::vector<CGAL_Point_3> box_corners;
			std::vector<std::pair<int, int> > box_edges;
			std::map<int, std::map<int, CGAL_Line_3> > lines;

			int N = init_support_planes_slices_and_lines(pt_min, pt_max, box_corners, box_edges, lines);
			decompose_primitives_in_subvolumes(pt_min, pt_max, box_corners, box_edges, lines, gx, gy, gz);

			// We also initialize the partition by creating an octree of the size of the system.
			partition = new Partition_Multiple(planes);
			Partition_Multiple* this_partition = dynamic_cast<Partition_Multiple*>(partition);

			KP_Stats::initialization_time = double(clock() - t_0) / CLOCKS_PER_SEC;

			std::cout << "** Initialized structure" << std::endl;

			// Part 2.
			// For each subvolume V of the main bounding box, identified by a triplet (i, j, k)

			int id_block = 0;
			int subvolumes = gx * gy * gz;
			int spaces = int(1 + floor(log10(subvolumes)));

			for (int i = 0; i < gx; ++i) {
				for (int j = 0; j < gy; ++j) {
					for (int k = 0; k < gz; ++k) {

						std::cout << "** Processes block : " << std::setw(spaces) << ++id_block << "/" << subvolumes << '\r' << std::flush;

						// Initializes a kinetic data structure
						init_multiple_kinetic_data_structure(i, j, k, N, lines);

						// Processes all events
						unstack();

						// Incrementally builds the partition,
						// by adding elements related to the processing of this kinetic data structure
						this_partition->incremental_build();

						// Deletes the kinetic data structure
						delete_multiple_kinetic_data_structure();

					}
				}
			}
			
			std::cout << std::endl;

			// Part 3.
			// Finalizes the partition by building facets on slicing planes and borders.

			this_partition->finalize(slices);
			export_partition();

			KP_Stats::destruction_time /= CLOCKS_PER_SEC;
			KP_Stats::running_time = double(clock() - t_0) / CLOCKS_PER_SEC;
			KP_Stats::p_initialization_time = 100 * KP_Stats::initialization_time / KP_Stats::running_time;
			KP_Stats::p_unstack_time = 100 * KP_Stats::unstack_time / KP_Stats::running_time;
			KP_Stats::p_partition_time = 100 * KP_Stats::partition_time / KP_Stats::running_time;
			KP_Stats::p_destruction_time = 100 * KP_Stats::destruction_time / KP_Stats::running_time;

		} catch (std::exception & e) {
			std::cout << e.what() << std::endl;
		}

		// delete Universe::params;

		// print_performances(perfs);
	}



	void Kinetic_Propagation_Multiple::decompose_primitives_in_subvolumes(const CGAL_Point_3 & pt_min,
		const CGAL_Point_3 & pt_max,
		const std::vector<CGAL_Point_3> & box_corners,
		const std::vector<std::pair<int, int> > & box_edges,
		const std::map<int, std::map<int, CGAL_Line_3> > & lines,
		const int gx, const int gy, const int gz)
	{
		// In this function, we initialize a vector of support planes which correspond to
		// the input primitives and the slicing planes. However, our goal is to decompose
		// the input primitives into subpolygons that will later be used to initialize each
		// subvolume of the main bounding box.

		// In practice, we do as if we were going to initialize a kinetic data structure 
		// with all the polygons. We initiliaze Support_Plane objects, Intersection_Lines,
		// and polygons that we decompose.

		int slices = gx + gy + gz - 3;
		
		/*std::vector<int> considered_indices (6 + slices);
		for (int i = 0; i < 6; ++i) considered_indices[i] = i;
		for (int i = 0; i < slices; ++i) considered_indices[i + 6] = int(planes.size()) - slices + i;*/

		init_supporting_planes(pt_min, pt_max, box_corners, box_edges, lines);

		size_t n = polygons.size();
		primitives_to_barycenters.reserve(n);
		primitives_to_directions.reserve(n);
		primitives_to_subvolumes.reserve(n);

		int spaces = int(1 + floor(log10(polygons.size())));

		for (size_t pol_id = 0; pol_id < polygons.size(); ++pol_id) {
			const std::vector<CGAL_Point_3> & P = polygons[pol_id];
			
			// std::cout << "** Initializes : " << std::setw(spaces) << pol_id << "/" << polygons.size() << '\r' << std::flush;

			CGAL_Point_2 P_barycenter;
			std::vector<std::pair<CGAL_Point_2, CGAL_Vector_2> > P_directions;
			std::set<std::tuple<int, int, int> > P_subvolumes;

			Support_Plane* SP = Universe::map_of_planes[polygons_to_planes[pol_id]];
			SP->project_and_decompose_with_respect_to_slices(P, gx, gy, gz, P_barycenter, P_directions, P_subvolumes);

			primitives_to_barycenters.push_back(P_barycenter);
			primitives_to_directions.push_back(P_directions);
			primitives_to_subvolumes.push_back(P_subvolumes);
		}

		// We delete the universe.

		for (size_t i = 0; i < Universe::map_of_planes.size(); i++) {
			delete Universe::map_of_planes[i];
		}
		Universe::map_of_planes.clear();

		Counters::id_planes = -1;
		Counters::id_objects = std::vector<int>();
	}



	void Kinetic_Propagation_Multiple::init_multiple_kinetic_data_structure(const int i, const int j, const int k, const int N,
		std::map<int, std::map<int, CGAL_Line_3> > & lines)
	{
		clock_t t_0 = clock();

		try {

			// Step 1.
			// We first delimit the bounding box of the current kinetic system.

			std::vector<int> bounding_planes_indices;
			CGAL_Point_3 pt_min, pt_max;
			std::vector<CGAL_Point_3> box_corners;
			std::vector<std::pair<int, int> > box_edges;
			init_current_bounding_box(i, j, k, N, bounding_planes_indices, pt_min, pt_max, box_corners, box_edges);

			// Step 2.
			// We identify the support planes whose primitives intersect the current bounding box.
			// The propagation will be performed on those planes.

			std::vector<Support_Plane_Basic_Info> inner_planes;
			std::vector<int> all_planes_indices;

			exhibit_considered_inner_planes(i, j, k, N, bounding_planes_indices, pt_min, pt_max, box_corners, box_edges, inner_planes, all_planes_indices);
			init_intersection_lines_partially(all_planes_indices, lines);

			// Step 3.
			// Finally, we initialize the kinetic data structure.

			init_current_support_planes(i, j, k, N, bounding_planes_indices, pt_min, pt_max, box_corners, box_edges, lines, inner_planes);
			init_queue_of_events();
			build_polygons(inner_planes);
			init_schedule();

		} catch (std::exception & except) {
			throw except;
		}

		KP_Stats::initialization_time += double(clock() - t_0) / CLOCKS_PER_SEC;

		// std::cout << "** Initialized structure" << std::endl;
	}



	int Kinetic_Propagation_Multiple::init_support_planes_slices_and_lines(CGAL_Point_3 & pt_min,
		CGAL_Point_3 & pt_max,
		std::vector<CGAL_Point_3> & box_corners,
		std::vector<std::pair<int, int> > & box_edges,
		std::map<int, std::map<int, CGAL_Line_3> > & lines)
	{
		int N = 0;

		try {
			// Builds an octree representing the point cloud
			const std::string & filename_point_cloud = Universe::params->path_point_cloud;
			if (!filename_point_cloud.empty() && Universe::params->stopping_condition != 0) {
				init_point_cloud_octree(filename_point_cloud);
			}

			// Initializes the plane equations
			pre_init_main_bounding_box();
			init_plane_equations();
			
			N = int(planes.size()) - 1;
			
			init_slicing_planes();

			int slices = Universe::params->grid_x + Universe::params->grid_y + Universe::params->grid_z - 3;
			init_colors(slices);

			init_main_bounding_box(pt_min, pt_max, box_corners, box_edges);
			init_intersection_lines_for_slices_and_facets(slices, lines);

			return N;

		} catch (std::exception & e) {
			throw e;
		}
	}



	void Kinetic_Propagation_Multiple::init_slicing_planes()
	{
		// For each direction of the main bounding box
		for (int p = 0 ; p < 3 ; ++p) {

			/*std::cout << "Axis : ";
			if (p == 0) std::cout << "X" << std::endl;
			if (p == 1) std::cout << "Y" << std::endl;
			if (p == 2) std::cout << "Z" << std::endl;*/
			
			// Evaluates the number of required slices for this direction
			int slices = 1;
			switch (p) {
			case 0: slices = Universe::params->grid_x; break;
			case 1: slices = Universe::params->grid_y; break;
			case 2: slices = Universe::params->grid_z; break;
			}
			if (slices == 1) continue;

			// Implicitely gets the orthogonal vector via P_ref
			const CGAL_Plane & P_ref = planes[2 * p];
			FT a = P_ref.a(), b = P_ref.b(), c = P_ref.c();

			// Gets the mininal and maximal abscissas along this direction
			FT d_inf = -planes[2 * p + 0].d(), d_sup = -planes[2 * p + 1].d();
			//std::cout << "[" << 2 * p + 0 << "]" << " " << d_inf << std::endl;
			FT d_step = (d_sup - d_inf) / slices;

			// Iterates between d_inf and d_sup
			// In total, we insert (grid_x - 1) + (grid_y - 1) + (grid_z - 1) new equations
			for (int s = 0 ; s < slices - 1 ; ++s) {
				FT d = d_inf + (s + 1) * d_step;
				CGAL_Plane P = CGAL_Plane(a, b, c, -d);
				//std::cout << "[" << planes.size() << "]" << " " << d << std::endl;
				planes.push_back(P);
			}

			//std::cout << "[" << 2 * p + 1 << "]" << " " << d_sup << std::endl;
		}
	}



	void Kinetic_Propagation_Multiple::init_intersection_lines_for_slices_and_facets(const int slices,
		std::map<int, std::map<int, CGAL_Line_3> > & lines)
	{
		int n = int(planes.size());

		// As this function is called, we have initialized a global kinetic data structure,
		// and all kinds of planes : bounding box facets, inner planes, and slicing planes.
		
		// In this function we compute some intersection lines between all these planes.
		// Couples of planes that are taken into account are : (F, F), (F, P), (F, S), 
		// (P, S), (S, S) where F is a facet, P is an inner plane, and S is a slicing plane.

		// Couples of planes (P, P') are currently discarded.

		int id_facet_min = 0, id_facet_max = 5;
		int id_primitive_min = 6, id_primitive_max = n - slices - 1;
		int id_slice_min = n - slices, id_slice_max = n - 1;

		int type_i = -1, type_j = -1;

		for (int i = 0; i < n; i++) {
			
			if (jin(id_facet_min, i, id_facet_max)) {
				type_i = 0; // the plane is a facet
			} else if (jin(id_primitive_min, i, id_primitive_max)) {
				type_i = 1; // the plane results from a primitive
			} else if (jin(id_slice_min, i, id_slice_max)) {
				type_i = 2; // the plane is an extra slice
			}

			const CGAL_Plane & P_i = planes[i];
			const CGAL_Vector_3 N_i = P_i.orthogonal_vector();

			for (int j = i + 1; j < n; j++) {
				if (jin(id_facet_min, j, id_facet_max)) {
					type_j = 0;
				} else if (jin(id_primitive_min, j, id_primitive_max)) {
					type_j = 1;
				} else if (jin(id_slice_min, j, id_slice_max)) {
					type_j = 2;
				}

				// The expected relationship between planes P_i and P_j (couples listed before)
				// can be expressed through a simple relationship.

				if (type_i == 0 || type_j == 2) {
					const CGAL_Plane & P_j = planes[j];
					const CGAL_Vector_3 N_j = P_j.orthogonal_vector();

					CGAL::cpp11::result_of<K::Intersect_3(CGAL_Plane, CGAL_Plane)>::type object = intersection(P_i, P_j);
					if (object) {
						if (const CGAL_Line_3* ptr_L_ij = boost::get<CGAL_Line_3>(&*object)) {
							CGAL_Line_3 L_ij = (*ptr_L_ij);
							lines[i][j] = L_ij;
							lines[j][i] = L_ij;
						}
					}
				}
			}
		}
	}



	void Kinetic_Propagation_Multiple::init_intersection_lines_partially(const std::vector<int> & indices,
		std::map<int, std::map<int, CGAL_Line_3> > & lines)
	{
		int n = int(indices.size());

		// We compute the intersection of each couple of planes (P_i, P_j).
		// If this intersection exists, and is a line, then it is added to the map of intersection lines.

		for (int i = 0; i < n; i++) {
			int ind_i = indices[i];
			const CGAL_Plane & P_i = planes[ind_i];
			const CGAL_Vector_3 N_i = P_i.orthogonal_vector();

			for (int j = i + 1; j < n; j++) {
				int ind_j = indices[j];
				
				bool L_ij_exists = false;
				if (lines.find(ind_i) != lines.end() && lines[ind_i].find(ind_j) != lines[ind_i].end()) {
					L_ij_exists = true;
				}

				if (L_ij_exists) {
					continue;
				} else {
					const CGAL_Plane & P_j = planes[ind_j];
					const CGAL_Vector_3 N_j = P_j.orthogonal_vector();

					CGAL::cpp11::result_of<K::Intersect_3(CGAL_Plane, CGAL_Plane)>::type object = intersection(P_i, P_j);
					if (object) {
						if (const CGAL_Line_3* ptr_L_ij = boost::get<CGAL_Line_3>(&*object)) {
							CGAL_Line_3 L_ij = (*ptr_L_ij);
							lines[ind_i][ind_j] = L_ij;
							lines[ind_j][ind_i] = L_ij;
						}
					}
				}
			}
		}
	}


	/*void Kinetic_Propagation_Multiple::project_primitives_to_their_respective_planes()
	{
		size_t nb_polygons = polygons.size();
		projected_primitives.reserve(nb_polygons);

		for (size_t pol_id = 0 ; pol_id < nb_polygons ; ++pol_id) {

			// For each input primitive, gets its 3D definition and its assigned plane
			const std::vector<CGAL_Point_3> & source = polygons[pol_id];
			const CGAL_Plane & assigned_plane = planes[polygons_to_planes[pol_id]];

			// Maps 3D points to 2D points included in the assigned plane
			std::vector<CGAL_Point_2> target;
			size_t size_polygon = source.size();
			target.reserve(size_polygon);

			for (size_t i = 0 ; i < size_polygon ; ++i) {
				target.push_back(assigned_plane.to_2d(source[i]));
			}

			// Adds the 2D polygon to the vector of projected primitives
			projected_primitives.push_back(target);
		}
	}*/



	void Kinetic_Propagation_Multiple::init_current_bounding_box(const int i, const int j, const int k,
		const int N,
		std::vector<int> & bounding_planes_indices,
		CGAL_Point_3 & pt_min, CGAL_Point_3 & pt_max,
		std::vector<CGAL_Point_3> & box_corners,
		std::vector<std::pair<int, int> > & box_edges)
	{
		int gx = Universe::params->grid_x;
		int gy = Universe::params->grid_y;
		int gz = Universe::params->grid_z;


		// Step 1.
		// Identifies the indices of the planes that delimit the bounding box (i, j, k).
		// Reminder about the structure of the vector 'planes' :

		// [0 .. 5] : main bounding box facets
		// [6 .. N - 1] : support planes of planar primitives

		// [N .. N + sx - 1] : slicing planes along the axis (Ox), sx = gx - 1
		// [N + sx .. N + sx + sy - 1] : slicing planes along the axis (Oy), sy = gy - 1
		// [N + sx + sy .. N + sx + sy + sz - 1] : slicing planes along the axis (Oz), sz = gz - 1

		bounding_planes_indices = std::vector<int>(6, 0);
		
		// Finds an index of bounding box facet, or slicing plane, along the axis (Ox)
		bounding_planes_indices[0] = (i == 0 ?      0 : N + i);
		bounding_planes_indices[1] = (i == gx - 1 ? 1 : N + i + 1);

		// Same for (Oy)
		bounding_planes_indices[2] = (j == 0 ?      2 : N + (gx - 1) + j);
		bounding_planes_indices[3] = (j == gy - 1 ? 3 : N + (gx - 1) + j + 1);

		// Same for (Oz)
		bounding_planes_indices[4] = (k == 0 ?      4 : N + (gx - 1) + (gy - 1) + k);
		bounding_planes_indices[5] = (k == gz - 1 ? 5 : N + (gx - 1) + (gy - 1) + k + 1);


		// Step 2.
		// Computes the corners of the bounding box.
		
		FT x_min = -planes[bounding_planes_indices[0]].d();
		FT x_max = -planes[bounding_planes_indices[1]].d();
		FT y_min = -planes[bounding_planes_indices[2]].d();
		FT y_max = -planes[bounding_planes_indices[3]].d();
		FT z_min = -planes[bounding_planes_indices[4]].d();
		FT z_max = -planes[bounding_planes_indices[5]].d();

		pt_min = CGAL_Point_3(x_min, y_min, z_min);
		pt_max = CGAL_Point_3(x_max, y_max, z_max);

		int vertices[8][3] = { {0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {1, 1, 0}, {0, 0, 1}, {1, 0, 1}, {0, 1, 1}, {1, 1, 1} };
		int edges[12][2] = { {0, 2}, {1, 3}, {4, 6}, {5, 7}, {0, 1}, {2, 3}, {4, 5}, {6, 7}, {0, 4}, {2, 6}, {1, 5}, {3, 7} };

		for (int i = 0; i < 8; i++) {
			FT x = (vertices[i][0] == 0 ? pt_min.x() : pt_max.x());
			FT y = (vertices[i][1] == 0 ? pt_min.y() : pt_max.y());
			FT z = (vertices[i][2] == 0 ? pt_min.z() : pt_max.z());
			box_corners.push_back(CGAL_Point_3(x, y, z));
		}

		for (int i = 0; i < 12; i++) {
			box_edges.push_back(std::make_pair(edges[i][0], edges[i][1]));
		}
	}



	void Kinetic_Propagation_Multiple::exhibit_considered_inner_planes(const int i, const int j, const int k,
		const int N,
		std::vector<int> & bounding_planes_indices,
		const CGAL_Point_3 & pt_min,
		const CGAL_Point_3 & pt_max,
		const std::vector<CGAL_Point_3> & box_corners,
		const std::vector<std::pair<int, int> > & box_edges,
		std::vector<Support_Plane_Basic_Info> & considered_planes,
		std::vector<int> & all_planes_indices)
	{
		std::list<Support_Plane_Basic_Info> res;
		std::copy(bounding_planes_indices.begin(), bounding_planes_indices.end(), std::back_inserter(all_planes_indices));

		for (int pl_id = 6 ; pl_id <= N ; ++pl_id) {
			const CGAL_Plane & P = planes[pl_id];

			// Let us check if one of the primitives assigned to the plane P
			// should be taken into account in the subvolume (i, j, k).

			std::list<int> PR;
			for (int pol_id : planes_to_polygons[pl_id]) {
				std::set<std::tuple<int, int, int> >::const_iterator it_begin = primitives_to_subvolumes[pol_id].begin();
				std::set<std::tuple<int, int, int> >::const_iterator it_end = primitives_to_subvolumes[pol_id].end();
				if (std::find(it_begin, it_end, std::make_tuple(i, j, k)) != it_end) {
					PR.push_back(pol_id);
				}
			}

			if (PR.empty()) continue;

			// If this plane must be taken into account,
			// then we compute its boundaries

			std::list<CGAL_Point_3> bounding_polygon_3d;
			std::vector<std::list<int> > bounding_facets_indices;

			Support_Plane::construct_bounding_polygon_of_support_plane(pt_min, pt_max, box_corners, box_edges,
				P, bounding_polygon_3d, bounding_facets_indices);

			res.push_back(Support_Plane_Basic_Info(pl_id, bounding_facets_indices, bounding_polygon_3d, PR));
			all_planes_indices.push_back(pl_id);
		}

		considered_planes = std::vector<Support_Plane_Basic_Info>(res.begin(), res.end());
	}



	void Kinetic_Propagation_Multiple::get_subset_of_3d_lines(const std::map<int, std::map<int, CGAL_Line_3> > & lines,
		std::map<int, std::map<int, CGAL_Line_3> > & subset)
	{
		// This function is called as all support planes associated to the current bounding box have just been created.
		// Contrary to the simple version of the kinetic partitioning algorithm,
		// the indices listed in the object 'lines' do not correspond to the indices of the support planes.
		// That's why we create a sub-table of lines adapted to the current partition.

		int n = int(Universe::map_of_planes.size());

		for (int i = 0 ; i < n ; ++i) {
			Support_Plane* SP_i = Universe::map_of_planes[i];
			int p_i = SP_i->real_id;

			for (int j = i + 1 ; j < n ; ++j) {
				Support_Plane* SP_j = Universe::map_of_planes[j];
				int p_j = SP_j->real_id;

				// Tests if there exists a line Int(SP_i, SP_j)
				auto it_l1 = lines.find(p_i);
				if (it_l1 != lines.end()) {
					auto it_l2 = it_l1->second.find(p_j);
					if (it_l2 != it_l1->second.end()) {
						const CGAL_Line_3 & L_ij = it_l2->second;
						subset[i][j] = L_ij;
						subset[j][i] = L_ij;
					}
				}
			}
		}
	}



	void Kinetic_Propagation_Multiple::init_current_support_planes(const int i, const int j, const int k,
		const int N,
		std::vector<int> & bounding_planes_indices,
		const CGAL_Point_3 & pt_min,
		const CGAL_Point_3 & pt_max,
		const std::vector<CGAL_Point_3> & box_corners,
		const std::vector<std::pair<int, int> > & box_edges,
		const std::map<int, std::map<int, CGAL_Line_3> > & lines,
		const std::vector<Support_Plane_Basic_Info> & considered_inner_planes)
	{
		int id_vertices[6][4] = { {0, 4, 6, 2}, {1, 5, 7, 3}, {1, 5, 4, 0}, {3, 7, 6, 2}, {1, 0, 2, 3}, {5, 4, 6, 7} };
		int id_facets[6][4] = { {2, 5, 3, 4}, {2, 5, 3, 4}, {1, 5, 0, 4}, {1, 5, 0, 4}, {2, 0, 3, 1}, {2, 0, 3, 1} };

		// When the multiple mode is turned on, and similarly with a simple kinetic propagation, 
		// There are two different types of support planes to build :
		// - The first category, from planes[0] to planes[5], corresponds to the facets of the current bounding box;
		// - The second category, corresponds to the planes whose assigned primitives intersect the current volume.


		// Step 1.
		// Over a first phase, all the Support_Plane objects are constructed using their plane equations.
		// We first process the bounding box facets, then the inner planes.

		Counters::id_objects = std::vector<int>(planes.size(), -1);

		for (int i = 0 ; i < 6 ; ++i) {
			int ind_i = bounding_planes_indices[i];
			new Support_Plane(ind_i, planes[ind_i], colors[ind_i]);
		}

		for (size_t i = 0 ; i < considered_inner_planes.size() ; ++i) {
			int ind_i = considered_inner_planes[i].id;
			new Support_Plane(ind_i, planes[ind_i], colors[ind_i]);
		}


		// Step 2.
		// Once all the Support_Planes exist, and their respective frames are initialized, 
		// we construct the Intersection_Lines and the bounding polygons.

		std::map<int, std::map<int, CGAL_Line_3> > subset_lines;
		get_subset_of_3d_lines(lines, subset_lines);

		size_t n = Universe::map_of_planes.size();

		for (size_t i = 0 ; i < n ; ++i) {
			Support_Plane* SP = Universe::map_of_planes[i];

			// As for the Intersection_Lines, their 3D versions can be obtained very easily
			std::map<int, std::map<int, CGAL_Line_3> >::const_iterator it_lines = subset_lines.find(i);
			const std::map<int, CGAL_Line_3> & L = it_lines->second;
			SP->init_intersection_lines(L);
			SP->init_polygon_set();

			// As for the bounding polygons, there are different ways to obtain them.
			// planes[0 .. 5] correspond to facets of the bounding box.
			// planes[6 ..] are planes associated to polygons, we run a specific algorithm to obtain BP.
			std::list<CGAL_Point_3> BP;
			std::vector<std::list<int> > BF;

			if (i < 6) {
				for (int j = 0; j < 4; j++) {
					BP.push_back(box_corners[id_vertices[i][j]]);
					BF.push_back(std::list<int>(1, id_facets[i][j]));
				}
			} else {
				BP = considered_inner_planes[i - 6].BP_3;
				BF = considered_inner_planes[i - 6].BF;
			}

			SP->init_bounding_polygon(BP, BF);
		}
	}



	void Kinetic_Propagation_Multiple::build_polygons(const std::vector<Support_Plane_Basic_Info> & considered_planes) const
	{
		// We loop on the inner support planes of the kinetic system
		// and build the polygons that propagate on these planes.

		size_t n = Universe::map_of_planes.size();

		for (size_t i = 6 ; i < n ; ++i) {
			Support_Plane* SP = Universe::map_of_planes[i];
			int id = SP->real_id;

			// We loop on the primitives that are assigned to SP and build them
			const std::list<int> & PR = considered_planes[i - 6].PR;
			for (int pol_id : PR) {
				const CGAL_Point_2 & B = primitives_to_barycenters[pol_id];
				const std::vector<std::pair<CGAL_Point_2, CGAL_Vector_2> > & dirs = primitives_to_directions[pol_id];
				SP->init_polygon(B, dirs);
			}
		}
	}



	void Kinetic_Propagation_Multiple::delete_multiple_kinetic_data_structure()
	{
		clock_t t_0 = clock();

		// Deletes the support planes and their objects
		for (size_t i = 0; i < Universe::map_of_planes.size(); i++) {
			delete Universe::map_of_planes[i];
		}
		Universe::map_of_planes.clear();

		// Deletes the queue
		delete Universe::event_queue;

		// Resets indices
		Counters::id_planes = -1;
		Counters::id_objects = std::vector<int>();

		// Forces the system to compute an appropriate tau at next execution
		Universe::params->D_is_set = false;

		KP_Stats::destruction_time += double(clock() - t_0);
	}
}