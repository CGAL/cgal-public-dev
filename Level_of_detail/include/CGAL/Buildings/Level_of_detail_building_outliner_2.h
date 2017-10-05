#ifndef CGAL_LEVEL_OF_DETAIL_BUILDING_OUTLINER_2_H
#define CGAL_LEVEL_OF_DETAIL_BUILDING_OUTLINER_2_H

// STL includes.
#include <map>
#include <vector>
#include <cassert>

// New CGAL includes.
#include <CGAL/Mylog/Mylog.h>
#include <CGAL/Level_of_detail_enum.h>

namespace CGAL {

	namespace LOD {

		template<class KernelTraits, class CDTInput>
		class Level_of_detail_building_outliner_2 {

		public:
			typedef KernelTraits Kernel;
			typedef CDTInput     CDT;

			typedef typename Kernel::FT FT;

			typedef typename CDT::Vertex_handle Vertex_handle;
			typedef typename CDT::Face_handle   Face_handle;
			typedef typename CDT::Edge 			Edge;

			// Extra.
			using Building  = CGAL::LOD::Building<FT, Vertex_handle, Face_handle>;
			using Buildings = std::map<int, Building>;

			using Building_iterator = typename Buildings::iterator;
			using Log = CGAL::LOD::Mylog;

			Level_of_detail_building_outliner_2() : m_save_info(false), m_max_outer_iters(1000000), m_max_inner_iters(1000) { }

			void save_info(const bool new_state) {
				m_save_info = new_state;
			}

			void set_max_outer_iterations(const size_t iters) {
				m_max_outer_iters = iters;
			}

			void set_max_inner_iterations(const size_t iters) {
				m_max_inner_iters = iters;
			}

			void find_boundaries(const CDT &cdt, Buildings &buildings) const {
				
				Log log;
				for (Building_iterator bit = buildings.begin(); bit != buildings.end(); ++bit) {
					
					if (m_save_info) log.out << "New building:" << std::endl;

					Building &building = (*bit).second;
					handle_building(cdt, building, log);

					if (m_save_info) log.skip_line();
				}
				if (m_save_info) log.save("tmp/outliner_internal_info");
			}

		private:
			bool m_save_info;

			size_t m_max_outer_iters;
			size_t m_max_inner_iters;

			// Work with each building separately.
			void handle_building(const CDT &cdt, Building &building, Log &log) const {

			 	const Face_handle   fh = find_triangle_to_start_traversal(cdt, building, log);
			 	const Vertex_handle vh = find_vertex_to_start_traversal(cdt, fh, log);

			 	get_boundary(cdt, fh, vh, building, log);
			}

			// Use statistics to guess the best boundary triangle to start traversal of the boundary.
			Face_handle find_triangle_to_start_traversal(const CDT &cdt, const Building &building, Log &log) const {
				
				Face_handle result;
				const size_t num_faces = building.faces.size();

				FT best_face_cost = FT(0);
				for (size_t i = 0; i < num_faces; ++i) {
					
					const Face_handle &fh = building.faces[i];
					const FT face_cost = compute_face_cost(cdt, fh);

					assert(face_cost != FT(0));
					if (face_cost >= best_face_cost) {

						best_face_cost = face_cost;
						result = fh;
					}
				}

				assert(best_face_cost != FT(0));
				if (m_save_info) log.out << "Start from face: ( " << cdt.triangle(result) << " )" << std::endl;

				return result;
			}

			FT compute_face_cost(const CDT &cdt, const Face_handle &fh) const {

				FT face_cost = FT(0);
				for (size_t i = 0; i < 3; ++i) {

					const Face_handle &fhn = fh->neighbor(i);
					const Edge edge = std::make_pair(fh, i);

					face_cost += compute_edge_cost(cdt, edge, fhn);
				}

				assert(face_cost != FT(0));
				return face_cost;
			}

			// Use statistics to guess the best vertex of the triangle to start traversal of the boundary.
			Vertex_handle find_vertex_to_start_traversal(const CDT &cdt, const Face_handle &fh, Log &log) const {

				Vertex_handle result;
				FT best_vertex_cost = FT(0);

				for (size_t i = 0; i < 3; ++i) {

					const Vertex_handle &vh = fh->vertex(i);
					const FT vertex_cost = compute_vertex_cost(cdt, fh, vh);
					
					assert(vertex_cost != FT(0));
					if (vertex_cost >= best_vertex_cost) {

						best_vertex_cost = vertex_cost;
						result = vh;
					}
				}

				assert(best_vertex_cost != FT(0));
				if (m_save_info) log.out << "Start from vertex: ( " << result->point() << " )" << std::endl;

				return result;
			}

			FT compute_vertex_cost(const CDT &cdt, const Face_handle &fh, const Vertex_handle &vh) const {

				FT vertex_cost = FT(0);

				vertex_cost += compute_label_cost(vh);
				vertex_cost += compute_neighbouring_edges_cost(cdt, fh, vh);

				assert(vertex_cost != FT(0));
				return vertex_cost;
			}

			FT compute_label_cost(const Vertex_handle &vh) const {

				const FT half = FT(1) / FT(2);

				const FT corner_cost  = FT(4);
				const FT linear_cost  = FT(2);
				const FT clutter_cost = FT(1);
				const FT default_cost = half;

				FT label_cost = FT(0);
				switch (vh->info().label) {

					case CGAL::LOD::Structured_label::CORNER:
						label_cost += corner_cost;

					case CGAL::LOD::Structured_label::LINEAR:
						label_cost += linear_cost;

					case CGAL::LOD::Structured_label::CLUTTER:
						label_cost += clutter_cost;

					default:
						label_cost += default_cost;
				}

				assert(label_cost != FT(0));
				return label_cost;
			}

			FT compute_neighbouring_edges_cost(const CDT &cdt, const Face_handle &fh, const Vertex_handle &vh) const {

				// Find edges.
				const int curr_index = fh->index(vh);

				const int prev_index =  fh->cw(curr_index);
				const int next_index = fh->ccw(curr_index);

				const Edge prev_edge = std::make_pair(fh, next_index);
				const Edge next_edge = std::make_pair(fh, prev_index);

				const Face_handle &fhn_prev = fh->neighbor(next_index);
				const Face_handle &fhn_next = fh->neighbor(prev_index);

				// Compute cost. 
				FT neighbouring_edges_cost = FT(0);

				neighbouring_edges_cost += compute_edge_cost(cdt, prev_edge, fhn_prev);
				neighbouring_edges_cost += compute_edge_cost(cdt, next_edge, fhn_next);

				assert(neighbouring_edges_cost != FT(0));
				return neighbouring_edges_cost;
			}

			FT compute_edge_cost(const CDT &cdt, const Edge &edge, const Face_handle &fh_neighbour) const {

				const FT half = FT(1) / FT(2);

				const FT outside_cost 	  = FT(4);
				const FT infinite_cost 	  = FT(4);
				const FT constrained_cost = FT(2);
				const FT inside_cost 	  = FT(1);
				const FT unknown_cost 	  = half;

				FT edge_cost = FT(0);

				if (cdt.is_constrained(edge)) 	   edge_cost += constrained_cost;
				if (cdt.is_infinite(fh_neighbour)) edge_cost += infinite_cost;

				if (fh_neighbour->info().in < half)  edge_cost += outside_cost;
				if (fh_neighbour->info().in == half) edge_cost += unknown_cost;
				if (fh_neighbour->info().in > half)  edge_cost += inside_cost;

				assert(edge_cost != FT(0));
				return edge_cost;
			}

			void get_boundary(const CDT &cdt, const Face_handle &fh, const Vertex_handle &vh, Building &building, Log &log) const {

				building.boundaries.clear();
				building.wedges.clear();

				// Temporarily, we use only one boundary here!
				building.boundaries.resize(1);
				building.wedges.resize(1);


				Vertex_handle curr_vh = vh;
				Face_handle   curr_fh = fh;

				// Can be removed later! -->
				if (m_save_info) {

					Edge tmp_edge;
					const Face_handle tmp_fh = get_next_face_and_edge(curr_fh, curr_vh, tmp_edge);

					log.out << "The first neighbouring face: ( ";

					if (cdt.is_infinite(tmp_fh)) log.out << "infinite";
					else log.out << cdt.triangle(tmp_fh);

					log.out << " )" << std::endl;
				} 
				// <--

				// Main loop.
				Edge edge;
				Face_handle curr_face_neigh;

				add_new_boundary_vertex(curr_vh, building);
				add_new_wedge_face(curr_fh, curr_vh, building);

				size_t iter = 0;
				do {

					curr_face_neigh = get_next_face_and_edge(curr_fh, curr_vh, edge);
					if (is_boundary_face(cdt, curr_face_neigh, edge)) {

						curr_vh = get_next_vertex_handle(curr_vh, curr_fh);
						if (curr_vh != vh) {

							add_new_boundary_vertex(curr_vh, building);
							add_new_wedge_face(curr_fh, curr_vh, building);
						}

					} else curr_fh = get_last_face(cdt, curr_face_neigh, curr_vh, building);
					
					++iter;
					assert(iter != m_max_outer_iters);

				} while (curr_vh != vh);
			}

			void add_new_boundary_vertex(const Vertex_handle &curr_vh, Building &building) const {
				building.boundaries[0].push_back(curr_vh);
			}

			void add_new_wedge_face(const Face_handle &curr_fh, const Vertex_handle &curr_vh, Building &building) const {
				building.wedges[0][curr_vh].push_back(curr_fh);
			}

			Vertex_handle get_next_vertex_handle(const Vertex_handle &curr_vh, const Face_handle &curr_fh) const {

				const int curr_index = curr_fh->index(curr_vh);
				const int next_index = curr_fh->ccw(curr_index);

				return curr_fh->vertex(next_index);
			}

			bool is_boundary_face(const CDT &cdt, const Face_handle &fh, const Edge &edge) const {
				
				const FT half = FT(1) / FT(2);

				if (cdt.is_infinite(fh))  							  return true;
				if (fh->info().in < half) 							  return true;
				if (fh->info().in > half && cdt.is_constrained(edge)) return true;

				assert(fh->info().in != half);
				return false;
			}

			Face_handle get_last_face(const CDT &cdt, const Face_handle &fh, const Vertex_handle &curr_vh, Building &building) const {

				Face_handle last_fh = fh;
				Face_handle curr_fh = fh;
				
				size_t iter = 0;
				Edge edge;

				do {

					last_fh = curr_fh;
					add_new_wedge_face(curr_fh, curr_vh, building);
					curr_fh = get_next_face_and_edge(curr_fh, curr_vh, edge);

					++iter;
					assert(iter != m_max_inner_iters);

				} while (!is_boundary_face(cdt, curr_fh, edge));

				return last_fh;
			}

			Face_handle get_next_face_and_edge(const Face_handle &curr_fh, const Vertex_handle &curr_vh, Edge &edge) const {

				const int curr_index = curr_fh->index(curr_vh);
				const int prev_index = curr_fh->cw(curr_index);

				edge = std::make_pair(curr_fh, prev_index);
				return curr_fh->neighbor(prev_index);
			}
		};
	}
}

#endif // CGAL_LEVEL_OF_DETAIL_BUILDING_OUTLINER_2_H