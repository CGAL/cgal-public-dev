#ifndef CGAL_LEVEL_OF_DETAIL_BUILDING_SPLITTER_2_H
#define CGAL_LEVEL_OF_DETAIL_BUILDING_SPLITTER_2_H

#if defined(WIN32) || defined(_WIN32) 
#define PSR "\\" 
#else 
#define PSR "/" 
#endif 

// STL includes.
#include <map>
#include <vector>
#include <cassert>
#include <iostream>

// CGAL includes.
#include <CGAL/IO/Color.h>
#include <CGAL/Random.h>
#include <CGAL/Unique_hash_map.h>
#include <CGAL/number_utils.h>

// New CGAL includes.
#include <CGAL/Mylog/Mylog.h>
#include <CGAL/Level_of_detail_enum.h>

namespace CGAL {

	namespace LOD {

		template<class KernelTraits, class CDTInput>
		class Level_of_detail_building_splitter_2 {

		public:
			typedef KernelTraits Kernel;
			typedef CDTInput     CDT;

			typedef typename Kernel::FT 	   FT;
			typedef typename Kernel::Point_2   Point_2;
			typedef typename Kernel::Segment_2 Segment_2;
			typedef typename Kernel::Line_2    Line_2;
			
			typedef CGAL::Color Color;

			typedef typename CDT::Face_handle 			Face_handle;
			typedef typename CDT::Finite_faces_iterator Face_iterator;
			typedef typename CDT::Vertex_handle 		Vertex_handle;
			typedef typename CDT::Edge 					Edge;

			// Extra.
			using Building  = CGAL::LOD::Building<FT, Vertex_handle, Face_handle>;
			using Buildings = std::map<int, Building>;

			using Log = CGAL::LOD::Mylog;
			using Segments = std::vector<Segment_2>;

			struct Building_data {
			
			public:
				int   index = -1;
				Color color = Color(192, 192, 192);
			};

			Level_of_detail_building_splitter_2() : m_silent(false) { }

			void make_silent(const bool new_state) {
				m_silent = new_state;
			}

			void use_custom_constraints(const bool new_state) {
				m_use_custom_constraints = new_state;
			}

			int split(CDT &cdt, Buildings &buildings, const Segments &segments) { 
				
				buildings.clear();

				size_t count = 0;
				CGAL::Unique_hash_map<Face_handle, int> F;

				Building_data bd; 
				for (Face_iterator fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit) {
				
					generate_new_building(bd);

					Face_handle fh = static_cast<Face_handle>(fit);
					flood(cdt, bd, segments, fh, buildings);

					F[fh] = count++;
				}

				const int number_of_buildings = static_cast<int>(buildings.size());
				assert(number_of_buildings >= 0);

				if (!m_silent) {
					Log log;
					log.save_cdt_ply(cdt, "tmp" + std::string(PSR) + "buildings", "bu");
				}

				// Remove later.
				/*
				log.clear();
				log.save_buildings_info(cdt, buildings, "tmp" + std::string(PSR) + "buildings_info"); */

				return number_of_buildings;
			}

		private:
			CGAL::Random m_rand;
			bool m_silent;
			bool m_use_custom_constraints;

			void generate_new_building(Building_data &bd) {

				++bd.index;
				bd.color = generate_random_color();
			}

			Color generate_random_color() {

				const FT r = m_rand.get_int(0, 256);
				const FT g = m_rand.get_int(0, 256);
				const FT b = m_rand.get_int(0, 256);

				return Color(CGAL::to_double(r), CGAL::to_double(g), CGAL::to_double(b));
			}

			void flood(const CDT &cdt, const Building_data &bd, const Segments &segments, Face_handle &fh, Buildings &buildings) const {
				
				// If this face is not valid due to some criteria, we do not handle this face and skip it.
				if (is_not_valid_face(cdt, fh)) return;


				// If this face is valid, we handle it below.
				update_face_and_buildings(bd, fh, buildings);
				for (size_t i = 0; i < 3; ++i) {

					Face_handle fhn = fh->neighbor(i);
					if (!is_constrained_edge(cdt, fh, fhn, segments)) 
						flood(cdt, bd, segments, fhn, buildings);
				}
			}

			bool is_not_valid_face(const CDT &cdt, const Face_handle &fh) const {

				// This face was already handled before.
				const int building_number = fh->info().bu;
				if (building_number >= 0) return true;


				// This is infinite face.
				if (cdt.is_infinite(fh)) return true;


				// This face is outside any building.
				const FT inside_probability = fh->info().in;
				const FT half = FT(1) / FT(2);
				
				assert(inside_probability != half); // if it is equal half, then most likely I use the wrong cdt as input! Comment it out if no need!
				if (!(inside_probability > half)) return true;


				// The face is valid.
				return false;
			}

			void update_face_and_buildings(const Building_data &bd, Face_handle &fh, Buildings &buildings) const {

				fh->info().bu 		= bd.index;
				fh->info().bu_color = bd.color;

				buildings[bd.index].faces.push_back(fh);
				buildings[bd.index].color = bd.color;
			}

			// This function can be improved!
			bool is_constrained_edge(const CDT &cdt, const Face_handle &fh, const Face_handle &fhn, const Segments &segments) const {

				const Vertex_handle vh = find_vertex_handle(fh, fhn);
				const int vertex_index = fh->index(vh); 			  // should be "i" from the function flood() above!

				const Edge edge = std::make_pair(fh, vertex_index);
				if (!m_use_custom_constraints) return cdt.is_constrained(edge);

				const Point_2 &p1 = fh->vertex((vertex_index + 1) % 3)->point();
				const Point_2 &p2 = fh->vertex((vertex_index + 2) % 3)->point();

				return check_for_constraint(p1, p2, segments);
			}

			Vertex_handle find_vertex_handle(const Face_handle &fh, const Face_handle &fhn) const {
				
				for (size_t i = 0; i < 3; ++i) {

					bool found = false;	
					for (size_t j = 0; j < 3; ++j) {

						if (fh->vertex(i) == fhn->vertex(j)) {
							found = true;
							break;
						}
					}
					if (!found) return fh->vertex(i);
				}

				assert(!"Wrong vertex handle!");
				return Vertex_handle();
			}

			bool check_for_constraint(const Point_2 &p1, const Point_2 &p2, const Segments &segments) const {
				
				assert(segments.size() > 0);
				for (size_t i = 0; i < segments.size(); ++i) {
					
					const Segment_2 &segment = segments[i];
					if (is_constrained(p1, p2, segment)) return true;
				}

				return false;
			}

			bool is_constrained(const Point_2 &p1, const Point_2 &p2, const Segment_2 &segment) const {

				Segment_2 edge(p1, p2);

				const Point_2 &source = segment.source();
				const Point_2 &target = segment.target();

				Line_2 line(source, target);
				
				const Point_2 pr1 = line.projection(p1);
				const Point_2 pr2 = line.projection(p2);

				const FT eps = FT(1);

				if (squared_distance(p1, pr1) > eps * eps) return false;
				if (squared_distance(p2, pr2) > eps * eps) return false;

				const FT tol = -FT(1) / FT(10);
				
				std::pair<FT, FT> bc = BC::compute_segment_coordinates_2(source, target, p1, Kernel());
				const bool state1 = bc.first > tol && bc.second > tol;

				bc = BC::compute_segment_coordinates_2(source, target, p2, Kernel());
				const bool state2 = bc.first > tol && bc.second > tol;

				bc = BC::compute_segment_coordinates_2(p1, p2, source, Kernel());
				const bool state3 = bc.first > tol && bc.second > tol;

				bc = BC::compute_segment_coordinates_2(p1, p2, target, Kernel());
				const bool state4 = bc.first > tol && bc.second > tol;

				if ( (state1 && state2) || (state3 && state4) ) return true;
				return false;
			}
		};
	}
}

#endif // CGAL_LEVEL_OF_DETAIL_BUILDING_SPLITTER_2_H