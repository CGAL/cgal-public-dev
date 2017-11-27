#ifndef CGAL_LEVEL_OF_DETAIL_BUILDING_SPLITTER_2_H
#define CGAL_LEVEL_OF_DETAIL_BUILDING_SPLITTER_2_H

#if defined(WIN32) || defined(_WIN32) 
#define PS "\\" 
#else 
#define PS "/" 
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

			typedef typename Kernel::FT FT;
			typedef CGAL::Color Color;

			typedef typename CDT::Face_handle 			Face_handle;
			typedef typename CDT::Finite_faces_iterator Face_iterator;
			typedef typename CDT::Vertex_handle 		Vertex_handle;
			typedef typename CDT::Edge 					Edge;

			// Extra.
			using Building  = CGAL::LOD::Building<FT, Vertex_handle, Face_handle>;
			using Buildings = std::map<int, Building>;

			using Log = CGAL::LOD::Mylog;

			struct Building_data {
			
			public:
				int   index = -1;
				Color color = Color(192, 192, 192);
			};

			Level_of_detail_building_splitter_2() { }

			int split(CDT &cdt, Buildings &buildings) { 
				
				buildings.clear();

				size_t count = 0;
				CGAL::Unique_hash_map<Face_handle, int> F;

				Building_data bd; 
				for (Face_iterator fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit) {
				
					generate_new_building(bd);

					Face_handle fh = static_cast<Face_handle>(fit);
					flood(cdt, bd, fh, buildings);

					F[fh] = count++;
				}

				const int number_of_buildings = static_cast<int>(buildings.size());
				assert(number_of_buildings >= 0);


				// Remove later.
				/*
				Log log;
				log.save_cdt_ply(cdt, "tmp" + std::string(PS) + "buildings", "bu");

				log.clear();
				log.save_buildings_info(cdt, buildings, "tmp" + std::string(PS) + "buildings_info"); */

				return number_of_buildings;
			}

		private:
			CGAL::Random m_rand;

			void generate_new_building(Building_data &bd) {

				++bd.index;
				bd.color = generate_random_color();
			}

			Color generate_random_color() {

				const FT r = m_rand.get_int(0, 256);
				const FT g = m_rand.get_int(0, 256);
				const FT b = m_rand.get_int(0, 256);

				return Color(r, g, b);
			}

			void flood(const CDT &cdt, const Building_data &bd, Face_handle &fh, Buildings &buildings) const {
				
				// If this face is not valid due to some criteria, we do not handle this face and skip it.
				if (is_not_valid_face(cdt, fh)) return;


				// If this face is valid, we handle it below.
				update_face_and_buildings(bd, fh, buildings);
				for (size_t i = 0; i < 3; ++i) {

					Face_handle fhn = fh->neighbor(i);
					if (!is_constrained_edge(cdt, fh, fhn)) flood(cdt, bd, fhn, buildings);
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
			bool is_constrained_edge(const CDT &cdt, const Face_handle &fh, const Face_handle &fhn) const {

				const Vertex_handle vh = find_vertex_handle(fh, fhn);
				const int vertex_index = fh->index(vh); 			  // should be "i" from the function flood() above!

				const Edge edge = std::make_pair(fh, vertex_index);
				return cdt.is_constrained(edge);
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
		};
	}
}

#endif // CGAL_LEVEL_OF_DETAIL_BUILDING_SPLITTER_2_H