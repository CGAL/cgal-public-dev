#ifndef CGAL_LEVEL_OF_DETAIL_RECONSTRUCTION_1_H
#define CGAL_LEVEL_OF_DETAIL_RECONSTRUCTION_1_H

// STL includes.
#include <vector>
#include <iostream>
#include <cassert>
#include <string>
#include <map>

// CGAL includes.
#include <CGAL/IO/Color.h>

// New CGAL includes.
#include <CGAL/Mylog/Mylog.h>
#include <CGAL/Level_of_detail_enum.h>

namespace CGAL {

	namespace LOD {

		template<class KernelTraits, class CDTInput, class BuildingsInput, class MeshOutput>
		class Level_of_detail_reconstruction_1 {

		public:
			typedef KernelTraits   Kernel;
			typedef CDTInput 	   CDT;
			typedef BuildingsInput Buildings;
			typedef MeshOutput 	   Mesh;

			// Extra.
			using Log = CGAL::LOD::Mylog;

			Level_of_detail_reconstruction_1() : m_save_as_ply(true) { }

			void save_as_ply(const bool new_state) {
				m_save_as_ply = new_state;
			}

			void reconstruct(const CDT &cdt, const Buildings &buildings, Mesh &) const {

				if (m_save_as_ply) save_as_ply(cdt, buildings);
				else {

					// use an incremental builder here!
				}
			}

		private:
			bool m_save_as_ply;

			// Save mesh as ply format.
			void save_as_ply(const CDT &, const Buildings &) const {

				Log log;

				add_header(log);
				add_vertices(log);
				add_faces(log);

				log.save("LOD1", ".ply");
			}

			void add_header(Log &log) const {

				log.out << 
				"ply\n"               					 << 
				"format ascii 1.0\n"     				 << 
				"element vertex "        				 << 5 << "\n" << 
				"property double x\n"    				 << 
				"property double y\n"    				 << 
				"property double z\n" 					 <<
				"element face " 						 << 2 << "\n" << 
				"property list uchar int vertex_index\n" <<
				"property uchar red\n" 					 <<
				"property uchar green\n" 				 <<
				"property uchar blue\n" 				 <<
				"end_header\n";
			}

			void add_vertices(Log &log) const {

				log.out << 0.0 << " " << 0.0 << " " << 0.0 << std::endl;
				log.out << 0.5 << " " << 1.0 << " " << 0.0 << std::endl;
				log.out << 1.0 << " " << 0.0 << " " << 0.0 << std::endl;
				log.out << 2.0 << " " << 1.0 << " " << 0.0 << std::endl;
				log.out << 1.0 << " " << 2.0 << " " << 0.0 << std::endl;
			}

			void add_faces(Log &log) const {

				log.out << 3 << " " << 0 << " " << 1 << " " << 2 << " " << 169 << " " << 0 << " " << 0 << std::endl;
				log.out << 4 << " " << 1 << " " << 3 << " " << 4 << " " << 2 << " " << 0 << " " << 169 << " " << 0 << std::endl;
			}
		};
	}
}

#endif // CGAL_LEVEL_OF_DETAIL_RECONSTRUCTION_1_H