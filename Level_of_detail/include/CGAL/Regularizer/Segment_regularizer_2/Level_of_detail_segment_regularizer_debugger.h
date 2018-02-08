#ifndef CGAL_LEVEL_OF_DETAIL_SEGMENT_REGULARIZER_DEBUGGER_H
#define CGAL_LEVEL_OF_DETAIL_SEGMENT_REGULARIZER_DEBUGGER_H

#if defined(WIN32) || defined(_WIN32) 
#define PS "\\"
#else 
#define PS "/" 
#endif

// STL includes.
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <cassert>
#include <iostream>

// CGAL includes.
#include <CGAL/Unique_hash_map.h>

// TODO:
// 1. Make all segments small boxes or fix the problem with polylines visualization in the polyhedron demo.

namespace CGAL {

	namespace LOD {

		class Level_of_detail_segment_regularizer_debugger {

		public:

            Level_of_detail_segment_regularizer_debugger() : m_prefix_path(std::string(std::getenv("LOD_LOG_PATH"))) { }

			void clear() {
				out.str(std::string());
			}

			template<typename SegmentRange, typename SegmentMap, typename Kernel>
			inline void export_segments(const SegmentRange &segments, const SegmentMap &segment_map, const std::string &filename) {
				this->print_segments<SegmentRange, SegmentMap, Kernel>(segments, segment_map, filename);
			}

			template<typename SegmentRange, typename SegmentMap, typename Kernel>
			void print_segments(const SegmentRange &segments, const SegmentMap &segment_map, const std::string &filename) {
                
				using Segment 		   = typename Kernel::Segment_2;
				using Segment_iterator = typename SegmentRange::const_iterator;

				clear();
				assert(segments.size() > 0);

				for (Segment_iterator sit = segments.begin(); sit != segments.end(); ++sit) {
					const Segment &segment = get(segment_map, *sit);

					out << "v " << segment.source() << " " << 0 << std::endl;
					out << "v " << segment.target() << " " << 0 << std::endl;
					out << "v " << segment.target() << " " << 0 << std::endl;
                    out << "v " << segment.source() << " " << 0 << std::endl;
				}

				for (size_t i = 0; i < segments.size() * 4; i += 4)
					out << "f " << i + 1 << " " << i + 2 << " " << i + 3 << " " << i + 4 << std::endl;

				save(filename, ".obj");
            }

			template<typename Values>
            void print_values(const Values &values, const std::string &name) const {

                std::cout << std::endl << name << " size: " << values.size() << std::endl;
                for (typename Values::const_iterator value = values.begin(); value != values.end(); ++value) std::cout << *value << std::endl;
                std::cout << std::endl;
            }

			template<typename Points>
			void print_sampled_segments(const Points &points, const std::string &filename) {
				
				assert(points.size() > 0);
				clear();

				for (typename Points::const_iterator point = points.begin(); point != points.end(); ++point)
					out << (*point).first << " " << 0 << std::endl;

				save(filename, ".xyz");
			}

			template<class Triangulation>
			void print_triangulation(const Triangulation &triangulation, const std::string &filename) {

				clear();

				typedef typename Triangulation::Vertex_handle Vertex_handle;
				CGAL::Unique_hash_map<Vertex_handle, int> V;

				int count = 0;
				for (typename Triangulation::Finite_vertices_iterator vit = triangulation.finite_vertices_begin(); vit != triangulation.finite_vertices_end(); ++vit) {
					
					out << "v " << (*vit) << " " << 0 << std::endl;
					V[vit] = count++;
				}

				for (typename Triangulation::Finite_faces_iterator fit = triangulation.finite_faces_begin(); fit != triangulation.finite_faces_end(); ++fit)
					out << "f " << V[(*fit).vertex(0)] + 1 << " " << V[(*fit).vertex(1)] + 1 << " " << V[(*fit).vertex(2)] + 1 << std::endl;

				save(filename, ".obj");
			}

        private:
            std::stringstream out;
			const std::string m_prefix_path;

			inline std::string data() const {
				return out.str();
			}

			void save(const std::string &filename, const std::string &extension = ".log", const std::string ending = ("logs" + std::string(PS) + "tmp" + std::string(PS))) const {
				const std::string default_path = m_prefix_path + ending;

				const std::string final_path = default_path + filename + extension;
				std::ofstream file(final_path.c_str(), std::ios_base::out);

				if (!file) std::cerr << std::endl << "ERROR: Error saving log file with the name " << filename << std::endl << std::endl;

				file << data() << std::endl;
				file.close();
			}
		};
	}
}

#endif // CGAL_LEVEL_OF_DETAIL_SEGMENT_REGULARIZER_LOG_H