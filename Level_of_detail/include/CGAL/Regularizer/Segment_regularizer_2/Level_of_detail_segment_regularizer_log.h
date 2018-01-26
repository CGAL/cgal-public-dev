#ifndef CGAL_LEVEL_OF_DETAIL_SEGMENT_REGULARIZER_LOG_H
#define CGAL_LEVEL_OF_DETAIL_SEGMENT_REGULARIZER_LOG_H

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
#include <iostream>

// TODO:
// 1. Make all segments small boxes or fix the problem with polylines visualization in the polyhedron demo.

namespace CGAL {

	namespace LOD {

		class Level_of_detail_segment_regularizer_log {

		public:

            Level_of_detail_segment_regularizer_log() : m_prefix_path(std::string(std::getenv("LOD_LOG_PATH"))) { }

			void clear() {
				out.str(std::string());
			}

			void save(const std::string &fileName, const std::string &extension = ".log", const std::string ending = ("logs" + std::string(PS) + "tmp" + std::string(PS))) const {
				const std::string default_path = m_prefix_path + ending;

				const std::string finalPath = default_path + fileName + extension;
				std::ofstream file(finalPath.c_str(), std::ios_base::out);

				if (!file) std::cerr << std::endl << "ERROR: Error saving log file with the name " << fileName << std::endl << std::endl;

				file << data() << std::endl;
				file.close();
			}

            template<typename Segments>
            void export_segments(const Segments &segments, const std::string &name) {
                
                clear();
				for (size_t i = 0; i < segments.size(); ++i) {

					out << "v " << segments[i].source() << " " << 0 << std::endl;
					out << "v " << segments[i].target() << " " << 0 << std::endl;
					out << "v " << segments[i].target() << " " << 0 << std::endl;
                    out << "v " << segments[i].source() << " " << 0 << std::endl;
				}

				for (size_t i = 0; i < segments.size() * 4; i += 4)
					out << "f " << i + 1 << " " << i + 2 << " " << i + 3 << " " << i + 4 << std::endl;

				save(name, ".obj");
            }

        private:
            std::stringstream out;
			const std::string m_prefix_path;

			std::string data() const {
				return out.str();
			}
		};
	}
}

#endif // CGAL_LEVEL_OF_DETAIL_SEGMENT_REGULARIZER_LOG_H