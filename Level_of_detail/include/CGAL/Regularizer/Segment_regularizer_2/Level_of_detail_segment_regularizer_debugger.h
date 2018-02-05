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

			template<typename Segments>
			void print_segments(const Segments &segments, const std::string &name) {
                
				assert(segments.size() > 0); 
                clear();

				for (typename Segments::const_iterator segment = segments.begin(); segment != segments.end(); ++segment) {

					out << "v " << segment->get().source() << " " << 0 << std::endl;
					out << "v " << segment->get().target() << " " << 0 << std::endl;
					out << "v " << segment->get().target() << " " << 0 << std::endl;
                    out << "v " << segment->get().source() << " " << 0 << std::endl;
				}

				for (size_t i = 0; i < segments.size() * 4; i += 4)
					out << "f " << i + 1 << " " << i + 2 << " " << i + 3 << " " << i + 4 << std::endl;

				save(name, ".obj");
            }

			template<typename Values>
            void print_values(const Values &values, const std::string &name) const {

                std::cout << std::endl << name << " size: " << values.size() << std::endl;
                for (typename Values::const_iterator value = values.begin(); value != values.end(); ++value) std::cout << *value << std::endl;
                std::cout << std::endl;
            }

        private:
            std::stringstream out;
			const std::string m_prefix_path;

			inline std::string data() const {
				return out.str();
			}

			void save(const std::string &fileName, const std::string &extension = ".log", const std::string ending = ("logs" + std::string(PS) + "tmp" + std::string(PS))) const {
				const std::string default_path = m_prefix_path + ending;

				const std::string finalPath = default_path + fileName + extension;
				std::ofstream file(finalPath.c_str(), std::ios_base::out);

				if (!file) std::cerr << std::endl << "ERROR: Error saving log file with the name " << fileName << std::endl << std::endl;

				file << data() << std::endl;
				file.close();
			}
		};
	}
}

#endif // CGAL_LEVEL_OF_DETAIL_SEGMENT_REGULARIZER_LOG_H