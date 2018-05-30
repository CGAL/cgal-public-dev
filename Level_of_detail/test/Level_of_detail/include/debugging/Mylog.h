#ifndef CGAL_LEVEL_OF_DETAIL_MYLOG_H
#define CGAL_LEVEL_OF_DETAIL_MYLOG_H

#if defined(WIN32) || defined(_WIN32) 
#define _SR_ "\\"
#define _NL_ "\r\n"
#else 
#define _SR_ "/" 
#define _NL_ "\n"
#endif

// STL includes.
#include <vector>
#include <string>
#include <sstream>
#include <fstream>

// CGAL includes.
#include <CGAL/Random.h>
#include <CGAL/IO/Color.h>

namespace CGAL {

	namespace Level_of_detail {

		class Mylog {

		private:
			using Colour = CGAL::Color;
			using Random = CGAL::Random;

			Random m_rand;

		public:
			void clear() {
				out.str(std::string());
			}

            template<class Elements, class Point_map>
            void save_points(const Elements &elements, const Point_map &point_map, const std::string &file_name) {

				clear();
                using Const_elements_iterator = typename Elements::const_iterator;

				for (Const_elements_iterator ce_it = elements.begin(); ce_it != elements.end(); ++ce_it) out << get(point_map, *ce_it) << std::endl;
				save(file_name, ".xyz");
            }

			template<class Vertex_range, class Point_map>
			void save_polygon(const Vertex_range &vertices, const Point_map &point_map, const std::string &file_name) {

				clear();
				using Const_vertices_iterator = typename Vertex_range::const_iterator;

				out << 
				"ply" 				   +  std::string(_NL_) + ""               			  << 
				"format ascii 1.0"     +  std::string(_NL_) + ""     			          << 
				"element vertex "      << vertices.size()  << "" + std::string(_NL_) + "" << 
				"property double x"    +  std::string(_NL_) + ""    			          << 
				"property double y"    +  std::string(_NL_) + ""    			          << 
				"property double z"    +  std::string(_NL_) + "" 				          <<
				"element face "        << 1                << "" + std::string(_NL_) + "" << 
				"property list uchar int vertex_indices"         + std::string(_NL_) + "" <<
				"property uchar red"   +  std::string(_NL_) + "" 				          <<
				"property uchar green" +  std::string(_NL_) + "" 				          <<
				"property uchar blue"  +  std::string(_NL_) + "" 				          <<
				"end_header"           +  std::string(_NL_) + "";

				size_t size = 0;
				for (Const_vertices_iterator cv_it = vertices.begin(); cv_it != vertices.end(); ++cv_it, ++size) 
					out << get(point_map, *cv_it) << std::endl;

				out << size << " ";
				for (size_t i = 0; i < size; ++i) out << i << " ";
				out << generate_random_colour() << std::endl;

				save(file_name, ".ply");
			}

        private:
            std::stringstream out;

			inline std::string data() const {
				return out.str();
			}

			void save(const std::string &file_name, const std::string &extension = ".log") const {

				const std::string final_path = file_name + extension;
				std::ofstream file(final_path.c_str(), std::ios_base::out);

				if (!file) std::cerr << std::endl << "ERROR: Error saving log file with the name " << file_name << std::endl << std::endl;

				file << data() << std::endl;
				file.close();
			}

			Colour generate_random_colour() {

				const int r = m_rand.get_int(0, 256);
				const int g = m_rand.get_int(0, 256);
				const int b = m_rand.get_int(0, 256);

				return Colour(r, g, b);
			}
        };

    } // Level_of_detail

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_MYLOG_H