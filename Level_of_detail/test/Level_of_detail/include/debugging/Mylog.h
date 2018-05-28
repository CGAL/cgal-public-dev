#ifndef CGAL_LEVEL_OF_DETAIL_MYLOG_H
#define CGAL_LEVEL_OF_DETAIL_MYLOG_H

// STL includes.
#include <vector>
#include <string>
#include <sstream>
#include <fstream>

namespace CGAL {

	namespace Level_of_detail {

		class Mylog {

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
        };

    } // Level_of_detail

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_MYLOG_H