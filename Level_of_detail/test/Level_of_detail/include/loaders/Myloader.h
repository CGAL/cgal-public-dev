#ifndef CGAL_LEVEL_OF_DETAIL_MYLOADER_H
#define CGAL_LEVEL_OF_DETAIL_MYLOADER_H

// STL includes.
#include <string>
#include <cassert>
#include <fstream>
#include <iostream>

// Boost includes.
#include <boost/tuple/tuple.hpp>

// CGAL includes.
#include <CGAL/array.h>
#include <CGAL/Point_set_3.h>

namespace CGAL {

	namespace Level_of_detail {

		template<class InputKernel>
		class Myloader {
		
		public:
			using Kernel = InputKernel;

			using FT     = typename Kernel::FT;
			using Point  = typename Kernel::Point_3;
			using Normal = typename Kernel::Vector_3;

			using Container = CGAL::Point_set_3<Point>;

			using Uchar = unsigned char;
			using Color = CGAL::cpp11::array<Uchar, 3>;
			using Label = int;

			using Point_map  = typename Container::Point_map;
			using Normal_map = typename Container::Vector_map;

			using Label_map  = typename Container:: template Property_map<Label>;
			using Color_map  = typename Container:: template Property_map<Color>;
			
			using Iterator = typename Container::const_iterator;

			void get_data(const std::string &file_path, Container &input) const {
            	
				// We use this to fix wrong order of the given input data in the file.
            	// std::ofstream saver((file_path + ".new").c_str(), std::ios_base::out); 
				// saver.precision(20);

				std::ifstream loader(file_path.c_str(), std::ios_base::in);
            	if (!loader) {

                	std::cerr << std::endl << std::endl << "ERROR: Error loading file with LOD data!" << std::endl << std::endl;
                	exit(EXIT_FAILURE);
            	}

				Color_map colors; Label_map labels;
				set_default_properties(input, colors, labels);

            	std::string stub;
            	std::getline(loader, stub);
            	std::getline(loader, stub);
            	std::getline(loader, stub);
            	std::getline(loader, stub);
            	std::getline(loader, stub);
            	std::getline(loader, stub);
            	std::getline(loader, stub);
            	std::getline(loader, stub);

            	size_t num_points;
            	loader >> stub >> stub >> num_points;

            	for (size_t i = 0; i < 13; ++i) std::getline(loader, stub);

				FT x, y, z, nx, ny, nz;
				int r, g, b, l;

            	for (size_t i = 0; i < num_points; ++i) {
            		loader >> x >> y >> z >> nx >> ny >> nz >> r >> g >> b >> stub >> l;
            		
					// We use this to fix wrong order of the given input data in the file.
					// saver << x << " " << y << " " << z << " " << nx << " " << ny << " " << nz << " " << r << " " << g << " " << b << " " << stub << " " << l << std::endl;

					const Iterator it = input.insert(Point(x, y, z), Normal(nx, ny, nz));

					 labels[*it] = l;
					 colors[*it] = {{ static_cast<Uchar>(r), static_cast<Uchar>(g), static_cast<Uchar>(b) }};
            	}
            	loader.close();
            	
				// We use this to fix wrong order of the given input data in the file.
            	// saver.close(); 
				// exit(0);
			}

		private:

			void set_default_properties(Container &input, Color_map &colors, Label_map &labels) const {

				bool success = false;
				input.add_normal_map();

				boost::tie(labels , success) = input. template add_property_map<Label>("label", -1);
				assert(success);

				boost::tie(colors , success) = input. template add_property_map<Color>("color", {{ 0, 0, 0 }});
				assert(success);
			}
		};
	
	} // Level_of_detail

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_MYLOADER_H