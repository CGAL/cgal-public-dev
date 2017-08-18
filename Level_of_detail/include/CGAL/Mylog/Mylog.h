#ifndef CGAL_MYLOG_H
#define CGAL_MYLOG_H

// STL includes.
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>

// CGAL includes.
#include <CGAL/array.h>
#include <CGAL/Point_set_3.h>

namespace CGAL {

	namespace LOD {

		class Mylog {

		public:
			std::string state() const {
				return "ok";
			}

			std::string data() const {
				return out.str();
			}

			void clear() {
				out.str(std::string());
			}

			template<typename T>
			void append(const T &token) {
				out << token;
			}

			bool save(const std::string &fileName) const {

				const std::string path = "/Users/danisimo/Documents/pipeline/logs/";
				const std::string finalPath = path + fileName + ".log";

				std::ofstream file(finalPath.c_str(), std::ios_base::out);

				if (!file) {
					std::cerr << "\nERROR: Error saving log file with the name " << fileName << "\n" << std::endl;
					return false;
				}

				file << data() << std::endl;
				file.close();

				return true;
			}

			template<class Traits, class Container>
			void save_ply(const Container &input, 
						  const std::string &fileName,
						  const bool withExtraProperties = false) {

				using Type  = unsigned char;
				using Color = CGAL::cpp11::array<Type, 3>;
				using Label = int;
				using Plane = typename Traits::Plane_3;
				using Index = int;

				using Color_map = typename Container:: template Property_map<Color>;
				using Label_map = typename Container:: template Property_map<Label>;
				using Plane_map = typename Container:: template Property_map<Plane>;
				using Index_map = typename Container:: template Property_map<Index>;

				typedef typename Container::const_iterator Iter;

				clear();

				Color_map colors;
				Label_map labels;
				Plane_map planes;
				Index_map indices;

				out << 
				"ply\n"                  << 
				"format ascii 1.0\n"     << 
				"element vertex "        << input.number_of_points() << "\n" << 
				"property double x\n"    << 
				"property double y\n"    << 
				"property double z\n"    <<
				"property double nx\n"   <<
				"property double ny\n"   <<
				"property double nz\n"   <<
				"property uchar red\n"   << 
				"property uchar green\n" <<
				"property uchar blue\n";

				boost::tie(colors,  boost::tuples::ignore) = input. template property_map<Color>("color");

				if (withExtraProperties) {

					out << 
					"property label\n"  <<
					"property index\n"  <<
					"property planea\n" <<
					"property planeb\n" <<
					"property planec\n" <<
					"property planed\n" <<
					"end_header\n";
					
					boost::tie(labels,  boost::tuples::ignore) = input. template property_map<Label>("label");
					boost::tie(indices, boost::tuples::ignore) = input. template property_map<Index>("index");
					boost::tie(planes,  boost::tuples::ignore) = input. template property_map<Plane>("plane");
				
				} else out << "end_header\n";

				for (Iter it = input.begin(); it != input.end(); ++it) {

					out.precision(10);

					out 
					<< input.point(*it)  << " " 
					<< input.normal(*it) << " "
					<< static_cast<int>(colors[*it][0]) << " " 
					<< static_cast<int>(colors[*it][1]) << " " 
					<< static_cast<int>(colors[*it][2]);

					if (withExtraProperties) out << " " <<  labels[*it] << " " << indices[*it] << " " <<  planes[*it];
					
					out << "\n";
				}
				save(fileName);
			}

			std::stringstream out;
		};
	}
}

#endif // CGAL_MYLOG_H