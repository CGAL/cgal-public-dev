#ifndef CGAL_LEVEL_OF_DETAIL_MYLOG_H
#define CGAL_LEVEL_OF_DETAIL_MYLOG_H

#if defined(WIN32) || defined(_WIN32) 
#define _NL_ "\r\n"
#else 
#define _NL_ "\n"
#endif

// CGAL includes.
#include <CGAL/Unique_hash_map.h>

// STL includes.
#include <vector>
#include <string>
#include <sstream>
#include <fstream>

namespace CGAL {

	namespace Level_of_detail {

		class Mylog {

		public:
			Mylog() {
				out.precision(20);
			}

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

			template<class Faces_range, class Point_map, class Colour_map>
			void save_faces(const Faces_range &faces_range, const Point_map &point_map, const Colour_map &colour_map, const std::string &file_name) {

				clear();
				using Const_faces_iterator = typename Faces_range::const_iterator;

				size_t num_faces    = 0;
				size_t num_vertices = 0;

				for (Const_faces_iterator cf_it = faces_range.begin(); cf_it != faces_range.end(); ++cf_it) {
					
					num_faces    += 1;
					num_vertices += (*cf_it).size();
				}

				out << 
				"ply" 				   +  std::string(_NL_) + ""               			  << 
				"format ascii 1.0"     +  std::string(_NL_) + ""     			          << 
				"element vertex "      << num_vertices     << "" + std::string(_NL_) + "" << 
				"property double x"    +  std::string(_NL_) + ""    			          << 
				"property double y"    +  std::string(_NL_) + ""    			          << 
				"property double z"    +  std::string(_NL_) + "" 				          <<
				"element face "        << num_faces        << "" + std::string(_NL_) + "" << 
				"property list uchar int vertex_indices"         + std::string(_NL_) + "" <<
				"property uchar red"   +  std::string(_NL_) + "" 				          <<
				"property uchar green" +  std::string(_NL_) + "" 				          <<
				"property uchar blue"  +  std::string(_NL_) + "" 				          <<
				"end_header"           +  std::string(_NL_) + "";

				for (Const_faces_iterator cf_it = faces_range.begin(); cf_it != faces_range.end(); ++cf_it) {
					const auto &vertices = *cf_it;

					for (auto cv_it = vertices.begin(); cv_it != vertices.end(); ++cv_it)
						out << get(point_map, *cv_it) << std::endl;
				}

				size_t count = 0;
				for (Const_faces_iterator cf_it = faces_range.begin(); cf_it != faces_range.end(); ++cf_it) {
					const auto &vertices = *cf_it;
					
					const size_t num_vertices = (*cf_it).size();
					out << num_vertices << " ";

					for (size_t i = 0; i < num_vertices; ++i) out << count++ << " ";
					out << get(colour_map, *cf_it) << std::endl;
				}

				save(file_name, ".ply");
			}

			template<class Regions_range, class Point_map, class Colour_map>
			void save_regions(const Regions_range &regions_range, const Point_map &point_map, const Colour_map &colour_map, const std::string &file_name) {

				clear();
				using Const_regions_iterator = typename Regions_range::const_iterator;

				size_t num_vertices = 0;
				for (Const_regions_iterator cr_it = regions_range.begin(); cr_it != regions_range.end(); ++cr_it) num_vertices += (*cr_it).size();

				out << 
				"ply" 				   +  std::string(_NL_) + ""               			  << 
				"format ascii 1.0"     +  std::string(_NL_) + ""     			          << 
				"element vertex "      << num_vertices     << "" + std::string(_NL_) + "" << 
				"property double x"    +  std::string(_NL_) + ""    			          << 
				"property double y"    +  std::string(_NL_) + ""    			          << 
				"property double z"    +  std::string(_NL_) + "" 				          <<
				"property uchar red"   +  std::string(_NL_) + "" 				          <<
				"property uchar green" +  std::string(_NL_) + "" 				          <<
				"property uchar blue"  +  std::string(_NL_) + "" 				          <<
				"end_header"           +  std::string(_NL_) + "";

				for (Const_regions_iterator cr_it = regions_range.begin(); cr_it != regions_range.end(); ++cr_it) {
					
					const auto &points  = *cr_it;
					const auto colour = get(colour_map, *cr_it);
					
					for (auto cp_it = points.begin(); cp_it != points.end(); ++cp_it) 
						out << get(point_map, *cp_it) << " " << colour << std::endl;
				}
				
				save(file_name, ".ply");
			}

			template<class Elements, class Segment_map_2>
			void save_segments(const Elements &elements, const Segment_map_2 &segment_map_2, const std::string &file_name) {
				
				clear();
				using Const_elements_iterator = typename Elements::const_iterator;

				size_t size = 0;
				for (Const_elements_iterator ce_it = elements.begin(); ce_it != elements.end(); ++ce_it, ++size) {
					const auto &segment = get(segment_map_2, *ce_it);

					out << "v " << segment.source() << " " << 0 << std::endl;
					out << "v " << segment.target() << " " << 0 << std::endl;
					out << "v " << segment.target() << " " << 0 << std::endl;
				}

				for (size_t i = 0; i < size * 3; i += 3)
					out << "f " << i + 1 << " " << i + 2 << " " << i + 3 << std::endl;
				
				save(file_name, ".obj");
			}

			template<class Triangulation, class Face_colour_map>
			void save_triangulation(const Triangulation &triangulation, const Face_colour_map &face_colour_map, const std::string &file_name) {

				clear();
				using Triangulation_vertex_handle 	  = typename Triangulation::Vertex_handle;
				using Triangulation_faces_iterator	  = typename Triangulation::Finite_faces_iterator;
				using Triangulation_vertices_iterator = typename Triangulation::Finite_vertices_iterator;
				
				const size_t num_faces    = triangulation.number_of_faces();
				const size_t num_vertices = triangulation.number_of_vertices();

				out << 
				"ply" 				   +  std::string(_NL_) + ""               			  << 
				"format ascii 1.0"     +  std::string(_NL_) + ""     			          << 
				"element vertex "      << num_vertices     << "" + std::string(_NL_) + "" << 
				"property double x"    +  std::string(_NL_) + ""    			          << 
				"property double y"    +  std::string(_NL_) + ""    			          << 
				"property double z"    +  std::string(_NL_) + "" 				          <<
				"element face "        << num_faces        << "" + std::string(_NL_) + "" << 
				"property list uchar int vertex_indices"         + std::string(_NL_) + "" <<
				"property uchar red"   +  std::string(_NL_) + "" 				          <<
				"property uchar green" +  std::string(_NL_) + "" 				          <<
				"property uchar blue"  +  std::string(_NL_) + "" 				          <<
				"end_header"           +  std::string(_NL_) + "";

				int count = 0;
				CGAL::Unique_hash_map<Triangulation_vertex_handle, int> vertices;

				for (Triangulation_vertices_iterator tv_it = triangulation.finite_vertices_begin(); tv_it != triangulation.finite_vertices_end(); ++tv_it) {
					const auto &point_2 = *tv_it;

					out << point_2 << " " << 0 << std::endl;
					vertices[tv_it] = count++;
				}

				for (Triangulation_faces_iterator tf_it = triangulation.finite_faces_begin(); tf_it != triangulation.finite_faces_end(); ++tf_it) {
					const auto &face_handle = *tf_it;

					out << "3 " << 
					vertices[face_handle.vertex(0)]   << " " << 
					vertices[face_handle.vertex(1)]   << " " << 
					vertices[face_handle.vertex(2)]   << " " <<
					get(face_colour_map, face_handle) << std::endl;
				}
				save(file_name, ".ply");
			}

			template<class Lod, class Ground_colour_map, class Wall_colour_map, class Roof_colour_map>
			void save_lod(
				const Lod &lod, 
				const Ground_colour_map &ground_colour_map,
				const Wall_colour_map &wall_colour_map,
				const Roof_colour_map &roof_colour_map,
				const std::string &file_name) {

				clear();
				
				const size_t num_faces    = lod.number_of_faces();
				const size_t num_vertices = lod.number_of_vertices();

				out << 
				"ply" 				   +  std::string(_NL_) + ""               			  << 
				"format ascii 1.0"     +  std::string(_NL_) + ""     			          << 
				"element vertex "      << num_vertices     << "" + std::string(_NL_) + "" << 
				"property double x"    +  std::string(_NL_) + ""    			          << 
				"property double y"    +  std::string(_NL_) + ""    			          << 
				"property double z"    +  std::string(_NL_) + "" 				          <<
				"element face "        << num_faces        << "" + std::string(_NL_) + "" << 
				"property list uchar int vertex_indices"         + std::string(_NL_) + "" <<
				"property uchar red"   +  std::string(_NL_) + "" 				          <<
				"property uchar green" +  std::string(_NL_) + "" 				          <<
				"property uchar blue"  +  std::string(_NL_) + "" 				          <<
				"end_header"           +  std::string(_NL_) + "";


				// Points.
				const int stub = 0;
				
				// Ground.
				for (auto ground_point = lod.ground_face().begin(); ground_point != lod.ground_face().end(); ++ground_point)
					out << *ground_point << std::endl;

				// Walls.
				for (auto wall_face = lod.wall_faces().begin(); wall_face != lod.wall_faces().end(); ++wall_face)
					for (auto wall_point = wall_face->begin(); wall_point != wall_face->end(); ++wall_point)
						out << *wall_point << std::endl;

				// Roofs.
				for (auto roof_face = lod.roof_faces().begin(); roof_face != lod.roof_faces().end(); ++roof_face)
					for (auto roof_point = roof_face->begin(); roof_point != roof_face->end(); ++roof_point)
						out << *roof_point << std::endl;


				// Faces.
				size_t count = 0;

				// Ground.
				out << lod.ground_face().size() << " ";
				for (size_t i = 0; i < lod.ground_face().size(); ++i) out << count++ << " ";
				out << get(ground_colour_map, stub) << std::endl;

				// Walls.
				for (auto wall_face = lod.wall_faces().begin(); wall_face != lod.wall_faces().end(); ++wall_face) {

					out << wall_face->size() << " ";
					for (size_t i = 0; i < wall_face->size(); ++i) out << count++ << " ";
					out << get(wall_colour_map, stub) << std::endl;
				}

				// Roofs.
				for (auto roof_face = lod.roof_faces().begin(); roof_face != lod.roof_faces().end(); ++roof_face) {

					out << roof_face->size() << " ";
					for (size_t i = 0; i < roof_face->size(); ++i) out << count++ << " ";
					out << get(roof_colour_map, stub) << std::endl;
				}

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
        };

    } // Level_of_detail

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_MYLOG_H