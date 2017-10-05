#ifndef CGAL_MYLOG_H
#define CGAL_MYLOG_H

// STL includes.
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <cmath>

// Boost includes.
#include <boost/tuple/tuple.hpp>

// CGAL includes.
#include <CGAL/array.h>
#include <CGAL/Unique_hash_map.h>
#include <CGAL/IO/Color.h>

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

			void skip_line() {
				out << std::endl;
			}

			template<typename T>
			void add_index(const T index) {
				out << "Index: " << index << std::endl;
			}

			bool save(const std::string &fileName, const std::string &extension = ".log", const std::string path = "/Users/danisimo/Documents/pipeline/logs/") const {

				const std::string finalPath = path + fileName + extension;
				std::ofstream file(finalPath.c_str(), std::ios_base::out);

				if (!file) {
					std::cerr << "\nERROR: Error saving log file with the name " << fileName << "\n" << std::endl;
					return false;
				}

				file << data() << std::endl;
				file.close();

				return true;
			}

			// Save mesh in a ply format.
			template<class Mesh, class Facet_colors>
			void save_mesh_as_ply(const Mesh &mesh, const Facet_colors &mesh_facet_colors, const std::string &filename) {

				clear();

				using Mesh_vertex_iterator = typename Mesh::Vertex_const_iterator;
				using Mesh_facet_iterator  = typename Mesh::Facet_const_iterator;
				using Mesh_facet_handle    = typename Mesh::Facet_const_handle;

				using Facet_halfedge_circulator = typename Mesh::Facet::Halfedge_around_facet_const_circulator;
				using Facet_vertex_handle 		= typename Mesh::Facet::Vertex_const_handle;


				// Add ply header.
				const size_t num_vertices = mesh.size_of_vertices();
				const size_t num_facets   = mesh.size_of_facets();

				out << 
				"ply\n"               					   << 
				"format ascii 1.0\n"     				   << 
				"element vertex "        				   << num_vertices << "\n" << 
				"property double x\n"    				   << 
				"property double y\n"    				   << 
				"property double z\n" 					   <<
				"element face " 						   << num_facets << "\n" << 
				"property list uchar int vertex_indices\n" <<
				"property uchar red\n" 					   <<
				"property uchar green\n" 				   <<
				"property uchar blue\n" 				   <<
				"end_header\n";


				// Add mesh vertices.
				int count = 0;
				CGAL::Unique_hash_map<Facet_vertex_handle, int> V;

				for (Mesh_vertex_iterator vit = mesh.vertices_begin(); vit != mesh.vertices_end(); ++vit) {
					
					V[static_cast<Facet_vertex_handle>(vit)] = count++;
					out << vit->point() << std::endl;
				}


				// Add mesh facets.
				for (Mesh_facet_iterator fit = mesh.facets_begin(); fit != mesh.facets_end(); ++fit) {
					out << (*fit).facet_degree() << " "; 

					Facet_halfedge_circulator hit = (*fit).facet_begin(), end = (*fit).facet_begin();
					do {

						out << V[(*hit).vertex()] << " ";
						++hit;

					} while (hit != end);
					out << mesh_facet_colors.at(static_cast<Mesh_facet_handle>(fit)) << std::endl;
				}


				// Save file.
				save(filename, ".ply");
			}

			template<class CDT, class Buildings>
			void save_buildings_info(const CDT &cdt, const Buildings &buildings, const std::string &filename) {

				typedef typename CDT::Vertex_handle 		   Vertex_handle;
				typedef typename CDT::Face_handle 			   Face_handle;
				typedef typename CDT::Finite_faces_iterator    Face_iterator;
				typedef typename CDT::Finite_vertices_iterator Vertex_iterator;
				typedef typename Buildings::const_iterator     Building_iterator;

				clear();
				
				CGAL::Unique_hash_map<Vertex_handle, int> V;
				CGAL::Unique_hash_map<Face_handle, int>   F;

				size_t count = 0;
				for (Face_iterator fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit) F[fit] = count++;

				count = 0;
				for (Vertex_iterator vit = cdt.finite_vertices_begin(); vit != cdt.finite_vertices_end(); ++vit) V[vit] = count++;

				count = 0;
				for (Building_iterator bit = buildings.begin(); bit != buildings.end(); ++bit, ++count) {
					
					const std::vector<Face_handle> &faces = (*bit).second.faces;
					const size_t num_faces = faces.size();

					const std::vector< std::vector<Vertex_handle> > &boundaries = (*bit).second.boundaries;
					const size_t num_boundaries = boundaries.size();

					out << "Building " << count << 
					" with color " << (*bit).second.color  << 
					" and height " << (*bit).second.height << 
					
					std::endl << "faces: ";
					for (size_t i = 0; i < num_faces; ++i) out << F[faces[i]] << " ";
					
					skip_line();
					for (size_t k = 0; k < num_boundaries; ++k) {

						out << "boundary: " << k << std::endl;
						const size_t num_vertices = boundaries[k].size();

						for (size_t i = 0; i < num_vertices; ++i) {

							out << V[boundaries[k][i]] << std::endl;

							if (!(*bit).second.wedges.empty()) {
								out << "with wedges: ";

								const std::vector<Face_handle> &wedges = (*bit).second.wedges[k].at(boundaries[k][i]);
								for (size_t j = 0; j < wedges.size(); ++j) out << F[wedges[j]] << " ";
								out << std::endl;
							}
						}
					}

					skip_line();
					skip_line();
				}
				save(filename);
			}

			template<class Points>
			void export_points(const std::string &name, Points &points) {

				for (typename Points::const_iterator it = points.begin(); it != points.end(); ++it)
					out << *it << " " << 0 << std::endl;

				save(name, ".xyz");
			}

			template<class Traits, class Container>
			void save_ply(const Container &input, 
						  const std::string &fileName,
						  const bool withExtraProperties = false) {

				using Type  = unsigned char;
				using Color = CGAL::cpp11::array<Type, 3>;
				
				using Label = int;
				using Types = int;
				using Index = int;

				// using Plane = typename Traits::Plane_3;

				using Color_map = typename Container:: template Property_map<Color>;
				using Label_map = typename Container:: template Property_map<Label>;
				using Types_map = typename Container:: template Property_map<Types>;
				using Index_map = typename Container:: template Property_map<Index>;

				// using Plane_map = typename Container:: template Property_map<Plane>;

				typedef typename Container::const_iterator Iter;

				clear();

				Color_map colors;
				Label_map labels;
				Types_map types;
				Index_map indices;

				// Plane_map planes;

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
					"property int label\n" <<
					"property int type\n"  <<
					"property int index\n" <<
					"end_header\n";
					
					boost::tie(labels,  boost::tuples::ignore) = input. template property_map<Label>("label");
					boost::tie(types ,  boost::tuples::ignore) = input. template property_map<Types>("types");
					boost::tie(indices, boost::tuples::ignore) = input. template property_map<Index>("index");

					// boost::tie(planes,  boost::tuples::ignore) = input. template property_map<Plane>("plane");
				
				} else out << "end_header\n";

				for (Iter it = input.begin(); it != input.end(); ++it) {

					// if (static_cast<int>(*it) % 3 == 0) out << "\n"; // remove if not needed

					out.precision(10);

					out 
					<< input.point(*it)  << " " 
					<< input.normal(*it) << " "
					<< static_cast<int>(colors[*it][0]) << " " 
					<< static_cast<int>(colors[*it][1]) << " " 
					<< static_cast<int>(colors[*it][2]);

					if (withExtraProperties) out << " " <<  labels[*it] << " " << types[*it] << " " << indices[*it];
					
					out << "\n";
				}
				save(fileName, ".log");
			}

			template<class CDT>
			void save_cdt_ply(const CDT &cdt, const std::string &filename, const std::string &color_type = "in") {

				clear();

				const auto number_of_vertices = cdt.number_of_vertices();
				const auto number_of_faces    = cdt.number_of_faces();

				out << 
				"ply\n"               					   << 
				"format ascii 1.0\n"     				   << 
				"element vertex "        				   << number_of_vertices << "\n" << 
				"property double x\n"    				   << 
				"property double y\n"    				   << 
				"property double z\n" 					   <<
				"element face " 						   << number_of_faces << "\n" << 
				"property list uchar int vertex_indices\n" <<
				"property uchar red\n" 					   <<
				"property uchar green\n" 				   <<
				"property uchar blue\n" 				   <<
				"end_header\n";

				typedef typename CDT::Vertex_handle Vertex_handle;
				CGAL::Unique_hash_map<Vertex_handle, int> V;

				int count = 0;
				for (typename CDT::Finite_vertices_iterator vit = cdt.finite_vertices_begin(); vit != cdt.finite_vertices_end(); ++vit) {
					
					out << (*vit) << " " << 0 << std::endl;
					V[vit] = count++;
				}

				for (typename CDT::Finite_faces_iterator fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit) {
					
					out << "3 " << V[(*fit).vertex(0)] << " " << V[(*fit).vertex(1)] << " " << V[(*fit).vertex(2)] << " ";
					
					if (color_type == "in") out << fit->info().in_color << std::endl;
					if (color_type == "bu") out << fit->info().bu_color << std::endl;
				}

				save(filename, ".ply");
			}

			template<class CDT>
			void save_cdt_obj(const CDT &cdt, const std::string &filename) {

				clear();

				typedef typename CDT::Vertex_handle Vertex_handle;
				CGAL::Unique_hash_map<Vertex_handle, int> V;

				int count = 0;
				for (typename CDT::Finite_vertices_iterator vit = cdt.finite_vertices_begin(); vit != cdt.finite_vertices_end(); ++vit) {
					
					out << "v " << (*vit) << " " << 0 << std::endl;
					V[vit] = count++;
				}

				for (typename CDT::Finite_faces_iterator fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit)
					out << "f " << V[(*fit).vertex(0)] + 1 << " " << V[(*fit).vertex(1)] + 1 << " " << V[(*fit).vertex(2)] + 1 << std::endl;

				save(filename, ".obj");
			}

			template<class CDT, class Container, class Segments>
			void save_visibility_eps(CDT &cdt, const Container &input, const Segments &segments, const std::string &fileName = "tmp/visibility", const std::string &color_type = "in") {

				clear();

		        // Compute bounding box.
		        double minbX, minbY, maxbX, maxbY;
		        bounding_box(cdt, minbX, minbY, maxbX, maxbY);

		        // Compute scale.
		        double scale = 1.0;
		        if (std::sqrt((maxbX - minbX) * (maxbX - minbX) + (maxbY - minbY) * (maxbY - minbY)) < 10.0 && scale == 1.0) scale *= 1000.0;

		        // Set header.
		        set_header(minbX * scale, minbY * scale, maxbX * scale, maxbY * scale);

		        // Start private namespace.
		        out << "0 dict begin gsave\n\n";

		        // Save mesh.
		        draw_mesh(cdt, scale, color_type);

		        // Save points.
		        draw_points(input, scale);

		        // Save segments.
		        draw_segments(segments, scale);

		        // Finish private namespace.
		        out << "grestore end\n\n";
		        out << "%%EOF\n";

		        save(fileName, ".eps");
			}

			template<class Point>
			void save_triangle_with_points_eps(const Point &a, const Point &b, const Point &c, const std::vector<Point> &samples, const std::string &fileName = "tmp/triangle") {

				clear();

		        // Compute bounding box.
		        double minbX, minbY, maxbX, maxbY;
		        bounding_box(a, b, c, minbX, minbY, maxbX, maxbY);

		        // Compute scale.
		        double scale = 1.0;
		        if (std::sqrt((maxbX - minbX) * (maxbX - minbX) + (maxbY - minbY) * (maxbY - minbY)) < 10.0 && scale == 1.0) scale *= 1000.0;

		        // Set header.
		        set_header(minbX * scale, minbY * scale, maxbX * scale, maxbY * scale);

		        // Start private namespace.
		        out << "0 dict begin gsave\n\n";

		        // Save mesh.
		        draw_triangle_with_points(a, b, c, samples, scale);

		        // Finish private namespace.
		        out << "grestore end\n\n";
		        out << "%%EOF\n";

		        save(fileName, ".eps");
			}

			template<class CDT>
			void save_visibility_eps(CDT &cdt, const std::string &fileName = "tmp/visibility", const std::string &color_type = "in") {

				clear();

		        // Compute bounding box.
		        double minbX, minbY, maxbX, maxbY;
		        bounding_box(cdt, minbX, minbY, maxbX, maxbY);

		        // Compute scale.
		        double scale = 1.0;
		        if (std::sqrt((maxbX - minbX) * (maxbX - minbX) + (maxbY - minbY) * (maxbY - minbY)) < 10.0 && scale == 1.0) scale *= 1000.0;

		        // Set header.
		        set_header(minbX * scale, minbY * scale, maxbX * scale, maxbY * scale);

		        // Start private namespace.
		        out << "0 dict begin gsave\n\n";

		        // Save mesh.
		        draw_mesh(cdt, scale, color_type);

		        // Finish private namespace.
		        out << "grestore end\n\n";
		        out << "%%EOF\n";

		        save(fileName, ".eps");
			}

			template<class Segments>
			void export_segments_as_obj(const std::string &name, const Segments &segments, const std::string &) {

				clear();
				for (size_t i = 0; i < segments.size(); ++i) {

					out << "v " << segments[i].source() << " " << 0 << std::endl;
					out << "v " << segments[i].target() << " " << 0 << std::endl;
					out << "v " << segments[i].target() << " " << 0 << std::endl;
				}

				for (size_t i = 0; i < segments.size() * 3; i += 3)
					out << "f " << i + 1 << " " << i + 2 << " " << i + 3 << std::endl;

				// save("segments", ".obj", default_path);
				save(name, ".obj");
			}

			template<class Projected_points>
			void export_projected_points_as_xyz(const std::string &name, const Projected_points &projected, const std::string &) {
				
				clear();
				for (typename Projected_points::const_iterator it = projected.begin(); it != projected.end(); ++it)
					out << (*it).second << " " << 0 << std::endl;

				// save(name, ".xyz", default_path);
				save(name, ".xyz");
			}

			template<class Planes_mapping, class Projected_points>
			void export_projected_points(const std::string &name, const Planes_mapping &planes, const Projected_points &projected) {

				clear();

				Projected_points tmp_points;
				auto plane_index = 0;

				using Plane_iterator = typename Planes_mapping::const_iterator;

				for (Plane_iterator it = planes.begin(); it != planes.end(); ++it, ++plane_index) {
					const auto num_points = (*it).second.size();

					tmp_points.clear();
					for (size_t i = 0; i < num_points; ++i) {
						
						const auto index = (*it).second[i];
						tmp_points[i] = projected.at(index);
					}

					export_projected_points(name + "_" + std::to_string(plane_index), tmp_points);
				}
			}

			std::stringstream out;

		private:
			template<class CDT>
			void bounding_box(const CDT &cdt, double &minbX, double &minbY, double &maxbX, double &maxbY) const {

				const double big_value = 100000.0;

		        minbX =  big_value, minbY =  big_value;
		        maxbX = -big_value, maxbY = -big_value;

		        for (typename CDT::Finite_vertices_iterator it = cdt.finite_vertices_begin(); it != cdt.finite_vertices_end(); ++it) {

		        	const double x = static_cast<double>((*it).point().x());
		        	const double y = static_cast<double>((*it).point().y());

		        	minbX = std::min(minbX, x);
		        	minbY = std::min(minbY, y);

		        	maxbX = std::max(maxbX, x);
		        	maxbY = std::max(maxbY, y);
		        }
			}

			template<class Point>
			void bounding_box(const Point &a, const Point &b, const Point &c, double &minbX, double &minbY, double &maxbX, double &maxbY) const {
				
				const double big_value = 100000.0;

		        minbX =  big_value, minbY =  big_value;
		        maxbX = -big_value, maxbY = -big_value;

		        minbX = std::min(minbX, a.x()); minbY = std::min(minbY, a.y());
		        maxbX = std::max(maxbX, a.x()); maxbY = std::max(maxbY, a.y());

		        minbX = std::min(minbX, b.x()); minbY = std::min(minbY, b.y());
		        maxbX = std::max(maxbX, b.x()); maxbY = std::max(maxbY, b.y());

		        minbX = std::min(minbX, c.x()); minbY = std::min(minbY, c.y());
		        maxbX = std::max(maxbX, c.x()); maxbY = std::max(maxbY, c.y());
			}

			void set_header(const double llx, const double lly, const double urx, const double ury) {
		        
		        out << "%!PS-Adobe-3.0 EPSF-3.0\n";
		        out << "%%BoundingBox: " << llx << " " << lly << " " << urx << " " << ury << "\n";
		        out << "%%Pages: 1\n";
		        out << "%%Creator: Dmitry Anisimov, danston@ymail.com\n";
		        out << "%%EndComments\n";
		        out << "%%EndProlog\n\n";
		        out << "%%Page: 1 1\n\n";
			}

			template<class CDT>
    		void draw_mesh(const CDT &cdt, const double scale, const std::string &color_type) {

				for (typename CDT::Finite_faces_iterator fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit) {

					out << (*(*fit).vertex(0)).point().x() * scale << " " << (*(*fit).vertex(0)).point().y() * scale << " moveto\n";
					out << (*(*fit).vertex(1)).point().x() * scale << " " << (*(*fit).vertex(1)).point().y() * scale << " lineto\n";
					out << (*(*fit).vertex(2)).point().x() * scale << " " << (*(*fit).vertex(2)).point().y() * scale << " lineto\n";
					out << (*(*fit).vertex(0)).point().x() * scale << " " << (*(*fit).vertex(0)).point().y() * scale << " lineto\n";

					out << "closepath\n\n";
					out << "gsave\n";

					/*
					const double visibility = static_cast<double>(fit->info().in);
					const double half = 0.5;

					if (visibility > half) out << "0.2 1 0.2 setrgbcolor\n";	  // INSIDE
					else if (visibility < half) out << "1 0.2 0.2 setrgbcolor\n"; // OUTSIDE
					else out << "1 0.8 0 setrgbcolor\n";						  // UNKNOWN */

					CGAL::Color fc;

					if (color_type == "in") fc = fit->info().in_color;
					if (color_type == "bu") fc = fit->info().bu_color;

					const double r = static_cast<double>(fc.red())   / 255.0;
					const double g = static_cast<double>(fc.green()) / 255.0;
					const double b = static_cast<double>(fc.blue())  / 255.0;

					out << r << " " << g << " " << b << " setrgbcolor\n";

					out << "fill\n";
					out << "grestore\n";
					out << "0 0 0 setrgbcolor\n";
					out << "0.5 setlinewidth\n";
		        	out << "stroke\n\n";
				}
    		}

    		template<class Point>
    		void draw_triangle_with_points(const Point &a, const Point &b, const Point &c, const std::vector<Point> &samples, const double scale) {

    			out << a.x() * scale << " " << a.y() * scale << " moveto\n";
    			out << b.x() * scale << " " << b.y() * scale << " lineto\n";
    			out << c.x() * scale << " " << c.y() * scale << " lineto\n";
    			out << a.x() * scale << " " << a.y() * scale << " lineto\n";

				out << "closepath\n\n";
				out << "0 0 0 setrgbcolor\n";
				out << "2 setlinewidth\n";
		        out << "stroke\n\n";

				for (size_t i = 0; i < samples.size(); ++i) draw_disc(samples[i], scale);
    		}

    		template<class Container>
    		void draw_points(const Container &input, const double scale) {
        		for (typename Container::const_iterator it = input.begin(); it != input.end(); ++it) draw_disc(input.point(*it), scale);
    		}

    		template<class Point>
    		void draw_disc(const Point &p, const double scale) {

		        out << "0 setgray\n";
		        out << "0 setlinewidth\n\n";
		        out << p.x() * scale << " " << p.y() * scale << " " << 10 << " 0 360 arc closepath\n\n";
		        out << "gsave\n";
		        out << "0 setgray fill\n";
		        out << "grestore\n";
		        out << "stroke\n\n";
    		}

    		template<class Segments>
    		void draw_segments(const Segments &segments, const double scale) {

    			if (segments.empty()) return;
    			for (size_t i = 0; i < segments.size(); ++i) {

    				for (size_t j = 0; j < segments[i].size() - 1; ++j) {
    					out << segments[i][j].x() * scale     << " " << segments[i][j].y() * scale     << " moveto\n";
    					out << segments[i][j + 1].x() * scale << " " << segments[i][j + 1].y() * scale << " lineto\n";
    				}
    			}

    			out << "0 0 0 setgray\n";
    			out << "2 setlinewidth\n";
    			out << "stroke\n\n";
    		}
		};
	}
}

#endif // CGAL_MYLOG_H