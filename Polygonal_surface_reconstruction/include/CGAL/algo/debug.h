#ifndef MY_DEBUG_H
#define MY_DEBUG_H

#include <CGAL/Point_set_with_segments.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/assertions.h>

#include <set>
#include <map>


//#define MY_DEBUG

#ifdef MY_DEBUG

namespace CGAL {

	/// \cond SKIP_IN_MANUA

	// for internal debug

	template <typename Kernel>
	class MyDebug
	{
	private:
		typedef typename Kernel::FT						FT;
		typedef typename Kernel::Point_3				Point;
		typedef typename Kernel::Point_2				Point2;
		typedef typename Kernel::Vector_3				Vector;
		typedef typename Kernel::Line_3					Line;
		typedef typename Kernel::Segment_3				Segment;
		typedef typename Kernel::Plane_3				Plane;
		typedef CGAL::Planar_segment<Kernel>			Planar_segment;
		typedef CGAL::Point_set_with_segments<Kernel>	Point_set;

		typedef CGAL::Surface_mesh<Point>				Mesh;
		typedef typename Mesh::Face_index				Face_descriptor;
		typedef typename Mesh::Edge_index				Edge_descriptor;
		typedef typename Mesh::Vertex_index				Vertex_descriptor;

		typedef typename Mesh::Halfedge_index			Halfedge_descriptor;

	public:
		static void save_face(Face_descriptor f, const Mesh& mesh, const std::string& file) {
			const Mesh::Property_map<Vertex_descriptor, Point>& coords = mesh.points();
			std::ofstream output(file.c_str());
			output << 1 << std::endl;
			output << mesh.degree(f) << " 3" << std::endl;

			Halfedge_around_face_circulator<Mesh> cir(mesh.halfedge(f), mesh), done(cir);
			do {
				Halfedge_descriptor hd = *cir;
				Vertex_descriptor vd = mesh.target(hd);
				output << coords[vd] << std::endl;
				++cir;
			} while (cir != done);
		}


		static void print_face(Face_descriptor f, const Mesh& mesh) {
			const Mesh::Property_map<Vertex_descriptor, Point>& coords = mesh.points();
			std::cout << "face degree: " << mesh.degree(f) << std::endl;
			Halfedge_around_face_circulator<Mesh> cir(mesh.halfedge(f), mesh), done(cir);
			do {
				Halfedge_descriptor hd = *cir;
				Vertex_descriptor vd = mesh.target(hd);
				std::cout << "\t" << coords[vd] << std::endl;
				++cir;
			} while (cir != done);
		}

		static void save_mesh(const Mesh& mesh, const std::string& file) {
			CGAL::write_off(std::ofstream(file.c_str()), mesh);
		}

		static void two_vertices_donot_coincide(Vertex_descriptor v0, Vertex_descriptor v1, const Mesh& mesh) {
			if (v0 == v1)
				std::cout << "Fatal error: two vertices are the same" << std::endl;

			const Mesh::Property_map<Vertex_descriptor, Point>& coords = mesh.points();
			const Point& p0 = coords[v0];
			const Point& p1 = coords[v1];
			FT sd = CGAL::squared_distance(p0, p1);
			if (sd < 1e-10) {
				std::cout << "Fatal error: two vertices coincide. The two vertices are: " << std::endl;
				std::cout << p0 << std::endl << p1 << std::endl;
			}
		}

		static void two_vertices_should_on_same_plane(Vertex_descriptor v0, Vertex_descriptor v1, const Mesh& mesh, const Plane* plane) {
			const Mesh::Property_map<Vertex_descriptor, Point>& coords = mesh.points();
			FT s_dist = std::sqrt(CGAL::squared_distance(*plane, coords[v0]));
			FT t_dist = std::sqrt(CGAL::squared_distance(*plane, coords[v1]));
			if (s_dist > 1e-5)
				std::cout << "vertex s doesn't lie on plane. Distance is " << s_dist << std::endl;
			if (t_dist > 1e-5)
				std::cout << "vertex t doesn't lie on plane. Distance is " << t_dist << std::endl;
		}
	};


} //namespace CGAL


#endif

#endif 
