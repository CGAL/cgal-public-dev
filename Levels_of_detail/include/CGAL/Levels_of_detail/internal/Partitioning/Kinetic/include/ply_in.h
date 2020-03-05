#pragma once
#include "defs_cgal.h"


namespace Skippy {

	class Ply_In 
	{
	public:
		Ply_In() {}

		~Ply_In() {}

	public:
		static void read(const std::string & filename, std::vector<std::vector<CGAL_Point_3> > & polygons);

		static void read(const std::string & filename, std::vector<CGAL_Point_3> & points);

		static void read(const std::string & filename, std::vector<CGAL_Point_3> & points, std::vector<CGAL_Vector_3> & normals);

		static void read(const std::string & filename, std::vector<CGAL_Inexact_Point_3> & points);

		static bool get_words(std::ifstream & file, std::vector<std::string> & words);

		static void get_number_of_vertices_and_facets(std::ifstream & file, int & V, int & F);

		static void get_vertices(std::ifstream & file, const int V, std::vector<CGAL_Point_3> & vertices);

		static void get_vertices(std::ifstream & file, const int V, 
			std::vector<CGAL_Point_3> & vertices, std::vector<CGAL_Vector_3> & normals);

		static void get_vertices(std::ifstream & file, const int V, std::vector<CGAL_Inexact_Point_3> & vertices);

		static void get_facets(std::ifstream & file, const int F, const std::vector<CGAL_Point_3> & vertices,
			std::vector<std::vector<CGAL_Point_3> > & polygons);
	};
}

