#include "../include/ply_in.h"
#include <string>
#include <iostream>
#include <fstream>


namespace Skippy {
	void Ply_In::read(const std::string & filename, std::vector<std::vector<CGAL_Point_3> > & polygons)
	{
		// Opens file
		std::ifstream file(filename);

		if (!file.is_open()) {
			throw std::ios_base::failure("Error : input file cannot be opened.");
		}

		// Declares the data we're interested in
		int V = -1, F = -1;
		std::vector<CGAL_Point_3> vertices;

		// Sequentially reads elements
		try {
			get_number_of_vertices_and_facets(file, V, F);
			get_vertices(file, V, vertices);
			get_facets(file, F, vertices, polygons);

		} catch (std::exception & except) {
			polygons.clear();
			file.close();
			throw except;
		}

		// Closes file
		file.close();
	}



	void Ply_In::read(const std::string & filename, std::vector<CGAL_Point_3> & points)
	{
		// Same as before, but we suppose that we read a file that contains points.

		std::ifstream file(filename);
		if (!file.is_open()) {
			throw std::ios_base::failure("Error : input file cannot be opened.");
		}

		// Declares the data we're interested in
		int V = -1, F = -1;

		// Sequentially reads elements
		try {
			get_number_of_vertices_and_facets(file, V, F);
			get_vertices(file, V, points);

		} catch (std::exception & except) {
			points.clear();
			file.close();
			throw except;
		}

		// Closes file
		file.close();
	}



	void Ply_In::read(const std::string & filename, std::vector<CGAL_Point_3> & points, std::vector<CGAL_Vector_3> & normals)
	{
		// Same as before, but we suppose that we read a file that contains points.

		std::ifstream file(filename);
		if (!file.is_open()) {
			throw std::ios_base::failure("Error : input file cannot be opened.");
		}

		// Declares the data we're interested in
		int V = -1, F = -1;

		// Sequentially reads elements
		try {
			get_number_of_vertices_and_facets(file, V, F);
			get_vertices(file, V, points, normals);

		} catch (std::exception & except) {
			points.clear();
			normals.clear();

			file.close();
			throw except;
		}

		// Closes file
		file.close();
	}



	void Ply_In::read(const std::string & filename, std::vector<CGAL_Inexact_Point_3> & points)
	{
		// Same as before, but we suppose that we read a file that contains points.

		std::ifstream file(filename);
		if (!file.is_open()) {
			throw std::ios_base::failure("Error : input file cannot be opened.");
		}

		// Declares the data we're interested in
		int V = -1, F = -1;

		// Sequentially reads elements
		try {
			get_number_of_vertices_and_facets(file, V, F);
			get_vertices(file, V, points);

		} catch (std::exception & except) {
			points.clear();
			file.close();
			throw except;
		}

		// Closes file
		file.close();
	}

	bool Ply_In::get_words(std::ifstream & file, std::vector<std::string> & words)
	{
		words.clear();
		std::string line;
		if (!std::getline(file, line)) return false;

		std::istringstream stream(line);
		std::copy(std::istream_iterator<std::string>(stream), std::istream_iterator<std::string>(), std::back_inserter(words));
		return true;
	}


	void Ply_In::get_number_of_vertices_and_facets(std::ifstream & file, int & V, int & F)
	{
		std::string line;
		std::vector<std::string> words;

		while (true) {
			if (!get_words(file, words))
				throw std::logic_error("Parsing error : unexpected end of file.");

			if (words.size() == 0) continue;

			if (words[0] == "end_header") {
				break;
			} else if (words[0] == "element") {
				if (words.size() != 3)
					throw std::logic_error("Parsing error : unspecified number of vertices or facets.");

				if (words[1] == "vertex") {
					V = atoi(words[2].c_str());
				} else if (words[1] == "face") {
					F = atoi(words[2].c_str());
				}
			}
		}

		if (V == -1 || F == -1)
			throw std::logic_error("Parsing error : unspecified number of vertices or facets.");
	}


	void Ply_In::get_vertices(std::ifstream & file, const int n, std::vector<CGAL_Point_3> & vertices)
	{
		std::string line;
		vertices.reserve(n);

		for (int i = 0; i < n; i++) {
			if (!std::getline(file, line)) {
				throw std::logic_error("Parsing error in PLY file : unexpected end of file.");
			}

			std::istringstream stream(line);
			CGAL_Point_3 pt;
			if (!(stream >> pt)) {
				throw std::logic_error("Parsing error in PLY file : missing vertices coordinates.");
			}

			vertices.push_back(pt);
		}
	}


	void Ply_In::get_vertices(std::ifstream & file, const int n, std::vector<CGAL_Point_3> & vertices, std::vector<CGAL_Vector_3> & normals)
	{
		std::string line;
		vertices.reserve(n);
		normals.reserve(n);

		for (int i = 0; i < n; i++) {
			if (!std::getline(file, line)) {
				throw std::logic_error("Parsing error in PLY file : unexpected end of file.");
			}

			std::istringstream stream(line);
			CGAL_Point_3 pt;
			if (!(stream >> pt)) {
				throw std::logic_error("Parsing error in PLY file : missing vertices coordinates.");
			}
			CGAL_Vector_3 pt_n;
			if (!(stream >> pt_n)) {
				throw std::logic_error("Parsing error in PLY file : missing normal coordinates.");
			}

			vertices.push_back(pt);
			normals.push_back(pt_n);
		}
	}


	void Ply_In::get_vertices(std::ifstream & file, const int n, std::vector<CGAL_Inexact_Point_3> & vertices)
	{
		std::string line;
		vertices.reserve(n);

		for (int i = 0; i < n; i++) {
			if (!std::getline(file, line)) {
				throw std::logic_error("Parsing error in PLY file : unexpected end of file.");
			}

			std::istringstream stream(line);
			CGAL_Inexact_Point_3 pt;
			if (!(stream >> pt)) {
				throw std::logic_error("Parsing error in PLY file : missing vertices coordinates.");
			}

			vertices.push_back(pt);
		}
	}


	
	void Ply_In::get_facets(std::ifstream & file, const int n, const std::vector<CGAL_Point_3> & vertices,
		std::vector<std::vector<CGAL_Point_3> > & polygons)
	{
		std::string line;
		polygons.reserve(n);

		for (int i = 0; i < n; i++) {
			if (!std::getline(file, line)) {
				throw std::logic_error("Parsing error in PLY file : unexpected end of file.");
			}

			std::istringstream stream(line);
			int p, q;
			if (!(stream >> p))
				throw std::logic_error("Parsing error in PLY file : facet size is missing.");

			std::vector<CGAL_Point_3> polygon;
			polygon.reserve(p);

			for (int j = 0; j < p; j++) {
				if (!(stream >> q))
					throw std::logic_error("Parsing error in PLY file : facet index is missing.");

				if (q < 0 || q >= vertices.size())
					throw std::logic_error("Parsing error in PLY file : facet index is invalid.");

				polygon.push_back(vertices[q]);
			}

			polygons.push_back(polygon);
		}
	}
}