#include "../include/polygon.h"

namespace Skippy 
{
	Polygon_Directions::Polygon_Directions(const CGAL_Point_2 & _O,
		const std::vector<std::pair<CGAL_Point_2, CGAL_Vector_2> > & _vertices)
		// const std::vector<Intersection_Line*> & _reference_lines)
		: O(_O)
	{
		// Copies vectors
		std::copy(_vertices.begin(), _vertices.end(), std::back_inserter(vertices));
		// std::copy(_reference_lines.begin(), _reference_lines.end(), std::back_inserter(reference_lines));
	}



	Polygon_Directions::~Polygon_Directions()
	{
		vertices.clear();
	}


	size_t Polygon_Directions::size() const
	{
		return vertices.size();
	}


	std::pair<CGAL_Point_2, CGAL_Vector_2> & Polygon_Directions::get_vertex(int i)
	{
		size_t n = vertices.size();
		return vertices[i % n];
	}


	/*Intersection_Line* Polygon_Directions::get_reference_line(int i)
	{
		return reference_lines[i];
	}*/


	const CGAL_Point_2 & Polygon_Directions::get_barycenter() const
	{
		return O;
	}
}