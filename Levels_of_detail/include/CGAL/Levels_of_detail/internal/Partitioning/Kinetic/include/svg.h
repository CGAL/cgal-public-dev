#pragma once
#include <iostream>
#include <fstream>
#include "defs.h"
#include "support_plane_objects.h"

namespace Skippy {

	namespace SVG {

		inline void markup_header(std::ostream & os, int rows, int cols)
		{
			os << "<?xml version=\"1.0\" encoding=\"utf-8\"?>" << std::endl;
			os << "<svg xmlns=\"http://www.w3.org/2000/svg\" version=\"1.1\" width=\"" << cols << "\" height=\"" << rows << "\">" << std::endl;
		}

		inline void markup_image(std::ostream & os, std::string path, int rows, int cols, double opacity)
		{
			os << "  <image xlink:href=\"" << path << "\" x=\"0\" y=\"0\" height=\"" << rows << "px\" width=\"" << cols << "px\" opacity=\"" << opacity << "\" />" << std::endl;
		}

		inline void markup_dashed_line(std::ostream & os, double db, double dw, double x1, double x2, double y1, double y2, CGAL_Color color)
		{
			os << "  <line stroke-dasharray=\"" << db << "," << dw << "\" x1=\"" << x1 << "\" y1=\"" << y1 << "\" x2=\"" << x2 << "\" y2=\"" << y2 << "\" style=\"stroke:rgb(" << int(color.r()) << "," << int(color.g()) << "," << int(color.b()) << ") \" />" << std::endl;
		}

		inline void markup_line(std::ostream & os, double x1, double x2, double y1, double y2, CGAL_Color color, std::string marker = "")
		{
			os << "  <line x1=\"" << x1 << "\" y1=\"" << y1 << "\" x2=\"" << x2 << "\" y2=\"" << y2 << "\" style=\"stroke:rgb(" << int(color.r()) << "," << int(color.g()) << "," << int(color.b()) << ")\"";
			if (!marker.empty()) {
				os << " marker-end=\"url(#" << marker << ")\"";
			}
			os << "/>" << std::endl;
		}

		inline void markup_point(std::ostream & os, double x, double y, double radius, std::string stroke_color, double stroke_width, std::string fill_color)
		{
			os << "  <circle cx=\"" << x << "\" cy=\"" << y << "\" r=\"" << radius << "\" stroke=\"" << stroke_color << "\" stroke-width=\"" << stroke_width << "\" fill=\"" << fill_color << "\" />" << std::endl;
		}

		inline void markup_polygon(std::ostream & os, std::vector<CGAL_Point_2> & P, CGAL_Color color, double stroke_width, double opacity)
		{
			os << "  <polygon points=\"";
			for (std::vector<CGAL_Point_2>::iterator it_p = P.begin(); it_p != P.end(); it_p++) {
				CGAL_Point_2 pt = (*it_p);
				os << pt.x() << "," << pt.y() << " ";
			}
			os << "\" style=\"fill:rgb(" << int(color.r()) << "," << int(color.g()) << "," << int(color.b()) << ");stroke:red;stroke-width:" << stroke_width << "\" fill-opacity=\"" << opacity << "\" />" << std::endl;
		}

		template<typename T>
		inline void markup_text(std::ostream & os, double x, double y, T text, int font_size, CGAL_Color color)
		{
			os << "  <text x=\"" << x << "\" y=\"" << y << "\" style=\"font-family: Times New Roman; font-size:" << font_size << "; stroke:rgb(" << int(color.r()) << "," << int(color.g()) << "," << int(color.b()) <<
				"); fill:rgb(" << int(color.r()) << "," << int(color.g()) << "," << int(color.b()) << ");\" >";
			os << text;
			os << "</text>" << std::endl;
		}

		inline void markup_footer(std::ostream & os)
		{
			os << "</svg>" << std::endl;
		}

		inline void markup_defs_header(std::ostream & os)
		{
			os << "  <defs>" << std::endl;
		}

		inline void markup_defs_footer(std::ostream & os)
		{
			os << "  </defs>" << std::endl;
		}

		inline void markup_marker(std::ostream & os, std::string id, CGAL_Color color)
		{
			os << "    <marker id=\"" << id << "\" markerWidth=\"10\" markerHeight=\"10\" refX=\"0\" refY=\"3\" orient=\"auto\" markerUnits=\"strokeWidth\">" << std::endl;
			os << "      <path d=\"M0, 0 L0, 6 L9, 3 z\" fill=\"rgb(" << int(color.r()) << "," << int(color.g()) << "," << int(color.b()) << ")\" />" << std::endl;
			os << "    </marker>" << std::endl;
		}
	}
}