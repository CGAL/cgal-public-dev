#pragma once
#include <iostream>
#include <fstream>
#include "defs.h"
#include "matrix.h"
#include "line_item.h"
#include "segment_ray.h"
#include "partition_elements.h"


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

	inline void markup_line(std::ostream & os, LineItem* l)
	{
		os << "  <line x1=\"" << l->x1 << "\" y1=\"" << l->y1 << "\" x2=\"" << l->x2 << "\" y2=\"" << l->y2 << "\" style=\"stroke:rgb(" << int(l->r) << "," << int(l->g) << "," << int(l->b) << ");stroke-width:1\" />" << std::endl;
	}

	inline void markup_line(std::ostream & os, LineItem* l, int r, int g, int b)
	{
		os << "  <line x1=\"" << l->x1 << "\" y1=\"" << l->y1 << "\" x2=\"" << l->x2 << "\" y2=\"" << l->y2 << "\" style=\"stroke:rgb(" << r << "," << g << "," << b << ");stroke-width:1\" />" << std::endl;
	}

	inline void markup_line(std::ostream & os, double x1, double x2, double y1, double y2, int r, int g, int b)
	{
		os << "  <line x1=\"" << x1 << "\" y1=\"" << y1 << "\" x2=\"" << x2 << "\" y2=\"" << y2 << "\" style=\"stroke:rgb(" << r << "," << g << "," << b << ");stroke-width:1\" />" << std::endl;
	}

	inline void markup_point(std::ostream & os, double x, double y, double radius, std::string stroke_color, double stroke_width, std::string fill_color)
	{
		os << "  <circle cx=\"" << x << "\" cy=\"" << y << "\" r=\"" << radius << "\" stroke=\"" << stroke_color << "\" stroke-width=\"" << stroke_width << "\" fill=\"" << fill_color << "\" />" << std::endl;
	}

	inline void markup_polygon(std::ostream & os, list<Point2d> & points, int r, int g, int b, double stroke_width)
	{
		os << "  <polygon points=\"";
		for (list<Point2d>::iterator it_v = points.begin() ; it_v != points.end() ; it_v++) {
			Point2d pt = (*it_v);
			os << pt.x << "," << pt.y << " ";
		}
		os << "\" style=\"fill:rgb(" << r << "," << g << "," << b << ");stroke:red;stroke-width:" << stroke_width << "\" />" << std::endl;
	}

	inline void markup_polygon(std::ostream & os, int rows, int cols, Face* f, int r, int g, int b, double stroke_width, double opacity)
	{
		os << "  <polygon points=\"";
		for (list<pair<Vertex *, Vertex_Type> >::iterator it_v = f->vertices.begin() ; it_v != f->vertices.end() ; it_v++) {
			Vertex* v = it_v->first;
			os << v->pt.x << "," << (rows - v->pt.y) << " ";
		}
		os << "\" style=\"fill:rgb(" << r << "," << g << "," << b << ");stroke:red;stroke-width:" << stroke_width << "\" fill-opacity=\"" << opacity << "\" />" << std::endl;
	}

	inline void markup_footer(std::ostream & os)
	{
		os << "</svg>" << std::endl;
	}


	inline void print_line_items(std::ostream & os, std::list<LineItem *> & lines, int r = -1, int g = -1, int b = -1)
	{
		if (r == -1 && g == -1 && b == -1) {
			for (std::list<LineItem *>::iterator it_l = lines.begin() ; it_l != lines.end() ; it_l++) {
				SVG::markup_line(os, (*it_l));
			}
		} else {
			for (std::list<LineItem *>::iterator it_l = lines.begin() ; it_l != lines.end() ; it_l++) {
				SVG::markup_line(os, (*it_l), r, g, b);
			}
		}
	}


	inline void print_line_items(std::ostream & os, std::list<LineItem *> & lines, int model_rows,
		double xa, double xb, double ya, double yb, int r = -1, int g = -1, int b = -1)
	{
		for (std::list<LineItem *>::iterator it_l = lines.begin() ; it_l != lines.end() ; it_l++) {
			LineItem* l = (*it_l);
			// Reverts line item to use the natural frame
			double x1 = l->x1, y1 = model_rows - l->y1, x2 = l->x2, y2 = model_rows - l->y2;


			double dx = x2 - x1, dy = y2 - y1;
			Vec2d v_12 = cv::normalize(Vec2d(dx, dy));
			if ((x1 < xa && x2 < xa) || (x1 > xb && x2 > xb) || (y1 < ya && y2 < ya) || (y1 > yb && y2 > yb)) continue;

			bool p1_located = false, p2_located = false;
			Point2d p1, p2;
			if (x1 >= xa && x1 <= xb && y1 >= ya && y1 <= yb) {
				p1 = Point2d(x1, y1);
				p1_located = true;
			} else {
				double t = FLT_MAX, u;
				if (fabs(dx) > 0) {
					u = (xa - x1) / v_12[0];
					if (u >= 0 && u < t) {
						double yu = y1 + u * v_12[1];
						if (yu >= ya && yu <= yb) t = u;
					}
					u = (xb - x1) / v_12[0];
					if (u >= 0 && u < t) {
						double yu = y1 + u * v_12[1];
						if (yu >= ya && yu <= yb) t = u;
					}
				}
				if (fabs(dy) > 0) {
					u = (ya - y1) / v_12[1];
					if (u >= 0 && u < t) {
						double xu = x1 + u * v_12[0];
						if (xu >= xa && xu <= xb) t = u;
					}
					u = (yb - y1) / v_12[1];
					if (u >= 0 && u < t) {
						double xu = x1 + u * v_12[0];
						if (xu >= xa && xu <= xb) t = u;
					}
				}
				if (t != FLT_MAX) {
					p1 = Point2d(x1 + t * v_12[0], y1 + t * v_12[1]);
					p1_located = true;
				}
			}
			if (!p1_located) continue;
			
			if (x2 >= xa && x2 <= xb && y2 >= ya && y2 <= yb) {
				p2 = Point2d(x2, y2);
				p2_located = true;
			} else {
				double t = FLT_MAX, u;
				if (fabs(dx) > 0) {
					u = -(xa - x2) / v_12[0];
					if (u >= 0 && u < t) {
						double yu = y2 - u * v_12[1];
						if (yu >= ya && yu <= yb) t = u;
					}
					u = -(xb - x2) / v_12[0];
					if (u >= 0 && u < t) {
						double yu = y2 - u * v_12[1];
						if (yu >= ya && yu <= yb) t = u;
					}
				}
				if (fabs(dy) > 0) {
					u = -(ya - y2) / v_12[1];
					if (u >= 0 && u < t) {
						double xu = x2 - u * v_12[0];
						if (xu >= xa && xu <= xb) t = u;
					}
					u = -(yb - y2) / v_12[1];
					if (u >= 0 && u < t) {
						double xu = x2 - u * v_12[0];
						if (xu >= xa && xu <= xb) t = u;
					}
				}
				if (t != FLT_MAX) {
					p2 = Point2d(x2 - t * v_12[0], y2 - t * v_12[1]);
					p2_located = true;
				}
			}
			if (!p2_located) continue;
			

			if (r == -1 && g == -1 && b == -1) {
				SVG::markup_line(os, p1.x - xa, p2.x - xa, yb - p1.y, yb - p2.y, l->r, l->g, l->b);
			} else {
				SVG::markup_line(os, p1.x - xa, p2.x - xa, yb - p1.y, yb - p2.y, 255, 0, 0);
			}
		}
	}


	inline void print_segments(std::ostream & os, std::vector<Segment *> & segments, int coords,
		double xa, double xb, double ya, double yb, int r, int g, int b)
	{
		for (uint i = 0; i < segments.size(); i++) {
			Segment* s = segments[i];
			double x1, y1, x2, y2;
			switch (coords) {
			case 0: x1 = s->end1.x; y1 = s->end1.y; x2 = s->end2.x; y2 = s->end2.y; break;
			case 1: x1 = s->interEnd1.x; y1 = s->interEnd1.y; x2 = s->interEnd2.x; y2 = s->interEnd2.y; break;
			case 2: x1 = s->finalEnd1.x; y1 = s->finalEnd1.y; x2 = s->finalEnd2.x; y2 = s->finalEnd2.y; break;
			}
			double dx = x2 - x1, dy = y2 - y1;
			Vec2d v_12 = cv::normalize(Vec2d(dx, dy));
			if ((x1 < xa && x2 < xa) || (x1 > xb && x2 > xb)
				|| (y1 < ya && y2 < ya) || (y1 > yb && y2 > yb)) continue;

			Point2d p1;
			if (x1 >= xa && x1 <= xb && y1 >= ya && y1 <= yb) {
				p1 = Point2d(x1, y1);
			} else {
				double t = FLT_MAX, u;
				if (fabs(dx) > 0) {
					u = (xa - x1) / v_12[0];
					if (u >= 0 && u < t) {
						double yu = y1 + u * v_12[1];
						if (yu >= ya && yu <= yb) t = u;
					}
					u = (xb - x1) / v_12[0];
					if (u >= 0 && u < t) {
						double yu = y1 + u * v_12[1];
						if (yu >= ya && yu <= yb) t = u;
					}
				}
				if (fabs(dy) > 0) {
					u = (ya - y1) / v_12[1];
					if (u >= 0 && u < t) {
						double xu = x1 + u * v_12[0];
						if (xu >= xa && xu <= xb) t = u;
					}
					u = (yb - y1) / v_12[1];
					if (u >= 0 && u < t) {
						double xu = x1 + u * v_12[0];
						if (xu >= xa && xu <= xb) t = u;
					}
				}
				p1 = Point2d(x1 + t * v_12[0], y1 + t * v_12[1]);
			}

			Point2d p2;
			if (x2 >= xa && x2 <= xb && y2 >= ya && y2 <= yb) {
				p2 = Point2d(x2, y2);
			} else {
				double t = FLT_MAX, u;
				if (fabs(dx) > 0) {
					u = -(xa - x2) / v_12[0];
					if (u >= 0 && u < t) {
						double yu = y2 - u * v_12[1];
						if (yu >= ya && yu <= yb) t = u;
					}
					u = -(xb - x2) / v_12[0];
					if (u >= 0 && u < t) {
						double yu = y2 - u * v_12[1];
						if (yu >= ya && yu <= yb) t = u;
					}
				}
				if (fabs(dy) > 0) {
					u = -(ya - y2) / v_12[1];
					if (u >= 0 && u < t) {
						double xu = x2 - u * v_12[0];
						if (xu >= xa && xu <= xb) t = u;
					}
					u = -(yb - y2) / v_12[1];
					if (u >= 0 && u < t) {
						double xu = x2 - u * v_12[0];
						if (xu >= xa && xu <= xb) t = u;
					}
				}
				p2 = Point2d(x2 - t * v_12[0], y2 - t * v_12[1]);
			}

			SVG::markup_line(os, p1.x - xa, p2.x - xa, yb - p1.y, yb - p2.y, r, g, b);
		}
	}
}