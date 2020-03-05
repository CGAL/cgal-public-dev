#pragma once
#include "defs.h"
#include "event.h"
#include "support_plane_object.h"
#include "segment_translation.h"
#include <list>
#include <map>
#include <vector>

namespace Skippy 
{
	class Intersection_Line;
	class Polygon_Vertex;
	class Polygon_Vertex_R;
	class Polygon_Vertex_S;

	class Polygon_Edge;
	class Polygon;
	class Polygon_Tree;
	class Polygon_Group;
	class Segment;
	class Planar_Segment;
	class Polygon_Segment;
	class Polygon_Segment_R;
	class Polygon_Segment_S;

	typedef std::pair<Polygon_Vertex_R*, bool> Companion;
	typedef std::pair<Intersection_Line*, Sign> Constraint;
}