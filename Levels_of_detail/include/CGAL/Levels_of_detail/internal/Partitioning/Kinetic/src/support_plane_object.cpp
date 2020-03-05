#include "../include/support_plane_object.h"
#include "../include/vars.h"
#include <iostream>


namespace Skippy {

	Support_Plane_Object::Support_Plane_Object(const int _id_plane)
		: id_object(++Counters::id_objects[_id_plane]),
		id_plane(_id_plane)
	{

	}


	Support_Plane_Object::~Support_Plane_Object()
	{

	}

}