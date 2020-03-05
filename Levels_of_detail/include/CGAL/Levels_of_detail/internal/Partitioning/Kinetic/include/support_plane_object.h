#pragma once


namespace Skippy {
	class Support_Plane_Object
	{
	protected:
		Support_Plane_Object(const int _id_plane);

	public:
		virtual ~Support_Plane_Object();

	public:
		const int id_object;
		const int id_plane;
	};
}