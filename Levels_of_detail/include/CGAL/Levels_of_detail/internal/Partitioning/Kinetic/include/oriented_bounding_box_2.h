#pragma once
#include "oriented_bounding_box.h"


namespace Skippy 
{
	class Oriented_Bounding_Box_2 : public Oriented_Bounding_Box
	{
	public:
		Oriented_Bounding_Box_2();

		Oriented_Bounding_Box_2(const std::string & point_cloud);

		Oriented_Bounding_Box_2(const std::vector<CGAL_Inexact_Point_3> & pts);

		~Oriented_Bounding_Box_2();

		CGAL_Point_3 transform(const CGAL_Point_3 & M) const;

		CGAL_Point_3 backtransform(const CGAL_Point_3 & M) const;

	protected:
		void search_best_2d_rotation_matrix(const std::vector<CGAL_Inexact_Point_3> & pts);

	public:
	};
}