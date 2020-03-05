#pragma once
#include "oriented_bounding_box.h"


namespace Skippy {

	class Oriented_Bounding_Box_3 : public Oriented_Bounding_Box
	{
	public:
		Oriented_Bounding_Box_3();

		Oriented_Bounding_Box_3(const std::string & point_cloud);

		Oriented_Bounding_Box_3(const std::vector<CGAL_Inexact_Point_3> & pts);
		
		~Oriented_Bounding_Box_3();

		CGAL_Point_3 transform(const CGAL_Point_3 & M) const;

		CGAL_Point_3 backtransform(const CGAL_Point_3 & M) const;

	protected:
		void search_best_3d_rotation_matrix(const std::vector<CGAL_Inexact_Point_3> & pts);

		void reinit_exact_translation_vector();

		void set_exact_translation_vector();

	protected:
		double m_pos[3];

	public:
		bool translation_vector_is_set;
		FT T[3];
	};

}