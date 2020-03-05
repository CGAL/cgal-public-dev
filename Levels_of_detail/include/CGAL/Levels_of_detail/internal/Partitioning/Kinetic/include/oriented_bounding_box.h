#pragma once
#include "defs_cgal.h"


namespace Skippy 
{
	class Oriented_Bounding_Box
	{
	public:
		Oriented_Bounding_Box();

		// Oriented_Bounding_Box(const std::string & point_cloud);

		// Oriented_Bounding_Box(const std::vector<CGAL_Inexact_Point_3> & pts);

		virtual ~Oriented_Bounding_Box();

		void get_bounding_box(const std::string & point_cloud, const std::vector<CGAL_Plane> & SP);

		virtual CGAL_Point_3 transform(const CGAL_Point_3 & M) const = 0;

		virtual CGAL_Point_3 backtransform(const CGAL_Point_3 & M) const = 0;

		CGAL_Plane transform(const CGAL_Plane & H) const;

		CGAL_Plane backtransform(const CGAL_Plane & H) const;

	protected:
		// void search_best_2d_rotation_matrix(const std::vector<CGAL_Inexact_Point_3> & pts);

		// void search_best_3d_rotation_matrix(const std::vector<CGAL_Inexact_Point_3> & pts);

		void reinit_exact_rotation_matrices();

		void set_exact_rotation_matrices();

	protected:
		double m_rot[3][3];

	public:
		bool transformation_matrices_are_set;
		FT R[3][3];
		FT Inv_R[3][3];
	};
}