#include "../include/oriented_bounding_box_2.h"
#include "../include/ply_in.h"
#include "../include/defs.h"
#include <CGAL/convex_hull_2.h>



namespace Skippy 
{
	Oriented_Bounding_Box_2::Oriented_Bounding_Box_2()
		: Oriented_Bounding_Box()
	{
		reinit_exact_rotation_matrices();
	}


	Oriented_Bounding_Box_2::Oriented_Bounding_Box_2(const std::string & point_cloud)
		: Oriented_Bounding_Box_2()
	{
		std::vector<CGAL_Inexact_Point_3> points;
		Ply_In::read(point_cloud, points);

		search_best_2d_rotation_matrix(points);
		set_exact_rotation_matrices();
	}


	Oriented_Bounding_Box_2::Oriented_Bounding_Box_2(const std::vector<CGAL_Inexact_Point_3> & pts)
		: Oriented_Bounding_Box_2()
	{
		search_best_2d_rotation_matrix(pts);
		set_exact_rotation_matrices();
	}



	Oriented_Bounding_Box_2::~Oriented_Bounding_Box_2()
	{
	}


	
	void Oriented_Bounding_Box_2::search_best_2d_rotation_matrix(const std::vector<CGAL_Inexact_Point_3> & pts)
	{
		// We are going to perform a rotation around the z-axis
		// We project all points onto a horizontal 2D plane

		std::vector<CGAL_Inexact_Point_2> projected;
		projected.reserve(pts.size());

		for (size_t i = 0 ; i < pts.size() ; ++i) {
			double x_i = pts[i].x(), y_i = pts[i].y();
			projected.push_back(CGAL_Inexact_Point_2(x_i, y_i));
		}

		std::vector<CGAL_Inexact_Point_2> hull_0, hull;
		CGAL::convex_hull_2(projected.begin(), projected.end(), std::back_inserter(hull_0));
		projected.clear();

		double x = 0, y = 0;
		for (size_t i = 0 ; i < hull_0.size() ; ++i) {
			double x_i = hull_0[i].x(), y_i = hull_0[i].y();
			x += x_i, y += y_i;
		}
		x /= hull_0.size();
		y /= hull_0.size();

		hull.reserve(hull_0.size());
		for (size_t i = 0 ; i < hull_0.size() ; ++i) {
			double x_i = hull_0[i].x(), y_i = hull_0[i].y();
			hull.push_back(CGAL_Inexact_Point_2(x_i - x, y_i - y));
		}

		double min_product_dims = FLT_MAX;
		int argmin_theta_d = -1;

		for (int theta_d = 0 ; theta_d <= 85 ; theta_d += 5) {

			double x_min = FLT_MAX, x_max = -FLT_MAX;
			double y_min = FLT_MAX, y_max = -FLT_MAX;

			for (size_t i = 0 ; i < hull.size() ; ++i) {
				double cos_theta = cos(theta_d * 3.1415926358978 / 180);
				double sin_theta = sin(theta_d * 3.1415926358978 / 180);
				double x_i = cos_theta * hull[i].x() - sin_theta * hull[i].y();
				double y_i = sin_theta * hull[i].x() + cos_theta * hull[i].y();
				if (x_i < x_min) x_min = x_i;
				if (x_i > x_max) x_max = x_i;
				if (y_i < y_min) y_min = y_i;
				if (y_i > y_max) y_max = y_i;
			}

			double product_dims = (x_max - x_min) * (y_max - y_min);

			if (product_dims < min_product_dims) {
				min_product_dims = product_dims;
				argmin_theta_d = theta_d;
			}
		}
		
		std::cout << "** Theta = " << argmin_theta_d << std::endl;

		double theta = argmin_theta_d * 3.1415926358978 / 180;
		m_rot[0][0] = cos(theta); m_rot[0][1] = -sin(theta); m_rot[0][2] = 0;
		m_rot[1][0] = sin(theta); m_rot[1][1] = cos(theta); m_rot[1][2] = 0;
		m_rot[2][0] = 0; m_rot[2][1] = 0; m_rot[2][2] = 1;
	}


	CGAL_Point_3 Oriented_Bounding_Box_2::transform(const CGAL_Point_3 & M) const
	{
		FT x = R[0][0] * M.x() + R[0][1] * M.y() + R[0][2] * M.z();
		FT y = R[1][0] * M.x() + R[1][1] * M.y() + R[1][2] * M.z();
		FT z = R[2][0] * M.x() + R[2][1] * M.y() + R[2][2] * M.z();

		return CGAL_Point_3(x, y, z);
	}


	CGAL_Point_3 Oriented_Bounding_Box_2::backtransform(const CGAL_Point_3 & M) const
	{
		FT x = Inv_R[0][0] * M.x() + Inv_R[0][1] * M.y() + Inv_R[0][2] * M.z();
		FT y = Inv_R[1][0] * M.x() + Inv_R[1][1] * M.y() + Inv_R[1][2] * M.z();
		FT z = Inv_R[2][0] * M.x() + Inv_R[2][1] * M.y() + Inv_R[2][2] * M.z();

		return CGAL_Point_3(x, y, z);	
	}
}