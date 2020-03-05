#include "../include/oriented_bounding_box_3.h"
#include "../include/ply_in.h"
#include "../include/defs.h"
#include <CGAL/convex_hull_3.h>
#include <CGAL/Polyhedron_3.h>


namespace Skippy 
{
	Oriented_Bounding_Box_3::Oriented_Bounding_Box_3()
		: Oriented_Bounding_Box()
	{
		reinit_exact_rotation_matrices();
		reinit_exact_translation_vector();
	}


	Oriented_Bounding_Box_3::Oriented_Bounding_Box_3(const std::string & point_cloud)
		: Oriented_Bounding_Box_3()
	{
		std::vector<CGAL_Inexact_Point_3> points;
		Ply_In::read(point_cloud, points);

		search_best_3d_rotation_matrix(points);
		set_exact_rotation_matrices();
		set_exact_translation_vector();
	}


	Oriented_Bounding_Box_3::Oriented_Bounding_Box_3(const std::vector<CGAL_Inexact_Point_3> & points)
		: Oriented_Bounding_Box_3()
	{
		search_best_3d_rotation_matrix(points);
		set_exact_rotation_matrices();
		set_exact_translation_vector();
	}


	Oriented_Bounding_Box_3::~Oriented_Bounding_Box_3()
	{
	}



	void Oriented_Bounding_Box_3::search_best_3d_rotation_matrix(const std::vector<CGAL_Inexact_Point_3> & pts)
	{
		typedef CGAL::Polyhedron_3<CGAL::Exact_predicates_inexact_constructions_kernel> CGAL_Inexact_Polyhedron_3;

		size_t n = pts.size();

		// Step 1.
		// Recenters points

		double x_0 = 0, y_0 = 0, z_0 = 0;
		for (size_t i = 0 ; i < n ; ++i) {
			x_0 += pts[i].x(), y_0 += pts[i].y(), z_0 += pts[i].z();
		}
		x_0 /= n, y_0 /= n, z_0 /= n;
		m_pos[0] = x_0, m_pos[1] = y_0, m_pos[2] = z_0;

		std::vector<CGAL_Inexact_Point_3> rc_pts;
		CGAL_Inexact_Polyhedron_3 polyhedron_ch_rc_pts;
		std::vector<CGAL_Inexact_Point_3> ch_rc_pts;

		rc_pts.reserve(n);
		for (size_t i = 0 ; i < n ; ++i) {
			const double &x = pts[i].x(), &y = pts[i].y(), &z = pts[i].z();
			rc_pts.push_back(CGAL_Inexact_Point_3(x - x_0, y - y_0, z - z_0));
		}

		CGAL::convex_hull_3(rc_pts.begin(), rc_pts.end(), polyhedron_ch_rc_pts);
		ch_rc_pts.reserve(polyhedron_ch_rc_pts.size_of_vertices());
		for (auto it_v = polyhedron_ch_rc_pts.vertices_begin() ; it_v != polyhedron_ch_rc_pts.vertices_end() ; ++it_v) {
			const CGAL_Inexact_Point_3 & pt = it_v->point();
			ch_rc_pts.push_back(pt);
		}

		// Step 2.
		// We consider couples (theta, phi) of rotation angles among the y-axis and the z-axis. 
		// We apply all transformations R(theta, phi) to rc_pts and obtain the couple (theta, phi)
		// that minimizes the volume of the bounding box for the rotated points.

		double min_volume = FLT_MAX;
		int argmin_volume_alpha = 0, argmin_volume_beta = 0, argmin_volume_gamma = 0;

		std::vector<double> t_cos(180), t_sin(180);
		for (int i = 0 ; i < 180 ; ++i) {
			double theta_i = PI * i / 180;
			t_cos[i] = cos(theta_i), t_sin[i] = sin(theta_i);
		}

		for (int alpha_d = 0 ; alpha_d < 90 ; alpha_d += 1) {
			double cos_alpha = t_cos[alpha_d], sin_alpha = t_sin[alpha_d];
			double rot_x[3][3] = { {1, 0, 0}, {0, cos_alpha, -sin_alpha}, {0, sin_alpha, cos_alpha} };

			for (int beta_d = 0; beta_d < 90; beta_d += 1) {
				double sin_beta = t_sin[beta_d], cos_beta = t_cos[beta_d];
				double rot_y[3][3] = { {cos_beta, 0, sin_beta}, {0, 1, 0}, {-sin_beta, 0, cos_beta} };

				double rot_xy[3][3] = { {0, 0, 0}, {0, 0, 0}, {0, 0, 0} };
				for (int i = 0; i < 3; ++i) {
					for (int j = 0; j < 3; ++j) {
						for (int k = 0; k < 3; ++k) rot_xy[i][j] += rot_x[i][k] * rot_y[k][j];
					}
				}

				for (int gamma_d = 0; gamma_d < 90; gamma_d += 1) {
					double sin_gamma = t_sin[gamma_d], cos_gamma = t_cos[gamma_d];
					double rot_z[3][3] = { {cos_gamma, -sin_gamma, 0}, {sin_gamma, cos_gamma, 0}, {0, 0, 1} };

					double rot_xyz[3][3] = { {0, 0, 0}, {0, 0, 0}, {0, 0, 0} };
					for (int i = 0; i < 3; ++i) {
						for (int j = 0; j < 3; ++j) {
							for (int k = 0; k < 3; ++k) rot_xyz[i][j] += rot_xy[i][k] * rot_z[k][j];
						}
					}

					double x_min = FLT_MAX, x_max = -x_min;
					double y_min = FLT_MAX, y_max = -y_min;
					double z_min = FLT_MAX, z_max = -z_min;
					for (size_t i = 0; i < ch_rc_pts.size(); ++i) {
						const CGAL_Inexact_Point_3 & pt = ch_rc_pts[i];
						const double &x_i = pt.x(), &y_i = pt.y(), &z_i = pt.z();
						double x = rot_xyz[0][0] * x_i + rot_xyz[0][1] * y_i + rot_xyz[0][2] * z_i;
						double y = rot_xyz[1][0] * x_i + rot_xyz[1][1] * y_i + rot_xyz[1][2] * z_i;
						double z = rot_xyz[2][0] * x_i + rot_xyz[2][1] * y_i + rot_xyz[2][2] * z_i;
						x_min = jmin(x, x_min); x_max = jmax(x, x_max);
						y_min = jmin(y, y_min); y_max = jmax(y, y_max);
						z_min = jmin(z, z_min); z_max = jmax(z, z_max);
					}

					double volume = (x_max - x_min) * (y_max - y_min) * (z_max - z_min);
					if (volume < min_volume) {
						min_volume = volume;
						argmin_volume_alpha = alpha_d;
						argmin_volume_beta = beta_d;
						argmin_volume_gamma = gamma_d;
						for (int i = 0; i < 3; ++i) {
							for (int j = 0; j < 3; ++j) {
								m_rot[i][j] = rot_xyz[i][j];
							}
						}
					}
				}
			}
		}

		std::cout << "** alpha = " << argmin_volume_alpha << ", beta = " << argmin_volume_beta << ", gamma = " << argmin_volume_gamma << std::endl;
	}


	void Oriented_Bounding_Box_3::reinit_exact_translation_vector()
	{
		m_pos[0] = 0, m_pos[1] = 0, m_pos[2] = 0;
		T[0] = FT(0), T[1] = FT(0), T[2] = FT(0);

		translation_vector_is_set = false;
	}


	void Oriented_Bounding_Box_3::set_exact_translation_vector()
	{
		std::stringstream stream;
		stream.precision(6);

		for (size_t i = 0 ; i < 3 ; ++i) {
			stream << m_pos[i] << " ";
		}

		for (size_t i = 0 ; i < 3 ; ++i) {
			stream >> T[i];
		}

		translation_vector_is_set = true;
	}


	CGAL_Point_3 Oriented_Bounding_Box_3::transform(const CGAL_Point_3 & M) const
	{
		CGAL_Point_3 M_t = CGAL_Point_3(M.x() - T[0], M.y() - T[1], M.z() - T[2]);

		FT x = R[0][0] * M_t.x() + R[0][1] * M_t.y() + R[0][2] * M_t.z();
		FT y = R[1][0] * M_t.x() + R[1][1] * M_t.y() + R[1][2] * M_t.z();
		FT z = R[2][0] * M_t.x() + R[2][1] * M_t.y() + R[2][2] * M_t.z();

		return CGAL_Point_3(x, y, z);
	}


	CGAL_Point_3 Oriented_Bounding_Box_3::backtransform(const CGAL_Point_3 & M) const
	{
		FT x = Inv_R[0][0] * M.x() + Inv_R[0][1] * M.y() + Inv_R[0][2] * M.z();
		FT y = Inv_R[1][0] * M.x() + Inv_R[1][1] * M.y() + Inv_R[1][2] * M.z();
		FT z = Inv_R[2][0] * M.x() + Inv_R[2][1] * M.y() + Inv_R[2][2] * M.z();

		return CGAL_Point_3(x + T[0], y + T[1], z + T[2]);	
	}
}