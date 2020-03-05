#include "../include/oriented_bounding_box.h"
#include "../include/ply_in.h"
#include "../include/defs.h"
#include <iostream>
#include <fstream>



namespace Skippy
{
	Oriented_Bounding_Box::Oriented_Bounding_Box()
	{
	}


	Oriented_Bounding_Box::~Oriented_Bounding_Box()
	{
	}



	void Oriented_Bounding_Box::get_bounding_box(const std::string & point_cloud, const std::vector<CGAL_Plane> & SP)
	{
		FT x_min = -SP[0].d(), x_max = -SP[1].d();
		FT y_min = -SP[2].d(), y_max = -SP[3].d();
		FT z_min = -SP[4].d(), z_max = -SP[5].d();

		std::cout << x_min << " " << x_max << " " << y_min << " " << y_max << " " << z_min << " " << z_max << std::endl;

		std::vector<CGAL_Point_3> P (8);
		P[0] = CGAL_Point_3(x_min, y_min, z_min);
		P[1] = CGAL_Point_3(x_min, y_max, z_min);
		P[2] = CGAL_Point_3(x_max, y_max, z_min);
		P[3] = CGAL_Point_3(x_max, y_min, z_min);
		P[4] = CGAL_Point_3(x_min, y_min, z_max);
		P[5] = CGAL_Point_3(x_min, y_max, z_max);
		P[6] = CGAL_Point_3(x_max, y_max, z_max);
		P[7] = CGAL_Point_3(x_max, y_min, z_max);

		std::vector<CGAL_Inexact_Point_3> points;
		Ply_In::read(point_cloud, points);

		int nb_points = 8 + int(points.size());

		std::string name = "bbox.ply";
		std::ofstream stream(name, std::ofstream::out);

		if (stream.is_open()) {
			stream << "ply" << std::endl
				<< "format ascii 1.0" << std::endl;
			stream << "element vertex " << nb_points << std::endl
				<< "property float x" << std::endl
				<< "property float y" << std::endl
				<< "property float z" << std::endl;
			stream << "element face 6" << std::endl
				<< "property list uchar int vertex_index" << std::endl
				<< "end_header" << std::endl;
			for (size_t i = 0 ; i < 8 ; ++i) {
				CGAL_Point_3 M = P[i];
				stream << M.x() << " " << M.y() << " " << M.z() << std::endl;
			}
			for (size_t i = 0 ; i < points.size() ; ++i) {
				CGAL_Point_3 M (points[i].x(), points[i].y(), points[i].z());
				M = transform(M);
				stream << M.x() << " " << M.y() << " " << M.z() << std::endl;
			}
			stream << "4 0 1 2 3" << std::endl
				<< "4 7 6 5 4" << std::endl
				<< "4 0 4 5 1" << std::endl
				<< "4 1 5 6 2" << std::endl
				<< "4 2 6 7 3" << std::endl
				<< "4 3 7 4 0" << std::endl;
			stream.close();
		}
	}


	void Oriented_Bounding_Box::reinit_exact_rotation_matrices()
	{
		for (size_t i = 0 ; i < 3 ; ++i) {
			for (size_t j = 0 ; j < 3 ; ++j) {
				if (i == j) {
					m_rot[i][j] = 1;
					R[i][j] = FT(1);
					Inv_R[i][j] = FT(1);
				} else {
					m_rot[i][j] = 0;
					R[i][j] = FT(0);
					Inv_R[i][j] = FT(0);
				}
			}
		}

		transformation_matrices_are_set = false;
	}


	void Oriented_Bounding_Box::set_exact_rotation_matrices()
	{
		// Usual frame to transformed frame
		std::stringstream stream;
		stream.precision(6);

		for (size_t i = 0 ; i < 3 ; ++i) {
			for (size_t j = 0 ; j < 3 ; ++j) {
				stream << m_rot[i][j] << " ";
			}
		}

		for (size_t i = 0 ; i < 3 ; ++i) {
			for (size_t j = 0 ; j < 3 ; ++j) {
				stream >> R[i][j];
			}
		}
		
		const FT &a = R[0][0], &b = R[0][1], &c = R[0][2],
			&d = R[1][0], &e = R[1][1], &f = R[1][2],
			&g = R[2][0], &h = R[2][1], &i = R[2][2];
		FT det_R = a * (e * i - f * h) - d * (b * i - c * h) + g * (b * f - c * e);

		// Computes Inv_R

		if (det_R == 0) {
			reinit_exact_rotation_matrices();
			return;
		}

		Inv_R[0][0] = (e * i - f * h) / det_R;
		Inv_R[0][1] = -(b * i - c * h) / det_R;
		Inv_R[0][2] = (b * f - c * e) / det_R;

		Inv_R[1][0] = -(d * i - f * g) / det_R;
		Inv_R[1][1] = (a * i - c * g) / det_R;
		Inv_R[1][2] = -(a * f - c * d) / det_R;

		Inv_R[2][0] = (d * h - e * g) / det_R;
		Inv_R[2][1] = -(a * h - b * g) / det_R;
		Inv_R[2][2] = (a * e - b * d) / det_R;

		for (size_t i = 0 ; i < 3 ; ++i) {
			for (size_t j = 0 ; j < 3 ; ++j) {
				FT sum_1 = 0, sum_2 = 0;
				for (size_t k = 0 ; k < 3 ; ++k) {
					sum_1 += R[i][k] * Inv_R[k][j];
					sum_2 += Inv_R[i][k] * R[k][j];
				}
				if (i == j) {
					assert(sum_1 == 1);
					assert(sum_2 == 1);
				} else {
					assert(sum_1 == 0);
					assert(sum_2 == 0);
				}
			}
		}

		transformation_matrices_are_set = true;
	}



	CGAL_Plane Oriented_Bounding_Box::transform(const CGAL_Plane & H) const
	{
		CGAL_Point_3 A = H.point();
		CGAL_Point_3 B = A + H.base1(), C = A + H.base2();
		CGAL_Point_3 a = transform(A), b = transform(B), c = transform(C);
		CGAL_Plane P(a, b, c);

		return P;
	}



	CGAL_Plane Oriented_Bounding_Box::backtransform(const CGAL_Plane & H) const
	{
		CGAL_Point_3 A = H.point();
		CGAL_Point_3 B = A + H.base1(), C = A + H.base2();
		CGAL_Point_3 a = backtransform(A), b = backtransform(B), c = backtransform(C);
		CGAL_Plane P(a, b, c);

		return P;
	}
}