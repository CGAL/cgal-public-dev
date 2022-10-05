// Copyright (c) 2022 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Hossam Saeed
//

#ifndef CGAL_POLYGON_MESH_PROCESSING_INTERPOLATED_CORRECTED_CURVATURE_MEASURES_H
#define CGAL_POLYGON_MESH_PROCESSING_INTERPOLATED_CORRECTED_CURVATURE_MEASURES_H

#include <CGAL/license/Polygon_mesh_processing/interpolated_corrected_curvature_measures.h>

#include <CGAL/assertions.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/Named_function_parameters.h>
#include <CGAL/property_map.h>
#include <CGAL/boost/graph/named_params_helper.h>
#include <Eigen/Eigenvalues>

#include <numeric>
#include <queue>
#include <unordered_set>
#include <functional>

#define EXPANDING_RADIUS_EPSILON 1e-6

namespace CGAL {

namespace Polygon_mesh_processing {

template<typename GT>
struct Principal_curvature {
    typename GT::FT min_curvature;
    typename GT::FT max_curvature;
    typename GT::Vector_3 min_direction;
    typename GT::Vector_3 max_direction;

    Principal_curvature() {
        min_curvature = 0;
        max_curvature = 0;
        min_direction = typename GT::Vector_3(0, 0, 0);
        max_direction = typename GT::Vector_3(0, 0, 0);
    }

    Principal_curvature(
        typename GT::FT min_curvature,
        typename GT::FT max_curvature,
        typename GT::Vector_3 min_direction,
        typename GT::Vector_3 max_direction) {
        this->min_curvature = min_curvature;
        this->max_curvature = max_curvature;
        this->min_direction = min_direction;
        this->max_direction = max_direction;
    }
};

namespace internal {

enum Curvature_measure_index {
    MU0_AREA_MEASURE, ///< corrected area density
    MU1_MEAN_CURVATURE_MEASURE, ///< corrected mean curvature density
    MU2_GAUSSIAN_CURVATURE_MEASURE ///< corrected gaussian curvature density
};

template<typename GT>
struct Vertex_curvatures {
    typename GT::FT area_measure = 0;
    typename GT::FT mean_curvature_measure = 0;
    typename GT::FT gaussian_curvature_measure = 0;
    std::array<typename GT::FT, 3 * 3> anisotropic_measure = { 0, 0, 0,
                                                               0, 0, 0,
                                                               0, 0, 0 };
};


template<typename PolygonMesh, typename GT>
typename GT::FT average_edge_length(const PolygonMesh& pmesh)
{
    const std::size_t n = edges(pmesh).size();
    if (n == 0)
        return 0;

    GT::FT avg_edge_length = 0;
    for (auto e : edges(pmesh))
        avg_edge_length += edge_length(e, pmesh);

    avg_edge_length /= n;
    return avg_edge_length;
}

template<typename GT>
Principal_curvature<GT> principal_curvature_from_anisotropic_measures(
    const std::array<typename GT::FT, 3 * 3> anisotropic_measure,
    const typename GT::FT v_mu0,
    const typename GT::Vector_3 u_GT
)
{
    Eigen::Matrix<typename GT::FT, 3, 3> v_muXY = Eigen::Matrix<typename GT::FT, 3, 3>::Zero();

    for (std::size_t ix = 0; ix < 3; ix++)
        for (std::size_t iy = 0; iy < 3; iy++)
            v_muXY(ix, iy) = anisotropic_measure[ix * 3 + iy];

    Eigen::Matrix<typename GT::FT, 3, 1> u(u_GT.x(), u_GT.y(), u_GT.z());
    const typename GT::FT K = 1000 * v_mu0;

    v_muXY = 0.5 * (v_muXY + v_muXY.transpose()) + K * u * u.transpose();

    Eigen::SelfAdjointEigenSolver<Eigen::Matrix <typename GT::FT, 3, 3>> eigensolver;

    eigensolver.computeDirect(v_muXY);

    if (eigensolver.info() != Eigen::Success)
        return Principal_curvature<GT>();

    const Eigen::Matrix<typename GT::FT, 3, 1> eig_vals = eigensolver.eigenvalues();
    const Eigen::Matrix<typename GT::FT, 3, 3> eig_vecs = eigensolver.eigenvectors();

    const typename GT::Vector_3 min_eig_vec(eig_vecs(0, 1), eig_vecs(1, 1), eig_vecs(2, 1));
    const typename GT::Vector_3 max_eig_vec(eig_vecs(0, 0), eig_vecs(1, 0), eig_vecs(2, 0));

    return Principal_curvature<GT>(
        (v_mu0 != 0.0) ? -eig_vals[1] / v_mu0 : 0.0,
        (v_mu0 != 0.0) ? -eig_vals[0] / v_mu0 : 0.0,
        min_eig_vec,
        max_eig_vec
        );
}


template<typename GT>
typename GT::FT interpolated_corrected_area_measure_face(const std::vector<typename GT::Vector_3>& u,
                                                         const std::vector<typename GT::Vector_3>& x = {})
{
    const std::size_t n = x.size();
    CGAL_precondition(u.size() == n);
    CGAL_precondition(n >= 3);

    typename GT::Construct_cross_product_vector_3 cross_product;

    // Triangle: use triangle formula
    if (n == 3)
    {
        const typename GT::Vector_3 um = (u[0] + u[1] + u[2]) / 3.0;
        return 0.5 * um * cross_product(x[1] - x[0], x[2] - x[0]);
    }
    // Quad: use bilinear interpolation formula
    else if (n == 4)
    {
        // for the formulas below, values of verices 2 & 3 are swapped (compared to paper) to correct order.
        // the indices in paper vs in here are: 00 = 0, 10 = 1, 11 = 2, 01 = 3
        return (1.0 / 36.0) * (
              (4 * u[0] + 2 * u[1] + 2 * u[3] + u[2]) * cross_product(x[1] - x[0], x[3] - x[0])
            + (2 * u[0] + 4 * u[1] + u[3] + 2 * u[2]) * cross_product(x[1] - x[0], x[2] - x[1])
            + (2 * u[0] + u[1] + 4 * u[3] + 2 * u[2]) * cross_product(x[2] - x[3], x[3] - x[0])
            + (u[0] + 2 * u[1] + 2 * u[3] + 4 * u[2]) * cross_product(x[2] - x[3], x[2] - x[1])
            );
    }
    // N-gon: split into n triangles by polygon center and use triangle formula for each
    else
    {
        typename GT::FT mu0 = 0;

        // getting center of points
        typename GT::Vector_3 xc =
            std::accumulate(x.begin(), x.end(), typename GT::Vector_3(0, 0, 0));
        xc /= n;

        // getting unit average normal of points
        typename GT::Vector_3 uc =
            std::accumulate(u.begin(), u.end(), typename GT::Vector_3(0, 0, 0));
        uc /= sqrt(uc * uc);

        // summing each triangle's measure after triangulation by barycenter split.
        for (std::size_t i = 0; i < n; i++)
        {
            mu0 += interpolated_corrected_area_measure_face<GT>(
                std::vector<typename GT::Vector_3> {u[i], u[(i + 1) % n], uc},
                std::vector<typename GT::Vector_3> {x[i], x[(i + 1) % n], xc}
            );
        }
        return mu0;
    }
}


template<typename GT>
typename GT::FT interpolated_corrected_mean_curvature_measure_face(const std::vector<typename GT::Vector_3>& u,
                                                                   const std::vector<typename GT::Vector_3>& x = {})
{
    const std::size_t n = x.size();
    CGAL_precondition(u.size() == n);
    CGAL_precondition(n >= 3);

    typename GT::Construct_cross_product_vector_3 cross_product;

    // Triangle: use triangle formula
    if (n == 3)
    {
        const typename GT::Vector_3 um = (u[0] + u[1] + u[2]) / 3.0;

        return 0.5 * um * (cross_product(u[2] - u[1], x[0])
            + cross_product(u[0] - u[2], x[1])
            + cross_product(u[1] - u[0], x[2]));
    }
    // Quad: use bilinear interpolation formula
    else if (n == 4)
    {
        // for the formulas below, values of verices 2 & 3 are swapped (compared to paper) to correct order.
        // the indices in paper vs in here are: 00 = 0, 10 = 1, 11 = 2, 01 = 3

        const typename GT::Vector_3 u02 = u[2] - u[0];
        const typename GT::Vector_3 u13 = u[3] - u[1];
        const typename GT::Vector_3 x0_cross = cross_product(u13, x[0]);
        const typename GT::Vector_3 x1_cross = -cross_product(u02, x[1]);
        const typename GT::Vector_3 x3_cross = cross_product(u02, x[3]);
        const typename GT::Vector_3 x2_cross = -cross_product(u13, x[2]);

        return (1.0 / 12.0) * (
              u[0] * (2 * x0_cross - cross_product((u[3] + u[2]), x[1]) + cross_product((u[1] + u[2]), x[3]) + x2_cross)
            + u[1] * (cross_product((u[3] + u[2]), x[0]) + 2 * x1_cross + x3_cross - cross_product((u[0] + u[3]), x[2]))
            + u[3] * (-cross_product((u[1] + u[2]), x[0]) + x1_cross + 2 * x3_cross + cross_product((u[0] + u[1]), x[2]))
            + u[2] * (x0_cross + cross_product((u[0] + u[3]), x[1]) - cross_product((u[0] + u[1]), x[3]) + 2 * x2_cross)
            );
    }
    // N-gon: split into n triangles by polygon center and use triangle formula for each
    else
    {
        typename GT::FT mu1 = 0;

        // getting center of points
        typename GT::Vector_3 xc =
            std::accumulate(x.begin(), x.end(), typename GT::Vector_3(0, 0, 0));
        xc /= n;

        // getting unit average normal of points
        typename GT::Vector_3 uc =
            std::accumulate(u.begin(), u.end(), typename GT::Vector_3(0, 0, 0));
        uc /= sqrt(uc * uc);

        // summing each triangle's measure after triangulation by barycenter split.
        for (std::size_t i = 0; i < n; i++)
        {
            mu1 += interpolated_corrected_mean_curvature_measure_face<GT>(
                std::vector<typename GT::Vector_3> {u[i], u[(i + 1) % n], uc},
                std::vector<typename GT::Vector_3> {x[i], x[(i + 1) % n], xc}
            );
        }
        return mu1;
    }
}


template<typename GT>
typename GT::FT interpolated_corrected_gaussian_curvature_measure_face(const std::vector<typename GT::Vector_3>& u,
                                                                       const std::vector<typename GT::Vector_3>& x = {})
{
    const std::size_t n = u.size();
    CGAL_precondition(n >= 3);

    typename GT::Construct_cross_product_vector_3 cross_product;

    // Triangle: use triangle formula
    if (n == 3)
    {
        return 0.5 * u[0] * cross_product(u[1], u[2]);
    }
    // Quad: use bilinear interpolation formula
    else if (n == 4)
    {
        // for the formulas below, values of verices 2 & 3 are swapped (compared to paper) to correct order.
        // the indices in paper vs in here are: 00 = 0, 10 = 1, 11 = 2, 01 = 3
        return (1.0 / 36.0) * (
              (4 * u[0] + 2 * u[1] + 2 * u[3] + u[2]) * cross_product(u[1] - u[0], u[3] - u[0])
            + (2 * u[0] + 4 * u[1] + u[3] + 2 * u[2]) * cross_product(u[1] - u[0], u[2] - u[1])
            + (2 * u[0] + u[1] + 4 * u[3] + 2 * u[2]) * cross_product(u[2] - u[3], u[3] - u[0])
            + (u[0] + 2 * u[1] + 2 * u[3] + 4 * u[2]) * cross_product(u[2] - u[3], u[2] - u[1])
            );
    }
    // N-gon: split into n triangles by polygon center and use triangle formula for each
    else
    {
        typename GT::FT mu2 = 0;

        // getting unit average normal of points
        typename GT::Vector_3 uc =
            std::accumulate(u.begin(), u.end(), typename GT::Vector_3(0, 0, 0));
        uc /= sqrt(uc * uc);

        // summing each triangle's measure after triangulation by barycenter split.
        for (std::size_t i = 0; i < n; i++)
        {
            mu2 += interpolated_corrected_gaussian_curvature_measure_face<GT>(
                std::vector<typename GT::Vector_3> {u[i], u[(i + 1) % n], uc}
            );
        }
        return mu2;
    }
}


template<typename GT>
std::array<typename GT::FT, 3 * 3> interpolated_corrected_anisotropic_measure_face(const std::vector<typename GT::Vector_3>& u,
                                                                                   const std::vector<typename GT::Vector_3>& x)
{
    const std::size_t n = x.size();
    CGAL_precondition(u.size() == n);
    CGAL_precondition(n >= 3);

    typename GT::Construct_cross_product_vector_3 cross_product;
    std::array<typename GT::FT, 3 * 3> muXY {0};

    // Triangle: use triangle formula
    if (n == 3)
    {
        const typename GT::Vector_3 u01 = u[1] - u[0];
        const typename GT::Vector_3 u02 = u[2] - u[0];
        const typename GT::Vector_3 x01 = x[1] - x[0];
        const typename GT::Vector_3 x02 = x[2] - x[0];
        const typename GT::Vector_3 um = (u[0] + u[1] + u[2]) / 3.0;

        for (std::size_t ix = 0; ix < 3; ix++)
        {
            typename GT::Vector_3 X;
            if (ix == 0)
                X = typename GT::Vector_3(1, 0, 0);
            if (ix == 1)
                X = typename GT::Vector_3(0, 1, 0);
            if (ix == 2)
                X = typename GT::Vector_3(0, 0, 1);

            for (std::size_t iy = 0; iy < 3; iy++)
                muXY[ix * 3 + iy] = 0.5 * um * (cross_product(u02[iy] * X, x01) - cross_product(u01[iy] * X, x02));
        }
    }
    // Quad: use bilinear interpolation formula
    else if (n == 4)
    {
        // for the formulas below, values of verices 2 & 3 are swapped (compared to paper) to correct order.
        // the indices in paper vs in here are: 00 = 0, 10 = 1, 11 = 2, 01 = 3
        for (std::size_t ix = 0; ix < 3; ix++)
        {
            typename GT::Vector_3 X;
            if (ix == 0)
                X = typename GT::Vector_3(1, 0, 0);
            if (ix == 1)
                X = typename GT::Vector_3(0, 1, 0);
            if (ix == 2)
                X = typename GT::Vector_3(0, 0, 1);

            const typename GT::Vector_3 u0xX = cross_product(u[0], X);
            const typename GT::Vector_3 u1xX = cross_product(u[1], X);
            const typename GT::Vector_3 u2xX = cross_product(u[2], X);
            const typename GT::Vector_3 u3xX = cross_product(u[3], X);

            for (std::size_t iy = 0; iy < 3; iy++)
                muXY[ix * 3 + iy] = (1.0 / 72.0) * (

                    u[0][iy] * (    u0xX * (    - x[0] - 11 * x[1] + 13 * x[3] -      x[2])
                                  + u1xX * ( -5 * x[0] -  7 * x[1] + 11 * x[3] +      x[2])
                                  + u3xX * (      x[0] -  7 * x[1] + 11 * x[3] -  5 * x[2])
                                  + u2xX * (    - x[0] -  5 * x[1] +  7 * x[3] -      x[2])
                      )
                  + u[1][iy] * (    u0xX * ( 13 * x[0] -      x[1] -  7 * x[3] -  5 * x[2])
                                  + u1xX * ( 17 * x[0] -  5 * x[1] -  5 * x[3] -  7 * x[2])
                                  + u3xX * (  5 * x[0] +      x[1] +      x[3] -  7 * x[2])
                                  + u2xX * (  7 * x[0] -      x[1] +  5 * x[3] - 11 * x[2])
                      )
                  + u[2][iy] * (    u0xX * (-11 * x[0] +  5 * x[1] -      x[3] +  7 * x[2])
                                  + u1xX * (- 7 * x[0] +      x[1] +      x[3] +  5 * x[2])
                                  + u3xX * (- 7 * x[0] -  5 * x[1] -  5 * x[3] + 17 * x[2])
                                  + u2xX * (- 5 * x[0] -  7 * x[1] -      x[3] + 13 * x[2])
                      )
                  + u[3][iy] * (    u0xX * (-     x[0] +  7 * x[1] -  5 * x[3] -      x[2])
                                  + u1xX * (- 5 * x[0] + 11 * x[1] -  7 * x[3] +      x[2])
                                  + u3xX * (      x[0] + 11 * x[1] -  7 * x[3] -  5 * x[2])
                                  + u2xX * (-     x[0] + 13 * x[1] - 11 * x[3] -      x[2])
                      )

                  );
        }
    }
    // N-gon: split into n triangles by polygon center and use triangle formula for each
    else
    {
        // getting center of points
        typename GT::Vector_3 xc =
            std::accumulate(x.begin(), x.end(), typename GT::Vector_3(0, 0, 0));
        xc /= n;

        // getting unit average normal of points
        typename GT::Vector_3 uc =
            std::accumulate(u.begin(), u.end(), typename GT::Vector_3(0, 0, 0));
        uc /= sqrt(uc * uc);

        // summing each triangle's measure after triangulation by barycenter split.
        for (std::size_t i = 0; i < n; i++)
        {
            std::array<typename GT::FT, 3 * 3> muXY_curr_triangle =
                interpolated_corrected_anisotropic_measure_face<GT>(
                    std::vector<typename GT::Vector_3> {u[i], u[(i + 1) % n], uc},
                    std::vector<typename GT::Vector_3> {x[i], x[(i + 1) % n], xc}
                );

            for (std::size_t ix = 0; ix < 3; ix++)
                for (std::size_t iy = 0; iy < 3; iy++)
                    muXY[ix * 3 + iy] += muXY_curr_triangle[ix * 3 + iy];
        }
    }
    return muXY;
}

//
//template<typename GT>
//typename GT::FT triangle_in_ball_ratio_1(const typename GT::Vector_3 x1,
//                                         const typename GT::Vector_3 x2,
//                                         const typename GT::Vector_3 x3,
//                                         const typename GT::FT r,
//                                         const typename GT::Vector_3 c,
//                                         const std::size_t res = 3)
//{
//    const typename GT::FT R = r * r;
//    const typename GT::FT acc = 1.0 / res;
//    std::size_t samples_in = 0;
//    for (GT::FT alpha = acc / 3; alpha < 1; alpha += acc)
//        for (GT::FT beta = acc / 3; beta < 1 - alpha; beta += acc)
//        {
//            if ((alpha * x1 + beta * x2 + (1 - alpha - beta) * x3 - c).squared_length() < R)
//                samples_in++;
//        }
//    return samples_in / (typename GT::FT)(res * (res + 1) / 2);
//}

template<typename GT>
typename GT::FT face_in_ball_ratio(const std::vector<typename GT::Vector_3>& x,
    const typename GT::FT r,
    const typename GT::Vector_3 c)
{
    const std::size_t n = x.size();

    // getting center of points
    typename GT::Vector_3 xm =
        std::accumulate(x.begin(), x.end(), typename GT::Vector_3(0, 0, 0));
    xm /= n;

    typename GT::FT d_min = (xm - c).squared_length();
    typename GT::FT d_max = d_min;

    for (const typename GT::Vector_3 xi : x)
    {
        const typename GT::FT d_sq = (xi - c).squared_length();
        d_max = (std::max)(d_sq, d_max);
        d_min = (std::min)(d_sq, d_min);
    }

    if (d_max <= r * r) return 1.0;
    else if (r * r <= d_min) return 0.0;

    d_max = sqrt(d_max);
    d_min = sqrt(d_min);

    return (r - d_min) / (d_max - d_min);
}

} // namespace internal

template<typename PolygonMesh, typename NamedParameters = parameters::Default_named_parameters>
class Interpolated_corrected_curvatures_computer
{
    typedef typename GetGeomTraits<PolygonMesh, NamedParameters>::type GT;

    typedef typename GT::FT FT;
    typedef typename GT::Point_3 Point_3;
    typedef typename GT::Vector_3 Vector_3;

    typedef typename GetVertexPointMap<PolygonMesh, NamedParameters>::const_type VPM;

    typedef dynamic_vertex_property_t<Vector_3> Vector_map_tag;
    typedef typename boost::property_map<PolygonMesh, Vector_map_tag>::const_type Default_vector_map;
    typedef typename internal_np::Lookup_named_param_def<internal_np::vertex_normal_map_t,
        NamedParameters,
        Default_vector_map>::type VNM;

    typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor halfedge_descriptor;
    typedef typename boost::graph_traits<PolygonMesh>::edge_descriptor edge_descriptor;
    typedef typename boost::graph_traits<PolygonMesh>::face_descriptor face_descriptor;
    typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor vertex_descriptor;

    typedef std::unordered_map<face_descriptor, FT> Face_scalar_measure_map_tag;
    typedef std::unordered_map<vertex_descriptor, FT> Vertex_scalar_measure_map_tag;
    typedef std::unordered_map<face_descriptor, std::array<FT, 3 * 3>> Face_anisotropic_measure_map_tag;
    typedef std::unordered_map<vertex_descriptor, Principal_curvature<GT>> Vertex_principal_curvature_map_tag;

private:
    const PolygonMesh& pmesh;
    VPM vpm;
    VNM vnm;
    FT ball_radius;

    bool mean_curvature_selected = false;
    bool gaussian_curvature_selected = false;
    bool principal_curvature_selected = false;

    bool area_measure_computed = false;
    bool mean_curvature_measure_computed = false;
    bool gaussian_curvature_measure_computed = false;
    bool anisotropic_measure_computed = false;

    bool area_measure_expanded = false;
    bool mean_curvature_computed = false;          // implicitly also means mean measure is expanded
    bool gaussian_curvature_computed = false;      // implicitly also means gaussiaan measure is expanded
    bool principal_curvature_computed = false;     // implicitly also means anisotropic measure is expanded

    Face_scalar_measure_map_tag mu0_map_init, mu1_map_init, mu2_map_init;
    Face_anisotropic_measure_map_tag muXY_map_init;
    Vertex_scalar_measure_map_tag mean_curvature_map_init, gaussian_curvature_map_init;
    Vertex_principal_curvature_map_tag principal_curvature_map_init;

    boost::associative_property_map<Face_scalar_measure_map_tag> mu0_map, mu1_map, mu2_map;
    boost::associative_property_map<Face_anisotropic_measure_map_tag> muXY_map;

public:

    boost::associative_property_map<Vertex_scalar_measure_map_tag> mean_curvature_map, gaussian_curvature_map;
    boost::associative_property_map<Vertex_principal_curvature_map_tag> principal_curvature_map;

    void set_named_params(const NamedParameters& np = parameters::default_values())
    {
        using parameters::choose_parameter;
        using parameters::get_parameter;
        using parameters::is_default_parameter;

        vpm = choose_parameter(get_parameter(np, CGAL::vertex_point),
            get_const_property_map(CGAL::vertex_point, pmesh));

        vnm = choose_parameter(get_parameter(np, internal_np::vertex_normal_map),
            get(Vector_map_tag(), pmesh));

        if (is_default_parameter<NamedParameters, internal_np::vertex_normal_map_t>::value)
            compute_vertex_normals(pmesh, vnm, np);

        ball_radius = choose_parameter(get_parameter(np, internal_np::ball_radius), 0);

        set_ball_radius(ball_radius);
    }

    void set_property_maps() {
        mu0_map = boost::associative_property_map<Face_scalar_measure_map_tag>(mu0_map_init);
        mu1_map = boost::associative_property_map<Face_scalar_measure_map_tag>(mu1_map_init);
        mu2_map = boost::associative_property_map<Face_scalar_measure_map_tag>(mu2_map_init);
        muXY_map = boost::associative_property_map<Face_anisotropic_measure_map_tag>(muXY_map_init);
        mean_curvature_map = boost::associative_property_map<Vertex_scalar_measure_map_tag>(mean_curvature_map_init);
        gaussian_curvature_map = boost::associative_property_map<Vertex_scalar_measure_map_tag>(gaussian_curvature_map_init);
        principal_curvature_map = boost::associative_property_map<Vertex_principal_curvature_map_tag>(principal_curvature_map_init);
    }

    Interpolated_corrected_curvatures_computer(const PolygonMesh& pmesh,
        const NamedParameters& np = parameters::default_values())
        :
        pmesh(pmesh)
    {
        set_named_params(np);
        set_property_maps();
    }

    Interpolated_corrected_curvatures_computer(const PolygonMesh& pmesh,
        bool mean_curvature_selected,
        bool gaussian_curvature_selected,
        bool principal_curvature_selected,
        const NamedParameters& np = parameters::default_values()
    ):
        pmesh(pmesh)
    {
        set_named_params(np);
        set_property_maps();
        this->mean_curvature_selected = mean_curvature_selected;
        this->gaussian_curvature_selected = gaussian_curvature_selected;
        this->principal_curvature_selected = principal_curvature_selected;
    }


    void set_ball_radius(const FT radius) {
        if (radius == 0)
            ball_radius = internal::average_edge_length<PolygonMesh, GT>(pmesh) * EXPANDING_RADIUS_EPSILON;

        area_measure_expanded = false;
        mean_curvature_computed = false;
        gaussian_curvature_computed = false;
        principal_curvature_computed = false;
    }


    /**
    * \ingroup PMP_corrected_curvatures_grp
    *
    * Computes the interpolated corrected curvature measure on each face of the mesh.
    *
    * @tparam PolygonMesh a model of `FaceListGraph`.
    * @tparam FaceMeasureMap a model of `WritablePropertyMap` with
    * `boost::graph_traits<PolygonMesh>::%face_descriptor` as key type and `GT::FT` as value type.
    * @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters".
    *
    * @param pmesh the polygon mesh.
    * @param fmm (face measure map) the property map used for storing the computed interpolated corrected measure.
    * @param mu_i an enum for choosing between computing
    *             the area measure, the mean curvature measure or the gaussian curvature measure.
    * @param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
    *
    * \cgalNamedParamsBegin
    *
    *   \cgalParamNBegin{vertex_point_map}
    *     \cgalParamDescription{a property map associating points to the vertices of `pmesh`}
    *     \cgalParamType{a class model of `ReadablePropertyMap` with
    *                    `boost::graph_traits<PolygonMesh>::%vertex_descriptor`
    *                    as key type and `%Point_3` as value type}
    *     \cgalParamDefault{`boost::get(CGAL::vertex_point, pmesh)`}
    *     \cgalParamExtra{If this parameter is omitted, an internal property map for
    *                     `CGAL::vertex_point_t` must be available in `PolygonMesh`.}
    *   \cgalParamNEnd
    *
    *   \cgalParamNBegin{vertex_normal_map}
    *     \cgalParamDescription{a property map associating normal vectors to the vertices of `pmesh`}
    *     \cgalParamType{a class model of `ReadablePropertyMap` with
    *                    `boost::graph_traits<PolygonMesh>::%vertex_descriptor`
    *                    as key type and `%Vector_3` as value type}
    *     \cgalParamDefault{`get(dynamic_vertex_property_t<GT::Vector_3>(), pmesh)`}
    *     \cgalParamExtra{If this parameter is omitted, vertex normals will be
    *                     computed using compute_vertex_normals()}
    *   \cgalParamNEnd
    *
    * \cgalNamedParamsEnd
    *
    * @see `interpolated_corrected_area_measure_face()`
    * @see `interpolated_corrected_mean_curvature_measure_face()`
    * @see `interpolated_corrected_gaussian_curvature_measure_face()`
    */
    void interpolated_corrected_measure_all_faces(Face_scalar_measure_map_tag fmm, const internal::Curvature_measure_index mu_i)
    {
        std::function
            <FT(const std::vector<Vector_3>&, const std::vector<Vector_3>&)>
            iccm_function;
        switch (mu_i)
        {
        case internal::MU0_AREA_MEASURE:
            iccm_function = &internal::interpolated_corrected_area_measure_face<GT>;
            break;
        case internal::MU1_MEAN_CURVATURE_MEASURE:
            iccm_function = &internal::interpolated_corrected_mean_curvature_measure_face<GT>;
            break;
        case internal::MU2_GAUSSIAN_CURVATURE_MEASURE:
            iccm_function = &internal::interpolated_corrected_gaussian_curvature_measure_face<GT>;
            break;
        }

        std::vector<Vector_3> x;
        std::vector<Vector_3> u;

        for (face_descriptor f : faces(pmesh))
        {

            for (vertex_descriptor v : vertices_around_face(halfedge(f, pmesh), pmesh))
            {
                Point_3 p = get(vpm, v);
                x.push_back(Vector_3(p.x(), p.y(), p.z()));
                u.push_back(get(vnm, v));
            }

            put(fmm, f, iccm_function(u, x));
            x.clear();
            u.clear();
        }

        switch (mu_i)
        {
        case internal::MU0_AREA_MEASURE:
            area_measure_computed = true;
            break;
        case internal::MU1_MEAN_CURVATURE_MEASURE:
            mean_curvature_measure_computed = true;
            break;
        case internal::MU2_GAUSSIAN_CURVATURE_MEASURE:
            gaussian_curvature_measure_computed = true;
            break;
        }
    }

    void interpolated_corrected_all_measures_all_faces()
    {
        bool mean_curvature_measure_needed = mean_curvature_selected && !mean_curvature_measure_computed;
        bool gaussian_curvature_measure_needed = gaussian_curvature_selected && !gaussian_curvature_measure_computed;
        bool anisotropic_measure_needed = principal_curvature_selected && !anisotropic_measure_computed;

        // if no computations are needed
        if (!(
            !area_measure_computed ||
            mean_curvature_measure_needed ||
            gaussian_curvature_measure_needed ||
            anisotropic_measure_needed
            )) return;

        std::vector<Vector_3> x;
        std::vector<Vector_3> u;

        for (face_descriptor f : faces(pmesh))
        {
            for (vertex_descriptor v : vertices_around_face(halfedge(f, pmesh), pmesh))
            {
                Point_3 p = get(vpm, v);
                x.push_back(Vector_3(p.x(), p.y(), p.z()));
                u.push_back(get(vnm, v));
            }
            if (!area_measure_computed)
                put(mu0_map, f, internal::interpolated_corrected_area_measure_face<GT>(u, x));

            if (mean_curvature_measure_needed)
                put(mu1_map, f, internal::interpolated_corrected_mean_curvature_measure_face<GT>(u, x));

            if (gaussian_curvature_measure_needed)
                put(mu2_map, f, internal::interpolated_corrected_gaussian_curvature_measure_face<GT>(u, x));

            if (anisotropic_measure_needed)
                put(muXY_map, f, internal::interpolated_corrected_anisotropic_measure_face<GT>(u, x));

            x.clear();
            u.clear();
        }
        if (!area_measure_computed)
            area_measure_computed = true;

        if (mean_curvature_measure_needed)
            mean_curvature_measure_computed = true;

        if (gaussian_curvature_measure_needed)
            gaussian_curvature_measure_computed = true;

        if (anisotropic_measure_needed)
            anisotropic_measure_computed = true;

        for (face_descriptor f : faces(pmesh))
        {
            std::cout << f << ": " << mu0_map[f] << ", " << mu1_map[f] << ", " << mu2_map[f] << "\n";
        }
    }


    /**
    * \ingroup PMP_corrected_curvatures_grp
    *
    * Computes the interpolated corrected anisotropic measure on each face of the mesh
    *
    * @tparam PolygonMesh a model of `FaceListGraph`.
    * @tparam FaceMeasureMap a model of `WritablePropertyMap` with
    * `boost::graph_traits<PolygonMesh>::%face_descriptor` as key type and `std::array<GT::FT, 3 * 3>` as value type.
    * @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters".
    *
    * @param pmesh the polygon mesh.
    * @param fmm (face measure map) the property map used for storing the computed interpolated corrected measure.
    * @param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
    *
    * \cgalNamedParamsBegin
    *
    *   \cgalParamNBegin{vertex_point_map}
    *     \cgalParamDescription{a property map associating points to the vertices of `pmesh`}
    *     \cgalParamType{a class model of `ReadablePropertyMap` with
    *                    `boost::graph_traits<PolygonMesh>::%vertex_descriptor`
    *                    as key type and `%Point_3` as value type}
    *     \cgalParamDefault{`boost::get(CGAL::vertex_point, pmesh)`}
    *     \cgalParamExtra{If this parameter is omitted, an internal property map for
    *                     `CGAL::vertex_point_t` must be available in `PolygonMesh`.}
    *   \cgalParamNEnd
    *
    *   \cgalParamNBegin{vertex_normal_map}
    *     \cgalParamDescription{a property map associating normal vectors to the vertices of `pmesh`}
    *     \cgalParamType{a class model of `ReadablePropertyMap` with
    *                    `boost::graph_traits<PolygonMesh>::%vertex_descriptor`
    *                    as key type and `%Vector_3` as value type}
    *     \cgalParamDefault{`get(dynamic_vertex_property_t<GT::Vector_3>(), pmesh)`}
    *     \cgalParamExtra{If this parameter is omitted, vertex normals will be
    *                     computed using compute_vertex_normals()}
    *   \cgalParamNEnd
    *
    * \cgalNamedParamsEnd
    *
    * @see `interpolated_corrected_anisotropic_measure_face()`
    * @see `interpolated_corrected_measure_mesh()`
    */
        void
        interpolated_corrected_anisotropic_measure_all_faces()
    {
        std::vector<Vector_3> x;
        std::vector<Vector_3> u;

        for (face_descriptor f : faces(pmesh))
        {

            for (vertex_descriptor v : vertices_around_face(halfedge(f, pmesh), pmesh))
            {
                Point_3 p = get(vpm, v);
                x.push_back(Vector_3(p.x(), p.y(), p.z()));
                u.push_back(get(vnm, v));
            }

            put(muXY_map, f, internal::interpolated_corrected_anisotropic_measure_face<GT>(u, x));
            x.clear();
            u.clear();
        }
        anisotropic_measure_computed = true;
    }

    /**
    * \ingroup PMP_corrected_curvatures_grp
    *
    * Expands given face area and curvature (mean or gaussian) measures around a vertex `v`.
    * Expansion is based on the inclusion ratio of each face in a ball of radius `r` around the vertex `v`.
    *
    * @tparam PolygonMesh a model of `FaceListGraph`.
    * @tparam FaceMeasureMap a model of `WritablePropertyMap` with
    * `boost::graph_traits<PolygonMesh>::%face_descriptor` as key type and `GT::FT` as value type.
    * @tparam VertexMeasureMap a model of `WritablePropertyMap` with
    * `boost::graph_traits<PolygonMesh>::%vertex_descriptor` as key type and `GT::FT` as value type.
    * @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters".
    *
    * @param pmesh the polygon mesh.
    * @param area_fmm (area face measure map) the property map storing the already computed area measure on each face.
    * @param curvature_fmm (curvature face measure map) the property map storing the already computed curvature measure on each face.
    * This curvature measure can be either the Mean Curvature or the Gaussian Curvature.
    * @param area_vmm (area vertex measure map) the property map provided to store the expanded area measure on each vertex.
    * @param curvature_vmm (curvature vertex measure map) the property map provided to store the expanded curvature measure on each vertex.
    * This curvature measure can be either the Mean Curvature or the Gaussian Curvature.
    * @param v (vertex) the vertex to expand the area and curvature measure around.
    * @param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
    *
    * \cgalNamedParamsBegin
    *
    *   \cgalParamNBegin{ball_radius}
    *     \cgalParamDescription{the radius of the ball around the vertex `v` to expand the area and curvature measure}
    *     \cgalParamType{`GT::FT`}
    *     \cgalParamDefault{`0.01`}
    *   \cgalParamNEnd
    *
    *   \cgalParamNBegin{vertex_point_map}
    *     \cgalParamDescription{a property map associating points to the vertices of `pmesh`}
    *     \cgalParamType{a class model of `ReadablePropertyMap` with
    *                    `boost::graph_traits<PolygonMesh>::%vertex_descriptor`
    *                    as key type and `%Point_3` as value type}
    *     \cgalParamDefault{`boost::get(CGAL::vertex_point, pmesh)`}
    *     \cgalParamExtra{If this parameter is omitted, an internal property map for
    *                     `CGAL::vertex_point_t` must be available in `PolygonMesh`.}
    *   \cgalParamNEnd
    *
    * \cgalNamedParamsEnd
    *
    * @see `expand_interpolated_corrected_anisotropic_measure_vertex()`
    * @see `face_in_ball_ratio()`
    */
    internal::Vertex_curvatures<GT> expand_interpolated_corrected_measure_vertex(vertex_descriptor v)
    {
        bool mean_curvature_measure_expansion_needed = mean_curvature_selected && !mean_curvature_computed;
        bool gaussian_curvature_measure_expansion_needed = gaussian_curvature_selected && !gaussian_curvature_computed;
        bool anisotropic_measure_expansion_needed = principal_curvature_selected && !principal_curvature_computed;

        std::queue<face_descriptor> bfs_queue;
        std::unordered_set<face_descriptor> bfs_visited;

        Point_3 vp = get(vpm, v);
        Vector_3 c = Vector_3(vp.x(), vp.y(), vp.z());

        internal::Vertex_curvatures<GT> vertex_curvatures;

        for (face_descriptor f : faces_around_target(halfedge(v, pmesh), pmesh)) {
            if (f != boost::graph_traits<PolygonMesh>::null_face())
            {
                bfs_queue.push(f);
                bfs_visited.insert(f);
            }
        }
        while (!bfs_queue.empty()) {
            face_descriptor fi = bfs_queue.front();
            bfs_queue.pop();

            // looping over vertices in face to get point coordinates
            std::vector<Vector_3> x;
            for (vertex_descriptor vi : vertices_around_face(halfedge(fi, pmesh), pmesh))
            {
                Point_3 pi = get(vpm, vi);
                x.push_back(Vector_3(pi.x(), pi.y(), pi.z()));
            }

            const FT f_ratio = internal::face_in_ball_ratio<GT>(x, ball_radius, c);

            if (f_ratio != 0.0)
            {
                if (!area_measure_expanded)
                    vertex_curvatures.area_measure += f_ratio * get(mu0_map, fi);

                if (mean_curvature_measure_expansion_needed)
                    vertex_curvatures.mean_curvature_measure += f_ratio * get(mu1_map, fi);

                if (gaussian_curvature_measure_expansion_needed)
                    vertex_curvatures.gaussian_curvature_measure += f_ratio * get(mu2_map, fi);

                if (anisotropic_measure_expansion_needed)
                {
                    const std::array<FT, 3 * 3> face_anisotropic_measure = get(muXY_map, fi);
                    for (std::size_t i = 0; i < 3 * 3; i++)
                        vertex_curvatures.anisotropic_measure[i] += f_ratio * face_anisotropic_measure[i];
                }

                for (face_descriptor fj : faces_around_face(halfedge(fi, pmesh), pmesh))
                {
                    if (bfs_visited.find(fj) == bfs_visited.end() && fj != boost::graph_traits<PolygonMesh>::null_face())
                    {
                        bfs_queue.push(fj);
                        bfs_visited.insert(fj);
                    }
                }
            }
        }

        std::cout << vertex_curvatures.area_measure << " " << vertex_curvatures.mean_curvature_measure << " " << vertex_curvatures.gaussian_curvature_measure << "\n";
        return vertex_curvatures;
    }


    void compute_selected_curvatures() {
        bool mean_curvature_needed = mean_curvature_selected && !mean_curvature_computed;
        bool gaussian_curvature_needed = gaussian_curvature_selected && !gaussian_curvature_computed;
        bool principal_curvature_needed = principal_curvature_selected && !principal_curvature_computed;

        if (!(mean_curvature_needed || gaussian_curvature_needed || principal_curvature_needed)) return;

        interpolated_corrected_all_measures_all_faces();

        for (vertex_descriptor v : vertices(pmesh))
        {
            const internal::Vertex_curvatures<GT> vertex_curvatures = expand_interpolated_corrected_measure_vertex(v);

            if (mean_curvature_needed) {
                vertex_curvatures.area_measure != 0 ?
                    put(mean_curvature_map, v, 0.5 * vertex_curvatures.mean_curvature_measure / vertex_curvatures.area_measure) :
                    put(mean_curvature_map, v, 0);
            }

            if (gaussian_curvature_needed) {
                vertex_curvatures.area_measure != 0 ?
                    put(gaussian_curvature_map, v, vertex_curvatures.gaussian_curvature_measure / vertex_curvatures.area_measure) :
                    put(gaussian_curvature_map, v, 0);
            }

            if (principal_curvature_needed) {
                const Vector_3  v_normal = get(vnm, v);
                const Principal_curvature<GT> principal_curvature = internal::principal_curvature_from_anisotropic_measures<GT>(
                        vertex_curvatures.anisotropic_measure,
                        vertex_curvatures.area_measure,
                        v_normal
                    );
                put(principal_curvature_map, v, principal_curvature);
            }
        }
        if (!area_measure_expanded)
            area_measure_expanded = true;

        if (mean_curvature_needed)
            mean_curvature_computed = true;

        if (gaussian_curvature_needed)
            gaussian_curvature_computed = true;

        if (principal_curvature_needed)
            principal_curvature_computed = true;
    }


    void compute_all_curvatures() {
        mean_curvature_selected = true;
        gaussian_curvature_selected = true;
        principal_curvature_selected = true;
        compute_selected_curvatures();
    }

    void compute_mean_curvature()
    {

    }


};

} // namespace Polygon_mesh_processing
} // namespace CGAL

#endif // CGAL_POLYGON_MESH_PROCESSING_INTERPOLATED_CORRECTED_CURVATURE_MEASURES_H