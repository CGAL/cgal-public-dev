// Copyright (c) 2023 GeometryFactory
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Robert Piel

#ifndef CGAL_SURFACE_MESH_APPROXIMATE_SHORTEST_PATH_SURFACE_MESH_APPROXIMATE_SHORTEST_PATH_H
#define CGAL_SURFACE_MESH_APPROXIMATE_SHORTEST_PATH_SURFACE_MESH_APPROXIMATE_SHORTEST_PATH_H

#include <CGAL/boost/graph/helpers.h>
#include <CGAL/boost/graph/iterator.h>

#include "CGAL/Surface_mesh_approximate_shortest_path/Enqueue_policies.h"
#include "CGAL/Surface_mesh_approximate_shortest_path/Skip_conditions.h"

#include <CGAL/Polygon_mesh_processing/locate.h>

namespace CGAL {

/*!
 * \brief The Face_values class is a struct that contains all the data we need to store per face.
 * \details It consists of the 'Kernel::FT' objects 'sigma' and 'd', respectively representing the geodesic distance from the source to
 * the current virtual source and the overall current distance from the source to the face target point (the face's barycenter unless otherwise specified, \sa 'add_target()').
 * The struct further contains a 'std::array<FT, 3>' storing the distance of the face's vertices to the current virtual source.
 */
template<class Kernel>
struct Face_values {
    typedef typename Kernel::FT FT;

    FT sigma;
    FT d;
    std::array<FT,3> d2verts;

    Face_values(FT _sigma=0., FT _d=std::numeric_limits<FT>::max(), std::array<FT,3> _d2verts = {-1., -1., -1.})
        : sigma(_sigma), d(_d), d2verts(_d2verts) {}

    friend std::ostream & operator <<(std::ostream& stream, const Face_values vals)
    {
        return ( stream << vals.sigma << "\t" << vals.d << "\t"
                       << vals.d2verts[0]  << "\t" << vals.d2verts[1] << "\t" << vals.d2verts[2] << std::endl);
    };
};

/*!
 * \brief The Unfolded_triangle_2 class is a struct containing two 'Point_2' objects, i.e. the tangent space points B and P as in the paper.
 */
template<class Kernel>
struct Unfolded_triangle_2 {
    typedef typename Kernel::Point_2 Point_2;

    //Point_2 A; // note that the first point is always (0,0)
    Point_2 B;
    Point_2 P;
};

class Surface_mesh_approximate_shortest_path_observer
{
public:
    int m_num_faces;
    int m_iter_counter;
    std::vector< std::vector<int> > m_update_counter;

public:
    Surface_mesh_approximate_shortest_path_observer() {};

    Surface_mesh_approximate_shortest_path_observer(int num_faces)
        : m_num_faces(num_faces)
    {
        m_iter_counter = 0;
    };

    void add_iteration_step()
    {
        m_iter_counter++;
        std::vector<int> iter_vector(m_num_faces);
        std::fill(iter_vector.begin(), iter_vector.end(), 0);
        m_update_counter.push_back(iter_vector);
    }

    void increment(int index)
    {
        m_update_counter[m_iter_counter-1][index]++;
    }

    std::vector<int>& get_update_counts_per_iteration(int iteration)
    {
        return m_update_counter[iteration];
    }

    std::vector<int>& get_cumulative_update_counts_after_iterations(int iterations)
    {
        std::vector<int> cumulative_update_counts = m_update_counter[0];
        for (int i = 1; i < iterations; i++)
        {
            std::transform(cumulative_update_counts.begin(), cumulative_update_counts.end(),
                           m_update_counter[i].begin(),
                           cumulative_update_counts.begin(),
                           std::plus<double>());
        }
        return cumulative_update_counts;
    }

    int number_of_total_face_updates()
    {
        std::vector<int> cumulative_update_counts = get_cumulative_update_counts_after_iterations(m_iter_counter-1); // try that without the -1, too

        int total_count = std::accumulate(cumulative_update_counts.begin(), cumulative_update_counts.end(), 0);

        return total_count;
    }
};

enum SourcesStatus
{
    NO_SOURCES_SET = 0,
    SOURCES_SET = 1,
    SOURCES_INITIALISED = 2
};

enum TargetsStatus
{
    NO_TARGETS_SET = 0,
    TARGETS_SET = 1,
};

enum PropagationStatus
{
    VALID_RESULTS = 1,
    INVALID_RESULTS = 0
};

/*!
 * \ingroup PkgSurfaceMeshApproximateShortestPathRef
 *
 * \brief Compute approximate shortest path distances between one or multiple source points and no, one or multiple target points on a surface mesh.
 *
 * \details This package provides an implementation of the paper "Geodesic Distance Computation via Virtual Source Propagation" by Trettner, Bommes and Kobbelt (2021).
 * It offers a simplified version of the algorithm behind the Surface_mesh_shortest_path package. Instead of geodesic windows, we only propagte virtual sources.
 * This leads to inexact reconstruction of the previous geodesic information and hence only gives us an approximate geodesic distance. In return, this simpler
 * algorithm is faster.
 *
 * \tparam Traits a model of 'SurfaceMeshApproximateShortestPathTraits'.
 * \tparam Skip_condition a model of SkipCondition, i.e. a functor that needs to return a 'SkipResult'.
 * \tparam Enqueue_policy a model of EnqueuePolicy, i.e. a functor that needs to return a 'EnqueueResult'.
 *
 * If no specific Skip_condition or Enqueue_policy are provided, the template arguments default to
 * Surface_mesh_approximate_shortest_path_3::Never_skip_condition and
 * Surface_mesh_approximate_shortest_path_3::Static_speed_limiter.
 */
template<class Traits,
         class Visibility_heuristic = typename Traits::Visibility_heuristic,
         class Skip_condition = Surface_mesh_approximate_shortest_path_3::Never_skip_condition,
         class Enqueue_policy = Surface_mesh_approximate_shortest_path_3::Static_speed_limiter<class Traits::Kernel>>
class Surface_mesh_approximate_shortest_path
{
public:
    typedef typename Traits::Kernel Kernel;
    typedef typename Kernel::FT FT;

    typedef typename Kernel::Point_2 Point_2;
    typedef typename Kernel::Point_3 Point_3;

    typedef typename Traits::Surface_mesh Surface_mesh;
    typedef boost::graph_traits<Surface_mesh> Graph_traits;

    typedef typename Graph_traits::vertex_descriptor vertex_descriptor;
    typedef typename Graph_traits::edge_descriptor edge_descriptor;
    typedef typename Graph_traits::halfedge_descriptor halfedge_descriptor;
    typedef typename Graph_traits::face_descriptor face_descriptor;

    typedef typename Surface_mesh::template Property_map<edge_descriptor, FT> Edge_property_map;

    typedef Unfolded_triangle_2<Kernel>     Unfolded_triangle;
    typedef Face_values<Kernel>             Face_values;

    typedef Polygon_mesh_processing::Barycentric_coordinates<FT> Barycentric_coordinates;
    typedef Polygon_mesh_processing::Face_location<Surface_mesh, FT> Face_location;

    typedef typename Surface_mesh::template Property_map<face_descriptor, Face_values> Face_values_map;
    typedef typename Surface_mesh::template Property_map<face_descriptor, face_descriptor> Shortest_path_tree;
    typedef Surface_mesh_approximate_shortest_path_observer Observer;

    /*!
     * \brief The IntersectionResult enum comprises all possible results from an intersection test of the
     * edge of the unfolded triangle and the line from the geodesic source to the heuristic point of the triangle.
     */
    enum IntersectionResult
    {
        EDGE_INTERSECTED = 0,
        LEFT_OF_EDGE = 1,
        RIGHT_OF_EDGE = 2
    };

    /*!
     * \brief The UpdateResult enum represents whether or not we have found a
     * shorter path and therefore whether or not we need to update our 'Face_values_map'.
     */
    enum UpdateResult
    {
        DO_NOT_UPDATE = 0,
        UPDATE = 1
    };

private:
    Surface_mesh& m_mesh;

    Edge_property_map m_edge_lengths;
    Face_values_map m_face_values;

    SourcesStatus m_sources_status;
    std::vector<Face_location> m_source_face_locations;
    std::vector<Point_3> m_source_physical_locations;

    TargetsStatus m_targets_status;
    std::vector<Face_location> m_target_face_locations;
    std::vector<Point_3> m_target_physical_locations;

    PropagationStatus m_propagation_status;
    typedef std::queue<halfedge_descriptor> Halfedge_queue;
    Halfedge_queue m_A, m_B;

    Visibility_heuristic m_Visibility_heuristic;
    Skip_condition m_Skip_condition;
    Enqueue_policy m_Enqueue_policy;

    bool m_build_shortest_path_tree;
    Shortest_path_tree m_shortest_path_tree;

    bool m_with_observer;

    Observer m_observer;

public:
    /*!
     * \brief Surface_mesh_approximate_shortest_path default constructor. The 'Visibility_condition' and the 'Skip_condition' and 'Enqueue_policy' are default constructed.
     */
    Surface_mesh_approximate_shortest_path(Surface_mesh& mesh,
                                           bool build_shortest_path_tree = false,
                                           bool with_observer = false)
        : m_mesh(mesh),
          m_A(), m_B(),
          m_sources_status(CGAL::NO_SOURCES_SET), m_targets_status(CGAL::NO_TARGETS_SET), m_propagation_status(CGAL::INVALID_RESULTS),
          m_Visibility_heuristic(), m_Skip_condition(), m_Enqueue_policy(),
        m_build_shortest_path_tree(build_shortest_path_tree), m_with_observer(with_observer)
    {
        init_property_maps(build_shortest_path_tree);
        if (m_with_observer) { init_observer(); }
    };

    /*!
     * \brief Surface_mesh_approximate_shortest_path constructor that allows for custom contruction of the 'Skip_condition' and 'Enqueue_policy'.
     * The 'Visibility_heuristic' is default constructed.
     */
    Surface_mesh_approximate_shortest_path(Surface_mesh& mesh,
                                           Skip_condition SkipCondition,
                                           Enqueue_policy EnqueuePolicy,
                                           bool build_shortest_path_tree = false,
                                           bool with_observer = false)
        : m_mesh(mesh),
          m_A(), m_B(),
          m_sources_status(CGAL::NO_SOURCES_SET), m_targets_status(CGAL::NO_TARGETS_SET), m_propagation_status(CGAL::INVALID_RESULTS),
          m_Visibility_heuristic(), m_Skip_condition(SkipCondition), m_Enqueue_policy(EnqueuePolicy),
          m_build_shortest_path_tree(build_shortest_path_tree), m_with_observer(with_observer)
    {
        init_property_maps(build_shortest_path_tree);
        if (m_with_observer) { init_observer(); }
    };

    /*!
     * \brief Surface_mesh_approximate_shortest_path constructor that allows for custom contruction of the 'Visibility_heuristic' and the 'Skip_condition' and 'Enqueue_policy'.
     */
    Surface_mesh_approximate_shortest_path(Surface_mesh& mesh,
                                           Visibility_heuristic VisibilityHeuristic,
                                           Skip_condition SkipCondition,
                                           Enqueue_policy EnqueuePolicy,
                                           bool build_shortest_path_tree = false,
                                           bool with_observer = false)
        : m_mesh(mesh),
          m_A(), m_B(),
          m_sources_status(CGAL::NO_SOURCES_SET), m_targets_status(CGAL::NO_TARGETS_SET), m_propagation_status(CGAL::INVALID_RESULTS),
          m_Visibility_heuristic(VisibilityHeuristic), m_Skip_condition(SkipCondition), m_Enqueue_policy(EnqueuePolicy),
          m_build_shortest_path_tree(build_shortest_path_tree), m_with_observer(with_observer)
    {
        init_property_maps(build_shortest_path_tree);
        if (m_with_observer) { init_observer(); }
    };

    void init_property_maps(bool build_shortest_path_tree)
    {
        init_propagation_property_maps();
        if (build_shortest_path_tree) { init_shortest_path_tree(); }
    }

    void init_propagation_property_maps()
    {
        bool created_edge_property_map, created_face_property_map;
        boost::tie(m_edge_lengths, created_edge_property_map) = m_mesh.template add_property_map<edge_descriptor, FT>("edge_lengths");
        assert(created_edge_property_map);

        boost::tie(m_face_values, created_face_property_map) = m_mesh.template add_property_map<face_descriptor, Face_values>("face_values");
        assert(created_face_property_map);
    }

    void init_shortest_path_tree()
    {
        bool created_shortest_path_tree;
        boost::tie(m_shortest_path_tree, created_shortest_path_tree) = m_mesh.template add_property_map<face_descriptor, face_descriptor>("geodesic_path_tree");
        assert(created_shortest_path_tree);
    }

    void init_observer()
    {
        m_observer = Observer(m_mesh.number_of_faces());
    }

    Edge_property_map& Get_edge_length_map() { return m_edge_lengths; };

    /// \name Constructions
    /// @{

    class Compute_squared_edge_length {
    public:
        Compute_squared_edge_length() {};

        FT operator() (Surface_mesh& tmesh, halfedge_descriptor h)
        {
            vertex_descriptor v1 = tmesh.source(h);
            vertex_descriptor v2 = tmesh.target(h);

            return operator() (tmesh, v1, v2);
        }

        FT operator() (Surface_mesh& tmesh, vertex_descriptor v1, vertex_descriptor v2)
        {
            Point_3 v1_point = tmesh.point(v1);
            Point_3 v2_point = tmesh.point(v2);

            double length = squared_distance(v1_point, v2_point);
            return length;
        }
    };

    /*!
     * \brief This function looks up the squared length of an edge in the internal 'edge_property_map' and only computes it if is has not been computed before.
     */
    class Find_edge_length_and_update_property_map
    {
    public:
        typedef FT result_type;

    public:
        Find_edge_length_and_update_property_map() {};

        result_type operator() (Surface_mesh& mesh, halfedge_descriptor h)
        {
            return Compute_squared_edge_length()(mesh, h);
        }

        result_type operator() (Surface_mesh& mesh, halfedge_descriptor h, Edge_property_map& edge_lengths)
        {
            edge_descriptor e = mesh.edge(h);
            FT edge_length = edge_lengths[e];

            if (edge_length == 0)
            {
                edge_length = Compute_squared_edge_length()(mesh, h);
                edge_lengths[e] = edge_length;
            }

            return edge_length;
        }
    };
    typedef Find_edge_length_and_update_property_map Find_edge_length;

    /*!
     * \brief This functor unfolds a triangular face in the tangent space of the face 'face(opposite(h, mesh), mesh)'.
     */
    class Unfold_triangle_3_along_halfedge
    {
    public:
        typedef Unfolded_triangle result_type;

        public:
            Unfold_triangle_3_along_halfedge() {}

            result_type operator() (Surface_mesh& tmesh, halfedge_descriptor h)
            {
                Unfolded_triangle unfolded_triangle;

                // edge length
                FT e0 = Compute_squared_edge_length()(tmesh, h);

                // second point
                unfolded_triangle.B = Point_2(sqrt(e0), 0.);

                // second edge length
                halfedge_descriptor nexth = tmesh.next(h);
                FT e1 = Compute_squared_edge_length()(tmesh, nexth);

                // third edge length
                FT e2 = Compute_squared_edge_length()(
                    tmesh, tmesh.target(nexth), tmesh.source(h));

                // third point
                FT Px = (e0 + (e2 - e1)) / (2*sqrt(e0));
                FT Py = sqrt(e2 - square(Px));
                unfolded_triangle.P = Point_2(Px, Py);

                return unfolded_triangle;
            }

            result_type operator() (Surface_mesh& tmesh, halfedge_descriptor h, Edge_property_map& edge_lengths)
            {
                Unfolded_triangle unfolded_triangle;

                // edge length
                FT e0 = Find_edge_length_and_update_property_map()(tmesh, h, edge_lengths);

                // second point
                unfolded_triangle.B = Point_2(sqrt(e0), 0.);

                // second edge length
                halfedge_descriptor next_halfedge = tmesh.next(h);
                FT e1 = Find_edge_length_and_update_property_map()(tmesh, next_halfedge, edge_lengths);

                // third edge length
                halfedge_descriptor prev_halfedge = tmesh.next(next_halfedge);
                FT e2 = Find_edge_length_and_update_property_map()(tmesh, prev_halfedge, edge_lengths);

                // third point
                FT Px = (e0 + (e2 - e1)) / (2*sqrt(e0));
                FT Py = sqrt(e2 - square(Px));
                unfolded_triangle.P = Point_2(Px, Py);

                return unfolded_triangle;
            }
        };
    typedef Unfold_triangle_3_along_halfedge unfold_triangle_3;

    /*!
     * \brief This functor reconstructs the geodesic source point in the tangent space of the face 'face(h, mesh)'.
     */
    class Reconstruct_source_point_in_triangle_tangent_space{
    public:
        typedef Point_2 result_type;

    public:
        Reconstruct_source_point_in_triangle_tangent_space() {}

        result_type operator() (Surface_mesh& tmesh,
                                halfedge_descriptor h,
                                Edge_property_map& edge_lengths,
                                Face_values_map& face_values)
        {
            FT e0 = Find_edge_length()(tmesh, h, edge_lengths);
            if (is_border(h, tmesh)) {
                std::cerr << "halfedge opposite to " << h << " is on border and hence there is no way to reconstruct the source" << std::endl;
            }

            // find the correct entries in the face_values_map
            face_descriptor face = tmesh.face(h);
            vertex_descriptor A = tmesh.target(h); // this is swapped because target(h) == source(opposite(h))
            vertex_descriptor B = tmesh.source(h);

            int A_loc = vertex_index_in_face(A, face, tmesh);
            int B_loc = vertex_index_in_face(B, face, tmesh);

            FT d2A = face_values[face].d2verts[A_loc];
            FT d2B = face_values[face].d2verts[B_loc];

            // first coordinate of the virtual geodesic source S
            FT Sx = (e0 + (d2A - d2B)) / (2.*sqrt(e0));
            FT d2ASx = d2A - square(Sx);
            FT Sy;

            if (-1e-10 < d2ASx && d2ASx < 0.)  // there should never be negative numbers unless d2A == Sx and numerical errors hit
            {
                //d2ASx = 0.;
                Sy = 0;
            }
            else
            {
                assert(d2ASx >= 0.);
                Sy = -sqrt(d2ASx);
            }

            // Source point in triangle tangent plane
            Point_2 S = {Sx, Sy};
            return S;
        }
    };
    typedef Reconstruct_source_point_in_triangle_tangent_space Reconstruct_source_point;

    /*!
     * \brief This functor constructs the centroid or barycenter of (the points of) an 'Unfolded_triangle_2'.
     */
    class Construct_triangle_centroid_2
    {
    public:
        typedef Point_2 result_type;

    public:
        Construct_triangle_centroid_2() {}

        result_type operator() (Point_2 B, Point_2 P)
        {
            FT cx = (B.x() + P.x()) / 3.;
            FT cy = P.y() / 3;

            return Point_2(cx, cy);
        }
    };

    /*!
     * \brief This functor constructs the centroid or barycenter of a triangle in 3D.
     */
    class Construct_triangle_centroid_3
    {
    public:
        typedef Point_3 result_type;

    public:
        Construct_triangle_centroid_3() {}

        result_type operator() (Surface_mesh& mesh, halfedge_descriptor halfedge)
        {
            // construct barycenter
            FT P_x(0.0);
            FT P_y(0.0);
            FT P_z(0.0);
            for (vertex_descriptor vd : mesh.vertices_around_face(halfedge))
            {
                Point_3 vert = mesh.point(vd);
                P_x += 0.3333333333333333 * vert.x();
                P_y += 0.3333333333333333 * vert.y();
                P_z += 0.3333333333333333 * vert.z();
            }
            return Point_3(P_x, P_y, P_z);
        }
    };

    /*!
     * \brief This function constructs the target points in the unfoldind tangent space.
     */
    Point_2 construct_target_point(halfedge_descriptor opposite_halfedge, Face_location target_face_location, Point_2 B, Point_2 P)
    {
        // get the barycentric coordinates of the target point (relative to some CGAL-intrinsic ordering)
        Barycentric_coordinates coords = target_face_location.second;
        // find the local vertex indices corresponding to the halfedge we unfolded along
        std::array<int, 3> local_vert_idx = get_local_vertex_indices(opposite_halfedge);

        // construct the target point in the tangent plant using the unfolded triangle point and the barycentric coordinates
        FT Tx = coords[local_vert_idx[1]] * B.x() + coords[local_vert_idx[2]] * P.x();
        FT Ty = coords[local_vert_idx[2]] * P.y();

        return Point_2(Tx, Ty);
    }

    /*!
     * \brief The Edge_intersection_test functor computes an 'IntersectionResult' based on
     *  the edge of the unfolded triangle (line A--B) and the visibility line S--Q.
     */
    class Edge_intersection_test
    {
    public:
        typedef typename Kernel::Left_turn_2    Left_turn_2;

    public:
        Edge_intersection_test() {};

        IntersectionResult operator() (Point_2 S, Point_2 Q, Point_2 A, Point_2 B)
        {
            bool left_of_edge = Left_turn_2()(S, A, Q);
            bool not_right_of_edge = Left_turn_2()(S, B, Q);

            IntersectionResult intersection_result;

            // check all the possible cases
            if (left_of_edge && !not_right_of_edge) {
                std::cerr << "Intersection test with edge failed. The source was identified to be both left and right of the edge (which should not be possible)";
            }
            else if (left_of_edge && not_right_of_edge) {
                // source is to the left of the edge
                intersection_result = LEFT_OF_EDGE;
            }
            else if (!left_of_edge && !not_right_of_edge) {
                // source is to the right of the edge
                intersection_result = RIGHT_OF_EDGE;
            }
            else if (!left_of_edge && not_right_of_edge) {
                // source is below the edge such that we get an intersection!
                intersection_result = EDGE_INTERSECTED;
            }

            return intersection_result;
        }
    };

    /// @}

    /// \name Addition and Removal of Source Points
    /// @{

    /*!
     * \brief Add a source point by specifying a integer vertex index
     */
    void add_source(int vertex_index)
    {
        vertex_descriptor vd(vertex_index);

        add_source(vd);
    }

    /*!
     * \brief Add a source point by specifying a 'vertex_descriptor'
     */
    void add_source(vertex_descriptor vd)
    {
        Point_3 source_point = m_mesh.point(vd);

        add_source(source_point);
    }

    /*!
     * \brief Add a source point by specifying a 'Point_3' in the physical domain.
     * \details Note that if the specfied point is not on the mesh, it will be projected onto the mesh.
     */
    void add_source(Point_3 source_point)
    {
        // locate source_point on mesh
        // NOTE: if the point is not on the mesh, it will be projected onto the mesh
        Face_location source_point_location = Polygon_mesh_processing::locate(source_point, m_mesh);

        add_source(source_point_location);
    };

    /*!
     * \brief Add a source point by specifying a 'Face_location'
     */
    void add_source(Face_location source_point)
    {
        m_source_face_locations.push_back(source_point);

        // construct physical point corresponding to target
        Point_3 physical_point = Polygon_mesh_processing::construct_point(source_point, m_mesh);
        m_source_physical_locations.push_back(physical_point);

        m_sources_status = CGAL::SOURCES_SET;
        m_propagation_status = CGAL::INVALID_RESULTS;
    };

    /*!
     * \brief Remove the 'i'th specified source point (0-based)
     */
    void remove_source(int i)
    {
        assert(i < m_source_face_locations.size());

        m_source_face_locations.erase(m_source_face_locations.begin() + i);
        m_source_physical_locations.erase(m_source_physical_locations.begin() + i);

        if (m_source_face_locations.empty()) { m_sources_status = CGAL::NO_SOURCES_SET; }
        m_propagation_status = CGAL::INVALID_RESULTS;
    }

    /*!
     * \brief Write the face locations and physical locations of all specified source points to standard out.
     */
    void print_source_information()
    {
        std::cout << "The following target points were set up:" << std::endl;
        for (int i = 0; i < m_source_face_locations.size(); i++)
        {
            Face_location source = m_source_face_locations[i];
            Point_3 physical_source = m_source_physical_locations[i];
            std::cout << "Source with index " << std::to_string(i)
                      << " in face " << std::to_string(source.first)
                      << " and barycentric coordinates [" << std::to_string(source.second[0])
                      << ", " << std::to_string(source.second[1]) << ", " << std::to_string(source.second[2]) << "]\n"
                      << "\t (Physical location: [" << std::to_string(physical_source[0])
                      << ", " << std::to_string(physical_source.second[1]) << ", "
                      << std::to_string(physical_source.second[2]) << "])" << std::endl;
        }
    }

    /// @}

    void initialize_sources()
    {
        if (m_sources_status != CGAL::SOURCES_SET)
        {
            std::cerr << "No sources set, the algorithm will terminate here." << std::endl;
        }

        for (int i = 0; i < m_source_face_locations.size(); i++)
        {
            // get physical and parametric coordinates of source_point
            Point_3 physical_source_point = m_source_physical_locations[i];
            face_descriptor face = m_source_face_locations[i].first;
            halfedge_descriptor h0 = m_mesh.halfedge(face);

            // set sigma and d for the source triangle
            m_face_values[face].sigma = 0.;
            Point_3 face_barycenter = Construct_triangle_centroid_3()(m_mesh, h0);
            m_face_values[face].d = sqrt(squared_distance(physical_source_point, face_barycenter));

            // set distances to face vertices
            for (halfedge_descriptor h : m_mesh.halfedges_around_face(h0))
            {
                Find_edge_length()(m_mesh, h, m_edge_lengths);
                vertex_descriptor v_idx = m_mesh.source(h);
                int local_v_idx = vertex_index_in_face(v_idx, face, m_mesh);

                FT dist = squared_distance(physical_source_point, m_mesh.point(v_idx));
                m_face_values[face].d2verts[local_v_idx] = dist;

                m_A.push(h);
            }
        }
    }

    bool is_source_face(face_descriptor face)
    {
        for (Face_location source_face_location : m_source_face_locations)
        {
            if (face == source_face_location.first) { return true; }
        }

        return false;
    }

    /// \name Addition and Removal of Target Points
    /// @{

    /*!
     * \brief Add a target point by specifying a integer vertex index
     */
    void add_target(int vertex_index)
    {
        vertex_descriptor vd(vertex_index);

        add_target(vd);
    }

    /*!
     * \brief Add a target point by specifying a 'vertex_descriptor'
     */
    void add_target(vertex_descriptor vd)
    {
        Point_3 target_point = m_mesh.point(vd);

        add_target(target_point);
    }

    /*!
     * \brief Add a target point by specifying a 'Point_3' in the physical domain.
     * \details Note that if the specfied point is not on the mesh, it will be projected onto the mesh.
     */
    void add_target(Point_3 target_point)
    {
        // locate source_point on mesh
        Face_location target_face_location = Polygon_mesh_processing::locate(target_point, m_mesh);

        add_target(target_face_location);
    };

    /*!
     * \brief Add a target point by specifying a 'Face_location'
     */
    void add_target(Face_location target_point)
    {
        m_target_face_locations.push_back(target_point);

        // construct physical point corresponding to target
        Point_3 physical_point = Polygon_mesh_processing::construct_point(target_point, m_mesh);
        m_target_physical_locations.push_back(physical_point);

        m_targets_status = CGAL::TARGETS_SET;
        m_propagation_status = CGAL::INVALID_RESULTS;
    }

    /*!
     * \brief Remove the 'i'th specified target point (0-based)
     */
    void remove_target(int i)
    {
        assert(i < m_target_face_locations.size());

        m_target_face_locations.erase(m_target_face_locations.begin() + i);
        m_target_physical_locations.erase(m_target_physical_locations.begin() + i);

        if (m_target_face_locations.empty()) { m_targets_status = CGAL::NO_TARGETS_SET; }
        m_propagation_status = CGAL::INVALID_RESULTS;
    }

    /*!
     * \brief Write the face locations and physical locations of all specified target points to standard out.
     */
    void print_target_information()
    {
        std::cout << "The following target points were set up:" << std::endl;
        for (int i = 0; i < m_target_face_locations.size(); i++)
        {
            Face_location target = m_target_face_locations[i];
            Point_3 physical_target = m_target_physical_locations[i];
            std::cout << "Target with index " << std::to_string(i)
                      << " in face " << std::to_string(target.first)
                      << " and barycentric coordinates [" << std::to_string(target.second[0])
                      << ", " << std::to_string(target.second[1]) << ", " << std::to_string(target.second[2]) << "]\n"
                      << "\t (Physical location: [" << std::to_string(physical_target[0])
                      << ", " << std::to_string(physical_target.second[1]) << ", "
                      << std::to_string(physical_target.second[2]) << "])" << std::endl;
        }
    }
    /// @}

    /*!
     * \brief Function to determine whether the propagated-to face contains a target point and if so, which target face.
     * \return 'std::pair<bool, int>' of a 'bool' ('true' if target face, 'false' otherwise) and the index of the target point in the target point vector.
     */
    std::pair<bool, int> belongs_to_target_face(halfedge_descriptor h)
    {
        std::pair<bool, int> is_in_target_face = { false, -1 };
        if (m_targets_status == CGAL::NO_TARGETS_SET) { return is_in_target_face; }

        face_descriptor face = m_mesh.face(h);
        for (int i = 0; i < m_target_face_locations.size(); i++)
        {
            if (m_target_face_locations[i].first == face)
            {
                is_in_target_face.first = true;
                is_in_target_face.second = i;
            }
        }

        return is_in_target_face;
    }

    /// \name Propagation
    /// @{

    /*!
     * \brief Function to evaluate the new geodesic distance given an 'IntersectionResult' and 'Unfolded_triangle' source and target points.
     */
    std::pair<FT, FT> get_new_dist(IntersectionResult intersection,
                                   halfedge_descriptor h,
                                   Point_2 C,
                                   Point_2 S)
    {
        face_descriptor face = m_mesh.face(h);
        FT e0 = sqrt(m_edge_lengths[m_mesh.edge(h)]);

        std::pair<FT, FT> prev_geodesic_dist;
        prev_geodesic_dist.first = m_face_values[face].sigma;
        prev_geodesic_dist.second = m_face_values[face].d;

        std::pair<FT, FT> new_geodesic_dist;
        if (intersection == EDGE_INTERSECTED)
        {
            new_geodesic_dist.first = prev_geodesic_dist.first;
            new_geodesic_dist.second = new_geodesic_dist.first + sqrt( squared_distance(C, S) );
        }
        else
        {
            if (intersection == RIGHT_OF_EDGE)
            {
                new_geodesic_dist.first = sqrt( square(S.x()-e0) + square(S.y()) ) + prev_geodesic_dist.first;
                new_geodesic_dist.second = new_geodesic_dist.first + sqrt( square(C.x()-e0) + square(C.y()) );
            }
            else
            {
                // left turn (source is to the left of A, t<0 case)
                new_geodesic_dist.first = sqrt(square(S.x()) + square(S.y())) + prev_geodesic_dist.first;
                new_geodesic_dist.second = new_geodesic_dist.first + sqrt(square(C.x()) + square(C.y()));
            }
        }

        return new_geodesic_dist;
    }

    /*!
     * \brief Function to get the local indices of the vertices in a triangle face. Needed for a consistent enumeration of vertices between adjacent cells.
     * \return 'std::array<int, 3>' of local vertex indices, i.e. permutations of {0, 1, 2}
     */
    std::array<int, 3> get_local_vertex_indices(halfedge_descriptor h)
    {
        face_descriptor face = m_mesh.face(h);
        vertex_descriptor A = m_mesh.source(h);
        vertex_descriptor B = m_mesh.target(h);
        vertex_descriptor P = m_mesh.target(m_mesh.next(h));

        int local_idx_A_in_face = vertex_index_in_face(A, face, m_mesh);
        int local_idx_B_in_face = vertex_index_in_face(B, face, m_mesh);
        int local_idx_P_in_face = vertex_index_in_face(P, face, m_mesh);

        std::array<int, 3> local_vertex_indices = {local_idx_A_in_face, local_idx_B_in_face, local_idx_P_in_face};

        return local_vertex_indices;
    }

    /*!
     * \brief Function to set the squared distances of the source to the vertices of an 'Unfolded_triangle' in the case of an intersection
     */
    void set_squared_vertex_distances_intersection(halfedge_descriptor oppo_h,
                                                   Point_2 P_coords,
                                                   Point_2 S_coords)
    {
        face_descriptor face = m_mesh.face(oppo_h);
        halfedge_descriptor h = m_mesh.opposite(oppo_h);
        face_descriptor prev_face = m_mesh.face(h);

        Face_values prev_face_values = m_face_values[prev_face];
        std::array<FT,3> prev_face_vertex_distances = prev_face_values.d2verts;

        vertex_descriptor A = m_mesh.target(oppo_h);
        int local_idx_A_in_face = vertex_index_in_face(A, face, m_mesh);
        int local_idx_A_in_prev_face = vertex_index_in_face(A, prev_face, m_mesh);
        m_face_values[face].d2verts[local_idx_A_in_face] = prev_face_vertex_distances[local_idx_A_in_prev_face];

        vertex_descriptor B = m_mesh.source(oppo_h);
        int local_idx_B_in_face = vertex_index_in_face(B, face, m_mesh);
        int local_idx_B_in_prev_face = vertex_index_in_face(B, prev_face, m_mesh);
        m_face_values[face].d2verts[local_idx_B_in_face] = prev_face_vertex_distances[local_idx_B_in_prev_face];

        vertex_descriptor P = m_mesh.target(m_mesh.next(oppo_h));
        int local_idx_P_in_face = vertex_index_in_face(P, face, m_mesh);
        m_face_values[face].d2verts[local_idx_P_in_face] = squared_distance(P_coords, S_coords);
    }

    /*!
     * \brief Function to set the squared distances of the source to the vertices of an 'Unfolded_triangle' in the case of a left turn
     */
    void set_squared_vertex_distances_left_turn(halfedge_descriptor h)
    {
        face_descriptor face = m_mesh.face(h);
        std::array<int, 3> local_vertex_indices = get_local_vertex_indices(h);

        m_face_values[face].d2verts[local_vertex_indices[0]] = 0.;
        m_face_values[face].d2verts[local_vertex_indices[1]] = m_edge_lengths[m_mesh.edge(h)];
        m_face_values[face].d2verts[local_vertex_indices[2]] = m_edge_lengths[m_mesh.edge(m_mesh.prev(h))];
    }

    /*!
     * \brief Function to set the squared distances of the source to the vertices of an 'Unfolded_triangle' in the case of a right turn
     */
    void set_squared_vertex_distances_right_turn(halfedge_descriptor h)
    {
        face_descriptor face = m_mesh.face(h);
        std::array<int, 3> local_vertex_indices = get_local_vertex_indices(h);

        m_face_values[face].d2verts[local_vertex_indices[0]] = m_edge_lengths[m_mesh.edge(h)];
        m_face_values[face].d2verts[local_vertex_indices[1]] = 0.;
        m_face_values[face].d2verts[local_vertex_indices[2]] = m_edge_lengths[m_mesh.edge(m_mesh.next(h))];
    }

    /*!
     * \brief Function to treat setting the squared distances of the source to
     * the vertices of an 'Unfolded_triangle' for different 'IntersectionResult's
     */
    void set_squared_vertex_distances(IntersectionResult intersection,
                                      halfedge_descriptor h,
                                      std::pair<FT, FT> new_geodesic_dist,
                                      Point_2 P,
                                      Point_2 S)
    {
        halfedge_descriptor oppo_h = m_mesh.opposite(h);
        face_descriptor next_face = m_mesh.face(oppo_h);
        m_face_values[next_face].sigma = new_geodesic_dist.first;
        m_face_values[next_face].d = new_geodesic_dist.second;

        if (intersection == EDGE_INTERSECTED)
        {
            set_squared_vertex_distances_intersection(oppo_h, P, S);
        }
        else
        {
            if (intersection == RIGHT_OF_EDGE)
            {
                set_squared_vertex_distances_right_turn(oppo_h);
            }
            else
            {
                // left turn (source is to the left of A, t<0 case)
                set_squared_vertex_distances_left_turn(oppo_h);
            }
        }
    }

    /*!
     * \brief Function to update the face values based on the 'IntersectionResult' and the unfolded geometry.
     * \return 'UpdateResult'
     */
    std::pair<UpdateResult, FT> update_face_values(IntersectionResult intersection,
                            halfedge_descriptor h,
                            Point_2 C,
                            Point_2 P,
                            Point_2 S)
    {
        std::pair<FT, FT> new_geodesic_dist = get_new_dist(intersection, h, C, S);
        face_descriptor face = m_mesh.face(m_mesh.opposite(h));

        if (new_geodesic_dist.second < m_face_values[face].d)
        {
            set_squared_vertex_distances(intersection, h, new_geodesic_dist, P, S);
            return { UPDATE, new_geodesic_dist.second };
        }

        return { DO_NOT_UPDATE, FT(-1.) };
    }

    /*!
     * \brief Function to interface the 'Enqueue_policy'. Determines whether new halfedges
     * are to be enqueued in the queues 'A' or 'B' as described in the paper.
     */
    void enqueue_new_halfedges(halfedge_descriptor h, FT new_geodesic_dist, int iter,
                               FT geodesic_dist_to_overall_target,
                               Point_3& prev_face_target_point,
                               Point_3& face_target_point,
                               Point_3& overall_target_point)
    {
        FT geodesic_radius(iter);

        //bool enqueue_in_A = m_Enqueue_policy();
        //EnqueueResult enqueue_result = m_Enqueue_policy(new_geodesic_dist, geodesic_radius);
        EnqueueResult enqueue_result = m_Enqueue_policy(
            new_geodesic_dist,
            geodesic_dist_to_overall_target,
            prev_face_target_point,
            face_target_point,
            overall_target_point);

        if (enqueue_result == CGAL::ENQUEUE_IN_A)
        {
            m_A.push(m_mesh.next(h));
            m_A.push(m_mesh.prev(h));
        }
        else if (enqueue_result == CGAL::ENQUEUE_IN_B)
        {
            m_B.push(m_mesh.next(h));
            m_B.push(m_mesh.prev(h));
        }

    }

    /*!
     * \brief Set an entry of the shortest path tree. For every face, we store its previous face.
     */
    void set_shortest_path_tree_entry(halfedge_descriptor h)
    {
        face_descriptor prev_face = m_mesh.face(h);
        face_descriptor face = m_mesh.face(m_mesh.opposite(h));
        m_shortest_path_tree[face] = prev_face;
    }

    /*!
     * \brief Function to carry out an entire propagation step.
     */
    void propagate_over_halfedge(halfedge_descriptor halfedge, int iter)
    {
        halfedge_descriptor opposite_halfedge = m_mesh.opposite(halfedge);
        if (m_mesh.is_border(opposite_halfedge)) { return; }
        FT e0 = sqrt(Find_edge_length()(m_mesh, halfedge, m_edge_lengths));

        // constructions
        auto unfolded_triangle = unfold_triangle_3()(m_mesh, opposite_halfedge, m_edge_lengths);
        Point_2 A(FT(0.), FT(0.));
        Point_2 B = unfolded_triangle.B;
        Point_2 P = unfolded_triangle.P;

        Point_2 C, Q;
        std::pair<bool, int> is_target_face = belongs_to_target_face(opposite_halfedge);
        if (is_target_face.first)
        {
            Face_location target_face_location = m_target_face_locations[is_target_face.second];
            C = construct_target_point(opposite_halfedge, target_face_location, B, P);
            Q = C;
        }
        else
        {
            C = Construct_triangle_centroid_2()(B, P);
            Q = m_Visibility_heuristic(m_mesh, m_edge_lengths, opposite_halfedge, P, C);
        }
        Point_2 S = Reconstruct_source_point()(m_mesh, halfedge, m_edge_lengths, m_face_values);

        // intersection test
        auto intersection = Edge_intersection_test()(S, Q, A, B);
        std::pair<UpdateResult, FT> update_result = update_face_values(intersection, halfedge, C, P, S);
        if (update_result.first == UPDATE)
        {
            Point_3 prev_face_target_point = Construct_triangle_centroid_3()(m_mesh, halfedge);
            Point_3 face_target_point = Construct_triangle_centroid_3()(m_mesh, opposite_halfedge);
            Point_3 overall_target_point = m_target_physical_locations[0];

            enqueue_new_halfedges(opposite_halfedge,
                                  update_result.second,
                                  iter,
                                  m_face_values[m_target_face_locations[0].first].d,
                                  prev_face_target_point,
                                  face_target_point,
                                  overall_target_point);

            if (m_build_shortest_path_tree) { set_shortest_path_tree_entry(halfedge); }
            int face_idx = m_mesh.face(opposite_halfedge);
            if (m_with_observer) { m_observer.increment(face_idx); }
        }
    }

    /*!
     * \brief Run the propagation algorithm. This function builds the 'Face_values_map' and hence the geodesic distances.
     */
    void propagate_geodesic_source()
    {
        if (m_propagation_status == CGAL::VALID_RESULTS)
        {
            return; // the results are still valid from an earlier computation (neither sources nor targets were changed)
        }

        // initialize the algorithm
        initialize_sources();

        // check for targets
        if (m_targets_status == CGAL::NO_TARGETS_SET)
        {
            std::cout << "No targets specified, will build the whole geodesic distance tree." << std::endl;
        }

        int iter = 1;
        while (!m_A.empty())
        {
            if (m_with_observer) { m_observer.add_iteration_step(); }

            // iterate over halfedges in queue and update face values
            while (!m_A.empty()) {
                halfedge_descriptor h = m_A.front();
                m_A.pop();

                if (m_Skip_condition() == CGAL::SKIP)
                {
                    m_B.push(h);
                }
                else
                {
                    propagate_over_halfedge(h, iter);
                }
            }

            std::swap( m_A, m_B );
            std::queue<halfedge_descriptor> empty;
            std::swap( m_B, empty );
            iter++;
        }

        m_propagation_status = CGAL::VALID_RESULTS;
    }

    /// @}

    /// \name Output Queries
    /// @{

    /*!
     * \brief Get the geodesic distances to all face target points (the barycenter unless specified otherwise).
     * \details All entries corresponding to faces that have not been visited by the propagation algorithm are set to -1.
     */
    std::vector<double> get_geodesic_distances()
    {
        std::vector<double> distances;

        get_geodesic_distances(distances);

        return distances;
    }

    /*!
     * \brief Get the geodesic distances to all face target points (the barycenter unless specified otherwise).
     * \details All entries corresponding to faces that have not been visited by the propagation algorithm are set to -1.
     */
    void get_geodesic_distances(std::vector<double>& distances)
    {
        assert(m_propagation_status == CGAL::VALID_RESULTS);

        distances.resize(m_mesh.num_faces());

        for (face_descriptor fd : faces(m_mesh))
        {
            double d = m_face_values[fd].d;

            // set the geodesic distances for all faces that have not been visited to -1
            if (d > 1e100) { d = -1.; }

            distances[fd.idx()] = d;
        }

    }

    /*!
     * \brief Get the geodesic distances to all specified target points.
     */
    std::vector<FT> get_geodesic_distance_to_targets()
    {
        assert(m_propagation_status == CGAL::VALID_RESULTS);

        std::vector<FT> target_distances;
        for (Face_location target_point : m_target_face_locations)
        {
            target_distances.push_back(m_face_values[target_point.first].d);
        }

        return target_distances;
    }

    /*!
     * \brief Returns the index of the closest specified target.
     */
    int get_index_of_closest_target()
    {
        assert(m_propagation_status == CGAL::VALID_RESULTS);

        int index = 0;

        for (int i = 1; i < m_target_face_locations.size(); i++)
        {
            if (m_face_values[i].d < m_face_values[i-1].d)
            {
                index = i;
            }
        }
    }

    /*!
     * \brief Returns the face location of the closest specified target.
     */
    Face_location closest_target_face_location()
    {
        int closest_target_index = get_index_of_closest_target();

        return m_target_face_locations[closest_target_index];
    }

    /*!
     * \brief Returns the physical location of the closest specified target.
     */
    Point_3 closest_target_point()
    {
        int closest_target_index = get_index_of_closest_target();

        return m_target_physical_locations[closest_target_index];
    }

    /*!
     * \brief Extracts the shortest path face sequence to a specified target point given its index.
     */
    std::vector<face_descriptor> extract_shortest_path_face_sequence(int target_index)
    {
        face_descriptor target_face_descriptor = m_target_face_locations[target_index].first;

        return extract_shortest_path_face_sequence(target_face_descriptor);
    }

    /*!
     * \brief Extracts the shortest path face sequence to any face as a 'std::vector<face_descriptor>'.
     */
    std::vector<face_descriptor> extract_shortest_path_face_sequence(face_descriptor target_face)
    {
        assert( "Shortest path tree was not build during propagation." && m_with_observer );

        std::vector<face_descriptor> face_location_sequence;
        face_location_sequence.push_back(target_face);

        face_descriptor prev_face(target_face);
        do {
            prev_face = m_shortest_path_tree[prev_face];
            face_location_sequence.push_back(prev_face);
        } while (!is_source_face(prev_face));

        return face_location_sequence;
    }

    /*!
     * \brief Extracts a 'std::vector<bool>' indicating whether the face with a given index belonds to the geodesic of the specified target.
     */
    std::vector<bool> extract_shortest_path_face_indicator_map(int target_index)
    {
        face_descriptor target_face_descriptor = m_target_face_locations[target_index].first;

        return extract_shortest_path_face_indicator_map(target_face_descriptor);
    }

    /*!
     * \brief Extracts a 'std::vector<bool>' indicating whether the face with a given index belonds to the geodesic starting at a given face target.
     */
    std::vector<bool> extract_shortest_path_face_indicator_map(face_descriptor target_face)
    {
        std::vector<face_descriptor> face_sequence = extract_shortest_path_face_sequence(target_face);

        std::vector<bool> face_indicator_map(m_mesh.number_of_faces());
        std::fill(face_indicator_map.begin(), face_indicator_map.end(), false);

        for (face_descriptor fd : face_sequence)
        {
            face_indicator_map[fd] = true;
            std::cout << fd << std::endl;
        }

        return face_indicator_map;
    }

    /// @}
};

}

template<class K, class SurfaceMesh>
class VTKWriter
{
public:
    typedef K Kernel;
    typedef SurfaceMesh Surface_mesh;
    typedef typename Kernel::Point_3 Point_3;

    typedef boost::graph_traits<Surface_mesh> Graph_traits;

    typedef typename Graph_traits::vertex_descriptor vertex_descriptor;
    typedef typename Graph_traits::halfedge_descriptor halfedge_descriptor;
    typedef typename Graph_traits::face_descriptor face_descriptor;

public:
    VTKWriter() {};

    void operator() (Surface_mesh& mesh, std::vector<double> face_data, std::string filename)
    {
        WriteMeshData(mesh, filename);
        WriteFaceData(mesh, face_data, filename);
    }

    void operator() (Surface_mesh& mesh, std::vector<bool> face_data, std::string filename)
    {
        WriteMeshData(mesh, filename);
        WriteFaceData(mesh, face_data, filename);
    }

    void WriteMeshData(Surface_mesh& mesh, std::string filename)
    {
        std::ofstream out(filename);

        // header for (legacy) vtk file
        out << "# vtk DataFile Version 2.0\n";
        out << "description: \n";
        out << "ASCII\n";

        // write geometry/topology
        out << "DATASET UNSTRUCTURED_GRID\n";
        out << "POINTS " << mesh.num_vertices() << " float\n";

        for (vertex_descriptor vd : vertices(mesh))
        {
            Point_3 v = mesh.point(vd);
            out << v.x() << " "
                << v.y() << " "
                << v.z() << std::endl;
        }

        // 4*_num_faces is the "size of the cell list"
        // see https://vtk.org/wp-content/uploads/2015/04/file-formats.pdf
        out << "CELLS " << mesh.num_faces() << " " << 4 * mesh.num_faces() << std::endl;
        for (face_descriptor fd : faces(mesh))
        {
            out << "3 ";
            halfedge_descriptor h0 = mesh.halfedge(fd);
            for (halfedge_descriptor hd : CGAL::halfedges_around_face(h0, mesh))
            {
                out << mesh.source(hd).idx() << " ";
            }
            out << std::endl;
        }
        // write cell types (5 = VTK_TRIANGLE)
        out << "CELL_TYPES " << mesh.num_faces() << std::endl;
        for (int face_num = 0; face_num < mesh.num_faces(); ++face_num)
        {
            out << "5" << std::endl;
        }

        out.close();
    }

    void WriteFaceData(Surface_mesh& mesh, std::vector<double> face_data, std::string filename)
    {
        std::ofstream out(filename, std::ios_base::app);

        out << "CELL_DATA " << mesh.num_faces() << std::endl;
        out << "SCALARS " << "cell_scalars " << "float " << "1" << std::endl;
        out << "LOOKUP_TABLE default" << std::endl;

        for (int i = 0; i < mesh.num_faces(); ++i)
        {
            out << face_data[i] << std::endl;
        }

        out.close();
    }

    void WriteFaceData(Surface_mesh& mesh, std::vector<bool> face_data, std::string filename)
    {
        std::ofstream out(filename, std::ios_base::app);

        out << "CELL_DATA " << mesh.num_faces() << std::endl;
        out << "SCALARS " << "cell_scalars " << "float " << "1" << std::endl;
        out << "LOOKUP_TABLE default" << std::endl;

        for (int i = 0; i < mesh.num_faces(); ++i)
        {
            if (face_data[i])
            {
                out << 1 << std::endl;
            }
            else
            {
                out << 0 << std::endl;
            }
        }

        out.close();
    }
};

#endif // CGAL_SURFACE_MESH_APPROXIMATE_SHORTEST_PATH_SURFACE_MESH_APPROXIMATE_SHORTEST_PATH_H
