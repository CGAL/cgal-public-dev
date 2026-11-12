// Copyright (c) 2026  GeometryFactory (France)
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Pranav Jain

/**
 * @file splat_surface_reconstruction.h
 * @brief Splat surface reconstruction for point clouds with normals.
 *
 * This header provides the public reconstruction entry point together with
 * an internal dense box-grid acceleration structure used for preprocessing,
 * local normal estimation, and splat-size estimation.
 */
#ifndef CGAL_SPLAT_SURFACE_RECONSTRUCTION_H
#define CGAL_SPLAT_SURFACE_RECONSTRUCTION_H

// #include <CGAL/license/Splat_surface_reconstruction_3.h>
#include <CGAL/property_map.h>
#include <CGAL/Kernel_traits.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/intersections.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/boost/graph/dijkstra_shortest_path.h>
#include <CGAL/number_utils.h>
#include <CGAL/Plane_3.h>
#include <CGAL/Projection_traits_3.h>
#include <CGAL/Kernel/global_functions_3.h>
#include <CGAL/Polygon_2.h>
#include <vector>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <optional>
#include <queue>
#include <memory>

#include <CGAL/internal/Splat_surface_reconstruction_3/geometry_utils.h>
#include <CGAL/internal/Splat_surface_reconstruction_3/splat_surface_approximation.h>
#include <CGAL/internal/Splat_surface_reconstruction_3/mesh_surface_approximation.h>


namespace CGAL {

  template <typename PointRange, typename NormalRange, typename TriangleMesh>
  struct Splat_surface_reconstruction_3_test_access;

  /**
   * @ingroup PkgSplatSurfaceReconstruction3Classes
   *
   * @brief Incrementally reconstructs a polygon mesh from splat candidates.
   *
   * The reconstruction grows a halfedge front from an initial seed. Candidate
   * vertices are generated from pairs of existing mesh vertices and accepted
   * according to proximity, projection, priority, and local topology tests.
   *
   * During growth, the class maintains an explicit halfedge connectivity graph.
   * Faces are created only after the growth phase by traversing the resulting
   * `next()` cycles.
   *
   * @tparam PointRange  Input point range type.
   * @tparam NormalRange Input normal range type.
   * @tparam TriangleMesh Output mesh type model of `HalfedgeListGraph`, `FaceListGraph`, and `MutableFaceGraph`
  */
  template <typename PointRange, typename NormalRange, typename TriangleMesh>
  class Splat_surface_reconstruction_3 {
    private:
      friend struct Splat_surface_reconstruction_3_test_access<PointRange, NormalRange, TriangleMesh>;

    public:
      using Point_3 = typename PointRange::value_type;
      using Kernel = typename Kernel_traits<Point_3>::Kernel;
      using Vector_3 = typename Kernel::Vector_3;
      using FT = typename Kernel::FT;
      using Grid = Box_grid<Kernel>;
      using Index = typename Grid::Index;

    private:
      using Vector_2 = typename Kernel::Vector_2;
      using Point_2  = typename Kernel::Point_2;

      using vertex_descriptor = typename boost::graph_traits<TriangleMesh>::vertex_descriptor;
      using edge_descriptor = typename boost::graph_traits<TriangleMesh>::edge_descriptor;
      using halfedge_descriptor = typename boost::graph_traits<TriangleMesh>::halfedge_descriptor;
      using face_descriptor = typename boost::graph_traits<TriangleMesh>::face_descriptor;

      using Surface_approximation = typename CGAL::internal::Splat_surface_reconstruction_3::Surface_approximation<Kernel>;
      using Splat_surface_approximation = typename CGAL::internal::Splat_surface_reconstruction_3::Splat_surface_approximation<Kernel>;
      using Mesh_surface_approximation = typename CGAL::internal::Splat_surface_reconstruction_3::Mesh_surface_approximation<Kernel, TriangleMesh>;
      using Initial_seed = typename Surface_approximation::Initial_seed;
      using Feature_segment = typename Surface_approximation::Feature_segment;

      struct Candidate {
        Point_3 position;
        Vector_3 normal;
        vertex_descriptor first;
        vertex_descriptor second;
        int priority = 0;

        Candidate(const Point_3& position,
                  const Vector_3& normal,
                  vertex_descriptor first,
                  vertex_descriptor second)
          : position(position),
            normal(normal),
            first(first),
            second(second) {}
      };


      struct FeatureVertexCandidate {
        Point_3 position;
        Vector_3 normal;
        std::size_t feature_segment_index;
        bool has_feature_segment;
        vertex_descriptor parent_vertex;
        bool has_parent_vertex;

        FeatureVertexCandidate(const Point_3& position,
          const Vector_3& normal,
          bool has_feature_segment = false,
          bool has_parent_vertex = false,
          std::size_t feature_segment_index = 0,
          vertex_descriptor parent_vertex = vertex_descriptor())
          : position(position)
          , normal(normal)
          , has_feature_segment(has_feature_segment)
          , has_parent_vertex(has_parent_vertex)
          , feature_segment_index(feature_segment_index)
          , parent_vertex(parent_vertex)
          {}
      };

    public:

      enum class Initialization_method
      {
        DEFAULT,
        FEATURE_AWARE
      };

      /**
       * @brief Constructs a reconstruction object and initializes mesh properties.
       *
       * Adds the per-vertex normal property map and per-edge property map,
       * initializes the spatial lookup structure, and creates the initial seed.
       *
       * @param grid Spatial grid containing input points, normals, and splat sizes.
       * @param output_mesh Mesh receiving the reconstructed surface.
      */
      Splat_surface_reconstruction_3(const std::vector<Point_3>& points,
        const std::vector<Vector_3>& normals,
        const FT d,
        TriangleMesh& output_mesh,
        Initialization_method init_method = Initialization_method::DEFAULT,
        FT feature_angle_minimum = FT(30.0))
      : mesh_(output_mesh),
        grid_(make_grid(points, d)),
        points_pm_(get(CGAL::vertex_point, mesh_)),
        surface_approximation_(std::make_unique<Splat_surface_approximation>(points, normals, grid_)) {


        // TODO add these writes back
        // grid_.write_splat_mesh_obj("splats_before.obj");
        // std::vector<FT> splat_sizes = grid_.estimate_individual_splat_sizes(); // estimate individual splat sizes based on local point distribution
        // std::cout<<"Estimated individual splat sizes for " << splat_sizes.size() << " points." << std::endl;
        // FT mean = FT(0.0);
        // for (const auto& s : splat_sizes) {
        //   mean += s;
        // }
        // mean /= splat_sizes.size();
        // std::cout << "Mean splat size: " << mean << std::endl;
        // grid_.write_splat_mesh_obj("splats_individual.obj");

        // grid_.write_splat_sizes_to_file("splat_sizes.txt");

        bool created = false;
        std::tie(normals_pm_, created) = mesh_.template add_property_map<vertex_descriptor, Vector_3>("v:normal", CGAL::NULL_VECTOR);

        bool created_w = false;
        std::tie(edge_weight_pm_, created_w) = mesh_.template add_property_map<edge_descriptor, FT>("e:weight", FT(1));

        if (created) {
          std::cout << "Created vertex normal property map." << std::endl;
        }

        mesh_vertices_per_cell_.resize(grid_.number_of_cells());

        if (init_method == Initialization_method::FEATURE_AWARE) {
          seed_from_features(feature_angle_minimum);
        }
        else {
          seed_from_grid();
        }
      }

      /**
       * @brief Constructs a reconstruction object and initializes mesh properties.
       *
       * Adds the per-vertex normal property map and per-edge property map,
       * initializes the spatial lookup structure, and creates the initial seed.
       *
       * @param grid Spatial grid containing input points, normals, and splat sizes.
       * @param output_mesh Mesh receiving the reconstructed surface.
      */
      Splat_surface_reconstruction_3(
        TriangleMesh& input_mesh,
        const FT d,
        TriangleMesh& output_mesh,
        Initialization_method init_method = Initialization_method::DEFAULT,
        FT feature_angle_minimum = FT(30.0))
      : mesh_(output_mesh),
        grid_(make_grid(input_mesh, d)),
        points_pm_(get(CGAL::vertex_point, mesh_)),
        surface_approximation_(std::make_unique<Mesh_surface_approximation>(input_mesh, grid_)) 
      {

        bool created = false;
        std::tie(normals_pm_, created) = mesh_.template add_property_map<vertex_descriptor, Vector_3>("v:normal", CGAL::NULL_VECTOR);

        bool created_w = false;
        std::tie(edge_weight_pm_, created_w) = mesh_.template add_property_map<edge_descriptor, FT>("e:weight", FT(1));

        if (created) {
          std::cout << "Created vertex normal property map." << std::endl;
        }

        mesh_vertices_per_cell_.resize(grid_.number_of_cells());

        if (init_method == Initialization_method::FEATURE_AWARE) {
          seed_from_features(feature_angle_minimum);
        }
        else {
          seed_from_grid();
        }
      }

      /**
       * @brief Executes the complete reconstruction pipeline.
       *
       * Candidates are processed by priority until no candidates remain. Accepted
       * vertices are inserted into the halfedge graph and new candidates are
       * generated around them. Once growth stops, open halfedge cycles are converted
       * into polygonal faces.
      */
      void run() {
        if (!seeded_) {
          std::cerr << "No initial seed was found. Mesh growth did not start.\n";
          return;
        }
        else {
          std::cout << "Initial seed vertices created. Starting mesh growth...\n";
        }

        std::size_t rejected_near = 0;
        std::size_t rejected_proj = 0;
        std::size_t accepted = 0;

        int iterationIndex = 0;

        // export the splats for visualization
        // export_halfedge_mesh_json(mesh_, "initial_splats.json");
        // CGAL_assertion(false);

        CGAL_assertion(CGAL::is_valid_halfedge_graph(mesh_, true));

        while (!candidate_queues_empty()) {
          iterationIndex++;
          Candidate cand(Point_3(0,0,0), CGAL::NULL_VECTOR, vertex_descriptor(), vertex_descriptor()); // Dummy initialization
          if (!pop_candidate_from_queue(cand)) {
            break;
          }
          size_t target_1 = 152;
          size_t target_2 = 535;
          bool is_target = (cand.first.idx() == target_1 || cand.second.idx () == target_1)
              && (cand.first.idx() == target_2 || cand.second.idx () == target_2);

          if (is_target) {
            std::cout << "Priority for candidate with parents " << cand.first << " and " << cand.second << ": " << cand.priority << std::endl;
          }

          if (has_vertex_near(cand.position)) {
            rejected_near++;
            //std::cout << "Candidate rejected due to nearby vertex: " << cand.first << ", " << cand.second << std::endl;
            if (is_target) {
              std::cout << "Candidate rejected due to nearby vertex: " << cand.first << ", " << cand.second << std::endl;
              //export_halfedge_mesh_json(mesh_, "candidate_rejected_near.json");
              //CGAL_assertion(false);
            }
            continue;
          }

          const Point_3 parent_1 = get(points_pm_, cand.first);
          const Point_3 parent_2 = get(points_pm_, cand.second);

          if (!projection_check(cand.position, parent_1, surface_approximation_->box_normal(cand.position))
            || !projection_check(cand.position, parent_2, surface_approximation_->box_normal(cand.position))
            || !projection_check(cand.position, parent_1, surface_approximation_->box_normal(parent_1))
            || !projection_check(cand.position, parent_2, surface_approximation_->box_normal(parent_1))
            || !projection_check(cand.position, parent_1, surface_approximation_->box_normal(parent_2))
            || !projection_check(cand.position, parent_2, surface_approximation_->box_normal(parent_2))) {
            rejected_proj++;
            if (is_target) {
              std::cout << "Candidate rejected due to projection check: " << cand.first << ", " << cand.second << std::endl;
            }
            continue;
          }

          int new_priority = compute_priority(cand);
          if (new_priority != cand.priority) {
            cand.priority = new_priority;
            push_candidate_in_queue(cand);
            if (is_target) {
              std::cout << "Candidate priority updated: " << cand.first << ", " << cand.second << " -> " << cand.priority << std::endl;
            }
            continue;
          }

          if (is_target) {
            std::cout << "Candidate accepted: " << cand.first << ", " << cand.second << std::endl;
          }

          vertex_descriptor nv = mesh_.add_vertex(cand.position);
          put(points_pm_, nv, cand.position);
          put(normals_pm_, nv, cand.normal);

          if (!CGAL::is_valid_halfedge_graph(mesh_, false)) {
            CGAL::is_valid_halfedge_graph(mesh_, true);
            export_halfedge_mesh_json(mesh_, "invalid_mesh.json");
          }

          CGAL_assertion(CGAL::is_valid_halfedge_graph(mesh_, false));


          accepted++;
          if (!build_graph(cand.first, cand.second, nv, cand.normal)) {
            remove_vertex(nv, mesh_);
            continue;
          }

          assert_valid_halfedges();

          insert_mesh_vertex(nv);

          push_candidates_from_vertex(nv);

          if (num_vertices(mesh_) % 1000 == 0) { // Print progress every 1000 vertices
            std::cout << num_vertices(mesh_) << " vertices in mesh." << "\n";
            export_halfedge_mesh_json(mesh_, ("progress_" + std::to_string(num_vertices(mesh_)) + ".json").c_str());
          }


          CGAL_assertion(CGAL::is_valid_halfedge_graph(mesh_, false));
        }

        //assert_valid_triangle_loops();

        CGAL_assertion(CGAL::is_valid_halfedge_graph(mesh_, false));

        export_halfedge_mesh_json(mesh_, "before_face_filling.json");

        std::cout<<"Setting Faces..." << std::endl;
        fill_faces_from_next_cycles();

        //assert_valid_triangle_faces();

        //assert_final_mesh_manifold();

        export_halfedge_mesh_json(mesh_, "final_mesh.json");

        std::cout << "Final mesh size: " << num_vertices(mesh_) << " vertices, " << num_edges(mesh_) << " edges, " << num_faces(mesh_) << " faces.\n";
        std::cout << "  Accepted: " << accepted << ", Rejected (near): " << rejected_near << ", Rejected (proj): " << rejected_proj << "\n";
      }

    private:

      Grid make_grid(const std::vector<Point_3>& points,
                const FT d)
      {
        const FT padding = FT(1.1) * d;

        const auto bbox = CGAL::bounding_box(points.begin(), points.end());

        return Grid(d,
                    FT(bbox.xmin()) - padding,
                    FT(bbox.ymin()) - padding,
                    FT(bbox.zmin()) - padding,
                    FT(bbox.xmax()) + padding,
                    FT(bbox.ymax()) + padding,
                    FT(bbox.zmax()) + padding);
      }

      Grid make_grid(const TriangleMesh& mesh,
                const FT d)
      {
        const FT padding = FT(1.1) * d;

        std::vector<Point_3> points;
        for (const vertex_descriptor v : vertices(mesh)) {
          points.push_back(mesh.point(v));
        }

        const auto bbox = CGAL::bounding_box(points.begin(), points.end());

        return Grid(d,
                    FT(bbox.xmin()) - padding,
                    FT(bbox.ymin()) - padding,
                    FT(bbox.zmin()) - padding,
                    FT(bbox.xmax()) + padding,
                    FT(bbox.ymax()) + padding,
                    FT(bbox.zmax()) + padding);
      }

      bool is_halfedge_inside(halfedge_descriptor h, vertex_descriptor v0, vertex_descriptor v1, vertex_descriptor v2) const {
        Vector_3 normal = normals_pm_[v0] + normals_pm_[v1] + normals_pm_[v2];
        normal = normal / std::sqrt(normal.squared_length());

        // Find the orientation of the triangle formed by the vertices with respect to the normal
        const auto& p0 = get(points_pm_, v0);
        const auto& p1 = get(points_pm_, v1);
        const auto& p2 = get(points_pm_, v2);

        return CGAL::orientation(p0, p1, p2, CGAL::ORIGIN + normal) == CGAL::POSITIVE;
      }

      // Asserts that all triangle faces in the mesh have valid halfedge connections
      void assert_valid_triangle_faces() const {
        size_t num_halfedges = halfedges(mesh_).size();

        std::vector<bool> halfedge_visited(num_halfedges, false);
        
        for (halfedge_descriptor h0 : halfedges(mesh_)) {
          const auto h0_target = target(h0, mesh_);
          const auto h0_source = source(h0, mesh_);

          for (halfedge_descriptor h1 : halfedges_around_source(h0_target, mesh_)) {

            // Not a halfedge that can complete the triangle
            if (source(h1, mesh_) != h0_target) {
              continue;
            }

            const auto h1_target = target(h1, mesh_);

            for (halfedge_descriptor h2 : halfedges_around_source(h1_target, mesh_)) {
              
              // Not a halfedge that can complete the triangle
              if (source(h2, mesh_) != h1_target) {
                continue;
              }

              const auto h2_target = target(h2, mesh_);

              if (!is_halfedge_inside(h0, h0_target, h1_target, h2_target) ||
                  !is_halfedge_inside(h1, h1_target, h2_target, h0_target) ||
                  !is_halfedge_inside(h2, h2_target, h0_target, h1_target)) {
                continue;
              }
              
              if (h2_target != h0_source) {
                continue;
              }

              bool loop_correct = (next(h0, mesh_) == h1 && next(h1, mesh_) == h2 && next(h2, mesh_) == h0);

              // Check that the edges are oriented correctly
              Point_3 p_a = points_pm_[h0_source];
              Point_3 p_b = points_pm_[h0_target];
              Point_3 p_c = points_pm_[h1_target];

              Vector_3 normal_a = normals_pm_[h0_source];
              Vector_3 normal_b = normals_pm_[h0_target];
              Vector_3 normal_c = normals_pm_[h1_target];
              Vector_3 final_normal = normal_a + normal_b + normal_c;
              final_normal = final_normal / std::sqrt(final_normal.squared_length());

              CGAL::Orientation orientation = orient3(p_a, p_b, p_c, final_normal);

              bool orientation_correct = (orientation == CGAL::LEFT_TURN);

              if (orientation_correct != loop_correct) {
                // export_halfedge_mesh_json(mesh_, "invalid_orientation.json");
                // std::cout << "Invalid orientation detected for triangle: h0=" << h0 << ", h1=" << h1 << ", h2=" << h2 << std::endl;
                // std::cout << h0 << ": " << h0_source << "->" << h0_target << " (next: " << next(h0, mesh_) << ")" << std::endl;
                // std::cout << h1 << ": " << source(h1, mesh_) << "->" << target(h1, mesh_) << " (next: " << next(h1, mesh_) << ")" << std::endl;
                // std::cout << h2 << ": " << source(h2, mesh_) << "->" << target(h2, mesh_) << " (next: " << next(h2, mesh_) << ")" << std::endl;
                // std::cout << "Orientation: " << orientation << std::endl;
                // std::cout << "Points: p_a=" << p_a << ", p_b=" << p_b << ", p_c=" << p_c << std::endl;
                // std::cout << "Final normal: " << final_normal << std::endl;
              }

              // If the orientation does not match the loop correctness, report an error.
              // A loop can have incorrect orientation but only if the loop is also incorrect.
              // However, if the loop is incorrect, this does not imply that the orientation is incorrect.
              //CGAL_assertion(implies(loop_correct, orientation_correct));
            }
          }
        }
        std::cout << "Finished checking all triangle faces." << std::endl;
      }

      bool implies(bool b1, bool b2) const {
        return !b1 || b2;
      }

      void assert_valid_triangle_loops() const {
        size_t num_halfedges = halfedges(mesh_).size();

        std::vector<bool> halfedge_visited(num_halfedges, false);
        
        for (halfedge_descriptor h0 : halfedges(mesh_)) {
          const auto h0_target = target(h0, mesh_);
          const auto h0_source = source(h0, mesh_);

          for (halfedge_descriptor h1 : halfedges_around_source(h0_target, mesh_)) {

            // Not a halfedge that can complete the triangle
            if (source(h1, mesh_) != h0_target) {
              continue;
            }

            const auto h1_target = target(h1, mesh_);

            for (halfedge_descriptor h2 : halfedges_around_source(h1_target, mesh_)) {
              
              // Not a halfedge that can complete the triangle
              if (source(h2, mesh_) != h1_target) {
                continue;
              }

              const auto h2_target = target(h2, mesh_);

              if (!is_halfedge_inside(h0, h0_target, h1_target, h2_target) ||
                  !is_halfedge_inside(h1, h1_target, h2_target, h0_target) ||
                  !is_halfedge_inside(h2, h2_target, h0_target, h1_target)) {
                continue;
              }
              
              if (h2_target != h0_source) {
                continue;
              }

              bool loop_correct = (next(h0, mesh_) == h1 && next(h1, mesh_) == h2 && next(h2, mesh_) == h0);

              // Check that the edges are oriented correctly
              Point_3 p_a = points_pm_[h0_source];
              Point_3 p_b = points_pm_[h0_target];
              Point_3 p_c = points_pm_[h1_target];

              Vector_3 normal_a = normals_pm_[h0_source];
              Vector_3 normal_b = normals_pm_[h0_target];
              Vector_3 normal_c = normals_pm_[h1_target];
              Vector_3 final_normal = normal_a + normal_b + normal_c;
              final_normal = final_normal / std::sqrt(final_normal.squared_length());

              CGAL::Orientation orientation = orient3(p_a, p_b, p_c, final_normal);

              bool orientation_correct = (orientation == CGAL::LEFT_TURN);

              // if (orientation_correct != loop_correct) {
              //   export_halfedge_mesh_json(mesh_, "invalid_orientation.json");
              //   std::cout << "Invalid orientation detected for triangle: h0=" << h0 << ", h1=" << h1 << ", h2=" << h2 << std::endl;
              //   std::cout << h0 << ": " << h0_source << "->" << h0_target << " (next: " << next(h0, mesh_) << ")" << std::endl;
              //   std::cout << h1 << ": " << source(h1, mesh_) << "->" << target(h1, mesh_) << " (next: " << next(h1, mesh_) << ")" << std::endl;
              //   std::cout << h2 << ": " << source(h2, mesh_) << "->" << target(h2, mesh_) << " (next: " << next(h2, mesh_) << ")" << std::endl;
              //   std::cout << "Orientation: " << orientation << std::endl;
              //   std::cout << "Points: p_a=" << p_a << ", p_b=" << p_b << ", p_c=" << p_c << std::endl;
              //   std::cout << "Final normal: " << final_normal << std::endl;
              // }

              // The loop being correct implies that the orientation should also be correct.
              // If the loop is not correct, the orientation may or may not be correct.
              //CGAL_assertion(implies(loop_correct, orientation_correct));
            }
          }
        }
        std::cout << "Finished checking all triangle loops." << std::endl;
      }

      void assert_final_mesh_manifold() const {
        // 1. Basic structural validity.
        CGAL_assertion(CGAL::is_valid_polygon_mesh(mesh_));

        // 2. The output should consist entirely of triangles.
        CGAL_assertion(CGAL::is_triangle_mesh(mesh_));

        // 3. No isolated vertices.
        for (const auto v : vertices(mesh_))
        {
          CGAL_assertion(halfedge(v, mesh_) != boost::graph_traits<TriangleMesh>
                                        ::null_halfedge());
        }

        // 4. No wire edges (edges with zero incident faces).
        for (const auto e : edges(mesh_))
        {
          const auto h = halfedge(e, mesh_);
          const bool has_face_1 = face(h, mesh_) !=
                                  boost::graph_traits<TriangleMesh>::null_face();
          const bool has_face_2 = face(opposite(h, mesh_), mesh_) !=
                                  boost::graph_traits<TriangleMesh>::null_face();

          CGAL_assertion(has_face_1 || has_face_2);
        }

        // All halfedges must be set correctly
        for (const auto h : halfedges(mesh_))
        {
          const auto h_next = next(h, mesh_);
          const auto h_prev = next(h_next, mesh_);
          CGAL_assertion(opposite(opposite(h, mesh_), mesh_) == h);

          const auto h_source = source(h, mesh_);
          const auto h_target = target(h, mesh_);
          const auto h_next_source = source(h_next, mesh_);

          CGAL_assertion(h_source == target(opposite(h, mesh_), mesh_));
          CGAL_assertion(h_target == source(opposite(h, mesh_), mesh_));

          CGAL_assertion(h_target == h_next_source);
        }

        std::size_t non_manifold_edges = 0;

        for (const auto h : halfedges(mesh_))
        {
          const auto a = source(h, mesh_);
          const auto b = target(h, mesh_);

          std::vector<std::size_t> incident_faces;

          // Find all faces containing the directed edge a -> b.
          for (const auto f : faces(mesh_))
          {
            auto hf = halfedge(f, mesh_);

            if (hf == TriangleMesh::null_halfedge())
              continue;

            const auto h0 = hf;
            std::size_t steps = 0;
            const std::size_t max_steps = num_halfedges(mesh_);

            do
            {
              if (steps++ >= max_steps)
              {
                std::cout << "Exceeded max steps in face cycle for face " << f.idx() << std::endl;
                CGAL_assertion(false);
                break;
              }

              if (source(hf, mesh_) == a && target(hf, mesh_) == b)
              {
                incident_faces.push_back(f.idx());
                break;
              }

              hf = next(hf, mesh_);

            } while (hf != h0);
          }

          if (incident_faces.size() > 2)
          {
            ++non_manifold_edges;

            std::cerr << "Directed edge "
                      << a.idx() << " -> " << b.idx()
                      << " belongs to " << incident_faces.size()
                      << " faces: ";

            for (const auto fi : incident_faces)
              std::cerr << fi << ' ';

            std::cerr << '\n';
          }
        }

        std::cerr << "Directed edges with more than two incident faces: "
                  << non_manifold_edges << '\n';

        assert(non_manifold_edges == 0 &&
              "An edge has more than two incident faces");
      }

      void assert_valid_halfedges() const {
        // Every halfedge should have exactly one previous halfedge
        std::vector<bool> has_prev(num_halfedges(mesh_), false);

        for (const auto h : halfedges(mesh_))
        {
          const auto h_next = next(h, mesh_);
          if (has_prev[h_next.idx()]) {
            std::cout << "found halfedge with two halfedges pointing to it: " << h_next << std::endl;
            export_halfedge_mesh_json(mesh_, "double_prev_halfedge.json");
          }
          CGAL_assertion(h_next != TriangleMesh::null_halfedge());
          CGAL_assertion(!has_prev[h_next.idx()]);
          has_prev[h_next.idx()] = true;
        }

        for (const auto h : halfedges(mesh_))
        {
          CGAL_assertion(has_prev[h.idx()] && "Halfedge has no previous halfedge");
        }
      }

      /**
       * @brief Inserts a mesh vertex into the spatial lookup structure.
       *
       * The vertex is stored in the grid cell containing its position.
      */
      void insert_mesh_vertex(vertex_descriptor vd) {
        int ix, iy, iz;

        const Point_3& p = get(points_pm_, vd);

        if (!grid_.to_grid_coordinates(p, ix, iy, iz))
            return;

        mesh_vertices_per_cell_[grid_.flat_index(ix, iy, iz)].push_back(vd);
      }

      /**
       * @brief Returns mesh vertices within a spatial neighborhood.
       *
       * @param p Query position.
       * @param radius Search radius.
       *
       * @return Mesh vertices in neighboring grid cells.
      */
      std::vector<vertex_descriptor> nearby_mesh_vertices(const Point_3& p, FT radius) const {
        std::vector<vertex_descriptor> out;

        int cx, cy, cz;

        if (!grid_.to_grid_coordinates(p, cx, cy, cz))
          return out;

        const int r = (std::max)(
          1,
          static_cast<int>(std::ceil(CGAL::to_double(radius) / CGAL::to_double(grid_.get_box_size()))));

        for (int dx = -r; dx <= r; ++dx)
          for (int dy = -r; dy <= r; ++dy)
            for (int dz = -r; dz <= r; ++dz) {
              if (!grid_.valid_coordinates(cx + dx, cy + dy, cz + dz))
                continue;

              const auto& cell = mesh_vertices_per_cell_[grid_.flat_index(cx + dx, cy + dy, cz + dz)];

              out.insert(out.end(), cell.begin(), cell.end());
            }

        return out;
      }

      /**
       * @brief Inserts a candidate into its priority queue.
       *
       * The candidate priority is clamped to the valid queue range before insertion.
      */
      void push_candidate_in_queue(Candidate cand) {
        CGAL_assertion(cand.priority >= 0);
        CGAL_assertion(cand.priority < NUM_PRIORITIES);

        candidate_queues_[cand.priority].push_back(cand);
      }

      /**
       * @brief Removes the highest-priority available candidate from the queues.
       *
       * @param cand Output candidate.
       *
       * @return `true` if a candidate was available.
      */
      bool pop_candidate_from_queue(Candidate& cand) {
        for (int p = NUM_PRIORITIES - 1; p >= 0; --p) {
          if (!candidate_queues_[p].empty()) {
            cand = candidate_queues_[p].front();
            candidate_queues_[p].pop_front();
            return true;
          }
        }

        return false;
      }

      /// @brief Returns whether every candidate queue is empty.
      bool candidate_queues_empty() const {
        for (const auto& q : candidate_queues_) {
          if (!q.empty())
            return false;
        }

        return true;
      }

      /**
       * @brief Initializes the mesh with the initial two-vertex seed.
       *
       * The seed creates a single undirected edge whose two halfedges form the
       * initial front. Candidate generation is then started from both seed vertices.
      */
      void seed_from_grid() {
        std::optional<Initial_seed> seed = surface_approximation_->get_initial_seed();
        if (!seed) {
          std::cerr << "No initial splat seed found.\n";
          return;
        }

        std::cout << "Adding initial seed vertices to mesh...\n";
        // v0_ = mesh_.add_vertex(seed->first);
        vertex_descriptor v0_ = add_vertex(mesh_);
        put(points_pm_, v0_, seed->first);
        put(normals_pm_, v0_, seed->normal);

        vertex_descriptor v1_ = add_vertex(mesh_);
        put(points_pm_, v1_, seed->second);
        put(normals_pm_, v1_, seed->normal);

        // Create the base edge between the two seed vertices.
        edge_descriptor e = add_edge(mesh_);
        halfedge_descriptor h = halfedge(e, mesh_);
        halfedge_descriptor ho = opposite(h, mesh_);

        set_target(h,  v1_, mesh_); // v0 -> v1
        set_target(ho, v0_, mesh_); // v1 -> v0
        set_halfedge(v0_, ho, mesh_);
        set_halfedge(v1_, h, mesh_);

        set_next(h, ho, mesh_);
        set_next(ho, h, mesh_);

        // Fully assert all properties of the initial vertices and edge
        CGAL_assertion(is_valid_edge_descriptor(e, mesh_));
        CGAL_assertion(is_valid_halfedge_descriptor(h, mesh_));
        CGAL_assertion(is_valid_halfedge_descriptor(ho, mesh_));
        CGAL_assertion(source(h, mesh_) == v0_);
        CGAL_assertion(source(ho, mesh_) == v1_);
        CGAL_assertion(target(h, mesh_) == v1_);
        CGAL_assertion(target(ho, mesh_) == v0_);
        CGAL_assertion(next(h, mesh_) == ho);
        CGAL_assertion(next(ho, mesh_) == h);
        CGAL_assertion(halfedge(v0_, mesh_) == ho);
        CGAL_assertion(halfedge(v1_, mesh_) == h);
        CGAL_assertion(is_valid_vertex_descriptor(v0_, mesh_));
        CGAL_assertion(is_valid_vertex_descriptor(v1_, mesh_));

        

        std::cout << "Found initial seed: " << seed->first << ", " << seed->second << ", " << seed->normal << std::endl;

        seeded_ = true;
        std::cout << "Mesh seeded with initial vertices. Starting growth...\n";

        insert_mesh_vertex(v0_);
        insert_mesh_vertex(v1_);

        push_candidates_from_vertex(v0_);
        push_candidates_from_vertex(v1_);
      }

      /**
       * @brief Generates new candidates around a newly inserted mesh vertex.
       *
       * Nearby mesh vertices are paired with the supplied vertex and candidate
       * splat intersections are generated for each pair.
       *
       * @param nv Newly inserted mesh vertex.
      */
      void push_candidates_from_vertex(vertex_descriptor nv) {
        const Point_3 p = get(points_pm_, nv);

        auto nearby = nearby_mesh_vertices(p, FT(8.0) * grid_.get_box_size());

        for (vertex_descriptor other_v : nearby) {

          if (other_v == nv) {
            continue;
          }

          std::vector<Candidate> candidates = generate_candidates_from_parents(nv, other_v);

          for (const Candidate& c : candidates) {
            push_candidate_in_queue(c);
          }
        }
      }

      /**
       * @brief Generates candidate vertices from a pair of parent vertices.
       *
       * A construction circle derived from the parent pair is intersected with
       * nearby splats. Each valid intersection becomes a reconstruction candidate.
       *
       * @param parent_a First parent vertex.
       * @param parent_b Second parent vertex.
       *
       * @return Candidate vertices generated from the parent pair.
      */
      std::vector<Candidate> generate_candidates_from_parents(vertex_descriptor parent_a,
                                                              vertex_descriptor parent_b) const {
        std::vector<Candidate> candidates;

        if (parent_a == parent_b) {
          return candidates;
        }

        const Point_3 pa = get(points_pm_, parent_a);
        const Point_3 pb = get(points_pm_, parent_b);

        const FT d = grid_.get_box_size();
        const FT d2 = CGAL::squared_distance(pa, pb);

        if (d2 < d * d - FT(1e-6)) {
          return candidates;
        }

        Point_3 circle_center;
        Vector_3 circle_normal;
        FT circle_radius;

        if (!CGAL::internal::Splat_surface_reconstruction_3::find_circle_between_points<Point_3, Vector_3, FT>(pa, pb, d, circle_center, circle_normal, circle_radius)) {
          return candidates;
        }

        // Now the circle is defined by circle_center, circle_normal, and circle_radius.

        // Find all intersections with the surface approximation.
        const std::vector<std::pair<Point_3, Vector_3>> hits =
          surface_approximation_->intersect_circle(circle_center, circle_normal, circle_radius);

        for (const auto& [p, normal] : hits) {
          // Changed: add this check
          // Don't add candidates that are outside the grid bounds
          if (!grid_.is_point_inside_grid(p))
            continue;

          // Changed: add this check
          // Don't add vertices that cannot find a proper box normal
          if (surface_approximation_->box_normal(p) == CGAL::NULL_VECTOR)
            continue;

          Candidate c(p, normal, parent_a, parent_b);
          c.priority = compute_priority(c);
          candidates.push_back(c);
        }

        return candidates;
      }

      /**
       * @brief Tests whether a mesh vertex already exists near a candidate position.
       *
       * @param p Candidate position.
       *
       * @return `true` when an existing mesh vertex is closer than the reconstruction
       *         tolerance.
      */
      bool has_vertex_near(const Point_3& p) const {
        const FT tol = grid_.get_box_size();
        const FT tol2 = tol * tol;
        const FT epsilon = FT(1e-6);

        auto nearby = nearby_mesh_vertices(p, FT(2.0) * grid_.get_box_size());
        for (vertex_descriptor vd : nearby) {
          if (CGAL::squared_distance(get(points_pm_, vd), p) < tol2 - epsilon) {
            return true;
          }
        }
        return false;
      }

      /**
       * @brief Creates an oriented halfedge between two vertices.
       *
       * The new edge is initialized with the requested orientation `va -> vb`.
       * The opposite halfedge is initialized with orientation `vb -> va`.
       *
       * Vertex halfedge pointers are initialized when the corresponding vertex does
       * not already have a halfedge.
       *
       * @param va Source vertex.
       * @param vb Target vertex.
       *
       * @return The halfedge oriented from `va` to `vb`.
      */
      halfedge_descriptor create_halfedge(vertex_descriptor va,
                                          vertex_descriptor vb)
      {
        edge_descriptor e = add_edge(mesh_);
        halfedge_descriptor h_vavb = halfedge(e, mesh_);
        halfedge_descriptor h_vbva = opposite(h_vavb, mesh_);

        // h_vavb : va -> vb
        set_target(h_vavb, vb, mesh_);
        set_target(h_vbva, va, mesh_);

        if (halfedge(va, mesh_) == boost::graph_traits<TriangleMesh>::null_halfedge()) {
          set_halfedge(va, h_vbva, mesh_);
        }
        if (halfedge(vb, mesh_) == boost::graph_traits<TriangleMesh>::null_halfedge()) {
          set_halfedge(vb, h_vavb, mesh_);
        }

        return h_vavb;
      }

      // DEBUG: (TODO: remove)

      template <typename Point>
      bool export_halfedge_mesh_json(const CGAL::Surface_mesh<Point>& mesh, const std::string& filename) const
      {
        std::cout << "Start exporting halfedge mesh to " << filename << std::endl;
        using Mesh = CGAL::Surface_mesh<Point>;

        std::ofstream out(filename);
        if (!out) {
          std::cerr << "Could not open " << filename << '\n';
          return false;
        }

        out << "{\n";

        // --------------------------------------------------
        // Vertices: descriptor ID and 3D position
        // --------------------------------------------------
        out << "  \"vertices\": [\n";

        bool first = true;
        for (const auto v : mesh.vertices()) {
          if (!first) out << ",\n";
          first = false;

          const Point& p = mesh.point(v);
          const Vector_3 normal = normals_pm_[v];
          const Vector_3 box_normal = surface_approximation_->box_normal(p);
          Vector_3 cross_product_normal = Vector_3(0, 0, 0);

          bool found_loop = false;

          // Find the longest attached loop
          size_t longest_loop_length = 0;
          halfedge_descriptor longest_loop_halfedge;
          for (const auto h : halfedges_around_target(v, mesh)) {
            size_t loop_length = 0;
            auto h_loop = h;
            bool found_loop_in_this_iteration = false;
            for (size_t i = 0; i < 100; i++) {
              h_loop = next(h_loop, mesh);
              ++loop_length;
              if (h_loop == h) {
                found_loop_in_this_iteration = true;
                found_loop = true;
                break;
              }
            }
            if (loop_length > longest_loop_length) {
              longest_loop_length = loop_length;
              longest_loop_halfedge = h;
            }
          }

          if (found_loop) {
            // Compute the normal based on the cross product of the edges of the longest loop
            Point_3 prev_v = points_pm_[source(longest_loop_halfedge, mesh)];
            Point_3 next_v = points_pm_[target(next(longest_loop_halfedge, mesh), mesh)];
            Vector_3 to_prev = prev_v - p;
            Vector_3 to_next = next_v - p;
            cross_product_normal = CGAL::cross_product(to_next, to_prev);
            if (cross_product_normal.squared_length() == 0) {
              cross_product_normal = Vector_3(0, 0, 0);
            }
            else {
              cross_product_normal = cross_product_normal / std::sqrt(cross_product_normal.squared_length());
            }
          }

          out << "    {\"id\": " << v.idx()
              << ", \"point\": ["
              << CGAL::to_double(p.x()) << ", "
              << CGAL::to_double(p.y()) << ", "
              << CGAL::to_double(p.z()) << "], "
              << "\"normal\": ["
              << CGAL::to_double(normal.x()) << ", "
              << CGAL::to_double(normal.y()) << ", "
              << CGAL::to_double(normal.z()) << "], "
              << "\"box_normal\": ["
              << CGAL::to_double(box_normal.x()) << ", "
              << CGAL::to_double(box_normal.y()) << ", "
              << CGAL::to_double(box_normal.z()) << "], "
              << "\"cross_product_normal\": ["
              << CGAL::to_double(cross_product_normal.x()) << ", "
              << CGAL::to_double(cross_product_normal.y()) << ", "
              << CGAL::to_double(cross_product_normal.z()) << "]}";
        }

        out << "\n  ],\n";

        // --------------------------------------------------
        // Faces: ordered vertex IDs around each face
        // --------------------------------------------------
        out << "  \"faces\": [\n";

        first = true;
        for (const auto f : mesh.faces()) {
          if (!first) out << ",\n";
          first = false;

          out << "    {\"id\": " << f.idx() << ", \"vertices\": [";

          const auto h0 = halfedge(f, mesh);
          auto h = h0;
          bool first_vertex = true;

          do {
            if (!first_vertex) out << ", ";
            first_vertex = false;

            out << source(h, mesh).idx();
            h = next(h, mesh);
          } while (h != h0);

          out << "]}";
        }

        out << "\n  ],\n";

        // --------------------------------------------------
        // Splats
        // --------------------------------------------------
        // out << "  \"splats\": [\n";

        // const std::vector<FT> splat_sizes = grid_.splat_sizes();
        // const std::vector<Vector_3> splat_normals = grid_.normals();
        // const std::vector<Point> splat_positions = grid_.points();

        // for (std::size_t i = 0; i < splat_sizes.size(); ++i) {
        //   if (i > 0) out << ",\n";

        //   out << "    {\"size\": " << CGAL::to_double(splat_sizes[i])
        //       << ", \"normal\": ["
        //       << CGAL::to_double(splat_normals[i].x()) << ", "
        //       << CGAL::to_double(splat_normals[i].y()) << ", "
        //       << CGAL::to_double(splat_normals[i].z()) << "], "
        //       << "\"position\": ["
        //       << CGAL::to_double(splat_positions[i].x()) << ", "
        //       << CGAL::to_double(splat_positions[i].y()) << ", "
        //       << CGAL::to_double(splat_positions[i].z()) << "]}";
        // }
        // out << "\n  ],\n";

        // --------------------------------------------------
        // Feature segments
        // --------------------------------------------------
        out << "  \"feature_segments\": [\n";

        for (std::size_t i = 0; i < feature_segments_.size(); ++i) {
          if (i > 0) out << ",\n";

          out << "    {\"point_1\": ["
              << CGAL::to_double(feature_segments_[i].point_1.x()) << ", "
              << CGAL::to_double(feature_segments_[i].point_1.y()) << ", "
              << CGAL::to_double(feature_segments_[i].point_1.z()) << "], "
              << "\"point_2\": ["
              << CGAL::to_double(feature_segments_[i].point_2.x()) << ", "
              << CGAL::to_double(feature_segments_[i].point_2.y()) << ", "
              << CGAL::to_double(feature_segments_[i].point_2.z()) << "]}";
        }
        out << "\n  ],\n";



        // --------------------------------------------------
        // Half-edges: complete directed connectivity
        // --------------------------------------------------
        out << "  \"halfedges\": [\n";

        first = true;
        for (const auto h : mesh.halfedges()) {
          if (!first) out << ",\n";
          first = false;

          const auto opp = opposite(h, mesh);
          const auto nxt = next(h, mesh);
          const auto f = face(h, mesh);

          // -1 represents the null face (a boundary half-edge).
          const long long face_id =
              (f == Mesh::null_face())
                  ? -1LL
                  : static_cast<long long>(f.idx());

          out << "    {\"id\": " << h.idx()
              << ", \"source\": " << source(h, mesh).idx()
              << ", \"target\": " << target(h, mesh).idx()
              << ", \"opposite\": " << opp.idx()
              << ", \"next\": " << nxt.idx()
              << ", \"face\": " << face_id
              << "}";
        }

        out << "\n  ]\n";
        out << "}\n";

        if (!out) {
          std::cerr << "Failed while writing " << filename << '\n';
          return false;
        }

        std::cout << "Exported mesh to " << filename << '\n';
        return true;
      }

      /**
       * @brief Inserts a new vertex into the local halfedge front.
       *
       * The two parent vertices define an existing front interval. The new vertex
       * is connected to both parents and the surrounding `next()` relationships are
       * rewired so that the old boundary interval is replaced by two new boundary
       * intervals.
       *
       * The local halfedge ordering around each parent is computed using the parent
       * box normals. The selected halfedges are first checked for compatibility
       * before any topology is modified.
       *
       * @param v0 First parent vertex.
       * @param v1 Second parent vertex.
       * @param nv Newly inserted vertex.
       * @param normal Normal associated with the new candidate.
       *
       * @return `true` if the local graph can be updated consistently, `false`
       *         otherwise.
      */
      bool build_graph(vertex_descriptor v0,
                       vertex_descriptor v1,
                       vertex_descriptor nv,
                       const Vector_3& normal) {

        const auto null_h = boost::graph_traits<TriangleMesh>::null_halfedge();
        const double two_pi = 2.0 * std::acos(-1.0);

        const Point_3 p0 = get(points_pm_, v0);
        const Point_3 p1 = get(points_pm_, v1);
        const Point_3 pn = get(points_pm_, nv);

        // The three vertices must be distinct.
        CGAL_assertion(v0 != v1);
        CGAL_assertion(v0 != nv);
        CGAL_assertion(v1 != nv);

        // The candidate triangle must not be degenerate.
        const Vector_3 e01 = p1 - p0;
        const Vector_3 e0n = pn - p0;
        CGAL_assertion(CGAL::cross_product(e01, e0n) != CGAL::NULL_VECTOR);

        // All three points should be distinct.
        CGAL_assertion(p0 != p1);
        CGAL_assertion(p0 != pn);
        CGAL_assertion(p1 != pn);

        auto normalize_vec = [](const Vector_3& n) -> Vector_3 {
          const double len2 = CGAL::to_double(n.squared_length());
          if (len2 <= 0.0) {
            return CGAL::NULL_VECTOR;
          }
          return n / std::sqrt(len2);
        };

        // Normalize the box normals before using them to construct tangent frames.
        Vector_3 n0 = normalize_vec(surface_approximation_->box_normal(p0));
        Vector_3 n1 = normalize_vec(surface_approximation_->box_normal(p1));
        //Vector_3 n = normalize_vec(n0+n1);
        Vector_3 n2 = normalize_vec(surface_approximation_->box_normal(pn));

        // // This normal can be very wrong sometimes
        // n0 = normals_pm_[v0];
        // n1 = normals_pm_[v1];

        //Vector_3 box_normal = normalize_vec(n0+n1+n2);

        // // Thus compute normal from the triangle (x impossible, v0 and v1 can be in any order)
        // Vector_3 n = normalize_vec(CGAL::cross_product(p1 - p0, pn - p0));
        // // Check that it agrees with the box normals.
        // if (n * box_normal < FT(0)) {
        //   n = -n;
        // }

        // n0 = n;
        // n1 = n;

        // CGAL_assertion(n0 != CGAL::NULL_VECTOR);
        // CGAL_assertion(n1 != CGAL::NULL_VECTOR);

        // if (n0 == CGAL::NULL_VECTOR || n1 == CGAL::NULL_VECTOR) {
        //   return false;
        // }

        // Construct an oriented tangent frame. The sign of v is adjusted so that
        // cross(u,v) agrees with the supplied normal.
        auto make_frame = [&](const Vector_3& n) -> std::pair<Vector_3, Vector_3> {
          std::vector<Vector_3> frame = CGAL::internal::Splat_surface_reconstruction_3::compute_local_tangent_frame<Point_3, Vector_3, FT>(n);
          Vector_3 u = frame[0];
          Vector_3 v = frame[1];

          if (u == CGAL::NULL_VECTOR || v == CGAL::NULL_VECTOR) {
            return {CGAL::NULL_VECTOR, CGAL::NULL_VECTOR};
          }

          if (CGAL::cross_product(u, v) * n < FT(0)) {
            v = -v;
          }

          CGAL_assertion(u != CGAL::NULL_VECTOR && v != CGAL::NULL_VECTOR);

          CGAL_assertion(u.squared_length() > FT(0));
          CGAL_assertion(v.squared_length() > FT(0));

          // Their cross product should agree with the normal
          CGAL_assertion(CGAL::cross_product(u, v) * n >= FT(0));

          // Check that u and v are perpendicular to the normal n.
          CGAL_assertion(
            std::abs(CGAL::to_double(u * n)) < 1e-8
          );
          CGAL_assertion(
            std::abs(CGAL::to_double(v * n)) < 1e-8
          );
          // Check that u and v are perpendicular to each other.
          CGAL_assertion(
            std::abs(CGAL::to_double(u * v)) < 1e-8
          );

          return {u, v};
        };

        // Return the polar angle of q around origin in the tangent frame.
        auto angle_of = [&](const Point_3& origin,
                            const Point_3& q,
                            const Vector_3& u,
                            const Vector_3& v) -> double {
          const Vector_3 d = q - origin;
          double a = std::atan2(CGAL::to_double(d * v), CGAL::to_double(d * u));
          if (a < 0.0) {
            a += two_pi;
          }
          if (a >= two_pi) {
            a = 0.0;
          }
          return a;
        };

        // Collect all incident halfedges around a parent and sort them by their
        // projected angle. These halfedges define the local front ordering.
        auto sorted_halfedges_around_vertex =
          [&](vertex_descriptor v, const Vector_3& local_normal)
            -> std::vector<std::pair<double, halfedge_descriptor>> {
          std::vector<std::pair<double, halfedge_descriptor>> out;
          const Point_3 pv = get(points_pm_, v);

          const auto [u, vaxis] = make_frame(local_normal);
          if (u == CGAL::NULL_VECTOR || vaxis == CGAL::NULL_VECTOR) {
            return out;
          }

          // If it has no outgoing halfedges, return an empty array
          if (out_degree(v, mesh_) == 0) {
            return out;
          }
          for (halfedge_descriptor h : halfedges_around_source(v, mesh_)) {
            CGAL_assertion(source(h, mesh_) == v);
            const vertex_descriptor t = target(h, mesh_);
            if (t == v) {
              continue;
            }
            out.emplace_back(angle_of(pv, get(points_pm_, t), u, vaxis), h);
          }

          std::sort(out.begin(), out.end(),
            [](const auto& a, const auto& b) { return a.first < b.first; });

          for (std::size_t i = 1; i < out.size(); ++i) {
            CGAL_assertion(out[i - 1].first <= out[i].first);
          }

          for (const auto& [angle, he] : out) {
            CGAL_assertion(he != null_h);
            CGAL_assertion(source(he, mesh_) == v);
            CGAL_assertion(angle >= 0.0 && angle < two_pi);
          }

          return out;
        };

        // Select the two front halfedges immediately clockwise and counter-clockwise
        // from the direction toward the new vertex.
        auto choose_closest_ccw_halfedges =
          [&](vertex_descriptor v,
              const Point_3& new_point,
              const Vector_3& local_normal)
            -> std::pair<halfedge_descriptor, halfedge_descriptor> {
          const Point_3 pv = get(points_pm_, v);
          const auto [u, vaxis] = make_frame(local_normal);
          if (u == CGAL::NULL_VECTOR || vaxis == CGAL::NULL_VECTOR) {
            return {null_h, null_h};
          }

          const double theta_new = angle_of(pv, new_point, u, vaxis);

          std::vector<std::pair<double, halfedge_descriptor>> sorted =
            sorted_halfedges_around_vertex(v, local_normal);

          if (sorted.empty()) {
            return {null_h, null_h};
          }

          auto it = std::lower_bound(
            sorted.begin(), sorted.end(), theta_new,
            [](const std::pair<double, halfedge_descriptor>& a, double theta) {
              return a.first < theta;
            });

          const auto& ccw = (it == sorted.end())   ? sorted.front() : *it;
          const auto& cw  = (it == sorted.begin()) ? sorted.back()  : *(it - 1);

          // const auto& ccw_next = (std::next(it) == sorted.end()) ? sorted.front() : *std::next(it);
          // const auto& cw_next = (it == sorted.begin() + 1) ? sorted.back() : *(it - 2);

          return {cw.second, ccw.second};
        };

        // Compute the local ordering using parent normals.
        auto [h1, h2] = choose_closest_ccw_halfedges(v0, pn, n0);
        auto [g2, g1] = choose_closest_ccw_halfedges(v1, pn, n1);

        // std::cout << "\n[BUILD_GRAPH ITERATION]" << std::endl;
        // std::cout << "v0: " << v0 << " v1: " << v1 << std::endl;
        // std::cout << "h1: " << h1 << " h2: " << h2 << " g1: " << g1 << " g2: " << g2 << std::endl;

        // g2 = h_v0v1;
        // h2 = h_v1v0;

        // // Then we can take h1 = opposite(next(h2))
        // h1 = opposite(next(h2, mesh_), mesh_);
        // // And g1 = previous(h2)

        // bool found_prev_h1 = false;
        // for (const auto& hi : mesh_.halfedges()) {
        //   if (next(hi, mesh_) == h2) {
        //     g1 = hi;
        //     found_prev_h1 = true;
        //     break;
        //   }
        // }
        
        // CGAL_assertion(found_prev_h1);

        // Ensure that the selected halfedges are oriented from the parent vertex to
        // the other vertex. If not, flip them to their opposite halfedge.

        // if (nv.idx() == 433) {
        //   std::cout << "nv.idx() == 433" << std::endl;
        //   std::cout << "v0: " << v0 << " v1: " << v1 << std::endl;
        //   std::cout << "h1: " << h1 << " h2: " << h2 << " g1: " << g1 << " g2: " << g2 << std::endl;
        //   export_halfedge_mesh_json(mesh_, "pre_433.json");
        // }

        const auto null_f = boost::graph_traits<TriangleMesh>::null_face();

        if (h1 != null_h) {
          if (target(h1, mesh_) != v0) h1 = opposite(h1, mesh_);
          CGAL_assertion(target(h1, mesh_) == v0);

          if (target(h2, mesh_) != v0) h2 = opposite(h2, mesh_);
          CGAL_assertion(target(h2, mesh_) == v0);
        
          CGAL_assertion(target(h1, mesh_) == v0);
          CGAL_assertion(target(h2, mesh_) == v0);
          CGAL_assertion(source(h1, mesh_) != v0);
          CGAL_assertion(source(h2, mesh_) != v0);

          CGAL_assertion(face(h1, mesh_) == null_f);
          CGAL_assertion(face(h2, mesh_) == null_f);
        }

        if (g1 != null_h) {
          if (target(g1, mesh_) != v1) g1 = opposite(g1, mesh_);
          CGAL_assertion(target(g1, mesh_) == v1);

          if (target(g2, mesh_) != v1) g2 = opposite(g2, mesh_);
          CGAL_assertion(target(g2, mesh_) == v1);
        
          CGAL_assertion(target(g1, mesh_) == v1);
          CGAL_assertion(target(g2, mesh_) == v1);
          CGAL_assertion(source(g1, mesh_) != v1);
          CGAL_assertion(source(g2, mesh_) != v1);

          CGAL_assertion(face(g1, mesh_) == null_f);
          CGAL_assertion(face(g2, mesh_) == null_f);
        }

        // std::cout << "after ensuring halfedge orientations:" << std::endl;
        // std::cout << "h1: " << h1 << " h2: " << h2 << " g1: " << g1 << " g2: " << g2 << std::endl;

        // If there exists an edge between v0 and v1, it should be represented by g2 and h2.
        if (halfedge(v0, v1, mesh_).second) {

          // std::cout << "halfedge between v0 and v1 exists." << std::endl;

          // std::cout << "actual: " << halfedge(v0, v1, mesh_).first << " " << halfedge(v1, v0, mesh_).first << std::endl;
          // std::cout << "g2: " << g2 << " h2: " << h2 << std::endl;
          // std::cout << "g1: " << g1 << " h1: " << h1 << std::endl;

          // It's either the pair g1, h1 or the pair g2, h2 that represents the edge between v0 and v1.
          if (g2 == halfedge(v0, v1, mesh_).first && h2 == halfedge(v1, v0, mesh_).first) {
            CGAL_assertion(g2 == opposite(h2, mesh_));
            CGAL_assertion(h2 == opposite(g2, mesh_));

            if (next(g1, mesh_) != h2) {
              std::cout << "Debugging halfedge mesh: next(h1, mesh_) != g2" << std::endl;
              std::cout << "h1: " << h1 << " next(h1, mesh_): " << next(h1, mesh_) << std::endl;
              std::cout << "g2: " << g2 << " next(g2, mesh_): " << next(g2, mesh_) << std::endl;
              std::cout << "g1: " << g1 << " next(g1, mesh_): " << next(g1, mesh_) << std::endl;
              std::cout << "h2: " << h2 << " next(h2, mesh_): " << next(h2, mesh_) << std::endl;
              export_halfedge_mesh_json(
                mesh_,
                "mesh_debug.json"
              );
            }

            // There must be continuity between g1 and h2
            CGAL_assertion(next(g1, mesh_) == h2);
          }
          else if (g1 == halfedge(v0, v1, mesh_).first && h1 == halfedge(v1, v0, mesh_).first) {
            CGAL_assertion(g1 == opposite(h1, mesh_));
            CGAL_assertion(h1 == opposite(g1, mesh_));

            if (next(h2, mesh_) != g1) {

              std::cout << "Debugging halfedge mesh: next(g2, mesh_) != h1" << std::endl;
              std::cout << "v0: " << v0 << " v1: " << v1 << std::endl;
              std::cout << "g2: " << g2 << " next(g2, mesh_): " << next(g2, mesh_) << std::endl;
              std::cout << "h1: " << h1 << " next(h1, mesh_): " << next(h1, mesh_) << std::endl;
              std::cout << "g1: " << g1 << " next(g1, mesh_): " << next(g1, mesh_) << std::endl;
              std::cout << "h2: " << h2 << " next(h2, mesh_): " << next(h2, mesh_) << std::endl;
              export_halfedge_mesh_json(
                mesh_,
                "mesh_debug.json"
              );
            }

            // There must be continuity between h2 and g1
            CGAL_assertion(next(h2, mesh_) == g1);
          }

          //   std::cout << "Debugging halfedge mesh: unexpected configuration" << std::endl;
          //   std::cout << "h1: " << h1 << " h2: " << h2 << " g1: " << g1 << " g2: " << g2 << std::endl;
          //   std::cout << "existing edge between v0 and v1: " << halfedge(v0, v1, mesh_).first << std::endl;
          //   std::cout << "existing edge between v1 and v0: " << halfedge(v1, v0, mesh_).first << std::endl;
          //   std::cout << "v0: " << v0 << std::endl;
          //   std::cout << "v1: " << v1 << std::endl;
          //   export_halfedge_mesh_json(
          //     mesh_,
          //     "mesh_debug.json"
          //   );
          //   CGAL_assertion(false);
          // }

          // std::cout << "halfedge between v0 and v1 exists." << std::endl;

          // h2 = halfedge(v1, v0, mesh_).first;
          // g2 = halfedge(v0, v1, mesh_).first;
          // h1 = opposite(next(h2, mesh_), mesh_);
          // g1 = prev(h2, mesh_);

          // // g2 should be the edge from v0 to v1
          // CGAL_assertion(source(g2, mesh_) == v0 && target(g2, mesh_) == v1);
          // // h2 should be the edge from v1 to v0
          // CGAL_assertion(source(h2, mesh_) == v1 && target(h2, mesh_) == v0);

          // CGAL_assertion(g2 == halfedge(v0, v1, mesh_).first);
          // CGAL_assertion(h2 == halfedge(v1, v0, mesh_).first);

          // // g2 should be opposite of h2
          // CGAL_assertion(g2 == opposite(h2, mesh_));
          // CGAL_assertion(h2 == opposite(g2, mesh_));

          // // Next of g1 should be h2
          // CGAL_assertion(next(g1, mesh_) == h2);
        }

        // Create the two new undirected edges connecting the new vertex to its parents.
        halfedge_descriptor h_v0nv = create_halfedge(v0, nv);
        halfedge_descriptor h_nvv0 = opposite(h_v0nv, mesh_);
        halfedge_descriptor h_v1nv = create_halfedge(v1, nv);
        halfedge_descriptor h_nvv1 = opposite(h_v1nv, mesh_);

        if (h1 == null_h) {
          h1 = h_nvv0;
          h2 = h_nvv0;
        }
        if (g1 == null_h) {
          g1 = h_nvv1;
          g2 = h_nvv1;
        }

        // Check that the halfedge cycles around each parent vertex don't have repeated vertices.

        // This ensures that the local front is a simple cycle.
        // We limit the number of steps to avoid infinite loops.
        // std::unordered_set<vertex_descriptor> visited;
        // halfedge_descriptor h = opposite(h1, mesh_);
        // bool reached_target = false;
        // for (std::size_t steps = 0; steps <= 2000; ++steps) {

        //   CGAL_assertion(h != null_h);
        //   CGAL_assertion(opposite(h, mesh_) != null_h);
        //   CGAL_assertion(next(h, mesh_) != null_h);

        //   if (!visited.insert(source(h, mesh_)).second)
        //     return false;

        //   if (h == g1) {
        //     reached_target = true;
        //     break;
        //   }

        //   h = next(h, mesh_);
        // }
        // CGAL_assertion(reached_target);

        // visited.clear();
        // h = opposite(g2, mesh_);
        // reached_target = false;
        // for (std::size_t steps = 0; steps <= 2000; ++steps) {

        //   CGAL_assertion(h != null_h);
        //   CGAL_assertion(opposite(h, mesh_) != null_h);
        //   CGAL_assertion(next(h, mesh_) != null_h);

        //   if (!visited.insert(source(h, mesh_)).second)
        //     return false;

        //   if (h == h2) {
        //     reached_target = true;
        //     break;
        //   }

        //   h = next(h, mesh_);
        // }

        // CGAL_assertion(reached_target);

        // First boundary update:
        // g1 -> v1 -> nv -> v0 -> opposite(h1)
        set_next(g1, h_v1nv, mesh_);
        set_next(h_v1nv, h_nvv0, mesh_);
        set_next(h_nvv0, opposite(h1, mesh_), mesh_);

        // Second boundary update:
        // h2 -> v0 -> nv -> v1 -> opposite(g2)
        set_next(h2, h_v0nv, mesh_);
        set_next(h_v0nv, h_nvv1, mesh_);
        set_next(h_nvv1, opposite(g2, mesh_), mesh_);

        // Verify all six new successor relationships.
        CGAL_assertion(next(g1, mesh_) == h_v1nv);
        CGAL_assertion(next(h_v1nv, mesh_) == h_nvv0);
        CGAL_assertion(next(h_nvv0, mesh_) == opposite(h1, mesh_));

        CGAL_assertion(next(h2, mesh_) == h_v0nv);
        CGAL_assertion(next(h_v0nv, mesh_) == h_nvv1);
        CGAL_assertion(next(h_nvv1, mesh_) == opposite(g2, mesh_));

        // And verify the source and target vertices of the new halfedges.
        CGAL_assertion(source(h_v0nv, mesh_) == v0);
        CGAL_assertion(target(h_v0nv, mesh_) == nv);

        CGAL_assertion(source(h_nvv0, mesh_) == nv);
        CGAL_assertion(target(h_nvv0, mesh_) == v0);

        CGAL_assertion(source(h_v1nv, mesh_) == v1);
        CGAL_assertion(target(h_v1nv, mesh_) == nv);

        CGAL_assertion(source(h_nvv1, mesh_) == nv);
        CGAL_assertion(target(h_nvv1, mesh_) == v1);

        // if (nv.idx() == 433) {
        //   std::cout << "nv.idx() == 433" << std::endl;
        //   std::cout << "v0: " << v0 << " v1: " << v1 << std::endl;
        //   std::cout << "h1: " << h1 << " h2: " << h2 << " g1: " << g1 << " g2: " << g2 << std::endl;
        //   export_halfedge_mesh_json(mesh_, "post_433.json");
        // }

        // std::cout << "g1: " << g1.idx() << ", h1: " << h1.idx() << ", g2: " << g2.idx() << ", h2: " << h2.idx() << std::endl;
        // std::cout << "h_v0nv: " << h_v0nv.idx() << ", h_nvv0: " << h_nvv0.idx() << ", h_v1nv: " << h_v1nv.idx() << ", h_nvv1: " << h_nvv1.idx() << std::endl;
        // std::cout << "h_v0nv: " << source(h_v0nv, mesh_).idx() << " -> " << target(h_v0nv, mesh_).idx() << ", \n";
        // std::cout << "h_nvv0: " << source(h_nvv0, mesh_).idx() << " -> " << target(h_nvv0, mesh_).idx() << ", \n";
        // std::cout << "h_v1nv: " << source(h_v1nv, mesh_).idx() << " -> " << target(h_v1nv, mesh_).idx() << ", \n";
        // std::cout << "h_nvv1: " << source(h_nvv1, mesh_).idx() << " -> " << target(h_nvv1, mesh_).idx() << std::endl;
        // std::cout << "next(h_v0nv): " << next(h_v0nv, mesh_).idx() << ", \n";
        // std::cout << "next(h_nvv0): " << next(h_nvv0, mesh_).idx() << ", \n";
        // std::cout << "next(h_v1nv): " << next(h_v1nv, mesh_).idx() << ", \n";
        // std::cout << "next(h_nvv1): " << next(h_nvv1, mesh_).idx() << std::endl;
        // std::cout << "opposite(h_v0nv): " << opposite(h_v0nv, mesh_).idx() << ", \n";
        // std::cout << "opposite(h_nvv0): " << opposite(h_nvv0, mesh_).idx() << ", \n";
        // std::cout << "opposite(h_v1nv): " << opposite(h_v1nv, mesh_).idx() << ", \n";
        // std::cout << "opposite(h_nvv1): " << opposite(h_nvv1, mesh_).idx() << std::endl;


        return true;
      }

      /**
       * @brief Converts closed `next()` cycles into polygonal faces.
       *
       * Every unassigned halfedge is used as the starting point of a boundary walk.
       * Closed cycles containing at least three halfedges are passed to the ear
       * clipping stage.
       *
       * Cycles that do not close within the safety limit are ignored.
      */
      void fill_faces_from_next_cycles()
      {
        const auto null_h = boost::graph_traits<TriangleMesh>::null_halfedge();
        const auto null_f = boost::graph_traits<TriangleMesh>::null_face();

        // Walk the current `next()` graph starting from an unassigned halfedge.
        for (halfedge_descriptor start : halfedges(mesh_)) {
          CGAL_assertion(start != null_h);

          // Only keep halfedges that are not yet assigned to a face.
          if (face(start, mesh_) != null_f) {
            continue;
          }

          std::vector<halfedge_descriptor> cycle;
          cycle.reserve(16);

          halfedge_descriptor h = start;
          const std::size_t max_steps = (std::min)((std::size_t)num_halfedges(mesh_), std::size_t(1000)); // If cycle is too long, it is considered a hole.
          bool closed = false;

          for (std::size_t i = 0; i < max_steps; ++i)
          {
            CGAL_assertion(h != null_h);

            if (face(h, mesh_) != null_f) {
              std::cout << "Encountered a halfedge that is already assigned to a face." << std::endl;
            }
            CGAL_assertion(face(h, mesh_) == null_f);

            cycle.push_back(h);

            h = next(h, mesh_);
            CGAL_assertion(h != null_h);

            if (h == start) {
              closed = true;
              break;
            }
          }

          CGAL_assertion(cycle.size() >= 3);

          // Only closed cycles with at least three halfedges can represent polygonal faces
          if (!closed || cycle.size() < 3) {
            continue;
          }

          // This assertion doesn't work as we cannot accurately determine whether an existing edge is part of another face or not.
          // // Verify that there are no preexisting halfedges between any non-consecutive vertices in the cycle.
          // for (std::size_t i = 0; i < cycle.size(); ++i) {
          //   const auto u = source(cycle[i], mesh_);

          //   for (std::size_t j = 0; j < cycle.size(); ++j) {
          //     // The first and last vertices are consecutive in a closed cycle.
          //     std::size_t i_prev = (i + cycle.size() - 1) % cycle.size();
          //     std::size_t i_next = (i + 1) % cycle.size();

          //     if (j == i || j == i_prev || j == i_next) {
          //       continue;
          //     }

          //     const auto v_prev = source(cycle[i_prev], mesh_);
          //     const auto v = source(cycle[j], mesh_);
          //     const auto v_next = source(cycle[i_next], mesh_);


          //     const auto [h_uv, found_uv] = halfedge(u, v, mesh_);
          //     const auto [h_vu, found_vu] = halfedge(v, u, mesh_);

          //     if (found_uv || found_vu) {
          //       Vector_3 from_prev = points_pm_[v_prev] - points_pm_[v];
          //       Vector_3 to_next = points_pm_[v] - points_pm_[v_next];
          //       Vector_3 to_u = points_pm_[u] - points_pm_[v];
          //       const Vector_3 normal = normals_pm_[v];

          //       to_next = to_next / std::sqrt(CGAL::scalar_product(to_next, to_next));
          //       to_u = to_u / std::sqrt(CGAL::scalar_product(to_u, to_u));
          //       from_prev = from_prev / std::sqrt(CGAL::scalar_product(from_prev, from_prev));

          //       Vector_3 right_next = CGAL::cross_product(to_next, normal);
          //       right_next = right_next / std::sqrt(CGAL::scalar_product(right_next, right_next));

          //       Vector_3 right_prev = CGAL::cross_product(from_prev, normal);
          //       right_prev = right_prev / std::sqrt(CGAL::scalar_product(right_prev, right_prev));

          //       const double d_next = CGAL::scalar_product(to_u, right_next);
          //       const double d_prev = CGAL::scalar_product(to_u, right_prev);


          //       if (d_next > 0 || d_prev > 0) {
          //         continue;
          //       }
                
          //       std::cout << "Found a preexisting halfedge between non-consecutive vertices with indices " << i << " and " << j << std::endl;
          //       std::cout << "Vertex u: " << u << std::endl;
          //       std::cout << "Vertex v: " << v << std::endl;
          //       std::cout << "Halfedge u->v: " << (found_uv ? "exists" : "does not exist") << std::endl;
          //       std::cout << "Halfedge v->u: " << (found_vu ? "exists" : "does not exist") << std::endl;
          //       if (found_uv) {
          //         std::cout << "Halfedge u->v descriptor: " << h_uv << std::endl;
          //         std::cout << "Halfedge u->v source: " << source(h_uv, mesh_) << ", target: " << target(h_uv, mesh_) << std::endl;
          //         std::cout << "Halfedge u->v face: " << face(h_uv, mesh_) << std::endl;

          //         std::cout << "All halfedges that go into vertex u:\n";
          //         for (auto h : halfedges_around_target(u, mesh_)) {
          //           std::cout << h << " " << source(h, mesh_) << "->" << target(h, mesh_) << ", next: " << next(h, mesh_) << "\n";
          //         }
          //         std::cout << std::endl;
          //         std::cout << "All halfedges that go into vertex v:\n";
          //         for (auto h : halfedges_around_target(v, mesh_)) {
          //           std::cout << h << " " << source(h, mesh_) << "->" << target(h, mesh_) << ", next: " << next(h, mesh_) << "\n";
          //         }
          //         std::cout << std::endl;
          //       }
          //       if (found_vu) {
          //         std::cout << "Halfedge v->u descriptor: " << h_vu << std::endl;
          //         std::cout << "Halfedge v->u source: " << source(h_vu, mesh_) << ", target: " << target(h_vu, mesh_) << std::endl;
          //         std::cout << "Halfedge v->u face: " << face(h_vu, mesh_) << std::endl;

          //         std::cout << "All halfedges that go into vertex v:\n";
          //         for (auto h : halfedges_around_target(v, mesh_)) {
          //           std::cout << h << " " << source(h, mesh_) << "->" << target(h, mesh_) << ", next: " << next(h, mesh_) << "\n";
          //         }
          //         std::cout << std::endl;
          //         std::cout << "All halfedges that go into vertex u:\n";
          //         for (auto h : halfedges_around_target(u, mesh_)) {
          //           std::cout << h << " " << source(h, mesh_) << "->" << target(h, mesh_) << ", next: " << next(h, mesh_) << "\n";
          //         }
          //         std::cout << std::endl;
          //       }
          //       std::cout << "In the cycle: ";
          //       for (std::size_t k = 0; k < cycle.size(); ++k) {
          //         std::cout << cycle[k] << " [" << source(cycle[k], mesh_) << "]"<< (k + 1 < cycle.size() ? " -> " : "");
          //       }
          //       std::cout << std::endl;

          //       export_halfedge_mesh_json(
          //         mesh_,
          //         "mesh_debug.json"
          //       );
          //     }

          //     CGAL_assertion(!found_uv && !found_vu);
          //   }
          // }

          // Ear clipping converts the polygon cycle into triangle faces while preserving
          // the existing halfedge connectivity.
          // TODO: should this have any reason to fail?
          bool success = ear_clip_and_add_faces(cycle);

          // Check that each halfedge in the cycle is now assigned to a non-null face
          if (success) {
            for (halfedge_descriptor h : cycle) {
              if (face(h, mesh_) == null_f) {
                std::cerr << "Halfedge " << h << " has null face!" << std::endl;
                export_halfedge_mesh_json(
                  mesh_,
                  "mesh_debug_null_face.json"
                );
              }
              CGAL_assertion(face(h, mesh_) != null_f);
            }
          }
        }
      }

      double triangle_angle(const Point_2& a, const Point_2& b, const Point_2& c) const {
        const double ux = CGAL::to_double(a.x() - b.x());
        const double uy = CGAL::to_double(a.y() - b.y());
        const double vx = CGAL::to_double(c.x() - b.x());
        const double vy = CGAL::to_double(c.y() - b.y());

        const double dot = ux * vx + uy * vy;
        const double nu = std::sqrt(ux * ux + uy * uy);
        const double nv = std::sqrt(vx * vx + vy * vy);
        if (nu <= 0.0 || nv <= 0.0)
          return 0.0;

        double cs = dot / (nu * nv);
        // cs = (std::max)(-1.0, (std::min)(1.0, cs));
        return std::acos(cs);
      }

      CGAL::Orientation orient3(const Point_3& a, const Point_3& b, const Point_3& c, const Vector_3& n) const {
        const Vector_3 ab = b - a;
        const Vector_3 ac = c - a;
        const Vector_3 cross = CGAL::cross_product(ab, ac);
        return (cross * n > 0) ? CGAL::LEFT_TURN : CGAL::RIGHT_TURN;
      }

      std::vector<std::pair<double, std::size_t>> find_ear_candidates(
          const std::vector<halfedge_descriptor>& cycle,
          const std::vector<vertex_descriptor>& verts,
          const std::vector<int>& idx)
      {
        // Verify inputs
        CGAL_assertion(cycle.size() >= 3);
        CGAL_assertion(verts.size() == cycle.size());
        CGAL_assertion(idx.size() == cycle.size());

        //bool target_found = false;

        // The cycle must be well-formed (properly connected).
        for (std::size_t i = 0; i < cycle.size(); ++i) {
          const halfedge_descriptor h = cycle[i];
          const halfedge_descriptor h_next = cycle[(i+1) % cycle.size()];

          CGAL_assertion(next(h, mesh_) == h_next);
          CGAL_assertion(target(h, mesh_) == source(h_next, mesh_));
          CGAL_assertion(source(h, mesh_) == verts[i]);
          CGAL_assertion(target(h, mesh_) == verts[(i+1) % cycle.size()]);

          // const vertex_descriptor v = verts[i];
          // if (v.idx() == 2068 && cycle.size() >= 6) {
          //   std::cout << "Boundary points: " << std::endl;
          //   for (std::size_t j = 0; j < cycle.size(); ++j) {
          //     const vertex_descriptor vj = verts[j];
          //     const Point_3 pj = points_pm_[vj];
          //     std::cout << "Point_3(" << pj.x() << ", " << pj.y() << ", " << pj.z() << ")," << std::endl;
          //   }
          //   std::cout << cycle[0] << std::endl;
          //   export_halfedge_mesh_json(mesh_, "boundary_points.json");
          //   CGAL_assertion(false);
          // }
        }
        

        std::vector<std::pair<double, std::size_t>> ear_candidates;

        // 2D orientation test.
        const auto orient2 = [](const Point_2& a, const Point_2& b, const Point_2& c) {
          return CGAL::orientation(a, b, c);
        };

        for (std::size_t i = 0; i < idx.size(); ++i)
        {
          // Candidate ear = (prev,current,next)
          const std::size_t i0 = idx[(i + idx.size() - 1) % idx.size()];
          const std::size_t i1 = idx[i];
          const std::size_t i2 = idx[(i + 1) % idx.size()];

          if (verts[i0] == verts[i2])
            continue;

          // If there already exists a connection between a and c, we cannot cut this ear 
          auto [he_ac, exists_ac] = halfedge(verts[i0], verts[i2], mesh_);
          if (exists_ac) {
            continue;
          }

          const Point_3& a3 = get(points_pm_, verts[i0]);
          const Point_3& b3 = get(points_pm_, verts[i1]);
          const Point_3& c3 = get(points_pm_, verts[i2]);

          // CHANGED: using normals_pm for angle computation vs box normal
          const Vector_3 b_box_normal = surface_approximation_->box_normal(b3);
          const Vector_3 b_normal = get(normals_pm_, verts[i1]);

          bool wrong_orientation_3d_new = orient3(a3, b3, c3, b_normal) != CGAL::LEFT_TURN;


          // Positive: left turn, negative: right turn in 3D projection.
          if (wrong_orientation_3d_new)
          {
            // This is fine since it can happen to concave polygons
            // However, if this happens without a preexisting connection between a and c, it is incorrect
            // Check if there is a preexisting connection between a and c
            //auto [he_ac, exists_ac] = halfedge(verts[i0], verts[i2], mesh_);
            // if (!exists_ac) {
            //   std::cout << "Halfedge between a and c exists: " << exists_ac << std::endl;
            //   std::cout << "Wrong orientation in 3D for triangle: " << i0 << ", " << i1 << ", " << i2 << std::endl;
            //   std::cout << "Triangle vertices (3D): a3=" << a3 << ", b3=" << b3 << ", c3=" << c3 << std::endl;
            //   std::cout << "Reference normal: " << b_normal << std::endl;
            //   std::cout << "Halfedge: " << halfedge(verts[i0], verts[i1], mesh_).first << std::endl;
            //   std::cout << "Orientation check: wrong_orientation_3d=" << wrong_orientation_3d << ", wrong_orientation_3d_new=" << wrong_orientation_3d_new << std::endl;
            //   // Wrong orientation in 3D, can happen on concave polygons
            //   export_halfedge_mesh_json(mesh_, "wrong_orientation.json");
            //   CGAL_assertion(false);
            // }
            continue;
          }

          // CHANGED: removed 2D intersection check

          const double ang = compute_triangle_angle(a3, b3, c3, b_normal);
          ear_candidates.emplace_back(ang, i); // Store candidate together with its interior angle.
        }

        return ear_candidates;
      }

      Point_2 project(const Point_3& p, const Point_3& origin, const Vector_3& u, const Vector_3& v) const {
        return Point_2(CGAL::to_double((p - origin) * u),
                       CGAL::to_double((p - origin) * v));
      }
      
      Vector_3 compute_cycle_normal(const std::vector<vertex_descriptor>& verts) const {
        Vector_3 n = CGAL::NULL_VECTOR;
        for (auto v: verts) {
          n = n + surface_approximation_->box_normal(get(points_pm_, v));
          //n = n + normals_pm_[v];
        }
        const double nlen2 = CGAL::to_double(n.squared_length());
        CGAL_assertion(nlen2 > 0.0);
        n = n / std::sqrt(nlen2);

        return n;
      }

      double compute_triangle_angle(const Point_3& a, const Point_3& b, const Point_3& c, const Vector_3& reference_normal) const {
        // First project the points onto the plane defined by the reference normal.
        const Vector_3 n = reference_normal;
        const Vector_3 u = (b - a) - ((b - a) * n) * n;
        const Vector_3 v = CGAL::cross_product(n, u);
        const Point_2 pa = project(a, a, u, v);
        const Point_2 pb = project(b, a, u, v);
        const Point_2 pc = project(c, a, u, v);

        // Compute the interior angle at vertex b in the projected 2D triangle.
        const double angle = triangle_angle(pa, pb, pc);
        return angle;
      }


      bool contains_subcycle(const std::vector<vertex_descriptor>& cycle)
      {
        std::vector<vertex_descriptor> path;
        path.reserve(cycle.size());

        for (const auto& v : cycle)
        {
          if (path.size() >= 2 && path[path.size() - 2] == v)
          {
            path.pop_back();
          }
          else if (std::find(path.begin(), path.end(), v) != path.end())
          {
            return true;
          }
          else
          {
            path.push_back(v);
          }
        }

        // The input itself is a cycle, so the first and last vertices
        // being equal is not a subcycle.
        return false;
      }

      /**
       * @brief Triangulates a polygonal halfedge cycle using ear clipping.
       *
       * The cycle is projected once onto a local tangent plane. Valid ears are
       * identified using 2D orientation tests, containment tests, and a 3D normal
       * consistency test. Candidate ears are ordered by their interior angle.
       *
       * After an ear is committed, the halfedge cycle is updated and the remaining
       * polygon is processed again. When only three halfedges remain, the final
       * triangle is created directly.
       *
       * @param cycle Halfedge cycle representing a polygon boundary.
       * @return True if the polygon was successfully triangulated, false otherwise.
      */
      bool ear_clip_and_add_faces(std::vector<halfedge_descriptor> cycle) {
        CGAL_assertion(cycle.size() >= 3);

        // 2D orientation test.
        const auto orient2 = [](const Point_2& a, const Point_2& b, const Point_2& c) {
          return CGAL::orientation(a, b, c);
        };

        // Returns true if p lies on or inside triangle (a,b,c).
        const auto point_in_triangle =
          [&](const Point_2& p,
              const Point_2& a,
              const Point_2& b,
              const Point_2& c) -> bool {
          const auto o1 = orient2(a, b, p);
          const auto o2 = orient2(b, c, p);
          const auto o3 = orient2(c, a, p);
          return (o1 == CGAL::LEFT_TURN &&
                  o2 == CGAL::LEFT_TURN &&
                  o3 == CGAL::LEFT_TURN);
        };

        std::vector<vertex_descriptor> verts;
        verts.reserve(cycle.size());
        for (halfedge_descriptor h : cycle) {
          verts.push_back(source(h, mesh_));
        }

        // --------------------------------------------------------------------------
        // Extract the current polygon vertices.
        //
        // A polygon cycle must contain each vertex only once. Repeated vertices mean
        // that the halfedge walk is actually a pinched/multi-cycle boundary and cannot
        // be safely triangulated as a simple polygon.
        // --------------------------------------------------------------------------
        // if (contains_subcycle(verts)) {
        //   std::cerr << "Warning: Ear clipping failed: cycle contains repeated vertices. Can result in holes.\n";
        //   for (halfedge_descriptor h : cycle) {
        //     vertex_descriptor v = source(h, mesh_);
        //     std::cerr << "  Vertex " << v << " at " << get(points_pm_, v) << "\n";
        //     std::cerr << "  Halfedge " << h << "\n";
        //   }
        //   export_halfedge_mesh_json(mesh_, "ear_clipping_failed.json");
        // }
        // CGAL_assertion(!contains_subcycle(verts));

        // Iteratively remove one ear until only one triangle remains.
        while (cycle.size() > 3) {
          verts.clear();

          verts.reserve(cycle.size());
          for (halfedge_descriptor h : cycle)
            verts.push_back(source(h, mesh_));

          std::vector<int> idx(cycle.size());
          std::iota(idx.begin(), idx.end(), 0);

          // ----------------------------------------------------------------------
          // Find every valid ear.
          // ----------------------------------------------------------------------

          std::vector<std::pair<double, std::size_t>> ear_candidates
            = find_ear_candidates(
                cycle,
                verts,
                idx
              );
          
          // There should always be at least one valid ear unless the polygon has 3 or fewer vertices.
          bool found_ear = !ear_candidates.empty();

          if (!found_ear && cycle.size() > 3) {
            // std::cerr << "No valid ear found for a polygon with more than 3 vertices.\n";
            // std::vector<Point_3> points_in_cycle;
            // for (halfedge_descriptor h : cycle) {
            //   points_in_cycle.push_back(mesh_.point(source(h, mesh_)));
            //   std::cout << "source in cycle: " << source(h, mesh_) << std::endl;
            // }
            // std::cout << "3d points: [";
            // for (std::size_t i = 0; i < points_in_cycle.size(); ++i) {
            //   std::cout << "[" << points_in_cycle[i].x() << ", " << points_in_cycle[i].y() << ", " << points_in_cycle[i].z() << "]";
            //   if (i + 1 < points_in_cycle.size()) {
            //     std::cout << ", ";
            //   }
            // }
            // std::cout << "]" << std::endl;
            
            // std::cout << "2d points: (projected on normal << " << n << ") [";
            // for (std::size_t i = 0; i < P.size(); ++i) {
            //   std::cout << "[" << P[i].x() << ", " << P[i].y() << "]";
            //   if (i + 1 < P.size()) {
            //     std::cout << ", ";
            //   }
            // }
            // std::cout << "]" << std::endl;

            // std::cout << "Find ear candidates output:" << std::endl;
            // find_ear_candidates(
            //     cycle,
            //     verts,
            //     P,
            //     n,
            //     origin,
            //     u,
            //     v,
            //     idx,
            //     true
            //   );


            // // int target_vertex = 349;
            // // for (halfedge_descriptor h : cycle) {
            // //   points_in_cycle.push_back(mesh_.point(source(h, mesh_)));
            // //   std::cout << "source in cycle: " << source(h, mesh_) << std::endl;
            // //   int target_vertex = 349;
            // //   if (source(h, mesh_).idx() == target_vertex) {
            // //     std::cout << "h: " << h << std::endl;
            // //     export_halfedge_mesh_json(mesh_, "no_ear_vertex_349.json");
            // //     std::cerr << "Found source vertex " << target_vertex << " in cycle." << std::endl;
            // //     CGAL_assertion(false);
            // //   }
            // // }


            // export_halfedge_mesh_json(mesh_, "no_ear.json");

            //CGAL_assertion(false);
            return false;
          }

          // There should always be at least one valid ear unless the polygon has 3 or fewer vertices.
          CGAL_assertion(found_ear || cycle.size() <= 3);

          // ----------------------------------------------------
          // commit_ear_step must:
          //  - add the ear triangle face,
          //  - update the mesh halfedge wiring,
          //  - erase the ear vertex from 'cycle'.
          // ----------------------------------------------------

          // Try ears from smallest angle to largest until one can be inserted topologically.
          std::sort(ear_candidates.begin(), ear_candidates.end());
          bool ear_committed = false;
          for (const auto& [angle, best_pos] : ear_candidates) {
            // commit_ear_step will now always return true, so no need to check its return value.
            // std::cout << "Trying to commit ear at position " << best_pos << " with angle " << angle << std::endl;
            // std::cout << "That is vertex " << verts[best_pos].idx() << std::endl;
            CGAL_assertion(commit_ear_step(cycle, best_pos));
            break;
          }
          ear_committed = true;

          // This code really shouldn't be reached
          if (!ear_committed) {
            if (cycle.size() == 4) {
              for (std::size_t i = 0; i < cycle.size(); ++i) {
                commit_ear_step(cycle, i);
              }
            }
            return true;
          }
        }

        CGAL_assertion(cycle.size() == 3);

        // Commit the last remaining triangle.
        if (cycle.size() == 3)
        {
          CGAL_assertion(commit_ear_step(cycle, 1));
        }

        CGAL_assertion(cycle.empty());

        return true;
      }

      /**
       * @brief Commits one ear removal to the halfedge mesh.
       *
       * For an ear
       *
       *   a -> b -> c
       *
       * a new diagonal `a -> c` is created, the triangle `a-b-c` is assigned a face,
       * and the remaining polygon boundary is rewired to replace
       *
       *   a -> b -> c
       *
       * by
       *
       *   a -> c.
       *
       * The supplied cycle is updated to remove the ear vertex.
       *
       * @param cycle Current polygon halfedge cycle.
       * @param ear_pos Position of the ear vertex in `cycle`.
       *
       * @return `true` if the ear was successfully committed.
      */
      bool commit_ear_step(std::vector<halfedge_descriptor>& cycle,
                           std::size_t ear_pos) {
        const auto null_h = boost::graph_traits<TriangleMesh>::null_halfedge();
        const auto null_f = boost::graph_traits<TriangleMesh>::null_face();
        const auto null_v = boost::graph_traits<TriangleMesh>::null_vertex();

        const std::size_t n = cycle.size();

        // Ensure that the cycle has at least three vertices and the ear position is valid.
        CGAL_assertion(n >= 3);
        CGAL_assertion(ear_pos < n);

        // The final cycle already consists of exactly three halfedges, so no new
        // diagonal is needed. Convert the cycle directly into a triangular face.
        if (n == 3) {

          // Here we must have a valid cycle of three triangles
          CGAL_assertion(cycle[0] != null_h && cycle[1] != null_h && cycle[2] != null_h);
          CGAL_assertion(target(cycle[0], mesh_) == source(cycle[1], mesh_));
          CGAL_assertion(target(cycle[1], mesh_) == source(cycle[2], mesh_));
          CGAL_assertion(target(cycle[2], mesh_) == source(cycle[0], mesh_));
          CGAL_assertion(face(cycle[0], mesh_) == null_f);
          CGAL_assertion(face(cycle[1], mesh_) == null_f);
          CGAL_assertion(face(cycle[2], mesh_) == null_f);
          CGAL_assertion(next(cycle[0], mesh_) == cycle[1]);
          CGAL_assertion(next(cycle[1], mesh_) == cycle[2]);
          CGAL_assertion(next(cycle[2], mesh_) == cycle[0]);

          face_descriptor f = add_face(mesh_);
          CGAL_assertion(f != null_f);

          set_halfedge(f, cycle[0], mesh_);
          for (halfedge_descriptor h : cycle) {
            set_face(h, f, mesh_);
          }

          // This shouldn't be necessary with the assertions above
          // set_next(cycle[0], cycle[1], mesh_);
          // set_next(cycle[1], cycle[2], mesh_);
          // set_next(cycle[2], cycle[0], mesh_);

          cycle.clear();
          return true;
        }

        const std::size_t i_ab     = (ear_pos + n - 1) % n;
        const std::size_t i_bc     = ear_pos;
        const std::size_t i_before = (ear_pos + n - 2) % n;
        const std::size_t i_after  = (ear_pos + 1) % n;

        // Identify the two halfedges entering and leaving the ear vertex.
        halfedge_descriptor h_ab = cycle[i_ab];
        halfedge_descriptor h_bc = cycle[i_bc];

        CGAL_assertion(h_ab != null_h && h_bc != null_h);
        if (h_ab == null_h || h_bc == null_h) {
          return false;
        }

        const vertex_descriptor a = source(h_ab, mesh_);
        const vertex_descriptor b = target(h_ab, mesh_);
        const vertex_descriptor c = target(h_bc, mesh_);

        CGAL_assertion(a != null_v && b != null_v && c != null_v);
        CGAL_assertion(source(h_bc, mesh_) == b);

        if (source(h_bc, mesh_) != b) {
          return false;
        }

        // Vector_3 normal = normals_pm_[b];

        // std::vector<Vector_3> frame = grid_.compute_local_tangent_frame(normal);
        // const Vector_3& u = frame[0];
        // const Vector_3& v = frame[1];

        // auto project_to_2d = [&](const Point_3& p) -> Point_2 {
        //   const Vector_3 d = p - points_pm_[b];
        //   return Point_2(d * u, d * v);
        // };

        // const bool found = halfedge(a, c, mesh_).second;

        // if (found) {

        //   // Check if the midpoint of the diagonal (a, c) lies inside the projected polygon.
        //   // If it does, this is a real issue because it indicates that the diagonal intersects the interior of the polygon.


        //   // Project all vertices of the loop onto the plane defined by the normal at vertex b.
        //   CGAL::Polygon_2<Kernel> projected_polygon;

        //   for (halfedge_descriptor h : cycle) {
        //     const Point_3 p = points_pm_[target(h, mesh_)];
        //     projected_polygon.push_back(project_to_2d(p));
        //   }

        //   if (projected_polygon.is_simple()) {
        //     Point_3 pa = points_pm_[a];
        //     Point_3 pc = points_pm_[c];

        //     Point_3 mid_ac = CGAL::midpoint(pa, pc);

        //     Point_2 sample_point = project_to_2d(points_pm_[b]);
        //     const auto side = projected_polygon.bounded_side(sample_point);

        //     // If it is inside, that indicates a wrong connection
        //     CGAL_assertion(side != CGAL::ON_BOUNDED_SIDE);

        //     // Definitely no need to add another edge
        //     return false;
        //   }

        //   // TODO: add non simple check?
        // }
        

        // Create the diagonal a -> c. Its opposite halfedge c -> a is used to close
        // the new triangular face.
        halfedge_descriptor h_ac = create_halfedge(a, c);
        halfedge_descriptor h_ca = opposite(h_ac, mesh_);
        if (h_ac.idx() == 744 || h_ca.idx() == 744) {
          std::cout << "creating edge 744 in ear committing" << std::endl;
          //CGAL_assertion(false);
        }

        // The two boundary halfedges of the ear must not already belong to faces.
        CGAL_assertion(face(h_ac, mesh_) == null_f && face(h_ca, mesh_) == null_f);
        CGAL_assertion(face(h_ab, mesh_) == null_f && face(h_bc, mesh_) == null_f);
        if (face(h_ac, mesh_) != null_f && face(h_ca, mesh_) != null_f) {
          std::cout << "both new halfedges already belong to a face" << std::endl;
          return false;
        }
        if (face(h_ab, mesh_) != null_f || face(h_bc, mesh_) != null_f) {
          std::cout << "one of the boundary halfedges already belongs to a face" << std::endl;
          return false;
        }
        
        // Construct the new triangular face:
        //     a -> b -> c -> a
        face_descriptor f = add_face(mesh_);
        CGAL_assertion(f != null_f);
        if (f == null_f) {
          return false;
        }

        if (face(h_ca, mesh_) == null_f) {
          set_halfedge(f, h_ab, mesh_);

          set_face(h_ab, f, mesh_);
          set_face(h_bc, f, mesh_);
          set_face(h_ca, f, mesh_);

          set_next(h_ab, h_bc, mesh_);
          set_next(h_bc, h_ca, mesh_);
          set_next(h_ca, h_ab, mesh_);
        }

        // Rewire the remaining boundary:
        // ... -> a -> b -> c -> ...
        // becomes
        // ... -> a -> c -> ...
        halfedge_descriptor h_before = cycle[i_before];
        halfedge_descriptor h_after  = cycle[i_after];
        CGAL_assertion(h_before != null_h && h_after != null_h);

        set_next(h_before, h_ac, mesh_);
        set_next(h_ac, h_after, mesh_);

        // Update the explicit cycle representation so that the removed ear vertex
        // is no longer part of the polygon being triangulated.
        cycle[i_ab] = h_ac;
        cycle.erase(cycle.begin() + static_cast<std::ptrdiff_t>(i_bc));

        // Verify that the cycle is still valid after removing the ear vertex.
        for (std::size_t i = 0; i < cycle.size(); ++i) {
          CGAL_assertion(next(cycle[i], mesh_) == cycle[(i + 1) % cycle.size()]);
          CGAL_assertion(target(cycle[i], mesh_) == source(cycle[(i + 1) % cycle.size()], mesh_));
          CGAL_assertion(face(cycle[i], mesh_) == null_f);
        }

        // Verify that the newly created face has the correct halfedges.
        CGAL_assertion(face(h_ab, mesh_) == f);
        CGAL_assertion(face(h_bc, mesh_) == f);
        CGAL_assertion(face(h_ca, mesh_) == f);
        CGAL_assertion(next(h_ab, mesh_) == h_bc);
        CGAL_assertion(next(h_bc, mesh_) == h_ca);
        CGAL_assertion(next(h_ca, mesh_) == h_ab);
        CGAL_assertion(source(h_ab, mesh_) == a);
        CGAL_assertion(source(h_bc, mesh_) == b);
        CGAL_assertion(source(h_ca, mesh_) == c);
        CGAL_assertion(target(h_ab, mesh_) == b);
        CGAL_assertion(target(h_bc, mesh_) == c);
        CGAL_assertion(target(h_ca, mesh_) == a);

        return true;
      }

      /**
       * @brief Checks whether a candidate is geometrically compatible with the
       *        surrounding mesh under a local planar projection.
       *
       * The parent-to-candidate edges are projected into the tangent plane defined
       * by the supplied normal. They are rejected when they intersect existing mesh
       * edges.
       *
       * @param cand Candidate reconstruction vertex.
       * @param n Projection normal.
       *
       * @return `true` when the candidate passes the local intersection test.
      */
      bool projection_check(
        const Point_3 point_1,
        const Point_3 point_2,
        const Vector_3& n) const
      {
        CGAL_assertion(n != CGAL::NULL_VECTOR);

        const Point_3 center = point_1 + (point_2 - point_1) / FT(2.0);

        auto local_frame = CGAL::internal::Splat_surface_reconstruction_3::compute_local_tangent_frame<Point_3, Vector_3, FT>(n);
        const Vector_3 u = local_frame[0];
        const Vector_3 v = local_frame[1];

        const Point_2 proj_1 = project(point_1, center, u, v);
        const Point_2 proj_2 = project(point_2, center, u, v);

        auto nearby = nearby_mesh_vertices(
          center, FT(3.0) * grid_.get_box_size()
        );

        const FT radius = FT(2.0) * grid_.get_box_size();
        const FT height = FT(1.0) * grid_.get_box_size();

        for (vertex_descriptor vtx : nearby) {
          const Point_3 v_point = get(points_pm_, vtx);

          if (!CGAL::internal::Splat_surface_reconstruction_3::point_in_cylinder<Point_3, Vector_3, FT>(v_point, center, n, radius, height)) {
            continue;
          }

          for (halfedge_descriptor h : halfedges_around_target(vtx, mesh_)) {
            const vertex_descriptor s = source(h, mesh_);
            const vertex_descriptor t = target(h, mesh_);

            const Point_3 q0 = get(points_pm_, s);
            const Point_3 q1 = get(points_pm_, t);

            const Point_2 E0 = project(q0, center, u, v);
            const Point_2 E1 = project(q1, center, u, v);

            if (CGAL::internal::Splat_surface_reconstruction_3::intersect_2d<Point_2, Kernel>(proj_1, proj_2, E0, E1)) {
              return false;
            }
          }
        }

        return true;
      }

      /**
       * @brief Searches the local halfedge neighborhood for another parent vertex.
       *
       * The search follows `next()` pointers for a bounded number of steps from
       * every incident halfedge.
       *
       * @param start_v Starting vertex.
       * @param target_v Vertex being searched for.
       *
       * @return `true` if the target vertex is encountered.
      */
      bool walk_finds_other_parent_from_vertex(vertex_descriptor start_v,
                                               vertex_descriptor target_v) const {

        for (halfedge_descriptor h0 : halfedges_around_source(start_v, mesh_)) {
          if (h0 == boost::graph_traits<TriangleMesh>::null_halfedge()) {
            return false;
          }

          halfedge_descriptor h = h0;
          for (std::size_t i = 0; i < window_size_; ++i) {
            const vertex_descriptor s = source(h, mesh_);
            const vertex_descriptor t = target(h, mesh_);

            if (s == target_v || t == target_v) {
              return true;
            }

            const halfedge_descriptor n = next(h, mesh_);
            if (n == boost::graph_traits<TriangleMesh>::null_halfedge()) {
              return false;
            }

            h = n;
          }
        }

        return false;
      }

      /**
       * @brief Determines whether a candidate connects two distinct boundary fronts.
       *
       * A candidate is considered to join two borders when neither parent can reach
       * the other through the local halfedge walks.
      */
      bool joins_two_borders(const Candidate& cand) const {
        const bool a_finds_b = walk_finds_other_parent_from_vertex(cand.first, cand.second);
        const bool b_finds_a = walk_finds_other_parent_from_vertex(cand.second, cand.first);

        return !a_finds_b && !b_finds_a;
      }

      /**
       * @brief Computes the processing priority of a reconstruction candidate.
       *
       * Higher priority is assigned to candidates adjacent to low-degree vertices,
       * followed by candidates that connect two separate boundary fronts.
       *
       * @param cand Candidate being prioritized.
       *
       * @return Integer priority in the range `[0, NUM_PRIORITIES-1]`.
      */
      int compute_priority(const Candidate& cand) const {

        // CHANGED: added this new priority rule
        // "At least one of the parent vertices has no edges incident to it"
        if (out_degree(cand.first, mesh_) == 0 || out_degree(cand.second, mesh_) == 0) {
          return 3; // highest priority
        }

        // "At least one of the parent vertices is a vertex with only one edge incident to it"
        if (out_degree(cand.first, mesh_) == 1 || out_degree(cand.second, mesh_) == 1) {
          return 2; // second highest priority
        }

        // "Inserting vc and its two edges joins two borders"
        if (joins_two_borders(cand)) {
          return 1; // medium priority
        }

        // "Inserting vc splits a border"
        return 0;   // lowest priority
      }














      // Feature awareness functions




      void seed_from_features(FT minimum_angle) {
        feature_vertices_per_cell_.clear();
        feature_vertices_per_cell_.resize(grid_.number_of_cells(), std::vector<Point_3>());
        feature_segment_indices_per_cell_.clear();
        feature_segment_indices_per_cell_.resize(grid_.number_of_cells(), std::vector<std::size_t>());
        feature_segment_is_used_.clear();
        feature_segment_is_used_.resize(feature_segments_.size(), false);

        feature_segments_ = surface_approximation_->find_feature_segments(minimum_angle);

        // Add all feature segments to their respective cells
        for (std::size_t i = 0; i < feature_segments_.size(); ++i) {
          const Feature_segment& seg = feature_segments_[i];

          // Add the feature segment to the grid cells
          // Add both ends of the segment
          const Point_3& p1 = seg.point_1;
          const Point_3& p2 = seg.point_2;

          const Point_3 p_min = Point_3(
            std::min(p1.x(), p2.x()),
            std::min(p1.y(), p2.y()),
            std::min(p1.z(), p2.z())
          );

          const Point_3 p_max = Point_3(
            std::max(p1.x(), p2.x()),
            std::max(p1.y(), p2.y()),
            std::max(p1.z(), p2.z())
          );

          int cell_min_x, cell_min_y, cell_min_z;
          int cell_max_x, cell_max_y, cell_max_z;

          grid_.to_grid_coordinates(p_min, cell_min_x, cell_min_y, cell_min_z);
          grid_.to_grid_coordinates(p_max, cell_max_x, cell_max_y, cell_max_z);

          for (int x = cell_min_x; x <= cell_max_x; ++x) {
            for (int y = cell_min_y; y <= cell_max_y; ++y) {
              for (int z = cell_min_z; z <= cell_max_z; ++z) {
                if (!CGAL::internal::Splat_surface_reconstruction_3::line_intersects_aabb<Point_3, Vector_3, FT>(
                      p1, p2,
                      grid_.cell_min_corner(x, y, z),
                      grid_.cell_max_corner(x, y, z))) {
                  continue;
                }

                std::size_t flat_index = grid_.flat_index(x, y, z);
                feature_segment_indices_per_cell_[flat_index].push_back(i);
              }
            }
          }
        }

        std::cout << "Found feature segments: " << feature_segments_.size() << std::endl;

        find_initial_feature_vertex_candidates(minimum_angle);
        std::cout << "Found initial feature vertex candidates: " << feature_vertex_candidate_queue_.size() << std::endl;

        // No feature segments have been used yet
        feature_segment_is_used_.assign(feature_segments_.size(), false);

        int loop_count = 0;
        while (true) {
          process_feature_vertex_candidates_from_queue();

          // Iterate through all feature segments, adding any unused segments as feature vertex candidates
          bool any_added = add_feature_vertices_from_unused_segments();

          ++loop_count;

          if (!any_added) {
            break;
          }
        }

        std::cout << "Total loops required for feature vertex generation: " << loop_count << std::endl;

        // close feature vertex gaps
        close_feature_vertex_gaps();



        std::cout << "Feature vertex generation finished" << std::endl;
        std::cout << "Total feature vertices: " << mesh_.number_of_vertices() << std::endl;
        std::cout << "Total feature edges: " << mesh_.number_of_edges() << std::endl;

        export_halfedge_mesh_json(mesh_, "feature_vertices.json");

        // insert all feature candidates from pairs of attached feature vertices
        for (const halfedge_descriptor& he : mesh_.halfedges()) {
          const vertex_descriptor& v = source(he, mesh_);
          const vertex_descriptor& u = target(he, mesh_);

          // We only need one call per pair of vertices, so we can skip if u < v to avoid duplicates
          if (u.idx() < v.idx()) {
            continue;
          }

          std::vector<Candidate> candidates = generate_candidates_from_parents(u, v);

          for (const Candidate& c : candidates) {
            push_candidate_in_queue(c);
          }
        }

        if (mesh_.number_of_vertices() == 0) {
          std::cout << "No feature vertices were generated. Seeding from grid instead." << std::endl;
          seed_from_grid();
        }
        else {
          seeded_ = true;
        }
      }


      void process_feature_vertex_candidates_from_queue() {
        
        const FT d = grid_.get_box_size();

        while (!feature_vertex_candidate_queue_.empty()) {
          FeatureVertexCandidate cand = feature_vertex_candidate_queue_.front();
          feature_vertex_candidate_queue_.pop_front();

          if (is_any_feature_vertex_too_close(cand.position, d)) {
            continue;
          }

          if (cand.has_parent_vertex && !projection_check(cand.position, points_pm_[cand.parent_vertex], surface_approximation_->box_normal(cand.position))) {
            continue;
          }

          // Now we create a feature vertex for this candidate.
          vertex_descriptor added_vertex = create_feature_vertex(cand);

          if (cand.has_feature_segment) {
            feature_segment_is_used_[cand.feature_segment_index] = true;
          }

          // And we need to add any new feature vertex candidates that arise from this insertion.
          add_new_feature_vertex_candidates(cand, added_vertex);
        }
      }

      bool add_feature_vertices_from_unused_segments() {

        for (std::size_t i = 0; i < feature_segments_.size(); ++i) {
          if (!feature_segment_is_used_[i]) {
            const Feature_segment& feature_segment = feature_segments_[i];
            Point_3 midpoint = CGAL::midpoint(feature_segment.point_1, feature_segment.point_2);
            FeatureVertexCandidate cand{
              midpoint,
              feature_segment.normal,
              true,  // has feature segment
              false, // no parent vertex
              i
            };
            feature_vertex_candidate_queue_.push_back(cand);

            // This segment can no longer be considered unused since we are adding a candidate for it.
            feature_segment_is_used_[i] = true;

            return true;
          }
        }
        return false;
      }

      void close_feature_vertex_gaps() {

        const FT d = grid_.get_box_size();

        // Iterate through all vertices
        for (vertex_descriptor v : vertices(mesh_)) {
          // Only keep it if it has valence = 1
          if (out_degree(v, mesh_) != 1) {
            continue;
          }

          // Now iterate over all vertices within distance 2d, also with valence 1
          int ix, iy, iz;
          if (!grid_.to_grid_coordinates(mesh_.point(v), ix, iy, iz))
            continue;

          // Find the closest vertex to v within this grid cell that has valence 1 and is not already connected to v
          vertex_descriptor closest_neighbor_v;
          FT closest_dist = 2 * d;
          bool found_closest = false;

          int reach = 2;
          for (int dx = -reach; dx <= reach; ++dx) {
            for (int dy = -reach; dy <= reach; ++dy) {
              for (int dz = -reach; dz <= reach; ++dz) {

                int neighbor_ix = ix + dx;
                int neighbor_iy = iy + dy;
                int neighbor_iz = iz + dz;

                if (!grid_.valid_coordinates(neighbor_ix, neighbor_iy, neighbor_iz)) {
                  continue;
                }

                std::size_t flat_index = grid_.flat_index(neighbor_ix, neighbor_iy, neighbor_iz);
                for (const vertex_descriptor& neighbor_v : mesh_vertices_per_cell_[flat_index]) {
                  if (neighbor_v == v || out_degree(neighbor_v, mesh_) != 1) {
                    continue;
                  }

                  // Don't add anything if the only attached edge is already between v and neighbor_v
                  if (halfedge(v, neighbor_v, mesh_).second) continue;

                  // Check that the distance between v and neighbor_v is less than 2d
                  const Point_3 p_v = mesh_.point(v);
                  const Point_3 p_neighbor_v = mesh_.point(neighbor_v);
                  const FT dist = std::sqrt(CGAL::squared_distance(p_v, p_neighbor_v));
                  if (dist >= 2 * d) {
                    continue;
                  }

                  if (dist < closest_dist) {
                    // We must make sure that the projection check passes for both vertices
                    if (!projection_check(p_v, p_neighbor_v, surface_approximation_->box_normal(p_v)) ||
                        !projection_check(p_v, p_neighbor_v, surface_approximation_->box_normal(p_neighbor_v))) {
                      continue;
                    }

                    closest_dist = dist;
                    closest_neighbor_v = neighbor_v;
                    found_closest = true;
                  }
                }
              }
            }
          }

          if (!found_closest) {
            continue;
          }

          // This is the current situation:
          // v0 -- v    n -- v1
          // where v0 and v1 are the neighbors of v and n respectively.

          const vertex_descriptor neighbor_v = closest_neighbor_v;
          const Point_3 p_v = mesh_.point(v);
          const Point_3 p_neighbor_v = mesh_.point(neighbor_v);
          const FT dist = std::sqrt(CGAL::squared_distance(p_v, p_neighbor_v));

          const halfedge_descriptor h_v_v0 = opposite(halfedge(v, mesh_), mesh_);
          CGAL_assertion(source(h_v_v0, mesh_) == v);
          const halfedge_descriptor h_v0_v = opposite(h_v_v0, mesh_);

          const halfedge_descriptor h_n_v1 = opposite(halfedge(neighbor_v, mesh_), mesh_);
          CGAL_assertion(source(h_n_v1, mesh_) == neighbor_v);
          const halfedge_descriptor h_v1_n = opposite(h_n_v1, mesh_);


          // Now we have found two vertices with valence 1 within distance 2d, 
          // so we can connect them by adding a new edge between them
          halfedge_descriptor h_v_n = create_halfedge(v, neighbor_v);
          halfedge_descriptor h_n_v = opposite(h_v_n, mesh_);
          if (h_v_n.idx() == 744 || h_n_v.idx() == 744) {
            std::cout << "creating edge 744 in the gap closing" << std::endl;
            //CGAL_assertion(false);
          }

          set_target(h_v_n, neighbor_v, mesh_); // v -> n
          set_target(h_n_v, v, mesh_); // n -> v
          set_halfedge(v, h_n_v, mesh_);
          set_halfedge(neighbor_v, h_v_n, mesh_);

          CGAL_assertion(source(h_v_n, mesh_) == v && target(h_v_n, mesh_) == neighbor_v);
          CGAL_assertion(source(h_n_v, mesh_) == neighbor_v && target(h_n_v, mesh_) == v);

          // Connect the new edge to the existing halfedges
          set_next(h_v0_v, h_v_n, mesh_);
          set_next(h_v_n, h_n_v1, mesh_);
          set_next(h_v1_n, h_n_v, mesh_);
          set_next(h_n_v, h_v_v0, mesh_);

          CGAL_assertion(CGAL::is_valid_halfedge_graph(mesh_, false));

          // Move the two vertices slightly towards each other to reduce the gap
          const FT factor = FT(0.33) * (dist - d) / dist;
          mesh_.point(v) =
              p_v + (CGAL::midpoint(p_v, p_neighbor_v) - p_v) * factor;

          mesh_.point(neighbor_v) =
              p_neighbor_v + (CGAL::midpoint(p_v, p_neighbor_v) - p_neighbor_v) * factor;
        }

      }


      void add_new_feature_vertex_candidates(const FeatureVertexCandidate& cand, vertex_descriptor parent_vertex_index) {
        int ix, iy, iz;

        if (!grid_.to_grid_coordinates(cand.position, ix, iy, iz))
          return;

        const FT sphere_radius = grid_.get_box_size();
        const Point_3& sphere_center = cand.position;

        // Iterate over all neighbouring cells to check for feature segments
        int reach = 2;
        for (int dx = -reach; dx <= reach; ++dx) {
          for (int dy = -reach; dy <= reach; ++dy) {
            for (int dz = -reach; dz <= reach; ++dz) {
              int neighbor_ix = ix + dx;
              int neighbor_iy = iy + dy;
              int neighbor_iz = iz + dz;

              if (!grid_.valid_coordinates(neighbor_ix, neighbor_iy, neighbor_iz)) {
                continue;
              }

              std::size_t flat_index = grid_.flat_index(neighbor_ix, neighbor_iy, neighbor_iz);
              for (const std::size_t& feature_segment_index : feature_segment_indices_per_cell_[flat_index]) {
                const Feature_segment& feature_segment = feature_segments_[feature_segment_index];
                std::vector<Point_3> intersection_points = 
                  CGAL::internal::Splat_surface_reconstruction_3::line_intersects_sphere<Point_3, Vector_3, FT>(
                    feature_segment.point_1,
                    feature_segment.point_2,
                    sphere_center,
                    sphere_radius);

                for (const Point_3& intersection_point : intersection_points)
                {
                  // We have found a new candidate
                  feature_vertex_candidate_queue_.push_back(
                    FeatureVertexCandidate{
                      intersection_point,
                      feature_segment.normal,
                      true, // does have an associated feature segment
                      true, // does have a parent vertex
                      feature_segment_index,
                      parent_vertex_index
                    });
                }
              }
            }
          }
        }
      }


      vertex_descriptor create_feature_vertex(const FeatureVertexCandidate& cand) {
        const auto null_h = boost::graph_traits<TriangleMesh>::null_halfedge();
        const double two_pi = 2.0 * std::acos(-1.0);

        auto make_frame = [&](const Vector_3& n) -> std::pair<Vector_3, Vector_3> {
          std::vector<Vector_3> frame = CGAL::internal::Splat_surface_reconstruction_3::compute_local_tangent_frame<Point_3, Vector_3, FT>(n);
          Vector_3 u = frame[0];
          Vector_3 v = frame[1];

          if (u == CGAL::NULL_VECTOR || v == CGAL::NULL_VECTOR) {
            return {CGAL::NULL_VECTOR, CGAL::NULL_VECTOR};
          }

          if (CGAL::cross_product(u, v) * n < FT(0)) {
            v = -v;
          }

          CGAL_assertion(u != CGAL::NULL_VECTOR && v != CGAL::NULL_VECTOR);

          CGAL_assertion(u.squared_length() > FT(0));
          CGAL_assertion(v.squared_length() > FT(0));

          // Their cross product should agree with the normal
          CGAL_assertion(CGAL::cross_product(u, v) * n >= FT(0));

          // Check that u and v are perpendicular to the normal n.
          CGAL_assertion(
            std::abs(CGAL::to_double(u * n)) < 1e-8
          );
          CGAL_assertion(
            std::abs(CGAL::to_double(v * n)) < 1e-8
          );
          // Check that u and v are perpendicular to each other.
          CGAL_assertion(
            std::abs(CGAL::to_double(u * v)) < 1e-8
          );

          return {u, v};
        };

        // Return the polar angle of q around origin in the tangent frame.
        auto angle_of = [&](const Point_3& origin,
                            const Point_3& q,
                            const Vector_3& u,
                            const Vector_3& v) -> double {
          const Vector_3 d = q - origin;
          double a = std::atan2(CGAL::to_double(d * v), CGAL::to_double(d * u));
          if (a < 0.0) {
            a += two_pi;
          }
          if (a >= two_pi) {
            a = 0.0;
          }
          return a;
        };

        // Collect all incident halfedges around a parent and sort them by their
        // projected angle. These halfedges define the local front ordering.
        auto sorted_halfedges_around_vertex =
          [&](vertex_descriptor v, const Vector_3& local_normal)
            -> std::vector<std::pair<double, halfedge_descriptor>> {
          std::vector<std::pair<double, halfedge_descriptor>> out;
          const Point_3 pv = get(points_pm_, v);

          const auto [u, vaxis] = make_frame(local_normal);
          if (u == CGAL::NULL_VECTOR || vaxis == CGAL::NULL_VECTOR) {
            return out;
          }

          for (halfedge_descriptor h : halfedges_around_source(v, mesh_)) {
            CGAL_assertion(source(h, mesh_) == v);
            const vertex_descriptor t = target(h, mesh_);
            if (t == v) {
              continue;
            }
            out.emplace_back(angle_of(pv, get(points_pm_, t), u, vaxis), h);
          }

          std::sort(out.begin(), out.end(),
            [](const auto& a, const auto& b) { return a.first < b.first; });

          for (std::size_t i = 1; i < out.size(); ++i) {
            CGAL_assertion(out[i - 1].first <= out[i].first);
          }

          for (const auto& [angle, he] : out) {
            CGAL_assertion(he != null_h);
            CGAL_assertion(source(he, mesh_) == v);
            CGAL_assertion(angle >= 0.0 && angle < two_pi);
          }

          return out;
        };

        // Select the two front halfedges immediately clockwise and counter-clockwise
        // from the direction toward the new vertex.
        auto choose_closest_ccw_halfedges =
          [&](vertex_descriptor v,
              const Point_3& new_point,
              const Vector_3& local_normal)
            -> std::pair<halfedge_descriptor, halfedge_descriptor> {
          const Point_3 pv = get(points_pm_, v);
          const auto [u, vaxis] = make_frame(local_normal);
          if (u == CGAL::NULL_VECTOR || vaxis == CGAL::NULL_VECTOR) {
            return {null_h, null_h};
          }

          const double theta_new = angle_of(pv, new_point, u, vaxis);

          std::vector<std::pair<double, halfedge_descriptor>> sorted =
            sorted_halfedges_around_vertex(v, local_normal);

          if (sorted.empty()) {
            return {null_h, null_h};
          }

          auto it = std::lower_bound(
            sorted.begin(), sorted.end(), theta_new,
            [](const std::pair<double, halfedge_descriptor>& a, double theta) {
              return a.first < theta;
            });

          const auto& ccw = (it == sorted.end())   ? sorted.front() : *it;
          const auto& cw  = (it == sorted.begin()) ? sorted.back()  : *(it - 1);

          return {cw.second, ccw.second};
        };

        insert_feature_vertex(cand.position);

        // Add this vertex to the mesh and maps
        vertex_descriptor v = add_vertex(mesh_);
        put(points_pm_, v, cand.position);
        put(normals_pm_, v, cand.normal);

        // Add the edges
        if (cand.has_parent_vertex) {
          bool has_previous_and_next_halfedges = false;
          halfedge_descriptor previous_halfedge, next_halfedge;
          std::vector<halfedge_descriptor> halfedges_out_of_parent;

          if (!is_isolated(cand.parent_vertex, mesh_)) {
            for (const auto h : halfedges_around_source(cand.parent_vertex, mesh_)) {
              halfedges_out_of_parent.push_back(h);
            }
          }

          // If there are no halfedges out of the parent, we cannot determine previous and next halfedges.
          if (halfedges_out_of_parent.empty()) {
            has_previous_and_next_halfedges = false;
          }
          // If there is only one halfedge out of the parent, we can set the previous and next halfedges accordingly.
          else if (halfedges_out_of_parent.size() == 1) {
            has_previous_and_next_halfedges = true;
            next_halfedge = halfedges_out_of_parent.front();
            previous_halfedge = opposite(next_halfedge, mesh_);
          }
          // If there are multiple halfedges out of the parent, we can choose the appropriate previous and next halfedges based on the local ordering.
          else {
            has_previous_and_next_halfedges = true;
            
            // Find the first clockwise and counterclockwise halfedges around the parent vertex.
            auto [cw, ccw] = choose_closest_ccw_halfedges(cand.parent_vertex, cand.position, normals_pm_[v]);

            previous_halfedge = ccw;
            next_halfedge = cw;

            // Make sure that they point the right way
            if (target(previous_halfedge, mesh_) != cand.parent_vertex) {
              previous_halfedge = opposite(previous_halfedge, mesh_);
            }
            if (source(next_halfedge, mesh_) != cand.parent_vertex) {
              next_halfedge = opposite(next_halfedge, mesh_);
            }
          }


          halfedge_descriptor h_vp = create_halfedge(v, cand.parent_vertex);
          halfedge_descriptor h_pv = opposite(h_vp, mesh_);

          if (has_previous_and_next_halfedges) {
            // ho is parent -> new vertex so it comes right after the previous halfedge around the parent vertex.
            set_next(previous_halfedge, h_pv, mesh_);

            // next_halfedge should come right after h around the parent vertex.
            set_next(h_vp, next_halfedge, mesh_);

            // next of h_pv should be h_vp (as it loops back around the new vertex)
            set_next(h_pv, h_vp, mesh_);
          }
          else {
            set_next(h_vp, h_pv, mesh_);
            set_next(h_pv, h_vp, mesh_);
          }
        }

        // Add the new vertex to the mesh too
        insert_mesh_vertex(v);

        if (!CGAL::is_valid_halfedge_graph(mesh_, false)) {
          CGAL::is_valid_halfedge_graph(mesh_, true);
          export_halfedge_mesh_json(mesh_, "mesh_debug.json");
        }

        CGAL_assertion(CGAL::is_valid_halfedge_graph(mesh_, false));

        return v;
      }


      void insert_feature_vertex(Point_3 p) {
        int ix, iy, iz;

        if (!grid_.to_grid_coordinates(p, ix, iy, iz))
          return;

        feature_vertices_per_cell_[grid_.flat_index(ix, iy, iz)].push_back(p);
      }

      void find_initial_feature_vertex_candidates(FT minimum_angle) {
        feature_vertex_candidate_queue_.clear();

        const std::vector<std::pair<Point_3, Vector_3>> initial_candidates = surface_approximation_->find_initial_feature_vertex_candidates(minimum_angle);

        for (const auto& [point, normal] : initial_candidates) {
          feature_vertex_candidate_queue_.push_back(
            FeatureVertexCandidate{
              point,
              normal,
              false, // does not have an associated feature segment
              false, // does not have a parent vertex
            });
        }
      }


      bool is_any_feature_vertex_too_close(const Point_3& point, FT min_distance) const {
        int ix, iy, iz;

        if (!grid_.to_grid_coordinates(point, ix, iy, iz))
          return false;

        int reach = 1;
        for (int dx = -reach; dx <= reach; ++dx) {
          for (int dy = -reach; dy <= reach; ++dy) {
            for (int dz = -reach; dz <= reach; ++dz) {
              int neighbor_ix = ix + dx;
              int neighbor_iy = iy + dy;
              int neighbor_iz = iz + dz;

              if (!grid_.valid_coordinates(neighbor_ix, neighbor_iy, neighbor_iz)) {
                continue;
              }

              std::size_t flat_index = grid_.flat_index(neighbor_ix, neighbor_iy, neighbor_iz);
              for (const Point_3& feature_vertex : feature_vertices_per_cell_[flat_index]) {
                if (CGAL::squared_distance(point, feature_vertex) < min_distance * min_distance - 1e-6) {
                  return true;
                }
              }
            }
          }
        }

        return false;
      }





    private:

      // Constructor to set the state for a test
      Splat_surface_reconstruction_3(const Grid& grid, TriangleMesh& mesh, std::vector<Vector_3> normals, int i)
        : grid_(grid),
          mesh_(mesh),
          points_pm_(get(CGAL::vertex_point, mesh))
      {
        std::tie(normals_pm_, std::ignore) =
          mesh_.template add_property_map<vertex_descriptor, Vector_3>(
            "v:normal", CGAL::NULL_VECTOR);

        for (vertex_descriptor v : vertices(mesh_))
        {
          normals_pm_[v] = normals[v.idx()];
        }

        std::tie(edge_weight_pm_, std::ignore) =
          mesh_.template add_property_map<edge_descriptor, FT>(
            "e:weight", FT(1));

        mesh_vertices_per_cell_.resize(grid_.number_of_cells());
      }

      Grid grid_;
      TriangleMesh& mesh_;
      std::unique_ptr<Surface_approximation> surface_approximation_;

      typename TriangleMesh:: template Property_map<vertex_descriptor,Point_3> points_pm_;
      typename TriangleMesh:: template Property_map<vertex_descriptor,Vector_3> normals_pm_;

      typename TriangleMesh::template Property_map<edge_descriptor, FT> edge_weight_pm_;

      vertex_descriptor v0_{boost::graph_traits<TriangleMesh>::null_vertex()};
      vertex_descriptor v1_{boost::graph_traits<TriangleMesh>::null_vertex()};
      bool seeded_ = false;

      static constexpr int NUM_PRIORITIES = 4; // 0 is lowest, NUM_PRIORITIES-1 is highest
      std::array<std::deque<Candidate>, NUM_PRIORITIES> candidate_queues_; // candidate_queues_[p] contains all candidates with priority p

      const std::size_t window_size_ = 80;

      std::vector<std::vector<vertex_descriptor>> mesh_vertices_per_cell_;

      // Variables for feature awareness
      std::deque<FeatureVertexCandidate> feature_vertex_candidate_queue_;
      std::vector<std::vector<Point_3>> feature_vertices_per_cell_;
      std::vector<Feature_segment> feature_segments_;
      std::vector<std::vector<std::size_t>> feature_segment_indices_per_cell_;

      std::vector<bool> feature_segment_is_used_;
  };

} // namespace CGAL

#endif // CGAL_SPLAT_SURFACE_RECONSTRUCTION_H
