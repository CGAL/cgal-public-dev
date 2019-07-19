#ifndef CGAL_SHAPE_REGULARIZATION_INTERNAL_GROUPING_SEGMENTS_ORDINATES_2
#define CGAL_SHAPE_REGULARIZATION_INTERNAL_GROUPING_SEGMENTS_ORDINATES_2

// #include <CGAL/license/Shape_regularization.h>

#include <vector>
#include <map>
#include <utility>

#include <CGAL/Shape_regularization/internal/Segment_data_2.h>

#include <eigen3/Eigen/SparseCore>

namespace CGAL {
namespace Regularization {
namespace internal {

  template<
    typename GeomTraits,
    typename InputRange>
  class Grouping_segments_ordinates_2 {
  public:
    using Traits = GeomTraits;
    using Input_range = InputRange;
    using FT = typename GeomTraits::FT;
    using Segment_data = typename internal::Segment_data_2<Traits>;
    using Sparse_matrix_FT  = Eigen::SparseMatrix<FT,  Eigen::RowMajor>;
    // using Sparse_matrix_int = Eigen::SparseMatrix<int, Eigen::RowMajor>;
    using FT_triplet  = Eigen::Triplet<FT>;
    using Int_triplet = Eigen::Triplet<int>;

    Grouping_segments_ordinates_2(
      const InputRange& input_range,
      const std::map <std::size_t, Segment_data> & segments) :
    m_input_range(input_range),
    m_segments(segments),
    theta_eps(FT(1) / FT(4)),
    y_eps(FT(1)),
    tolerance(FT(1) / FT(1000000)) {

      CGAL_precondition(input_range.size() > 0);
      CGAL_precondition(segments.size() > 0);

    }

    void make_groups(const std::map <std::pair<std::size_t, std::size_t>, FT> & t_ijs,
                     const FT mu_ij, const std::vector<FT> & orientations,
                     const std::map<FT, std::vector<std::size_t>> & parallel_groups_by_angles,
                     std::map<FT, std::vector<std::size_t>> & collinear_groups_by_angles) { 

      CGAL_precondition(t_ijs.size() > 0);
      CGAL_precondition(orientations.size() > 0);
      // parallel_groups_by_angles.clear();
      collinear_groups_by_angles.clear();
      
      const std::size_t n = m_input_range.size();
      std::map<std::size_t, std::vector<std::size_t>> collinear_groups;
      std::vector<int> segments_to_groups_hashmap(n, -1);

      build_eigen_matrix(t_ijs);

      std::vector<std::size_t> vec;
      std::size_t g = 0, p = 0;
      int g_i, g_j;

      for (std::size_t k = 0; k < m_targets.outerSize(); ++k) {
        for (typename Sparse_matrix_FT::InnerIterator it_tar(m_targets, k); it_tar; ++it_tar) {

          const std::size_t i = it_tar.row();
          const std::size_t j = it_tar.col();

          if (CGAL::abs(orientations[n + p]) < tolerance) { 

            // Then segments i and j belong to the same group of parallel segments.
            // For the moment, these groups are materialized by integers.
            if (segments_to_groups_hashmap[i] == -1 && segments_to_groups_hashmap[j] == -1) {

              // Segments i and j are not assigned to any group of parallel segments
              // So we create one with them.
              segments_to_groups_hashmap[i] = g;
              segments_to_groups_hashmap[j] = g;

              collinear_groups[g].push_back(i);
              collinear_groups[g].push_back(j);

              // nodes_to_groups[m_input_segments[i]->parallel_node].push_back(g);
              ++g;
            }
            else if (segments_to_groups_hashmap[i] == -1 && segments_to_groups_hashmap[j] != -1) {

              // Assigns segment i to the group of the segment j.
              g_j = segments_to_groups_hashmap[j];
              segments_to_groups_hashmap[i] = g_j;
              collinear_groups[g_j].push_back(i);
            }
            else if (segments_to_groups_hashmap[i] != -1 && segments_to_groups_hashmap[j] == -1) {

              // Assigns segment j to the group of the segment i.
              g_i = segments_to_groups_hashmap[i];
              segments_to_groups_hashmap[j] = g_i;
              collinear_groups[g_i].push_back(j);
            }
            else {
              g_i = segments_to_groups_hashmap[i];
              g_j = segments_to_groups_hashmap[j];
              if (g_i != g_j) {
                     
                // Segments i and j have been assigned to different groups, but in fact
                // they belong to the same group. That's why we merge them.
                for (const auto gr : collinear_groups[g_j]) {
                  segments_to_groups_hashmap[gr] = g_i;
                  collinear_groups[g_i].push_back(gr);
                }
                collinear_groups[g_j].clear();

              }
            }
          }
          ++p;
        }
      }
      CGAL_postcondition(collinear_groups.size() > 0);

      std::map<int, FT> ordinates;

      // build_map_of_ordinates(orientations, collinear_groups, 
                          // segments_to_groups_hashmap, ordinates);

      // Try to assign segments whose orientation has not been optimized thanks to the regularization process, to an existing group.
      // assign_segments_to_groups(parallel_groups, segments_to_groups_hashmap, angles);

      // build_parallel_groups_angle_map(segments_to_groups_hashmap, angles, parallel_groups_by_angles);
      // */

    }

  private:
    const Input_range& m_input_range;
    const std::map <std::size_t, Segment_data> & m_segments;
    Sparse_matrix_FT  m_targets;
    const FT theta_eps;
    const FT y_eps;
    const FT tolerance;

    void build_eigen_matrix(const std::map <std::pair<std::size_t, std::size_t>, FT> & t_ijs) {

      std::vector<FT_triplet> vec_targets;
      FT_triplet t_ij_triplet;
      for (const auto& it : t_ijs) {
        t_ij_triplet = FT_triplet(it.first.first, it.first.second, it.second);
        vec_targets.push_back(t_ij_triplet);
      }
      CGAL_postcondition(vec_targets.size() == t_ijs.size());

      const std::size_t n = m_input_range.size();
      m_targets.resize(n, n);
      m_targets.setFromTriplets(vec_targets.begin(), vec_targets.end());
      m_targets.makeCompressed();
      CGAL_postcondition(m_targets.nonZeros() == t_ijs.size());

    }

    void build_map_of_ordinates(const std::vector<FT> & orientations,
                             std::vector<std::vector<std::size_t>> & collinear_groups,
                             std::vector<int> & segments_to_groups_hashmap,
                             std::map<int, FT> & ordinates) {
      /*
      for (std::size_t i = 0; i < segments_to_groups_hashmap.size(); ++i) {
        const int g_i = segments_to_groups_hashmap[i];

        if (g_i != -1 && (ordinates.find(g_i) == ordinates.end())) {

          FT y = m_segments[i].m_reference_coordinates + orientations[i];
          
          // Check if the angle that seems to be associated to this group of segments is not too close to another value.
          int g_j = -1;
          for (auto it_m = ordinates.begin(); it_m != ordinates.end(); ++it_m) {
            if (CGAL::abs(it_m->second - y) < y_eps) {

            }
              g_j = it_angle->first;
          }

          if (g_j == -1) 
            ordinates[g_i] = theta;
          else {                       
            // Merge groups.
            for (const auto gr : collinear_groups[g_i]) {
              segments_to_groups_hashmap[gr] = g_j;
              collinear_groups[g_j].push_back(gr);
            }
            collinear_groups[g_i].clear();
          }

        }
      }
      */
    }

    void assign_segments_to_groups(std::vector<std::vector<std::size_t>> & parallel_groups,
                                   std::vector<int> & segments_to_groups_hashmap,
                                   std::map<int, FT> & angles) {

      for (std::size_t i = 0; i < segments_to_groups_hashmap.size(); ++i) {
        int g_i = segments_to_groups_hashmap[i];
        if (g_i == -1) {
          const FT alpha = m_segments[i].m_orientation;
          int g_j = -1;
          for (auto it_angle = angles.begin(); it_angle != angles.end(); ++it_angle) {
            const FT alpha_j = it_angle->second;
            for (int k = -1; k <= 1; ++k) {
              if (CGAL::abs(alpha_j - alpha + static_cast<FT>(k) * FT(180)) < theta_eps) {
                g_j = it_angle->first;
                break;
              }
            }
            if (g_j != -1) 
              break;
          }
          if (g_j == -1) {                  
            g_i = angles.rbegin()->first + 1;
            angles[g_i] = alpha;
          } else g_i = g_j;
          segments_to_groups_hashmap[i] = g_i; // Segmentation fault (core dumped) for the example with 3 segments and bounds = 25
          parallel_groups[g_i].push_back(i);
        }
      }

    }

    void build_parallel_groups_angle_map(const std::vector<int> & segments_to_groups_hashmap,
                                         std::map<int, FT> & angles, 
                                         std::map<FT, std::vector<std::size_t>> & parallel_groups_by_angles) {

      for (auto it_angle = angles.begin(); it_angle != angles.end(); ++it_angle) {
        const FT angle = angles[it_angle->first];
        if (parallel_groups_by_angles.find(angle) == parallel_groups_by_angles.end()) 
          parallel_groups_by_angles[angle] = std::vector<std::size_t>();
      }

      for (std::size_t i = 0; i < segments_to_groups_hashmap.size(); ++i) {
        // If segment s_i is included in a group of parallel segments,
        // then it should be assigned to a leaf of the regularization tree.
        const FT angle = angles[segments_to_groups_hashmap[i]];
        if (parallel_groups_by_angles.find(angle) != parallel_groups_by_angles.end()) 
          parallel_groups_by_angles[angle].push_back(i);
      }

    }

  };

} // namespace internal
} // namespace Regularization
} // namespace CGAL

#endif // CGAL_SHAPE_REGULARIZATION_INTERNAL_GROUPING_SEGMENTS_ORDINATES_2