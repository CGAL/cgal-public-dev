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

    void make_groups(const std::map <std::pair<std::size_t, std::size_t>, std::pair<FT, std::size_t>> & t_ijs,
                     const std::map <std::size_t, Segment_data> & temp_segments,
                     const std::vector<FT> & orientations,
                     std::map<FT, std::vector<std::size_t>> & collinear_groups_by_ordinates) { 

      CGAL_precondition(t_ijs.size() > 0);
      CGAL_precondition(orientations.size() > 0);
      // parallel_groups_by_angles.clear();
      collinear_groups_by_ordinates.clear();
      
      const std::size_t n = m_input_range.size();
      std::map<std::size_t, std::vector<std::size_t>> collinear_groups;
      std::map<std::size_t, int> segments_to_groups_hashmap;

      for (const auto & it : temp_segments) {
        std::size_t seg_index = it.second.m_index;
        segments_to_groups_hashmap[seg_index] = -1;
      }
      // std::vector<int> segments_to_groups_hashmap(n, -1);

      // build_eigen_matrix(t_ijs);


      // std::vector<std::size_t> vec;
      // std::size_t g = 0, p = 0;
      // int g_i, g_j;

      std::size_t g = 0;
      int g_i, g_j;

      for (const auto & tar_it : t_ijs) {

        const std::size_t i = tar_it.first.first;
        const std::size_t j = tar_it.first.second;
        const std::size_t p = tar_it.second.second;

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

      }

      // CGAL_postcondition(collinear_groups.size() > 0);
      // if (collinear_groups.size() == 0) {
      //   collinear_groups[0];
      // }

      std::map<int, FT> ordinates;

      build_map_of_ordinates(orientations, temp_segments, collinear_groups, 
                          segments_to_groups_hashmap, ordinates);

      // Try to assign segments whose orientation has not been optimized thanks to the regularization process, to an existing group.
      assign_segments_to_groups(temp_segments, collinear_groups, segments_to_groups_hashmap, ordinates);

      build_collinear_groups_ordinate_map(segments_to_groups_hashmap, ordinates, collinear_groups_by_ordinates);
      // */

    }

  private:
    const Input_range& m_input_range;
    const std::map <std::size_t, Segment_data> & m_segments;
    Sparse_matrix_FT  m_targets;
    const FT theta_eps;
    const FT y_eps;
    const FT tolerance;

    void build_eigen_matrix(const std::map <std::pair<std::size_t, std::size_t>, std::pair<FT, std::size_t>> & t_ijs) {

      if(m_targets.nonZeros() > 0)
        m_targets.resize(0, 0);

      std::vector<FT_triplet> vec_targets;
      FT_triplet t_ij_triplet;
      for (const auto& it : t_ijs) {
        t_ij_triplet = FT_triplet(it.first.first, it.first.second, it.second.first);
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
                             const std::map <std::size_t, Segment_data> & temp_segments,
                             std::map<std::size_t, std::vector<std::size_t>> & collinear_groups,
                             std::map<std::size_t, int> & segments_to_groups_hashmap,
                             std::map<int, FT> & ordinates) {

      for (const auto & sm_i : segments_to_groups_hashmap) {
        const int g_i = sm_i.second;

        if (g_i != -1 && (ordinates.find(g_i) == ordinates.end())) {
          const std::size_t seg_index = sm_i.first;
          const Segment_data & seg_data = temp_segments.at(seg_index);
          const FT y = seg_data.m_reference_coordinates.y() + orientations[seg_index];
          
          // Check if the angle that seems to be associated to this group of segments is not too close to another value.
          int g_j = -1;
          for (const auto & it_m : ordinates) {
            if (CGAL::abs(it_m.second - y) < y_eps) {
              g_j = it_m.first;
            }
          }
          if (g_j == -1) 
            ordinates[g_i] = y;
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

    }

    void assign_segments_to_groups(const std::map <std::size_t, Segment_data> & temp_segments,
                                   std::map<std::size_t, std::vector<std::size_t>> & collinear_groups,
                                   std::map<std::size_t, int> & segments_to_groups_hashmap,
                                   std::map<int, FT> & ordinates) {

      for (const auto & sm_i : segments_to_groups_hashmap) {
        int g_i = sm_i.second;

        if (g_i == -1) {
          const std::size_t seg_index = sm_i.first;
          const Segment_data & seg_data = temp_segments.at(seg_index);
          const FT y = seg_data.m_reference_coordinates.y();

          int g_j = -1;
          for (const auto & it_m : ordinates) {
            const FT y_j = it_m.second;
            if (CGAL::abs(y_j - y) < y_eps) {
              g_j = it_m.first;
            }
            if (g_j != -1) 
              break;
          }
          if (g_j == -1) {                  
            g_i = ordinates.rbegin()->first + 1;
            ordinates[g_i] = y;
          } 
          else 
            g_i = g_j;

          segments_to_groups_hashmap[seg_index] = g_i; 
          collinear_groups[g_i].push_back(seg_index);
        }
      }

    }

    void build_collinear_groups_ordinate_map(const std::map<std::size_t, int> & segments_to_groups_hashmap,
                                         const std::map<int, FT> & ordinates, 
                                         std::map<FT, std::vector<std::size_t>> & collinear_groups_by_ordinates) {

      for (const auto & it_m : ordinates) {
        const FT y = it_m.second;
        if (collinear_groups_by_ordinates.find(y) == collinear_groups_by_ordinates.end()) 
          collinear_groups_by_ordinates[y] = std::vector<std::size_t>();
      }

      for (const auto & sm_i : segments_to_groups_hashmap) {
        const FT y = ordinates.at(sm_i.second);
        if (collinear_groups_by_ordinates.find(y) != collinear_groups_by_ordinates.end()) 
          collinear_groups_by_ordinates[y].push_back(sm_i.first);
      }

    }

  };

} // namespace internal
} // namespace Regularization
} // namespace CGAL

#endif // CGAL_SHAPE_REGULARIZATION_INTERNAL_GROUPING_SEGMENTS_ORDINATES_2