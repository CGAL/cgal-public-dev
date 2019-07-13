#ifndef CGAL_SHAPE_REGULARIZATION_INTERNAL_GROUPING_SEGMENTS_2
#define CGAL_SHAPE_REGULARIZATION_INTERNAL_GROUPING_SEGMENTS_2

// #include <CGAL/license/Shape_regularization.h>

#include <vector>
#include <map>
#include <utility>

#include <CGAL/Shape_regularization/internal/Segment_data_2.h>
// #include <CGAL/Shape_regularization/internal/utils.h>

#include <eigen3/Eigen/SparseCore>

namespace CGAL {
namespace Regularization {
namespace internal {

  template<
    typename GeomTraits>
  class Grouping_segments_2 {
  public:
    using Traits = GeomTraits;
    using FT = typename GeomTraits::FT;
    // using Segment = typename GeomTraits::Segment_2;
    using Segment_data = typename internal::Segment_data_2<Traits>;
    using Sparse_matrix_FT  = Eigen::SparseMatrix<FT,  Eigen::RowMajor>;
    using Sparse_matrix_int = Eigen::SparseMatrix<int, Eigen::RowMajor>;
    using FT_triplet  = Eigen::Triplet<FT>;
    using Int_triplet = Eigen::Triplet<int>;

    Grouping_segments_2(
      const std::vector<Segment_data> & segments) :
    m_segments(segments),
    theta_eps(FT(1) / FT(4)),
    tolerance(FT(1) / FT(1000000)) {

      CGAL_precondition(segments.size() > 0);

    }

    void make_groups(const std::map <std::pair<std::size_t, std::size_t>, FT> & t_ijs,
                     const std::map <std::pair<std::size_t, std::size_t>, FT> & r_ijs,
                     const FT mu_ij, const std::vector<FT> & orientations ,
                     std::map<FT, std::vector<std::size_t>> & parallel_groups_by_angles) { 

      CGAL_precondition(t_ijs.size() > 0);
      CGAL_precondition(r_ijs.size() > 0);
      CGAL_precondition(orientations.size() > 0);
      parallel_groups_by_angles.clear();
      
      const std::size_t n = m_segments.size();
      std::vector<std::vector<std::size_t>> parallel_groups;
      std::vector<int> segments_to_groups_hashmap(n, -1);
      parallel_groups.resize(n);

      build_eigen_matrices(t_ijs, r_ijs);

      std::vector<std::size_t> vec;
      std::size_t g = 0, p = 0;
      int g_i, g_j;
      for (std::size_t k = 0; k < m_targets.outerSize(); ++k) {
        typename Sparse_matrix_FT::InnerIterator  it_targets(m_targets, k);
        typename Sparse_matrix_int::InnerIterator it_relations(m_relations, k);

        while (it_targets && it_relations) {
          const std::size_t i = it_targets.row();
          const std::size_t j = it_targets.col();
          const std::size_t r = it_relations.value();

          if (CGAL::abs(orientations[n + p]) < tolerance) { 
            if (segments_to_groups_hashmap[i] == -1 && segments_to_groups_hashmap[j] == -1) {
              switch (r) {
                case 0:
                  // Then segments i and j belong to the same group of parallel segments.
                  // We should create a group of segments, that is initialized with these two individuals.
                  segments_to_groups_hashmap[i] = g;
                  segments_to_groups_hashmap[j] = g;
                  parallel_groups[g].push_back(i);
                  parallel_groups[g].push_back(j);
                  ++g;
                  break;
                case 1:
                  // The segments i and j are orthogonal.
                  // We create two different groups of parallel segments.
                  segments_to_groups_hashmap[i] = g;
                  parallel_groups[g].push_back(i);
                  segments_to_groups_hashmap[j] = ++g;
                  parallel_groups[g].push_back(j);
                  ++g;
                  break;
              }
            }
            else if (segments_to_groups_hashmap[i] == -1 && segments_to_groups_hashmap[j] != -1) {
              switch (r) {
                case 0:
                  // Then segment i is parallel to j, and can be assigned to the same group.
                  g_j = segments_to_groups_hashmap[j];
                  segments_to_groups_hashmap[i] = g_j;
                  parallel_groups[g_j].push_back(i);
                  break;
                case 1:
                  // Then segment i is orthogonal to j, and we should initialize a new group with this segment.
                  segments_to_groups_hashmap[i] = g;
                  parallel_groups[g].push_back(i);
                  ++g;
                  break;
              }
            }
            else if (segments_to_groups_hashmap[i] != -1 && segments_to_groups_hashmap[j] == -1) {
              // Symmetrical situation to before.
              switch (r) {
                case 0:
                  g_i = segments_to_groups_hashmap[i];
                  segments_to_groups_hashmap[j] = g_i;
                  parallel_groups[g_i].push_back(j);
                  break;
                case 1:
                  segments_to_groups_hashmap[j] = g;
                  parallel_groups[g].push_back(j);
                  ++g;
                  break;
              }
            }
            else {
              g_i = segments_to_groups_hashmap[i];
              g_j = segments_to_groups_hashmap[j];
              if (g_i != g_j) {
                if (r == 0) {                       
                  // Segments i and j have been assigned to different groups, but in fact
                  // they are parallel and belong to the same group. That's why we merge them.
                  for (const auto gr : parallel_groups[g_j]) {
                    segments_to_groups_hashmap[gr] = g_i;
                    parallel_groups[g_i].push_back(gr);
                  }
                  parallel_groups[g_j].clear();
                } else if (r == 1) {
                  // We do nothing here.
                }
              }
            }
          }
          ++p;
          ++it_targets;
          ++it_relations;
        }
      }
      CGAL_postcondition(parallel_groups.size() > 0);
      CGAL_postcondition(g <= n);
      parallel_groups.resize(g);

      std::map<int, FT> angles;

      for (std::size_t i = 0; i < segments_to_groups_hashmap.size(); ++i) {

        g_i = segments_to_groups_hashmap[i];
        if (g_i != -1 && (angles.find(g_i) == angles.end())) {
          FT theta = m_segments[i].m_orientation + orientations[i];
          if (theta < FT(0)) 
            theta += FT(180);
          else if (theta > FT(180)) 
            theta -= FT(180);
          
          // Check if the angle that seems to be associated to this group of segments is not too close to another value.
          g_j = -1;
          for (typename std::map<int, FT>::iterator it_angle = angles.begin(); it_angle != angles.end(); ++it_angle) {
            if (CGAL::abs(it_angle->second - theta) < theta_eps) 
              g_j = it_angle->first;
          }

          if (g_j == -1) 
            angles[g_i] = theta;
          else {                       
            // Merge groups.
            for (const auto gr : parallel_groups[g_i]) {
              segments_to_groups_hashmap[gr] = g_j;
              parallel_groups[g_j].push_back(gr);
            }
            parallel_groups[g_i].clear();
          }
        }
      }

      // Try to assign segments whose orientation has not been optimized thanks to the regularization process, to an existing group.
      for (std::size_t i = 0; i < segments_to_groups_hashmap.size(); ++i) {
        g_i = segments_to_groups_hashmap[i];
        if (g_i == -1) {
          const FT alpha = m_segments[i].m_orientation;
          g_j = -1;
          for (typename std::map<int, FT>::iterator it_angle = angles.begin(); it_angle != angles.end(); ++it_angle) {
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
        }
      }

      for (std::size_t i = 0; i < segments_to_groups_hashmap.size(); ++i) {
        // If segment s_i is included in a group of parallel segments,
        // then it should be assigned to a leaf of the regularization tree.
        const FT angle = angles[segments_to_groups_hashmap[i]];
        parallel_groups_by_angles[angle].push_back(i);
      }

    }

  private:
    const std::vector<Segment_data> & m_segments;
    Sparse_matrix_FT  m_targets;
    Sparse_matrix_int m_relations;
    const FT theta_eps;
    const FT tolerance;

    void build_eigen_matrices(const std::map <std::pair<std::size_t, std::size_t>, FT> & t_ijs,
                              const std::map <std::pair<std::size_t, std::size_t>, FT> & r_ijs) {

      std::vector<FT_triplet> vec_targets;
      FT_triplet t_ij_triplet;
      for (const auto& it : t_ijs) {
        t_ij_triplet = FT_triplet(it.first.first, it.first.second, it.second);
        vec_targets.push_back(t_ij_triplet);
      }
      CGAL_postcondition(vec_targets.size() == t_ijs.size());

      std::vector<Int_triplet> vec_relations;
      Int_triplet r_ij_triplet;
      for (const auto& it : r_ijs) {
        r_ij_triplet = Int_triplet(it.first.first, it.first.second, it.second);
        vec_relations.push_back(r_ij_triplet);
      }
      CGAL_postcondition(vec_relations.size() == r_ijs.size());

      const std::size_t n = m_segments.size();
      m_targets.resize(n, n);
      m_targets.setFromTriplets(vec_targets.begin(), vec_targets.end());
      m_targets.makeCompressed();
      CGAL_postcondition(m_targets.nonZeros() == t_ijs.size());

      m_relations.resize(n, n);
      m_relations.setFromTriplets(vec_relations.begin(), vec_relations.end());
      m_relations.makeCompressed();
      CGAL_postcondition(m_relations.nonZeros() == r_ijs.size());

    }

  };

} // namespace internal
} // namespace Regularization
} // namespace CGAL

#endif // CGAL_SHAPE_REGULARIZATION_INTERNAL_GROUPING_SEGMENTS_2