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
                     const FT mu_ij, const std::vector<FT> & orientations,
                     std::vector<std::vector<std::size_t>> & result) { 

      CGAL_precondition(t_ijs.size() > 0);
      CGAL_precondition(r_ijs.size() > 0);
      CGAL_precondition(orientations.size() > 0);
      result.clear();

      const std::size_t n = m_segments.size();
      result.resize(n);
      build_eigen_matrices(t_ijs, r_ijs);
      std::vector<int> segments_to_groups_hashmap(n, -1);
      // std::map<int, std::vector<int>> groups_to_segments;

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
                  result[g].push_back(i);
                  result[g].push_back(j);
                  ++g;
                  break;
                case 1:
                  // The segments i and j are orthogonal.
                  // We create two different groups of parallel segments.
                  segments_to_groups_hashmap[i] = g;
                  result[g].push_back(i);
                  segments_to_groups_hashmap[j] = ++g;
                  result[g].push_back(j);
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
                  result[g_j].push_back(i);
                  break;
                case 1:
                  // Then segment i is orthogonal to j, and we should initialize a new group with this segment.
                  segments_to_groups_hashmap[i] = g;
                  result[g].push_back(i);
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
                  result[g_i].push_back(j);
                  break;
                case 1:
                  segments_to_groups_hashmap[j] = g;
                  result[g].push_back(j);
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
                  for (const auto gr : result[g_j]) {
                    segments_to_groups_hashmap[gr] = g_i;
                    result[g_i].push_back(gr);
                  }
                  result[g_j].clear();
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
      CGAL_postcondition(result.size() > 0);
      CGAL_postcondition(g <= n);
      result.resize(g);


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



    /*
    using List_element  = std::list<int>;
    void build_tree() {
      const int n = static_cast<int>(m_input_range.size());
      std::vector<int> segments_to_groups(n, -1); //using Segments_to_groups = std::vector<int>;
      std::map<int, List_element> groups_to_segments; // using Groups_to_segments = std::map<int, List_element>;

      // Categorize segments.
      int g = 0, p = 0;
      for (int k = 0; k < targets_matrix.outerSize(); ++k) {
        Targets_iterator     it_targets(  targets_matrix, k);
        Relations_iterator it_relations(relations_matrix, k);

        while (it_targets && it_relations) {
          const int i = it_targets.row();
          const int j = it_targets.col();
          const int r = it_relations.value();

          if (CGAL::abs(m_orientations[n + p]) < m_tolerance) {
            // case-->
            if (segments_to_groups[i] == -1 && segments_to_groups[j] == -1) {
              if (r == 0) { 
                // Then segments i and j belong to the same group of parallel segments.
                // We should create a group of segments, that is initialized with these two individuals.
                segments_to_groups[i] = segments_to_groups[j] = g;
                groups_to_segments[g].push_back(i);
                groups_to_segments[g].push_back(j);
                ++g;
              } else if (r == 1) {              
                // The segments i and j are orthogonal.
                // We create two different groups of parallel segments.
                segments_to_groups[i] = g;
                groups_to_segments[g].push_back(i);
                segments_to_groups[j] = ++g;
                groups_to_segments[g].push_back(j);
                ++g;
              }
            }
            // case--> 
            else if (segments_to_groups[i] == -1 && segments_to_groups[j] != -1) {
              if (r == 0) {
                // Then segment i is parallel to j, and can be assigned to the same group.
                const int g_j = segments_to_groups[j];
                segments_to_groups[i] = g_j;
                groups_to_segments[g_j].push_back(i);
              } else if (r == 1) {               
                // Then segment i is orthogonal to j, and we should initialize a new group with this segment.
                segments_to_groups[i] = g;
                groups_to_segments[g].push_back(i);
                ++g;
              }
            }
            // case-->
            else if (segments_to_groups[i] != -1 && segments_to_groups[j] == -1) {
              // Symmetrical situation to before.
              if (r == 0) {
                const int g_i = segments_to_groups[i];
                segments_to_groups[j] = g_i;
                groups_to_segments[g_i].push_back(j);
              } else if (r == 1) {
                segments_to_groups[j] = g;
                groups_to_segments[g].push_back(j);
                ++g;
              }
            }
            // case-->
            else {
              const int g_i = segments_to_groups[i];
              const int g_j = segments_to_groups[j];
              if (g_i != g_j) {
                if (r == 0) {                       
                  // Segments i and j have been assigned to different groups, but in fact
                  // they are parallel and belong to the same group. That's why we merge them.
                  for (List_iterator it_list = groups_to_segments[g_j].begin(); it_list != groups_to_segments[g_j].end(); ++it_list) {
                    segments_to_groups[*it_list] = g_i;
                    groups_to_segments[g_i].push_back(*it_list);
                  }
                  groups_to_segments[g_j].clear();
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

      // Prepare for construction of the regularization tree.
      Angles angles;

      for (size_t i = 0; i < segments_to_groups.size(); ++i) {
        const int g_i = segments_to_groups[i];

        if (g_i != -1) {
          if (angles.find(g_i) == angles.end()) {
            Vector v_i = compute_direction(i);
            FT theta = compute_orientation(v_i) + m_orientations[i];

            if (theta < FT(0)) 
              theta += FT(180);
            else if (theta > FT(180)) 
              theta -= FT(180);
            
            // Check if the angle that seems to be associated to this group of segments is not too close to another value.
            int g_j = -1;
            for (Angles_iterator it_angle = angles.begin(); it_angle != angles.end(); ++it_angle)
              if (CGAL::abs(it_angle->second - theta) < theta_eps) 
                g_j = it_angle->first;

            if (g_j == -1) 
              angles[g_i] = theta;
            else {                       
              // Merge groups.
              for (List_iterator it_list = groups_to_segments[g_i].begin(); it_list != groups_to_segments[g_i].end(); ++it_list) {    
                segments_to_groups[*it_list] = g_j;
                groups_to_segments[g_j].push_back(*it_list);
              }
              groups_to_segments[g_i].clear();
            }
          }
        }
      }

      // Try to assign segments whose orientation has not been optimized thanks to the regularization process, to an existing group.
      for (size_t i = 0; i < segments_to_groups.size(); ++i) {
        int g_i = segments_to_groups[i];
        if (g_i == -1) {
          Vector v_i = compute_direction(i);
          const FT alpha = compute_orientation(v_i);
          int g_j = -1;

          for (Angles_iterator it_angle = angles.begin(); it_angle != angles.end(); ++it_angle) {
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

      // Build regularization tree.

      for (Angles_iterator it_angle = angles.begin(); it_angle != angles.end(); ++it_angle) 
        create_parallel_node(angles[it_angle->first]);

      for (size_t i = 0; i < segments_to_groups.size(); ++i) {
          
        // If segment s_i is included in a group of parallel segments,
        // then it should be assigned to a leaf of the regularization tree.
        assign_to_parallel_node(angles[segments_to_groups[i]], m_input_range[i]);
      }


    }
    */

  };

} // namespace internal
} // namespace Regularization
} // namespace CGAL

#endif // CGAL_SHAPE_REGULARIZATION_INTERNAL_GROUPING_SEGMENTS_2