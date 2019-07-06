#ifndef CGAL_SHAPE_REGULARIZATION_INTERNAL_TREE_H
#define CGAL_SHAPE_REGULARIZATION_INTERNAL_TREE_H

// #include <CGAL/license/Shape_regularization.h>

#include <map>
#include <list>
#include <vector>
#include <cassert>

// Eigen includes.
#include <eigen3/Eigen/SparseCore>

namespace CGAL {
namespace Regularization {
namespace internal {

  template<
    typename GeomTraits, 
    typename InputRange>
  class Tree {
  public:
    using Traits = GeomTraits;
    using Input_range = InputRange;
    using FT = typename GeomTraits::FT;
    using Point = typename GeomTraits::Point_2;
    using Segment = typename GeomTraits::Segment_2;
    using Segments = std::vector<Segment>;
    using Orientations = std::vector<FT>;
    using Vector  = typename GeomTraits::Vector_2;

    using List_element  = std::list<int>;
    using List_iterator = typename List_element::const_iterator;

    using Segments_to_groups = std::vector<int>;
    using Groups_to_segments = std::map<int, List_element>;
    using Parallel_segments          = std::map<FT, Segments>;
    using Parallel_segments_iterator       = typename Parallel_segments::iterator;

    using Targets_matrix   = Eigen::SparseMatrix<FT,  Eigen::RowMajor>;
    using Relations_matrix = Eigen::SparseMatrix<int, Eigen::RowMajor>;
    using FT_triplet  = Eigen::Triplet<FT>;
    using Int_triplet = Eigen::Triplet<int>;
    using Mus       = std::vector<FT_triplet>;
    using Targets   = std::vector<FT_triplet>;
    using Relations = std::vector<Int_triplet>;

    using Targets_iterator   = typename Targets_matrix::InnerIterator;
    using Relations_iterator = typename Relations_matrix::InnerIterator;

    using Angles          = std::map<int, FT>;
    using Angles_iterator = typename Angles::const_iterator;

    Tree(
      InputRange& input_range,
      const std::map <std::pair<std::size_t, std::size_t>, FT> t_ijs,
      const std::map <std::pair<std::size_t, std::size_t>, FT> r_ijs,
      const FT mu_ij,
      const Orientations &orientations) :
    m_input_range(input_range),
    m_t_ijs(t_ijs),
    m_r_ijs(r_ijs),
    m_mu_ij(mu_ij),
    m_orientations(orientations) {

      CGAL_precondition(input_range.size() > 0);

      build_tree();

    }

    void apply_new_orientations() {
      for (Parallel_segments_iterator it_ps = m_parallel_segments.begin(); it_ps != m_parallel_segments.end(); ++it_ps) {
        const FT theta = it_ps->first;
        const Segments &subtree = it_ps->second;

        // Each group of parallel segments has a normal vector that we compute with alpha.
        const FT x = static_cast<FT>(cos(CGAL::to_double(theta * static_cast<FT>(CGAL_PI) / FT(180))));
        const FT y = static_cast<FT>(sin(CGAL::to_double(theta * static_cast<FT>(CGAL_PI) / FT(180))));

        Vector v_dir = Vector(x, y);
        const Vector v_ort = Vector(-v_dir.y(), v_dir.x());
        
        const FT a = v_ort.x();
        const FT b = v_ort.y();

        // Rotate segments with precision.
        for (int i = 0; i < subtree.size(); ++i) {
          Segment segment_pointer = subtree[i];

          // Compute equation of the supporting line of the rotated segment.
          const Point &barycentre = compute_barycentre(segment_pointer);
          const FT c = -a * barycentre.x() - b * barycentre.y();

          int j = find_segment(subtree[i]);

          if(j >= 0) {
            Vector v_j = compute_direction(j);
            set_orientation(j, theta - compute_orientation(v_j), a, b, c, v_dir);
          }
          
        }

      }
    }

  private:
    Input_range& m_input_range;
    std::map <std::pair<std::size_t, std::size_t>, FT> m_t_ijs;
    std::map <std::pair<std::size_t, std::size_t>, FT> m_r_ijs;
    const FT m_mu_ij;
    const Orientations     &m_orientations;
    const FT   m_tolerance = FT(1) / FT(1000000);
    Parallel_segments m_parallel_segments;
    // Mus_matrix       m_mus_matrix;
    // Targets_matrix   m_targets_matrix;
    // Relations_matrix m_relations_matrix;

    Point compute_barycentre(const Segment& m_segment) {
      const FT half = FT(1) / FT(2);
      const Point &source = m_segment.source();
      const Point &target = m_segment.target();
      const FT x = half * (source.x() + target.x());
      const FT y = half * (source.y() + target.y());
      return Point(x, y);
    }

    Vector compute_direction(const int i) {
      Vector v = m_input_range[i].to_vector(); 
      if (v.y() < FT(0) || (v.y() == FT(0) && v.x() < FT(0))) 
        v = -v;
      return v;
    }
    
    FT compute_orientation(Vector v) {
      const FT atan = static_cast<FT>(std::atan2(CGAL::to_double(v.y()), CGAL::to_double(v.x())));
      FT orientation = atan * FT(180) / static_cast<FT>(CGAL_PI);
      if (orientation < FT(0)) 
        orientation += FT(180);
      return orientation;
    }

    int find_segment(Segment seg) {
      for (int i = 0; i < m_input_range.size(); i++) {
        if (seg == m_input_range[i])
          return i;
      }
      return -1;
    }

    void set_orientation(int i, const FT new_orientation, const FT a, const FT b, const FT c, const Vector &direction) {
      FT m_orientation = new_orientation;
      FT m_a = a;
      FT m_b = b;
      FT m_c = c;
      Vector m_direction = direction;
      if (m_direction.y() < FT(0) || (m_direction.y() == FT(0) && m_direction.x() < FT(0))) 
        m_direction = -m_direction;
      FT x1, y1, x2, y2;
      Point m_barycentre = compute_barycentre(m_input_range[i]);
      FT m_length = static_cast<FT>(CGAL::sqrt(CGAL::to_double(m_input_range[i].squared_length())));
      if (CGAL::abs(m_direction.x()) > CGAL::abs(m_direction.y())) { 
        x1 = m_barycentre.x() - m_length * m_direction.x() / FT(2);
        x2 = m_barycentre.x() + m_length * m_direction.x() / FT(2);

        y1 = (-m_c - m_a * x1) / m_b;
        y2 = (-m_c - m_a * x2) / m_b;
      } else {
        y1 = m_barycentre.y() - m_length * m_direction.y() / FT(2);
        y2 = m_barycentre.y() + m_length * m_direction.y() / FT(2);

        x1 = (-m_c - m_b * y1) / m_a;
        x2 = (-m_c - m_b * y2) / m_a;
      }
      const Point source = Point(x1, y1);
      const Point target = Point(x2, y2);

      m_input_range[i] = Segment(source, target);

    }

    void build_tree() {
      const int n = static_cast<int>(m_input_range.size());
      Segments_to_groups segments_to_groups(n, -1);
      Groups_to_segments groups_to_segments;

      // Mus mus;
      Targets targets;
      Relations relations;
      for (typename std::map<std::pair<std::size_t, std::size_t>, FT>::iterator it = m_t_ijs.begin(); it!=m_t_ijs.end(); ++it) {
        targets.push_back(FT_triplet(it->first.first, it->first.second, it->second));
      }
      for (typename std::map<std::pair<std::size_t, std::size_t>, FT>::iterator it = m_r_ijs.begin(); it!=m_r_ijs.end(); ++it) {
        relations.push_back(Int_triplet(it->first.first, it->first.second, it->second));
        // mus.push_back(FT_triplet(it->first.first, it->first.second, m_mu_ij));
      }

      Targets_matrix targets_matrix;
      Relations_matrix relations_matrix;
      // Mus_matrix mus_matrix;

      const size_t global_size = m_input_range.size();

      targets_matrix.resize(global_size, global_size);
      targets_matrix.setFromTriplets(targets.begin(), targets.end());
      targets_matrix.makeCompressed();

      relations_matrix.resize(global_size, global_size);
      relations_matrix.setFromTriplets(relations.begin(), relations.end());
      relations_matrix.makeCompressed();

    /*  mus_matrix.resize(global_size, global_size);
      mus_matrix.setFromTriplets(mus.begin(), mus.end());
      mus_matrix.makeCompressed(); */

      const FT theta_eps = FT(1) / FT(4);

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

    void create_parallel_node(const FT angle) {
      // Segment parallel_segment;
      if (m_parallel_segments.find(angle) == m_parallel_segments.end()) 
        m_parallel_segments[angle] = Segments();
      // parallel_segments  = Parallel_segments();
    }

    void assign_to_parallel_node(const FT angle, Segment segment_pointer) {
      if (m_parallel_segments.find(angle) != m_parallel_segments.end())
        m_parallel_segments[angle].push_back(segment_pointer);
    }

  };

} // namespace internal
} // namespace Regularization
} // namespace CGAL

#endif // CGAL_SHAPE_REGULARIZATION_INTERNAL_TREE_H
