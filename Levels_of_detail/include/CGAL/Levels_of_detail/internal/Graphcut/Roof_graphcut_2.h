// Copyright (c) 2019 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is a part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY, AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
//
//
// Author(s)     : Dmitry Anisimov, Simon Giraudot, Pierre Alliez, Florent Lafarge, and Andreas Fabri
//

#ifndef CGAL_LEVELS_OF_DETAIL_INTERNAL_ROOF_GRAPHCUT_2_H
#define CGAL_LEVELS_OF_DETAIL_INTERNAL_ROOF_GRAPHCUT_2_H

// CGAL includes.
#include <CGAL/assertions.h>

#define CGAL_DO_NOT_USE_BOYKOV_KOLMOGOROV_MAXFLOW_SOFTWARE
#include <CGAL/internal/Surface_mesh_segmentation/Alpha_expansion_graph_cut.h>

// LOD includes.
#include <CGAL/Levels_of_detail/enum.h>

// Internal includes.
#include <CGAL/Levels_of_detail/internal/utils.h>

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

  template<typename GeomTraits>
  class Roof_graphcut_2 {

  public:
    using Traits = GeomTraits;
    
		using Partition_2 = internal::Partition_2<Traits>;

    using FT = typename Traits::FT;

    using Face = typename Partition_2::Face;
    using Edge = typename Partition_2::Edge;

		using Size_pair = std::pair<std::size_t, std::size_t>;
		using Alpha_expansion = CGAL::internal::Alpha_expansion_graph_cut_boost;

    Roof_graphcut_2(
      const FT graphcut_beta,
      const std::size_t num_roofs) : 
    m_beta(graphcut_beta),
    m_num_roofs(num_roofs) { 

      CGAL_assertion(m_num_roofs != std::size_t(-1));
    }

    void apply(Partition_2& partition) const {
      
      if (partition.empty()) return;
      auto& pfaces = partition.faces;
      auto& pedges = partition.edges;

			compute_weights(pfaces);
			compute_weights(pedges);

      std::vector<Size_pair> edges;
      std::vector<double> edge_weights;
			set_graph_edges(pedges, edges, edge_weights);

      std::vector< std::vector<double> > cost_matrix;
			set_cost_matrix(pfaces, cost_matrix);

      std::vector<std::size_t> labels;
			set_initial_labels(pfaces, labels);

      compute_graphcut(edges, edge_weights, cost_matrix, labels);
			apply_new_labels(labels, pfaces);
    }

  private:
    const FT m_beta;
    const std::size_t m_num_roofs;

    template<typename Object>
		void compute_weights(
			std::vector<Object>& objects) const {

			FT sum = FT(0);
			for (auto& object : objects) {
				object.compute_weight();
				sum += object.weight;
			}
			CGAL_assertion(sum > FT(0));
			for (auto& object : objects)
				object.weight /= sum;
		}

    void set_graph_edges(
      const std::vector<Edge>& pedges, 
      std::vector<Size_pair>& edges,
      std::vector<double>& edge_weights) const {

			edges.clear();
			edge_weights.clear();
			for (const auto& pedge : pedges) {
				
				const FT edge_weight = pedge.weight;
				const auto& neighbors = pedge.neighbors;
				
				const int idx1 = neighbors.first;
				const int idx2 = neighbors.second;

				// Boundary edges.
				if (idx1 < 0 && idx2 >= 0)
					continue;
				if (idx2 < 0 && idx1 >= 0)
					continue;

				// Internal edges.
				CGAL_assertion(idx1 >= 0);
				const std::size_t id1 = static_cast<std::size_t>(idx1);
				CGAL_assertion(idx2 >= 0);
				const std::size_t id2 = static_cast<std::size_t>(idx2);

				CGAL_assertion(edge_weight >= 0.0);
				edges.push_back(std::make_pair(id1, id2));
				edge_weights.push_back(get_graph_edge_cost(edge_weight));
			}
		}

    double get_graph_edge_cost(const FT edge_weight) const {
			return CGAL::to_double(m_beta * edge_weight);
		}

    void set_cost_matrix(
      const std::vector<Face>& pfaces,  
      std::vector< std::vector<double> >& cost_matrix) const {

			cost_matrix.clear();
			cost_matrix.resize(m_num_roofs + 1);
      for (auto& vec : cost_matrix)
        vec.resize(pfaces.size());

			for (std::size_t i = 0; i < pfaces.size(); ++i) {
        const auto& pface = pfaces[i];
        const std::size_t label = pface.label;
				const FT face_weight = pface.weight;
        
				if (label == std::size_t(-1)) {
					for (std::size_t k = 0; k < m_num_roofs; ++k)
						cost_matrix[k][i] = get_graph_face_cost(FT(0), face_weight);
					cost_matrix[m_num_roofs][i] = get_graph_face_cost(FT(1), face_weight);
					continue;
				}

				const auto& probabilities = pface.probabilities;
				CGAL_precondition(face_weight >= FT(0));
				for (std::size_t k = 0; k < m_num_roofs; ++k) {
					const FT probability = probabilities[k];
					cost_matrix[k][i] = get_graph_face_cost(probability, face_weight);
				}
				cost_matrix[m_num_roofs][i] = get_graph_face_cost(FT(0), face_weight);
			}
		}

		double get_graph_face_cost(
      const FT face_prob, const FT face_weight) const {
			
			const double weight = CGAL::to_double(face_weight);
			const double value = (1.0 - CGAL::to_double(face_prob));
      return weight * value;
		}

    void set_initial_labels(
      const std::vector<Face>& pfaces,  
      std::vector<std::size_t>& labels) const {

			labels.clear();
			labels.resize(pfaces.size());

			for (std::size_t i = 0; i < pfaces.size(); ++i) {
				const auto& pface = pfaces[i];
				const std::size_t label = pface.label;

				if (label == std::size_t(-1)) labels[i] = m_num_roofs;
				else labels[i] = label;
			}
			
			/* std::cout << "labels are set" << std::endl; */
		}

    void compute_graphcut(
      const std::vector<Size_pair>& edges,
      const std::vector<double>& edge_weights,
      const std::vector< std::vector<double> >& cost_matrix,
      std::vector<std::size_t>& labels) const {

      Alpha_expansion graphcut;
      graphcut(edges, edge_weights, cost_matrix, labels);
			
			/* std::cout << "gc computed" << std::endl; */
    }

		void apply_new_labels(
      const std::vector<std::size_t>& labels,
			std::vector<Face>& pfaces) const {

			for (std::size_t i = 0; i < labels.size(); ++i) {
				auto& pface = pfaces[i];

				if (labels[i] == m_num_roofs) pface.label = std::size_t(-1);
				else pface.label = labels[i];
			}
		}
  };

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_INTERNAL_ROOF_GRAPHCUT_2_H
