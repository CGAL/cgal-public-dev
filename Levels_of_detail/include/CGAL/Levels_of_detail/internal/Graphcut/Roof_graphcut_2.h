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

			compute_face_weights(pfaces);
			compute_edge_weights(pfaces, pedges);

      std::vector<Size_pair> edges;
      std::vector<double> edge_weights;
			set_graph_edges(pfaces, pedges, edges, edge_weights);

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

		void compute_face_weights(
			std::vector<Face>& pfaces) const {

			FT sum = FT(0);
			for (auto& pface : pfaces) {
				if (pface.visibility == Visibility_label::OUTSIDE)
					continue;
				
				pface.compute_weight();
				sum += pface.weight;
			}
			
			std::size_t count = 0;
			CGAL_assertion(sum > FT(0));
			for (auto& pface : pfaces) {
				if (pface.visibility == Visibility_label::OUTSIDE)
					continue;
			
				pface.weight /= sum;
				pface.index = count; 
				++count;
			}
		}

		void compute_edge_weights(
			const std::vector<Face>& pfaces,
			std::vector<Edge>& pedges) const {

			FT sum = FT(0);
			for (auto& pedge : pedges) {
				const auto& neighbors = pedge.neighbors;
				
				const int idx1 = neighbors.first;
				const int idx2 = neighbors.second;
					
				if (idx1 < 0 && idx2 >= 0)
					continue;
				if (idx2 < 0 && idx1 >= 0)
					continue;	

				CGAL_assertion(idx1 >= 0);
				const std::size_t id1 = static_cast<std::size_t>(idx1);
				CGAL_assertion(idx2 >= 0);
				const std::size_t id2 = static_cast<std::size_t>(idx2);

				if (pfaces[id1].visibility != pfaces[id2].visibility)
					continue;

				if (
					pfaces[id1].visibility == Visibility_label::OUTSIDE &&
					pfaces[id2].visibility == Visibility_label::OUTSIDE)
					continue;

				pedge.compute_weight();
				sum += pedge.weight;
			}

			CGAL_assertion(sum > FT(0));
			for (auto& pedge : pedges) {
				const auto& neighbors = pedge.neighbors;
				
				const int idx1 = neighbors.first;
				const int idx2 = neighbors.second;
					
				if (idx1 < 0 && idx2 >= 0)
					continue;
				if (idx2 < 0 && idx1 >= 0)
					continue;	

				CGAL_assertion(idx1 >= 0);
				const std::size_t id1 = static_cast<std::size_t>(idx1);
				CGAL_assertion(idx2 >= 0);
				const std::size_t id2 = static_cast<std::size_t>(idx2);	

				if (pfaces[id1].visibility != pfaces[id2].visibility)
					continue;

				if (
					pfaces[id1].visibility == Visibility_label::OUTSIDE &&
					pfaces[id2].visibility == Visibility_label::OUTSIDE)
					continue;

				pedge.weight /= sum;
			}
		}

    void set_graph_edges(
			const std::vector<Face>& pfaces,
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

				if (pfaces[id1].visibility != pfaces[id2].visibility)
					continue;

				if (
					pfaces[id1].visibility == Visibility_label::OUTSIDE &&
					pfaces[id2].visibility == Visibility_label::OUTSIDE)
					continue;

				CGAL_assertion(edge_weight >= 0.0);
				edges.push_back(std::make_pair(
					pfaces[id1].index, pfaces[id2].index));
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
			cost_matrix.resize(m_num_roofs);

			std::size_t count = 0;
			for (const auto& pface : pfaces) {
				if (pface.visibility == Visibility_label::OUTSIDE)
					continue;
				++count;
			}

			for (auto& vec : cost_matrix)
				vec.resize(count);	

			count = 0;
			for (const auto& pface : pfaces) {
				if (pface.visibility == Visibility_label::OUTSIDE)
					continue;

				const FT face_weight = pface.weight;
				CGAL_precondition(face_weight >= FT(0));
				const auto& probabilities = pface.probabilities;
				
				for (std::size_t k = 0; k < m_num_roofs; ++k) {
					const FT probability = probabilities[k];
					cost_matrix[k][count] = 
						get_graph_face_cost(probability, face_weight);
				}
				++count;
			}
		}

		double get_graph_face_cost(
      const FT face_prob, const FT face_weight) const {
			
			const double weight = CGAL::to_double(face_weight);
			const double value  = (1.0 - CGAL::to_double(face_prob));
      return weight * value;
		}

    void set_initial_labels(
      const std::vector<Face>& pfaces,  
      std::vector<std::size_t>& labels) const {

			labels.clear();
			for (const auto& pface : pfaces) {
				if (pface.visibility == Visibility_label::OUTSIDE) 
					continue;
				labels.push_back(pface.label);
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

			std::size_t count = 0;
			for (auto& pface : pfaces) {
				if (pface.visibility == Visibility_label::OUTSIDE)
					continue;
				pface.label = labels[count]; ++count;
			}
		}
  };

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_INTERNAL_ROOF_GRAPHCUT_2_H
