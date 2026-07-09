#include <boost/graph/graph_traits.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>


template <typename FaceGraph, typename OutputIterator>
typedef boost::graph_traits<FaceGraph>::vertex_descriptor vertex_descriptor;
void min_length_homotopy_basis(const FaceGraph& mesh,
												vertex_descriptor basepoint,
												OutputIterator basis_out) {
													
}
