/*!
\ingroup PkgTopologicalInvariantConcepts
\cgalConcept

\cgalRefines Boost::Graph

The concept GraphOnSurface is a refinement of HalfedgeGraph and adds the notion of a graph embedding.

##Notation##
<dl>
<dt>`G`</dt> 	  <dd>A type that is a model of `GraphOnSurface`.</dd>
<dt>`g`</dt> 	  <dd>An object of type `G`.</dd>
<dt>`e`</dt> 	  <dd>An edge descriptor.</dd>
<dt>`p`</dt>	  <dd>A path</dd>
</dl>

##Associated Types##
Type											|	Description
------------------------------------------------|-----------------------------------------------------------------
boost::graph_traits<G>::path_descriptor			| a path

##Valid Expressions##
Expression							|		Returns										|		Description											
------------------------------------|---------------------------------------------------|-----------------------------------------------------------
path(g, e)							| path_descriptor									| the path of a specified edge
in_graph(g, path_descriptor)		| bool												| return true if the path is in the graph

*/
class GraphOnSurface{}; 
