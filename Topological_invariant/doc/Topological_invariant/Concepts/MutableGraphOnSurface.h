/*!
\ingroup PkgTopologicalInvariantConcepts
\cgalConcept

\cgalRefines Boost::GraphOnSurface

The concept MutableGraphOnSurface is a refinement of GraphOnSurface and adds the requirements for operations to add vertices and edges, and to update the incidence information between vertices and halfedges.

##Notation##
<dl>
<dt>`G`</dt> 	        <dd>A type that is a model of `MutableGraphOnSurface`.</dd>
<dt>`g`</dt> 	        <dd>An object of type `G`.</dd>
<dt>`e`</dt> 	        <dd>An edge descriptor.</dd>
<dt>`v, v1, v2`</dt>    <dd>A vertex descriptor.</dd>
<dt>`p`</dt>	        <dd>A path</dd>
</dl>

##Valid Expressions##
Expression							|		Returns										|		Description											
------------------------------------|---------------------------------------------------|--------------------------------------------
add_vertex(g)    		        	| vertex_descriptor									| Adds a new vertex to the graph
add_edge(g, v1, v2, p)              | edge_descriptor                                   | Adds a new edge beteween v1 and v2 with the path p
remove_vertex(g, v)                 | void                                              | remove a vertex
remove_edge(g, e)                   | void                                              | remove a edge
*/
class MutableGraphOnSurface{}; 
