/*!
\ingroup PkgTopologicalInvariantConcepts
\cgalConcept

The concept GeneralizedMapGraph is a refinement of the concept of HalfedgeGraph Graph. this concept add a signature for each edge.

##notation##

G : A type that is a model of GeneralizedMapGraph.

g : An object of type G. 

he : An halfedge descriptor. 

d : A dart descriptor

p : a path descriptor a model of Path

e : a embeded graph a model of GraphOnSurface

##valide Expressions##
Expression		|		Returns		|		Description											|
----------------|-------------------|-----------------------------------------------------------|
is_flip(he)    	| bool				| ...                                                       |
flip(d)         | void              | ...                                                       |
signature(d)    | bool              | ...                                                       |
halfedge(d)     | he                | ...                                                       |
Order(d)        | ?                 | ...                                                       |
Order(he)       | ?                 | ...                                                       |
create_path(g)  | p                 | ...                                                       |                      
erase_path(p)   | void              | ...                                                       |
create_graph(g) | e                 | ...                                                       |
erase_graph(e)  | void              | ...                                                       |
dart(he)        | d                 | ...                                                       |

*/
class GeneralizedHalfedgeGraph{}; 
