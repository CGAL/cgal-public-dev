CGAL_BEGIN_NAMESPACE

/////////////////
// TDS members //
/////////////////

template < class Vb, class Cb >
template < class Cell_handle > 
struct Triangulation_data_structure_3<Vb,Cb>::
less_Cell_handle
{
	bool operator()(const Cell_handle &s1, const Cell_handle &s2) const
	{
		Vertex_handle v1[] = {
				s1->vertex(0),
				s1->vertex(1),
				s1->vertex(2),
				s1->vertex(3) };

		Vertex_handle v2[] = {
				s2->vertex(0),
				s2->vertex(1),
				s2->vertex(2),
				s2->vertex(3) };

		std::sort(v1, v1+4);
		std::sort(v2, v2+4);

		int i;
		for(i=0; i<4; i++)
			if (v1[i] != v2[i])
				return v1[i] < v2[i];

		// equivalent cells
		return false;
	}
};

template < class Vb, class Cb >
template < class Vertex_handle > 
struct Triangulation_data_structure_3<Vb,Cb>::
less_Vertex_handle
{
	bool operator()(const Vertex_handle &s1, const Vertex_handle &s2) const
	{
		return s1 < s2;
	}
};

template < class Vb, class Cb >
template < class Edge, class Vertex_handle  >
struct Triangulation_data_structure_3<Vb,Cb>::
CUnoriented_edge
{
	std::pair<Vertex_handle, Vertex_handle> handles;

	CUnoriented_edge(Edge edge) 
	{
		Cell_handle cell = edge.first;
		int v1_index = edge.second;
		int v2_index = edge.third;

		handles.first = cell->vertex(v1_index);
		handles.second = cell->vertex(v2_index);

		if (handles.first < handles.second)
			std::swap(handles.first, handles.second);		
	}

	bool operator<(const CUnoriented_edge &another_edge) const
	{
		// different verteces are allocated on different places
		return handles < another_edge.handles;
	}	
};


template < class Vb, class Cb >
typename Triangulation_data_structure_3<Vb,Cb>::Vertex_handle
Triangulation_data_structure_3<Vb,Cb>::
facet_vertex(Facet facet, int index)
{
	Cell_handle cell = facet.first;
	return cell->vertex(index);
}

template < class Vb, class Cb >
void
Triangulation_data_structure_3<Vb,Cb>::
collect_vertices_and_edges_from_link(Vertex_handle v,
                                              Vertex_handle_set& vertices,
                                              Unoriented_edge_set& edges)
{
	std::list<Cell_handle> cells;		
	incident_cells(v, std::back_inserter(cells));
	typename std::list<Cell_handle>::iterator it;
	for (it = cells.begin(); it != cells.end(); it++)
	{
		Cell_handle cell = *it;
		int index = cell->index(v);
		Facet facet(cell, index);
		Vertex_handle v0 = facet_vertex(facet, 0);
		Vertex_handle v1 = facet_vertex(facet, 1);
		Vertex_handle v2 = facet_vertex(facet, 2);
		vertices.insert(v0);
		vertices.insert(v1);
		vertices.insert(v2);
		edges.insert(Unoriented_edge(Edge(cell, (index+1)%4, (index+2)%4)));
		edges.insert(Unoriented_edge(Edge(cell, (index+1)%4, (index+3)%4)));
		edges.insert(Unoriented_edge(Edge(cell, (index+2)%4, (index+3)%4)));
	}
}
    
template < class Vb, class Cb >
bool
Triangulation_data_structure_3<Vb,Cb>::
is_top_collapsible(const Edge& edge)
    {
        //m_timer_top_tests.start();
        bool test = do_is_top_collapsible(edge);
        //m_timer_top_tests.stop();
        //m_nb_top_tests_computed++;
        return test;
    }
    
template < class Vb, class Cb >
bool
Triangulation_data_structure_3<Vb,Cb>::
do_is_top_collapsible(const Edge& edge)
{	
	// edge pq
	Cell_handle cell = edge.first;
	Vertex_handle source = source_vertex(edge);
	Vertex_handle target = target_vertex(edge);

	// get revolving vertices and edges
	Vertex_handle_set revolving_vertices;
	Unoriented_edge_set revolving_edges; 
	get_revolving_vertices(edge, revolving_vertices);
	get_revolving_uedges(edge, revolving_edges);

	// get vertices and edges in link of source, then target
	Vertex_handle_set svertices, tvertices;
	Unoriented_edge_set sedges, tedges;
	collect_vertices_and_edges_from_link(source, svertices, sedges);
	collect_vertices_and_edges_from_link(target, tvertices, tedges);
	 
	// compute l(a) inter l(b) for vertices
	Vertex_handle_set vinter;
	typename Vertex_handle_set::iterator vit;
	for (vit = svertices.begin(); vit != svertices.end(); vit++)
	    if (tvertices.find(*vit) != tvertices.end())
		vinter.insert(*vit);
	
	// compare the two sets in size then element by element
	if (vinter.size() != revolving_vertices.size())
	    return false;
	for (vit = vinter.begin(); vit != vinter.end(); vit++)
	    if (revolving_vertices.find(*vit) == revolving_vertices.end())
		return false;

	// compute l(a) inter l(b) for edges
	Unoriented_edge_set einter;
	typename Unoriented_edge_set::iterator eit;
	for (eit = sedges.begin(); eit != sedges.end(); eit++)
	    if (tedges.find(*eit) != tedges.end())
		einter.insert(*eit);

	// compare the two sets in size then element by element
	if (einter.size() != revolving_edges.size())
	    return false;
	for (eit = einter.begin(); eit != einter.end(); eit++)
	    if (revolving_edges.find(*eit) == revolving_edges.end())
		return false;
	
	// finally
	return true;
}

template < class Vb, class Cb >
void
Triangulation_data_structure_3<Vb,Cb>::
get_revolving_vertices(const Edge& edge, Vertex_handle_set& vertices)
{
	Cell_handle seed_cell = edge.first;
	Vertex_handle s = source_vertex(edge);
	Vertex_handle t = target_vertex(edge);
        
	Cell_handle cell = seed_cell;
	Vertex_handle v = get_any_other_vertex(seed_cell, s, t);
	do {
		vertices.insert(v);
		int index = cell->index(v);
		v = get_remaining_vertex(cell, s, t, v);
		cell = cell->neighbor(index);
	} while (cell != seed_cell);
}
    
template < class Vb, class Cb >
void
Triangulation_data_structure_3<Vb,Cb>::
get_revolving_uedges(const Edge& edge, Unoriented_edge_set& uedges)
{
	Vertex_handle s = source_vertex(edge);
	Vertex_handle t = target_vertex(edge);
        
	Cell_circulator cell = incident_cells(edge); 
	Cell_circulator end = cell;
	
	CGAL_For_all(cell, end)
        {
		Vertex_handle u = get_any_other_vertex(cell, s, t);
		Vertex_handle v = get_remaining_vertex(cell, s, t, u);
		Edge e(cell, cell->index(u), cell->index(v));
		uedges.insert(Unoriented_edge(e));
	}
}

template < class Vb, class Cb >
typename Triangulation_data_structure_3<Vb,Cb>::Vertex_handle
Triangulation_data_structure_3<Vb,Cb>::
source_vertex(const Edge& edge)
{
	return edge.first->vertex(edge.second);
}

template < class Vb, class Cb >
typename Triangulation_data_structure_3<Vb,Cb>::Vertex_handle
Triangulation_data_structure_3<Vb,Cb>::
target_vertex(const Edge& edge)
{
	return edge.first->vertex(edge.third);
}

template < class Vb, class Cb >
typename Triangulation_data_structure_3<Vb,Cb>::Edge
Triangulation_data_structure_3<Vb,Cb>::
get_twin_edge(const Edge& edge)
{
	return Edge(edge.first, edge.third, edge.second);
}

template < class Vb, class Cb >
typename Triangulation_data_structure_3<Vb,Cb>::Vertex_handle
Triangulation_data_structure_3<Vb,Cb>::
get_any_other_vertex(
		Cell_handle cell,
		Vertex_handle va,
		Vertex_handle vb)
{
	int i;
	for(i=0;i<4;i++)
		if(cell->vertex(i) != va && 
			cell->vertex(i) != vb)
			return cell->vertex(i);
	assert(false); // never come here
	return Vertex_handle();
}

// return the vertex of cell different from va, vb and vc
template < class Vb, class Cb >
typename Triangulation_data_structure_3<Vb,Cb>::Vertex_handle
Triangulation_data_structure_3<Vb,Cb>::
get_remaining_vertex(Cell_handle cell,
	Vertex_handle va,
	Vertex_handle vb,
	Vertex_handle vc)
{
	int i;
	for(i=0;i<4;i++)
		if(cell->vertex(i) != va && 
			cell->vertex(i) != vb && 
			cell->vertex(i) != vc)
			return cell->vertex(i);
	assert(false); // never come here
	return Vertex_handle();
}
/*
template < class Vb, class Cb>
bool
My_triangulation_data_structure_3<Vb,Cb>::
is_simplex( Cell_handle c ) const
{
  switch(dimension()) {
    case 3 : return is_cell(c);
    case 2 : return is_facet(c, 3);
    case 1 : return is_edge(c, 0, 1);
    case 0 : return is_vertex(c->vertex(0));
    case -1 : return c == cells().begin();
  }
  return false;
}

template < class Vb, class Cb>
bool
My_triangulation_data_structure_3<Vb,Cb>::
is_vertex(Vertex_handle v) const
{
    Vertex_iterator vit = vertices_begin();
    while (vit != vertices_end() && v != vit)
        ++vit;
    return v == vit;
}

template < class Vb, class Cb>
bool
My_triangulation_data_structure_3<Vb,Cb>::
is_edge(Vertex_handle u, Vertex_handle v,
	Cell_handle &c, int &i, int &j) const
  // returns false when dimension <1 or when indices wrong
{
    CGAL_triangulation_expensive_precondition( is_vertex(u) && is_vertex(v) );

    if (u==v)
        return false;

    std::vector<Cell_handle> cells;
    cells.reserve(64);
    incident_cells(u, std::back_inserter(cells));

    for (typename std::vector<Cell_handle>::iterator cit = cells.begin();
	      cit != cells.end(); ++cit)
        if ((*cit)->has_vertex(v, j)) {
	    c = *cit;
	    i = c->index(u);
	    return true;
	}
    return false;
}

template < class Vb, class Cb>
bool
My_triangulation_data_structure_3<Vb,Cb>::
is_edge(Vertex_handle u, Vertex_handle v) const
{
    Cell_handle c;
    int i, j;
    return is_edge(u, v, c, i, j);
}

template < class Vb, class Cb>
bool
My_triangulation_data_structure_3<Vb,Cb>::
is_edge(Cell_handle c, int i, int j) const
  // returns false when dimension <1
{
  if ( i==j ) return false;
  if ( (i<0) || (j<0) ) return false;
  if ( (dimension() == 1) && ((i>1) || (j>1)) ) return false;
  if ( (dimension() == 2) && ((i>2) || (j>2)) ) return false;
  if ((i>3) || (j>3)) return false;

  for(Cell_iterator cit = cells().begin(); cit != cells_end(); ++cit)
    if (c == cit)
/	return true;
  return false;
}

template < class Vb, class Cb>
bool
My_triangulation_data_structure_3<Vb,Cb>::
is_facet(Vertex_handle u, Vertex_handle v,
	 Vertex_handle w,
	 Cell_handle & c, int & i, int & j, int & k) const
  // returns false when dimension <2 or when indices wrong
{
    CGAL_triangulation_expensive_precondition( is_vertex(u) &&
					       is_vertex(v) &&
					       is_vertex(w) );

    if ( u==v || u==w || v==w )
	return false;
    if (dimension() < 2)
	return false;

    std::vector<Cell_handle> cells;
    cells.reserve(64);
    incident_cells(u, std::back_inserter(cells));

    for (typename std::vector<Cell_handle>::iterator cit = cells.begin();
	      cit != cells.end(); ++cit)
        if ((*cit)->has_vertex(v, j) && (*cit)->has_vertex(w, k)) {
	    c = *cit;
	    i = c->index(u);
	    return true;
	}
    return false;
}

template < class Vb, class Cb>
bool
My_triangulation_data_structure_3<Vb,Cb>::
is_facet(Cell_handle c, int i) const
  // returns false when dimension <2
{
    CGAL_triangulation_precondition(i>=0 && i<4);
    if ( (dimension() == 2) && (i!=3) )
        return false;

    Cell_iterator cit = cells().begin(); // needs to work in dim 2.
    while (cit != cells_end() && c != cit)
        ++cit;
    return c == cit;
}

template < class Vb, class Cb>
bool
My_triangulation_data_structure_3<Vb,Cb>::
is_cell( Cell_handle c ) const
  // returns false when dimension <3
{
    if (dimension() < 3)
        return false;

    Cell_iterator cit = cells_begin();
    while (cit != cells_end() && c != cit)
        ++cit;
    return c == cit;
}

template < class Vb, class Cb>
bool
My_triangulation_data_structure_3<Vb,Cb>::
is_cell(Vertex_handle u, Vertex_handle v,
	Vertex_handle w, Vertex_handle t,
	Cell_handle & c, int & i, int & j, int & k, int & l) const
  // returns false when dimension <3
{
    CGAL_triangulation_expensive_precondition( is_vertex(u) &&
					       is_vertex(v) &&
					       is_vertex(w) &&
					       is_vertex(t) );

    if ( u==v || u==w || u==t || v==w || v==t || w==t )
        return false;

    std::vector<Cell_handle> cells;
    cells.reserve(64);
    incident_cells(u, std::back_inserter(cells));

    for (typename std::vector<Cell_handle>::iterator cit = cells.begin();
	      cit != cells.end(); ++cit)
        if ((*cit)->has_vertex(v, j) && (*cit)->has_vertex(w, k) &&
	    (*cit)->has_vertex(t, l)) {
	    c = *cit;
	    i = c->index(u);
	    return true;
	}
    return false;
}

template < class Vb, class Cb>
bool
My_triangulation_data_structure_3<Vb,Cb>::
is_cell(Vertex_handle u, Vertex_handle v,
	Vertex_handle w, Vertex_handle t)
    const
  // returns false when dimension <3
{
    Cell_handle c;
    int i, j, k, l;
    return is_cell(u, v, w, t, c, i, j, k, l);
}

template < class Vb, class Cb>
inline
bool
My_triangulation_data_structure_3<Vb,Cb>::
has_vertex(Cell_handle c, int i, Vertex_handle v, int & j) const
  // computes the index j of the vertex in the cell c giving the query
  // facet (c,i)
  // j has no meaning if false is returned
{
  CGAL_triangulation_precondition( dimension() == 3 );
  return ( c->has_vertex(v,j) && (j != i) );
}

template < class Vb, class Cb>
inline
bool
My_triangulation_data_structure_3<Vb,Cb>::
has_vertex(Cell_handle c, int i, Vertex_handle v) const
  // checks whether the query facet (c,i) has vertex v
{
  CGAL_triangulation_precondition( dimension() == 3 );
  int j;
  return ( c->has_vertex(v,j) && (j != i) );
}

template < class Vb, class Cb>
inline
bool
My_triangulation_data_structure_3<Vb,Cb>::
has_vertex(const Facet & f, Vertex_handle v, int & j) const
{
  return has_vertex(f.first, f.second, v, j);
}

template < class Vb, class Cb>
inline
bool
My_triangulation_data_structure_3<Vb,Cb>::
has_vertex(const dFacet & f, Vertex_handle v) const
{
  return has_vertex(f.first, f.second, v);
}

template < class Vb, class Cb>
bool
My_triangulation_data_structure_3<Vb,Cb>::
are_equal(Cell_handle c, int i, Cell_handle n, int j) const
  // tests whether facets c,i and n,j, have the same 3 vertices
  // the triangulation is supposed to be valid, the orientation of the
  // facets is not checked here
  // the neighbor relations between c and  n are not tested either,
  // which allows to use this method before setting these relations
  // (see remove in Delaunay_3)
  //   if ( c->neighbor(i) != n ) return false;
  //   if ( n->neighbor(j) != c ) return false;
{
  CGAL_triangulation_precondition( dimension() == 3 );

  if ( (c==n) && (i==j) ) return true;

  int j1,j2,j3;
  return( n->has_vertex( c->vertex((i+1)&3), j1 ) &&
	  n->has_vertex( c->vertex((i+2)&3), j2 ) &&
	  n->has_vertex( c->vertex((i+3)&3), j3 ) &&
	  ( j1+j2+j3+j == 6 ) );
}

template < class Vb, class Cb>
bool
My_triangulation_data_structure_3<Vb,Cb>::
are_equal(const Facet & f, const Facet & g) const
{
  return are_equal(f.first, f.second, g.first, g.second);
}

template < class Vb, class Cb>
bool
My_triangulation_data_structure_3<Vb,Cb>::
are_equal(const Facet & f, Cell_handle n, int j) const
{
  return are_equal(f.first, f.second, n, j);
}
*/

///////////////////
// COLLAPSE EDGE //
///////////////////

template < class Vb, class Cb >
bool
Triangulation_data_structure_3<Vb,Cb>::
collapse_edge(const Edge& edge)
{
	// edge source-target: source is removed
	Cell_handle seed_cell = edge.first;
	Vertex_handle source = seed_cell->vertex(edge.second);
	Vertex_handle target = seed_cell->vertex(edge.third);
	assert(source != target);

	// collect all cells incident to source - as this vertex will disappear
	// we will have to update source to target later in these cells.
	std::list<Cell_handle> scells;
	incident_cells(source, std::back_inserter(scells)); // precondition: dimension()==3

	// stores all cells incident to edge pq (they will be deleted)
	std::list<Cell_handle> rcells; 

	// revolve on cells around edge pq and update
	// incidence between mate cells.
	Cell_handle cell = seed_cell;
	Vertex_handle u = any_other_vertex(cell,target,source);
	
	do
	{
		// invariant: there is a circular list of cells around the edge (source,target)

		// add cell to revolving cells
		rcells.push_back(cell);

		// get mate cells
		Vertex_handle v = remaining_vertex(cell,target,source,u);
		Cell_handle cell_suv = cell->neighbor(cell->index(target));
		Cell_handle cell_tuv = cell->neighbor(cell->index(source));
		Vertex_handle mt = mirror_vertex(cell,cell->index(target));
		Vertex_handle ms = mirror_vertex(cell,cell->index(source));

		// update incidence relationships
		int index_mt = cell_suv->index(mt);
		int index_ms = cell_tuv->index(ms);

		assert(!cell_suv->has_neighbor(cell_tuv));
		assert(!cell_tuv->has_neighbor(cell_suv));

		if(cell_suv->has_neighbor(cell_tuv) || cell_tuv->has_neighbor(cell_suv))
		{
			std::cerr << "ERROR: non-manifold edge collapse - must DELETE the triangulation as it is invalid" << std::endl; 
			return false;
		}

		cell_suv->set_neighbor(index_mt, cell_tuv);
		cell_tuv->set_neighbor(index_ms, cell_suv);

		// set incident cell to ensure that vertices do not point
		// to a deleted cell (among scells).
		u->set_cell(cell_tuv);
		target->set_cell(cell_tuv);

		//std::cout << cell->vertex(0)->id() << " "
		//					<< cell->vertex(1)->id() << " "
		//					<< cell->vertex(2)->id() << " "
		//					<< cell->vertex(3)->id() << std::endl;

		// revolves to neighboring cell
		cell = cell->neighbor(cell->index(u));
		u = v;
	}
	while(cell != seed_cell);

	// update cells incident to source
	typename std::list<Cell_handle>::iterator it;
	for(it = scells.begin();
		it != scells.end();
		it++)
	{
		Cell_handle cell = *it;
		cell->set_vertex(cell->index(source),target);
	}

	// delete cells and vertex source
	delete_cells(rcells.begin(),rcells.end());
	delete_vertex(source);

	std::cerr << "Edge collapsed!" << std::endl;

	return true;
}

/////////////////////////////
// TRIANGULATION_3 members //
/////////////////////////////

template < class GT, class Tds >
bool
Triangulation_3<GT,Tds>::
is_collapsible(const Edge& edge)
{
//	if (is_edge_dummy(edge))        
//	{
//	    std::cerr << "one dummy edge tested against is_collapsible" << std::endl;
//	    return false;
//	}

	if (is_infinite(edge))          return false;
//	if (is_edge_pinned(edge))       return false;
//	if (!is_geom_collapsible(edge)) return false;            
	if (!_tds.is_top_collapsible(edge))  return false;

	return true;
}

CGAL_END_NAMESPACE

