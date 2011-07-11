CGAL_BEGIN_NAMESPACE

/////////
// TDS //
/////////

template < class Vb, class Cb >
void
Triangulation_data_structure_3<Vb,Cb>::
collect_vertices_and_edges_from_link(Vertex_handle v,
                                              Vertex_handle_set& vertices,
                                              Unoriented_edge_set& edges)
    {
        std::list<Cell_handle> cells;		
        Dt::incident_cells(v, std::back_inserter(cells));
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
        m_timer_top_tests.start();
        bool test = do_is_top_collapsible(edge);
        m_timer_top_tests.stop();
        m_nb_top_tests_computed++;
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
        /*
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
        */
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
get_revolving_uedges(const Edge& edge, 
                              Unoriented_edge_set& uedges)
    {
        Vertex_handle s = source_vertex(edge);
        Vertex_handle t = target_vertex(edge);
        
        Cell_circulator cell = Dt::incident_cells(edge); 
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
any_other_vertex(
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
remaining_vertex(Cell_handle cell,
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
has_vertex(const Facet & f, Vertex_handle v) const
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

//////////////
// COLLAPSE //
//////////////

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

/////////////////////
// TRIANGULATION_3 //
/////////////////////

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
	if (!is_top_collapsible(edge))  return false;

	return true;
}

/*
template < class Kernel, class TDS >
class DT3 : public CGAL::Delaunay_triangulation_3<Kernel, TDS>
{
public:
    typedef DT3<Kernel,TDS> Dt;
    
    typedef typename Kernel::FT FT;
    typedef typename Kernel::Point_3 Point;
    typedef typename Kernel::Plane_3 Plane;
    typedef typename Kernel::Vector_3 Vector;
    typedef typename Kernel::Segment_3 Segment;
    typedef typename Kernel::Triangle_3 Triangle;
    typedef typename Kernel::Tetrahedron_3 Tetrahedron;
    typedef typename CGAL::Bbox_3 Bbox;
    
    typedef std::list<Point> Point_list;
    typedef typename Point_list::iterator Point_list_iterator;
    typedef typename Point_list::const_iterator Point_list_const_iterator;
    
    typedef typename Dt::Vertex                   Vertex;
    typedef typename Dt::Vertex_handle            Vertex_handle;
    typedef typename Dt::Vertex_iterator          Vertex_iterator;
    typedef typename Dt::Finite_vertices_iterator Finite_vertices_iterator;
    
    typedef typename Dt::Edge                  Edge;
    typedef typename Dt::Edge_iterator         Edge_iterator;
    typedef typename Dt::Finite_edges_iterator Finite_edges_iterator;
    
    typedef typename Dt::Facet                   Facet;
    typedef typename Dt::Facet_iterator          Facet_iterator;
    typedef typename Dt::Facet_circulator        Facet_circulator;
    typedef typename Dt::Finite_facets_iterator  Finite_facets_iterator;
    
    typedef typename Dt::Cell                  Cell;
    typedef typename Dt::Cell_handle           Cell_handle;
    typedef typename Dt::Cell_iterator         Cell_iterator;
    typedef typename Dt::Cell_circulator       Cell_circulator;
    typedef typename Dt::Finite_cells_iterator Finite_cells_iterator;
    
    typedef std::set< Cell_handle, less_Cell_handle<Cell_handle> > Cell_handle_set;
    typedef std::set< Vertex_handle, less_Vertex_handle<Vertex_handle> > Vertex_handle_set;    
    
    typedef typename Render_transparent<Kernel>::Projected_facet Proj_facet;
    
    typedef CSample<Kernel> Sample;
    typedef std::list<Sample*> Sample_list;
    typedef typename Sample_list::iterator Sample_list_iterator;
    typedef typename Sample_list::const_iterator Sample_list_const_iterator;
    
    typedef CSFacet<Facet, Vertex_handle, Cell_handle, Sample, Kernel> SFacet;
    typedef std::list<SFacet> SFacet_list;
    typedef std::set<SFacet> SFacet_set;
    
    typedef CUnoriented_edge<Edge, Vertex_handle> Unoriented_edge;
    typedef std::set<Unoriented_edge> Unoriented_edge_set;
    
    typedef CPEdge<Edge, Vertex_handle> PEdge;
    
    typedef CPQueue<PEdge> PQueue;
    
    typedef CTransport<Kernel> Transport;
    
    typedef CCost<Kernel> Cost;
   
   //------------------//
    // COLLAPSIBLE TEST //
    //------------------//
    
    bool is_collapsible(const Edge& edge)
    {
        if (is_edge_dummy(edge))        
        {
            std::cerr << "one dummy edge tested against is_collapsible" << std::endl;
            return false;
        }
        
        if (is_infinite(edge))          return false;
        if (is_edge_pinned(edge))       return false;
        if (!is_geom_collapsible(edge)) return false;            
        if (!is_top_collapsible(edge))  return false;
        
        return true;
    }
    
    // check if edge has both vertices pinned
    bool is_edge_dummy(const Edge& edge)
    {
        return (source_vertex(edge)->id() == target_vertex(edge)->id());
    }
    
    // check if edge has one vertex pinned
    bool is_edge_pinned(const Edge& edge)
    {
        Vertex_handle v0 = source_vertex(edge);        
        if (v0->pinned()) return true;        
        Vertex_handle v1 = target_vertex(edge);        
        if (v1->pinned()) return true;        
        return false;
    }
    
    bool is_facet_pinned(const Facet& facet)
    {
        if (facet_vertex(facet,0)->pinned()) return true;
        if (facet_vertex(facet,1)->pinned()) return true;
        if (facet_vertex(facet,2)->pinned()) return true;
        return false;
    }

    bool is_geom_collapsible(const Edge& edge)
    {
        m_timer_geom_tests.start();
        m_nb_geom_tests_computed++;
        bool test = do_is_geom_collapsible(edge);
        m_timer_geom_tests.stop();
        return test;
    }
    
    bool do_is_geom_collapsible(const Edge& edge)
    {
        std::list<Triangle> triangles;
        Vertex_handle s = source_vertex(edge);
        Vertex_handle t = target_vertex(edge);
        add_kernel_triangles_around(s, t, triangles);
        return is_visible(t->point(), triangles.begin(), triangles.end());
    }
    
    // triangles of link of s which do not go through t
    void add_kernel_triangles_around(Vertex_handle s, Vertex_handle t, 
                                     std::list<Triangle>& triangles)
    {
        assert(s != Vertex_handle());
        
        std::list<Cell_handle> cells;
        Dt::incident_cells(s, std::back_inserter(cells));
        
        typename std::list<Cell_handle>::const_iterator it;
        for (it = cells.begin(); it != cells.end(); it++)
        {
            Cell_handle cell = *it;
            int index = cell->index(s);
            
            Facet facet(cell, index);
            
            // skip infinite facets 
            if (Dt::is_infinite(facet))
                continue;
            
            // FIXME: make option
            //if (is_facet_pinned(facet)) continue;
            
            Vertex_handle va = facet_vertex(facet,0);
            Vertex_handle vb = facet_vertex(facet,1);
            Vertex_handle vc = facet_vertex(facet,2);
            if (index%2 == 0) std::swap(vb, vc);
            
            if (va->id() == t->id()) continue;
            if (vb->id() == t->id()) continue;
            if (vc->id() == t->id()) continue;
            
            const Point& pa = va->point();
            const Point& pb = vb->point();
            const Point& pc = vc->point();
            
            // prevents adding degenerate triangles
            if (pa == pb || pb == pc || pc == pa) continue;
            
            triangles.push_back(Triangle(pa, pb, pc));
        }
    }
    
    // check if a point query is inside the kernel of a set of triangles
    template < class Iterator > // value_type = Triangle
    bool is_visible(const Point& query, Iterator begin, Iterator end)
    {
        for (Iterator it = begin; it != end; it++)
        {
            const Triangle& triangle = *it;
            const Point& a = triangle[0];
            const Point& b = triangle[1];
            const Point& c = triangle[2];
            if (CGAL::orientation(a, b, c, query) == CGAL::NEGATIVE)
                return false; 
        }
        return true;
    }
};
*/

CGAL_END_NAMESPACE

