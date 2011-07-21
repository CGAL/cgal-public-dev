// the functions from this file will be copied to the TDS and Tri_3

#include<cmath>

CGAL_BEGIN_NAMESPACE

#define TDS_FUNC(type) \
	template < class Vb, class Cb >\
	type\
	Triangulation_data_structure_3<Vb,Cb>::

#define TRI_FUNC(type) \
	template < class Vb, class Cb >\
	type\
	Triangulation_3<Vb,Cb>::

/////////////////
// TDS structs //
/////////////////

template < class Vb, class Cb >
template < class Cell_handle > 
struct Triangulation_data_structure_3<Vb,Cb>::
less_Cell_handle
{    
    void get_id(const Cell_handle& cell, Vertex_handle& a0, Vertex_handle& a1, Vertex_handle& a2, Vertex_handle& a3) const
    {
        a0 = cell->vertex(0);
        a1 = cell->vertex(1);
        a2 = cell->vertex(2);
        a3 = cell->vertex(3);
    }
    
    bool operator()(const Cell_handle& a, const Cell_handle& b) const
    {
        Vertex_handle a0, a1, a2, a3;
        get_id(a, a0, a1, a2, a3);
        
        Vertex_handle b0, b1, b2, b3;
        get_id(b, b0, b1, b2, b3);
        
        if (a0 < b0) return true;
        if (a0 > b0) return false;
        
        if (a1 < b1) return true;
        if (a1 > b1) return false;
        
        if (a2 < b2) return true;
        if (a2 > b2) return false;
        
        if (a3 < b3) return true;
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
template < class Edge >
struct Triangulation_data_structure_3<Vb,Cb>::
less_Edge
{
    void get_id(const Edge& edge, Vertex_handle& a, Vertex_handle& b) const
    {
        a = edge.first->vertex(edge.second);
        b = edge.first->vertex(edge.third );
        if ( !(a < b) ) std::swap(a, b);
    }
    
    bool operator()(const Edge& a, const Edge& b) const
    {
        std::pair<Vertex_handle, Vertex_handle> a_handle, b_handle;

        get_id(a, a_handle.first, a_handle.second);
        get_id(b, b_handle.first, b_handle.second);
        
        return a_handle < b_handle;
    }
};

/////////////////
// TDS members //
/////////////////

template < class Vb, class Cb >
typename Triangulation_data_structure_3<Vb,Cb>::Vertex_handle
Triangulation_data_structure_3<Vb,Cb>::
facet_vertex(Facet facet, int index)
{
	Cell_handle cell = facet.first;
	return cell->vertex(index);
}

template < class Vb, class Cb >
typename Triangulation_data_structure_3<Vb,Cb>::Vertex_handle
Triangulation_data_structure_3<Vb,Cb>::
get_source_vertex(const Edge& edge) const
{
	return edge.first->vertex(edge.second);
}

template < class Vb, class Cb >
typename Triangulation_data_structure_3<Vb,Cb>::Vertex_handle
Triangulation_data_structure_3<Vb,Cb>::
get_target_vertex(const Edge& edge) const
{
	return edge.first->vertex(edge.third);
}
/*
template < class Vb, class Cb >
typename Triangulation_data_structure_3<Vb,Cb>::Edge
Triangulation_data_structure_3<Vb,Cb>::
get_twin_edge(const Edge& edge)
{
	return Edge(edge.first, edge.third, edge.second);
}
*/
template < class Vb, class Cb >
typename Triangulation_data_structure_3<Vb,Cb>::Vertex_handle
Triangulation_data_structure_3<Vb,Cb>::
get_any_other_vertex(
		Cell_handle cell,
		Vertex_handle va,
		Vertex_handle vb) const
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
	Vertex_handle vc) const
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

template < class Vb, class Cb >
typename Triangulation_data_structure_3<Vb,Cb>::Vertex_handle
Triangulation_data_structure_3<Vb,Cb>::
get_vertex_from_facet(const Facet& facet, int index) const
{
        return facet.first->vertex( (facet.second + (index%3) + 1) % 4 );
}

TDS_FUNC(void) get_vertices_from_facet(const Facet& facet, 
                                Vertex_handle& va,
                                Vertex_handle& vb,
                                Vertex_handle& vc) const
{
	va = get_vertex_from_facet(facet, 0);
	vb = get_vertex_from_facet(facet, 1);
	vc = get_vertex_from_facet(facet, 2);
	
	if (facet.second % 2 == 0)
		std::swap(vb, vc);
}

// GET LINK //
    
TDS_FUNC(void)
get_vertices_from_edge_link(const Edge& edge, Vertex_handle_set& vertices) const
{
	Vertex_handle s = get_source_vertex(edge);
	Vertex_handle t = get_target_vertex(edge);
	Cell_circulator ccirc = incident_cells(edge); 
	Cell_circulator cend  = ccirc;

	CGAL_For_all(ccirc, cend)
	{
		Cell_handle cell = ccirc;
		Vertex_handle u = get_any_other_vertex(cell, s, t);
		Vertex_handle v = get_remaining_vertex(cell, s, t, u);
		vertices.insert(u);
		vertices.insert(v);
	}
}
    
TDS_FUNC(void)
get_edges_from_edge_link(const Edge& edge, Edge_set& edges) const
{
	Vertex_handle s = get_source_vertex(edge);
	Vertex_handle t = get_target_vertex(edge); 
	Cell_circulator ccirc = incident_cells(edge); 
	Cell_circulator cend  = ccirc;
	CGAL_For_all(ccirc, cend)
	{
		Cell_handle cell = ccirc;
		Vertex_handle u = get_any_other_vertex(cell, s, t);
		Vertex_handle v = get_remaining_vertex(cell, s, t, u);
		Edge e(cell, cell->index(u), cell->index(v));
		edges.insert(e);
	}
}
    
TDS_FUNC(void)
get_vertices_from_vertex_link(Vertex_handle vertex, Vertex_handle_set& vertices) const
{
	std::list<Cell_handle> cells;		
	incident_cells(vertex, std::back_inserter(cells));
	
	typename std::list<Cell_handle>::const_iterator it;
        for (it = cells.begin(); it != cells.end(); it++)
        {
		Cell_handle cell = *it;
		Vertex_handle v0, v1, v2;
		Facet facet(cell, cell->index(vertex));
		get_vertices_from_facet(facet, v0, v1, v2);
		vertices.insert(v0);
		vertices.insert(v1);
		vertices.insert(v2);
	}
}
        
TDS_FUNC(void)
get_edges_from_vertex_link(Vertex_handle vertex, Edge_set& edges) const
{
	std::list<Cell_handle> cells;		
	incident_cells(vertex, std::back_inserter(cells));
	typename std::list<Cell_handle>::const_iterator it;
	for (it = cells.begin(); it != cells.end(); it++)
	{
		Cell_handle cell = *it;
		Vertex_handle v0, v1, v2;
		Facet facet(cell, cell->index(vertex));
		get_vertices_from_facet(facet, v0, v1, v2);
		int i0 = cell->index(v0);
		int i1 = cell->index(v1);
		int i2 = cell->index(v2);
		edges.insert( Edge(cell, i0, i1) );
		edges.insert( Edge(cell, i1, i2) );
		edges.insert( Edge(cell, i2, i0) );
	}
}
/*
template < class Vb, class Cb >
typename Triangulation_data_structure_3<Vb,Cb>::Facet
Triangulation_data_structure_3<Vb,Cb>::
get_twin_facet(const Facet& facet) const
{        
	Cell_handle c  = facet.first;
	Cell_handle nc = c->neighbor(facet.second);
	return Facet(nc, nc->index(c));
}
*/ 
TDS_FUNC(void)
get_facets_from_link(Vertex_handle vertex, Facet_list& hull, bool outward = false) const
{
	std::list<Cell_handle> cells;
	incident_cells(vertex, std::back_inserter(cells));
	
	typename std::list<Cell_handle>::const_iterator it;
	for (it = cells.begin(); it != cells.end(); ++it)
	{
		Cell_handle cell = *it;

		Facet facet(cell, cell->index(vertex));
		//if (outward) facet = get_twin_facet(facet);
		hull.push_back(facet);
	}
}

TDS_FUNC(bool)
check_link_test(const Edge& edge, int verbose = 0) const
{
	return 	check_link_test_for_vertices(edge, verbose)	
		&& check_link_test_for_edges(edge, verbose);
}
    
TDS_FUNC(bool)
check_link_test_for_vertices(const Edge& edge, int verbose = 0) const
{
	Cell_handle cell = edge.first;
	Vertex_handle source = get_source_vertex(edge);
	Vertex_handle target = get_target_vertex(edge);

	Vertex_handle_set svertices;
	get_vertices_from_vertex_link(source, svertices);

	Vertex_handle_set tvertices;
	get_vertices_from_vertex_link(target, tvertices);

        Vertex_handle_set ivertices;
	typename Vertex_handle_set::const_iterator vit;
	for (vit = svertices.begin(); vit != svertices.end(); ++vit)
		if (tvertices.find(*vit) != tvertices.end())
			ivertices.insert(*vit);
                
	Vertex_handle_set evertices;
	get_vertices_from_edge_link(edge, evertices);

	if (ivertices.size() != evertices.size())
		return false;
 
	for (vit = ivertices.begin(); vit != ivertices.end(); ++vit)
		if (evertices.find(*vit) == evertices.end())
			return false;
        
	return true;
}
    
TDS_FUNC(bool)
check_link_test_for_edges(const Edge& edge, int verbose = 0) const
{
	Cell_handle cell = edge.first;

	Vertex_handle source, target;
	source = get_source_vertex(edge);
	target = get_target_vertex(edge);

	Edge_set sedges, tedges;
	get_edges_from_vertex_link(source, sedges);
	get_edges_from_vertex_link(target, tedges);

#ifdef PESHO_DEBUG
	print_edges(sedges, "sedges");
	print_edges(tedges, "tedges");
#endif        

	Edge_set iedges;
	typename Edge_set::const_iterator eit;
	for (eit = sedges.begin(); eit != sedges.end(); ++eit)
		if (tedges.find(*eit) != tedges.end())
			iedges.insert(*eit);
        
	Edge_set eedges;
	get_edges_from_edge_link(edge, eedges);
        
	if (iedges.size() != eedges.size()) {
		//std::cout << "#1";
		return false;
	}

	for (eit = iedges.begin(); eit != iedges.end(); ++eit)
		if (eedges.find(*eit) == eedges.end()) {
			//std::cout << "#2";
			return false;
		}

        return true;        
}

///////////////////
// COLLAPSE EDGE //
///////////////////

TDS_FUNC(bool)
collapse_edge(Edge& edge)
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
	//std::cout << "scells.size() = " << scells.size() << std::endl;

	// stores all cells incident to edge pq (they will be deleted)
	std::list<Cell_handle> rcells; 

	// revolve on cells around edge pq and update
	// incidence between mate cells.
	Cell_handle cell = seed_cell;
	Vertex_handle u = get_any_other_vertex(cell,target,source);
	
	do {
		// add cell to revolving cells
		rcells.push_back(cell);

		// get mate cells
		Vertex_handle v = get_remaining_vertex(cell,target,source,u);
		Cell_handle cell_suv = cell->neighbor(cell->index(target));
		Cell_handle cell_tuv = cell->neighbor(cell->index(source));
		Vertex_handle mt = mirror_vertex(cell,cell->index(target));
		Vertex_handle ms = mirror_vertex(cell,cell->index(source));

		// update incidence relationships
		int index_mt = cell_suv->index(mt);
		int index_ms = cell_tuv->index(ms);

		assert(!cell_suv->has_neighbor(cell_tuv));
		assert(!cell_tuv->has_neighbor(cell_suv));

		if(cell_suv->has_neighbor(cell_tuv) || cell_tuv->has_neighbor(cell_suv)) {
			std::cerr << "ERROR: non-manifold edge collapse - must DELETE the triangulation as it is invalid" << std::endl; 
//			return false;
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
	} while(cell != seed_cell);

	// update cells incident to source
	typename std::list<Cell_handle>::iterator it;
	for(it = scells.begin(); it != scells.end(); it++) {
		Cell_handle cell = *it;
		cell->set_vertex(cell->index(source),target);
	}

	// delete cells and vertex source
	delete_cells(rcells.begin(),rcells.end());
	delete_vertex(source);

	//std::cerr << "Edge collapsed!" << std::endl;

	return true;
}

/////////////////////////////
// TRIANGULATION_3 members //
/////////////////////////////

TRI_FUNC(bool)
check_kernel_test(const Edge& edge) const
{
        Vertex_handle s = _tds.get_source_vertex(edge);
        Vertex_handle t = _tds.get_target_vertex(edge);

	CGAL_assertion( !infinite_vertex(s) );
	CGAL_assertion( !infinite_vertex(t) );
	        
        Facet_list hull;
        _tds.get_facets_from_link(s, hull);
        _tds.get_facets_from_link(t, hull); // not needed

        return is_in_kernel(t, hull.begin(), hull.end());
}
  
TRI_FUNC(bool)
check_kernel_test(const Edge& edge, const Point& point) const
{
	Vertex_handle s = _tds.get_source_vertex(edge);
        Vertex_handle t = _tds.get_target_vertex(edge);

	CGAL_assertion( !infinite_vertex(s) );
	CGAL_assertion( !infinite_vertex(t) );
	        
        Facet_list hull;
        _tds.get_facets_from_link(s, hull);
        _tds.get_facets_from_link(t, hull);

        return is_in_kernel(t, hull.begin(), hull.end());
}  

template < class Vb, class Cb >
template < class Iterator > // value_type = Face
bool
Triangulation_3<Vb,Cb>::
is_in_kernel(Vertex_handle query, Iterator begin, Iterator end) const
{
	for (Iterator it = begin; it != end; ++it) {
		Facet facet = *it;
		if (is_infinite(facet))
			continue;
		
		Vertex_handle va, vb, vc;
		_tds.get_vertices_from_facet(facet, va, vb, vc);
	
		if (CGAL::orientation(query->point(), va->point(), vb->point(), vc->point()) == CGAL::POSITIVE)
			return false; 
	}

	return true;
}

///////////
// DEBUG //
///////////

template < class Vb, class Cb >
template < class Cont >
void
Triangulation_data_structure_3<Vb,Cb>::
print_cells(const Cont S, std::string note) const
{
	std::cout << "Cells, " << note << std::endl;

	typename Cont::iterator cit;
	for (cit = S.begin(); cit != S.end(); cit++) {
		Cell_handle c = *cit;

		Vertex_handle v0 = c->vertex(0);
		Vertex_handle v1 = c->vertex(1);
		Vertex_handle v2 = c->vertex(2);
		Vertex_handle v3 = c->vertex(3);

		std::cout << "[  (" << v0->point() << ")   " << std::endl;
		std::cout << "   (" << v1->point() << ")   " << std::endl;
		std::cout << "   (" << v2->point() << ")   " << std::endl;
		std::cout << "   (" << v3->point() << ")  ]" << std::endl;
	}

	std::cout << std::endl;
}	

template < class Vb, class Cb >
template < class Cont >
void
Triangulation_data_structure_3<Vb,Cb>::
print_vertices(const Cont S, std::string note) const
{
	std::cout << "Vertices, " << note << std::endl;

	typename Cont::iterator vit;
	for (vit = S.begin(); vit != S.end(); vit++) {
		Vertex_handle v = *vit;
		std::cout << "(" << v->point() << ") ";
	}

	std::cout << std::endl;
}	

template < class Vb, class Cb >
template < class Cont >
void
Triangulation_data_structure_3<Vb,Cb>::
print_edges(const Cont S, std::string note) const
{
	std::cout << "Edges, " << note << std::endl;

	typename Cont::iterator eit;
	for (eit = S.begin(); eit != S.end(); eit++) {
		Edge e = *eit;
		std::cout << "(" << get_source_vertex(e)->point() << ") -> (" << get_target_vertex(e)->point() << ") " << std::endl;
	}

	std::cout << std::endl;
}	

CGAL_END_NAMESPACE

