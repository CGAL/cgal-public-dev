// the functions from this file will be copied to the TDS and Tri_3

CGAL_BEGIN_NAMESPACE

/////////////////
// TDS members //
/////////////////

// TODO: check how to identify the vertex without its address

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

		std::cout << "less_Cell_handle: " << v1[0]->point() << " " << v1[1]->point() << " " << v1[2]->point() << " " << v1[3]->point() << std::endl;

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

		// compares the addresses of the vertices
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

	//print_cells< std::list<Cell_handle> >(cells, "Incident to v");

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

	// DEBUG	
	print_vertices<Vertex_handle_set>(svertices, "svertices");
	print_vertices<Vertex_handle_set>(tvertices, "tvertices");
	print_vertices<Vertex_handle_set>(vinter, "vinter");
	print_vertices<Vertex_handle_set>(revolving_vertices, "revolving_vertices");

	// compare the two sets in size then element by element
	if (vinter.size() != revolving_vertices.size())
	    return false;
	for (vit = vinter.begin(); vit != vinter.end(); vit++)
	    if (revolving_vertices.find(*vit) == revolving_vertices.end())
		return false;
/*
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

// check if edge has both vertices pinned
template < class Vb, class Cb >
bool
Triangulation_data_structure_3<Vb,Cb>::
is_edge_dummy(const Edge& edge)
{
	return (source_vertex(edge) == target_vertex(edge));
}

///////////////////
// COLLAPSE EDGE //
///////////////////

template < class Vb, class Cb >
bool
Triangulation_data_structure_3<Vb,Cb>::
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

	// stores all cells incident to edge pq (they will be deleted)
	std::list<Cell_handle> rcells; 

	// revolve on cells around edge pq and update
	// incidence between mate cells.
	Cell_handle cell = seed_cell;
	Vertex_handle u = any_other_vertex(cell,target,source);
	
	do
	{
		// invariant: here is a circular list of cells around the edge (source,target)

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
is_geom_collapsible(const Edge& edge)
{
	std::list<Triangle> triangles;
	
	Vertex_handle s = source_vertex(edge);
	Vertex_handle t = target_vertex(edge);
	add_kernel_triangles_around(s, t, triangles);
	
	return is_visible(t->point(), triangles.begin(), triangles.end());
}

template < class GT, class Tds >
bool
Triangulation_3<GT,Tds>::
is_geom_collapsible(const Edge& edge, const Point& point)
{
	// TODO
	return true;	
}
    
// triangles of link of s which do not go through t
template < class GT, class Tds >
void
Triangulation_3<GT,Tds>::
add_kernel_triangles_around(Vertex_handle s, Vertex_handle t, 
			     std::list<Triangle>& triangles)
{
	assert(s != Vertex_handle());

	std::list<Cell_handle> cells;
	incident_cells(s, std::back_inserter(cells));

	typename std::list<Cell_handle>::const_iterator it;
	for (it = cells.begin(); it != cells.end(); it++)
	{
		Cell_handle cell = *it;
		int index = cell->index(s);

		Facet facet(cell, index);

		// skip infinite facets 
		if (is_infinite(facet))
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
template < class GT, class Tds >
template < class Iterator > // value_type = Triangle
bool
Triangulation_3<GT,Tds>::
is_visible(const Point& query, Iterator begin, Iterator end)
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

///////////
// DEBUG //
///////////

template < class Vb, class Cb >
template < class Cont >
void
Triangulation_data_structure_3<Vb,Cb>::
print_cells(Cont S, std::string note)
{
	return;
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
print_vertices(Cont S, std::string note)
{
	return;
	std::cout << "Vertices, " << note << std::endl;

	typename Cont::iterator vit;
	for (vit = S.begin(); vit != S.end(); vit++) {
		Vertex_handle v = *vit;
		std::cout << "(" << v->point() << ") ";
	}

	std::cout << std::endl;
}	

CGAL_END_NAMESPACE

